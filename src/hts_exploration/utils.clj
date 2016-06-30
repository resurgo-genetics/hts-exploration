(ns hts-exploration.utils
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :as math]
            iota
            [clojure.java.jdbc :as jdbc]
            edu.bc.utils.probs-stats
            gibbs-sampler
            [smith-waterman :as sw]
            [edu.bc.bio.seq-utils :refer (reverse-compliment markov-step generate-rand-seq)])
  (:use [clojure.contrib.core :only [dissoc-in]]
        hts-exploration.globals
        hts-exploration.db-queries
        hts-exploration.hts-utils.file
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.sequtils.files
        edu.bc.utils.fold-ops))



(defn partition-into [n coll]
  (partition-all (/ (count coll) n) coll))

(defn hamming-dist [s t]
  (->> (map #(if (= %1 %2) 0 1) s t)
       (reduce +)))
  
(defn normalized-dist [st1 st2]
  (/ (bpdist st1 st2 :bpdist true)
     (count st1)))

(defn relative-dist [bpdist len] (/ bpdist len))


(defn get-db-seq
  "more generic than get-seqs"

  [query]
  (->> (sql-query query)
       (mapv (fn [m] (subs (m :sequence) (m :usable_start) (m :usable_stop))))))

(defn struct->vector [st] (map #(if (= % \.) 0 1) st))

(defn standard-chars? [s] (nil? (re-find #"[^ACGTU]" s)))

(defn qual-good? [thr s] (every? #(>= (- (int %) 33) thr) (seq s)))

(defn qual-remove
  "Removes sequences which have a base with lower than 20 quality
  = (- (int 5) 33)"

  [data & {:keys [thr] :or {thr 20}}]
  (r/remove (complement #(qual-good? thr (last %))) data))

(defn parasite? [s]
  (let [parasite (str "(" (str/join ")|(" parasite) ")")]
    (re-find (re-pattern parasite) s)))

(defn parasite-rev? [s]
  (let [parasite (str "(" (str/join ")|(" parasite-rev) ")")]
    (re-find (re-pattern parasite) s)))

(defn parasite-remove
  "remove sequencing containing a parasite"

  [data]
  (r/remove #(-> % second parasite?) data))

(defn calc-dist
  "Calculates the base pair distance of sequences in the dataset from
  a reference sequence"
  
  [n ref-struct data]
  (let [tmpfiles (repeatedly n #(fs/tempfile))
        write-fasta (fn [tmp lines] (io/with-out-writer tmp
                                     (println ">ref-struct")
                                     (println ref-struct)
                                     (doseq [[nm s st] lines]
                                       (println (str ">" nm))
                                       (println (str s "\n" st)))) 
                      tmp)]
    (->> (partition-into n data)
         (map #(write-fasta %1 %2) tmpfiles)
         (pmap bpdist-fasta)
         (apply merge))))

(defn dist-filter
  "Keeps sequences which fold into a structure similar to the
  reference structure. Returns a reducible"

  [n thr ref-struct data]
  (let [dist-map (calc-dist n ref-struct data)
        len (count ref-struct)]
    (r/filter #(< (relative-dist (dist-map (first %)) len) thr) data)))

(defn add-structure
  "Adds structure to a fasta type entry which has both name and
  seq. Will execute n RNAfold operations on fasta files. Returns a
  vector of name, seq, structure."

  [n data]
  (let [tmpfiles (repeatedly n #(fs/tempfile))]
    (as-> (into [] data) d
          (partition-into n d)
          (map #(write-fasta %1 %2) tmpfiles d)
          (pmap fold-fasta d)
          (apply concat d)
          (mapv (fn [x st] (conj x st)) data d))))

(defn round [n x]
  (/ (int (* x (Math/pow 10 n)))
     (Math/pow 10 n)))

(defn mutant-neighbor
  "Creates a list of n-mutant neighbors of type :RNA or :DNA. :RNA is
  used by default."

  [n inseq & type]
  (let [b (case (first type)
            :RNA [\A \C \G \U]
            :DNA [\A \C \G \T]
            [\A \C \G \U]);rna default
        fun (fn fun [n init inseq]
              (for [i (range init (count inseq))
                    mut b
                    :let [cur-base (str/get inseq i)
                          cur-seq (str (str/take i inseq) mut (subs inseq (inc i)))]
                    :when (not= cur-base mut)]
                (if (= n 1)
                  cur-seq
                  (fun (dec n) i cur-seq))))]
    (->> (fun n 0 inseq)
         flatten
         distinct
         (remove #(= inseq %)))))

(defn distinct-hits [hits]
  (-> (into {} hits)
      clojure.set/map-invert 
      clojure.set/map-invert))

(defn str-partition [n s] (map #(apply str %) (partition n 1 s)))

(defn str-partition-all [n s] (map #(apply str %) (partition-all n 1 s)))

(defn str-take-last [n s] (apply str (take-last n s)))

(defn str-insert-at [n insertion s]
  (str (str/take n s) insertion (subs s n)))

(defn str-remove-at [n start s]
  (str (str/take start s) (str/drop (+ start n) s)))

(defn str-replace-at [n replacement s]
  (str (str/take n s) replacement (subs s (inc n))))

(defn str-re-pos [re s]
  (loop [m (re-matcher re s)
         res (sorted-map)]
    (if (.find m)
      (recur m (assoc res (.start m) (.group m)))
      res)))

(defn seq-get-middle 
  "Cuts the primers off the sequence"
  [s]
  (subs s 15 (- (count s) 15)))

(defn rand-exponential [lambda] (/ (- (Math/log (- 1 (rand)))) lambda))
(defn rand-geometric [p] (math/floor
                          (/ (Math/log (- 1 (rand)))
                             (Math/log (- 1 p)))))

(defn kmer-freqn [n inseqs]
  (->> (map #(freqn n %) inseqs)
       (apply merge-with +)
       (reduce-kv (fn [M k v]
                    (assoc M k (double (/ v (count inseqs)))))
                  {})))

(defn kmer-probs
  "Takes a length k and a set of sequences and computes the mean
  probability of seeing the kmer in a sequence."
  
  [k inseqs]
  (->> (map #(probs k %) inseqs)
       (apply merge-with +)
       (reduce (fn [M [kmer v]]
                 (assoc M kmer (double (/ v (count inseqs)))))
               {})))

(defn kmer-distance [pdist qdist]
  (math/sqrt
   (sum
    (map (fn [i] (sqr
                 (- (get pdist i 0)
                    (get qdist i 0))))
         (clojure.set/union
          (set (keys pdist))
          (set (keys qdist)))))))

(defn melt-temp
  "Finds the melting temp of a short seq of length len using the 2-4
  rule. Only returns those seq with temp >= 60. "

  [s]
  (let [s (str/upper-case s)
        m {\G 4 \C 4 \T 2 \A 2}
        scorefn (fn [s] (reduce + (map m s)))
        calc-temp (fn [len s]
                    (->> (str-partition len s)
                         (map-indexed (fn [i s] [(inc i) (scorefn s) len s]))
                         (filter #(>= (second %) 60))))]
    (->> (for [len (range 19 25)] (calc-temp len s))
         (remove empty?)
         first)))

(defn design-primer
  "Takes a sequence and designs primers for it so that it can be
  created by annealing the 2 primers together. These primers need to
  overlap ~20 bases in the middle. "
  
  [s]
  (let [mid (quot (count s) 2)
        near-mid? (fn [flen] (<  (- mid 10) flen (+ mid 10)))
        overlap (last
                 (take-while (fn [x]
                               (let [pos (first x)]
                                 (< pos mid)))
                             (melt-temp s)))
        flen (+ (overlap 0) (overlap 2))];forward primer length
    (vec
     (concat [:seq s :overlap overlap]
             (if (near-mid? flen)
               [:f (str/take flen s) :r (reverse-compliment (str/drop (overlap 0) s))]
               [:f (str/take (+ mid 10) s) :r (reverse-compliment (str/drop (- mid 10) s))])))))

(defn const-start
  "attempts to find the start of the constant region in a sequence"
  
  [s]
  (if (> (count s) 20)
    (let [dists (map-indexed (fn [i x] [i (levenshtein prime3-const x)])
                            (str-partition 20 s))
         m (apply min (map second dists))]
     (-> (drop-while #(not= (second %) m) dists)
         ffirst inc))
    1))

(defn xseq-start
  "attempts to find the start of the constant region in a sequence"
  
  [xseq s]
  (let [xlen (count xseq)
        num-windows (cond (> xlen 20) 7
                          (> xlen 15) 6
                          :else 5)
        windowsize (- xlen num-windows)]
    (if (> (count s) windowsize)
     (let [dists (map-indexed (fn [i x] [i (levenshtein xseq x)])
                              (str-partition windowsize s))
           m (apply min (map second dists))]
       (-> (drop-while #(not= (second %) m) dists)
           ffirst inc))
     1)))

(def const-start (partial xseq-start prime3-const))

(defn rand-sequence-set
  "Generate a set of n random sequences of length len."
  ([n len]
   (repeatedly n #(generate-rand-seq len (probs 1 "ACGT"))))
  ([n len probs]
   (repeatedly n #(generate-rand-seq len probs))))

(defn rand-background [n]
  (map #(str prime5-const
             %
             prime3-const
             cds)
       (rand-sequence n 30)))

(defn get-seqs
  "Get sequences from the database by round number. Returns a map of
  the sequences' variable+const region (k) and the the IDs which also have
  the same sequence (v).  "
  
  [round]
  (->> (round-all-usable-seqs round) sql-query
       (map (juxt :id #(->> % :hit_seq (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT")) :cs))
       (filter #(if-let [cs (% 2)] (<= 28 cs 32)))
       (filter second)
       (reduce (fn [M [id [hit_seq sq] cs]]
                 (let [k sq
                       cval (get M k [])]
                   (assoc M k (conj cval id))))
               {}))) 

(defn get-seqs-variable
  "Get sequences from the database by round number. Returns a map of
  the sequences' variable region (k) and the the IDs which also have
  the same sequence (v).  "
  
  [round]
  (->> (round-all-usable-seqs round) sql-query
       (map (juxt :id #(->> % :hit_seq (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT")) :cs))
       (filter #(if-let [cs (% 2)] (<= 28 cs 32)))
       (filter second)
       (reduce (fn [M [id [hit_seq sq] cs]]
                 (let [k (str/take cs sq)
                       cval (get M k [])]
                   (assoc M k (conj cval id))))
               {})))

(defn parse-ps [f]
  (with-open [infile (clojure.java.io/reader f)]
    (let [ex-file (line-seq infile)
          sq (->> (drop-while (complement #(= % "/sequence { (\\")) ex-file)
                  second
                  (str/replace-re #"U" "T")
                  (str/butlast 1)
                  (str/drop 15)
                  (str/butlast 15))
          bprobs (drop-while (complement #(= % "%start of base pair probability data")) ex-file)
          ]
      [sq
       (fs/basename f)
       (->> (filter #(re-find #"ubox" %) bprobs)
            (map #(vec (str/split #" " %)))
            (reduce (fn [M x]
                      (let [i (read-string (x 0))
                            j (read-string (x 1))
                            p (read-string (x 2))]
                        (assoc M [i j] (* p p))))
                    {})
            )])))

(defn parse-lunp
  "reads a lunp file and returns a hash-map of the probabilities that
  the interval [i-x+1 .. i] is unpaired. The NA's are also read in and
  they can be found by comparing to (symbol NA)."

    [f]
    (with-open [infile (clojure.java.io/reader f)]
      (let [ex-file (->> infile
                         line-seq
                         (drop 2)
                         (map #(str/split #"\t" %))
                         (map (partial map read-string) ))]; strings -> numbers
        (reduce (fn [M [i & line]]
                  ;; make hmap of current line then add it to existing hmap M
                  (into M (map-indexed (fn [j pij]
                                         [[(- i j) i] pij]); [i-x+1 .. i] prob
                                       line)))
                {} ex-file))))

(defn ensemble-dist
  "Calculate the distance between 2 ensembles x and y. "
  
  [x y]
  (->> (set (concat (keys x) (keys y))) ;S
       (mapv (fn [k]
               (let [diff (- (get x k 0)
                             (get y k 0))] ;(PAij-PBij)
                 (* diff diff))));sq
       (reduce +)));sum

(defn pairwise 
  "takes a collection and applies the function f to the elements in a
  pairwise fashion. f takes 2 args."

  [f coll]
  (loop [x coll
         Vx []]
    (if (seq x)
      (recur (rest x)
             (concat Vx
                     (reduce (fn [Vy y]
                               (conj Vy (f (first x) y)));inner loop
                             [] (vec (rest x)))))
      Vx)))

(defn add-letters
  "Takes a set of letters and creates all the combinations of it with n letters"
  
  [n xs]
  (letfn [(helper [n i cur]
            (if (< i n)
              (let [new (for [c cur x xs] (str c x))]
                (helper n (inc i) new))
              cur))]
    (helper n 1 xs)))

(defn freqs-id-seq
  "Takes id-seq pairs and returns a map where k=sequence and the value
  is a list of ids which map to that particular sequence."
  [id-seq]
  (r/fold (fn ([] {})
            ([& m] (apply merge-with (comp vec concat) m)))
          (fn ([] {})
            ([M [id s]] (assoc M s (conj (get M s []) id))))
          id-seq))

(defn intersect-maps 
  "Finds the intersection of maps according to keys and concatenates
  all the values associated with the keys. This is different than
  merge because only intersecting keys are retained in a new map."
  [& maps]
  (reduce (fn [M x]
            (assoc M x (map #(get % x) maps)))
          {}
          (->> (map keys maps)
               (map set)
               (apply clojure.set/intersection))))

(defn vfold
  "Fold function f over a collection or set of collections (c1, ...,
   cn) producing a collection (concrete type of vector).  Arity of f
   must be equal to the number of collections being folded with
   parameters matching the order of the given collections.  Folding
   here uses the reducer lib fold at base, and thus uses work stealing
   deque f/j to mostly relieve the partition problem.  In signatures
   with N given, N provides the packet granularity, or if (< N 1),
   granularity is determined automatically (see below) as in the
   base (non N case) signature.
   While vfold is based on r/fold, it abstracts over the combiner,
   reducer, work packet granularity, and transforming multiple
   collections for processing by f by chunking the _transpose_ of the
   collection of collections.
   As indicated above, vfold's fold combiner is monoidal on vectors:
   it constructs a new vector from an l and r vector, and has identity
   [] (empty vector).  In line with this, vfold's reducer builds up
   new vectors from elements by conjing (f a1, ... an) onto a prior
   result or the combiner identity, [], as initial result.
   Packet granularity is determined automatically (base case or N=0)
   or supplied with N > 1 in signatures with N.  Automatic
   determination tries to balance significant work chunks (keep thread
   overhead low), with chunks that are small enough to have multiple
   instances per worker queue (supporting stealing by those workers
   that finish ahead of others).
  "
  ([f coll]
     (let [cores (.. Runtime getRuntime availableProcessors)
           workers (int (math/floor (* 3/4 cores)))
           base (* 8 cores)
           n (max 2 (int (math/floor (/ (count coll) (* 2 workers)))))
           n (min base n)]
       #_(println :>>N n)
       (vfold f n coll)))
  ([f n coll]
     (assert (integer? n)
             (str "VFOLD: Folding granularity N " n " must be an integer."))
     (if (< n 1)
       (vfold f coll)
       (r/fold n
               (fn
                 ([] [])
                 ([l r] (apply conj l r)))
               (fn[v x] (conj v (f x)))
               (vec coll))))
  ([f n coll & colls]
     (vfold (fn[v] (apply f v)) n (apply transpose coll colls))))

(defn euclidean-dist
  "calculate the euclidean distance between 2 vectors."
  
  [v1 v2]
  (->> (map (comp sqr -) v1 v2); sq difference
       sum
       Math/sqrt))

(defn select-vals
  "get a collection of values associated with keyseq from a
  map (m). Values are returned in the same order as the keyseq."
  [m keyseq]
  (map m keyseq))

(def concatv (comp vec concat))

(defn rotations
  "Returns a lazy seq of all rotations of a seq"
  [x]
  (if (seq x)
    (map
     (fn [n _]
       (lazy-cat (drop n x) (take n x)))
     (iterate inc 0) x)
    (list nil)))

(defn deep-merge-with
  "Copied verbatim from the defunct clojure-contrib"
  [f & maps]
  (apply
   (fn m [& maps]
     (if (every? map? maps)
       (apply merge-with m maps)
       (apply f maps)))
   maps))

(defn deep-probs [M]
  (if (not-any? map? (vals M))            ; no maps left? no recuring
    (probs M)                             ; base case
    (zipmap (keys M)
            (map deep-probs (vals M))))) ; recursive call
