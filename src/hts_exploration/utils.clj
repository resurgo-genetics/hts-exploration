(ns hts-exploration.utils
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :as math]
            [clojure.java.jdbc :as jdbc]
            edu.bc.utils.probs-stats
            gibbs-sampler
            [smith-waterman :as sw])
  (:use [clojure.contrib.core :only [dissoc-in]]
        hts-exploration.globals
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.sequtils.files
        [edu.bc.bio.seq-utils :only (markov-step)]
        edu.bc.utils.fold-ops))

(def mysql-ds
  {:classname "com.mysql.jdbc.Driver"
   :subprotocol "mysql"
   :subname "//127.0.0.1:3306/s15selexdata"
   :user "root"
   :password "rna314rulz"})

(defn sql-query [stmt & {:keys [f p] :or {f identity p false}}]
  (let [q (partial jdbc/query mysql-ds)]
    (when p (println stmt))
    (cond (string? stmt) (q [stmt])
          (vector? stmt) (q stmt)
          :else "Invalid query")))

(defn read-aln
  "read clustalW file and makes it unblocked"

  [aln & {:keys [info]
          :or {info :data}}]
  (let [seq-lines (second (join-sto-fasta-lines aln ""))
        sl (reduce (fn [v [nm [_ sq]]]
                     (let [sq (.toUpperCase sq)]
                       (if (re-find #"[^ACGTU-]" sq)
                         v
                         (case info
                           :both (conj v [nm sq])
                           :name (conj v nm)
                           :data (conj v sq)))))
                   [] seq-lines)]
    (->> sl rest )))

(defn hts-name-parts
  "Breaks the name line from a fastq file into the components and
  stores it in a map"

  [nm]
  (let [[x y] (str/split #" " nm)
        xx (interleave [:inst-name :run-id
                        :flowcell-id :flowcell-lane
                        :flowcell-tile-num
                        :xcoord :ycoord]
                       (str/split #":" x))
        yy (when-not (nil? y)
             (interleave [:pair :filter?
                          :control-bits :index-seq]
                         (str/split #":" y)))]
    (reduce (fn [m [k v]]
              (assoc m k v))
            {} (partition-all 2 (concat xx yy)))))

(defn column-align
  "prints out the variable region from reads with the high scoring
  motif in caps and everything else in lowercase. Also aligns the
  sequences around the high scoring motif"

  [s prob-vec & {:keys [range]}]
  (let [s (if range (apply subs s range) s)
        likely-seq (second (gibbs-sampler/score s prob-vec))
        [s1 s2] (str/split (re-pattern likely-seq) 2 s)
        [c1 c2] (map count [s1 s2])
        lflank (take (- 40 c1) (repeat "-"))  ;left flank
        rflank (take (- 80 (count lflank) (count s)) (repeat "-"))] ;right
    (->> [lflank
          (str/lower-case (if (empty? s1) "" s1))
          likely-seq
          (str/lower-case (if (empty? s2) "" s2))
          rflank]
         (map #(apply str %))
         (str/join ""))))


(defn partition-into [n coll]
  (partition-all (/ (count coll) n) coll))

(defn hamming-dist [s t]
  (->> (map #(if (= %1 %2) 0 1) s t)
       (reduce +)))
  
(defn normalized-dist [st1 st2]
  (/ (bpdist st1 st2 :bpdist true)
     (count st1)))

(defn relative-dist [bpdist len] (/ bpdist len))

(defn rpartition-all [n coll]
  (loop [p coll
         V []]
    (if (seq p)
      (recur (drop n p)
             (conj V (r/take n p)))
      V)))

(defn write-fasta
  "Write out a fasta given sequences. If inseq is not a pair then
  numerical names are given for each sequence."

  [outfile inseqs & {:keys [info]
                     :or {info :both}}]
  (cond
   (= info :both)
   (io/with-out-writer outfile
     (doseq [[nm s] inseqs]
       (println (str ">" nm))
       (println s)))
   (= info :data)
   (write-fasta outfile
                (->> (interleave (range (count inseqs)) inseqs)
                     (partition-all 2))))
  outfile)

(defn read-fasta-csv [f]
  (if (seq? f)
    (map #(str/split #"," 3 %) f)
    (r/map #(str/split #"," 3 %) f)))

(defn get-usable
  "Returns a vector of seqs which are in the forward direction"
  
  [f]
  (->> (iota/vec f)
       read-fasta-csv
       (r/reduce (fn [V [nm inseq]]
                   (let [sq (re-find forward-re inseq)]
                     (if sq
                       (conj V [nm sq])
                       V)))
                 [])))

(defn get-db-seq [query]
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

(defn mutant-neighbor [n inseq & type]
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

(defn str-regex-index [re s]
  (let [[start x] (first (str-re-pos re s))]
    [start (+ start (count x))]))

(defn starting-gaps [s]
  (if-let [gaps (re-find #"^\-+" s)] (count gaps) 0))

(defn- -align
  [alignfn inseqs]
  (->> inseqs
       (r/map alignfn)
       (r/fold concat
               (fn ([] [])
                 ([V [_ s1 s2]]
                    (let [d (starting-gaps s1)]
                      (conj V [(str/drop d s1) (str/drop d s2)])))))))

(defn align-to-const [inseqs]
  (let [const prime3-const
        get-str (fn [s]
                  (let [s (->> s (str-take-last 52) (str/butlast 15))
                        rtstr (str-take-last 15 s)
                        re #"AATGTC"]
                    (if-let [start-codon (re-find re rtstr)]
                      (str/butlast  (- 15 (.indexOf rtstr start-codon)) s)
                      s)))
        to-const (fn [x] (first (sw/sw const x :gap-open-weight -2
                                      :gap-ext-weight -1 :anchor-right true)))]
    (->> inseqs
         (r/map get-str)
         (-align to-const))))

(defn align-to-primer [inseqs]
  (let [primers (str prime5-const cds)
        to-primer (fn [x] (first (sw/sw primers x :gap-open-weight -2)))]
    (->> inseqs
         (r/map #(str (str/take 15 %) (str-take-last 15 %)))
         (-align to-primer))))

(defn count-subs

  [s1 s2]
  (->> (map (fn [x1 x2]
              (if (or (= x1 \-)
                      (= x2 \-))
                0
                (if (= x1 x2) 0 1)))
            s1 s2)
       (reduce +)))

(defn count-indel [s1 s2]
  (let [helper (fn [x] (map count (re-seq #"\-+" x)))
        ins-len (helper s1)
        del-len (helper s2)]
    [(count ins-len) ;insertion
     ins-len
     (count del-len) ;deletion
     del-len]))

(defn rand-exponential [lambda] (/ (- (Math/log (- 1 (rand)))) lambda))
(defn rand-geometric [p] (clojure.contrib.math/floor
                          (/ (Math/log (- 1 (rand)))
                             (Math/log (- 1 p)))))

(defn count-mutations
  "Counts the number of subs,ins, and dels when comparing each set of
  aligned sequences. Returns a collection of lists containing
  substition rate, insertion rate, insertion lengths, deletion rate,
  deletion lengths for each aligned seq pair."
  
  [aligned-seqs]
  (map cons (map #(apply count-subs %) aligned-seqs)
       (map #(apply count-indel %) aligned-seqs)))

(defn get-mutation-rates
  "Combines the mutation rates for each individual pair together to
  get a mean mutation rate. Returns a vector of mutation rates and
  indel lengths as [[sub insertion deletion] [insertion deletion]]."

  [counts]
  (let [[subs ins ins-lens dels del-lens] (transpose counts)
        ins-lens (probs 1 (flatten ins-lens))
        del-lens (probs 1 (flatten del-lens))
        overall (->> (transpose [subs ins dels]) (map #(apply + %)) mean)
        totals (map #(apply + %) [subs ins dels])
        ]
    [(mapv mean [subs ins dels])
     [ins-lens del-lens]
     overall
     (->> (map #(double (/ % (apply + totals))) totals)
          (interleave [:sub :ins :del])
          (apply hash-map))]))

(defn mutation-matrix [aligned-seqs]
  (let [mutations (fn [s1 s2]
                    (->> (map (fn [x1 x2]
                                (when (and (not= x1 x2)
                                           (not= x2 \-))
                                  [x1 x2])) 
                              s1 s2)
                         (remove nil?)))
        init-transitions (as-> (repeat {\A 0.1 \C 0.1 \G 0.1 \T 0.1}) x
                               (interleave [\A \C \G \T \-] x)
                               (apply assoc {} x)
                               (reduce #(dissoc-in %1 %2);remove self transitions
                                       x (map #(repeat 2 %) [\A \C \G \T])))]
    (->> (mapcat #(apply mutations %) aligned-seqs)
         (reduce (fn [M [wt mut]] 
                   (assoc-in M [wt mut] (inc (get-in M [wt mut]))))
                 init-transitions)
         (reduce (fn [M [k m]] (assoc M k (probs m))) {}))))

(defn mut-locations
  "Generates locations for any mutations based on the assumption that
  mutations follow a Poisson distribution. Here the positions are
  determined by using the exponential distribution to predict the
  number of nucleotides until the next mutation. Returns a collection
  of positions which will not exceed the length."

  [rate len]
  (->> (repeatedly #(rand-exponential rate))
       (map math/floor)
       (reductions +)
       (take-while #(< % len) )))

(defn rand-mutations [s rate] (-> (mut-locations rate (count s)) reverse vec))

(defmulti ^{:doc "Makes mutations of a type to the a sequence"
            :arglists '([s loc mutation-rates mutant-type])}

  add-mutation (fn [s loc mutation-rates mutant-type] mutant-type))

(defmethod add-mutation :sub [s loc Mrates mutant-type]
  (let [{tmatrix :tmatrix} Mrates
        replacement (->> (str/get s loc) tmatrix markov-step)]
    (str-replace-at loc replacement s)))

(defmethod add-mutation :ins [s loc Mrates mutant-type]
  (let [{:keys [lrates tmatrix]} Mrates
        ins-len (-> lrates first markov-step)]
    (as-> ins-len n
          (repeatedly n #(markov-step (tmatrix \-)))
          (apply str n)
          (str/lower-case n)
          (str-insert-at loc n s))))

(defmethod add-mutation :del [s loc Mrates mutant-type]
  (-> (Mrates :lrates) second markov-step ;del length
      (str-remove-at loc s)))

(defn simulate [start-seq mutation-rates]
  (let [[srate irate drate] (mutation-rates :mrates)
        mutate (fn [type s loc] 
                 (add-mutation s loc mutation-rates type))
        add-del (partial mutate :del)
        add-subs (partial mutate :sub)
        add-ins (partial mutate :ins)]
    (as-> start-seq st
          (reduce add-del st (rand-mutations st drate))
          (reduce add-subs st (rand-mutations st srate))
          (reduce add-ins st (rand-mutations st irate)))))

(defn simulate2 [start-seq mutation-rates]
  (let [mrate (mutation-rates :overall)
        totals (mutation-rates :totals)
        mutate (fn [[s last-pos] loc type] 
                 (let [loc (if (>= loc (count s))
                             (dec (count s)) loc)
                       ;(prn :s s :loc loc :type type)
                       newstr (add-mutation s loc mutation-rates type)]
                                        
                   [newstr loc]))
        locs (->> (rand-mutations start-seq mrate)
                  (map (fn [loc] [loc (markov-step totals)]))
                  #_(group-by second))
        add-mutation (fn [st type] (apply mutate st type))]
    (first (reduce #(apply mutate %1 %2) [start-seq -1] locs))))

(defn kmer-freqn [n inseqs]
  (->> (map #(freqn n %) inseqs)
       (apply merge-with +)
       (reduce-kv (fn [M k v]
                    (assoc M k (double (/ v (count inseqs)))))
                  {})))
