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
            [smith-waterman :as sw])
  (:use [clojure.contrib.core :only [dissoc-in]]
        hts-exploration.globals
        hts-exploration.db-queries
        hts-exploration.hts-utils.file
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.sequtils.files
        [edu.bc.bio.seq-utils :only (markov-step)]
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

(defn str-partition [n s] (map #(apply str %) (partition n 1 s)))

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

(defn rand-exponential [lambda] (/ (- (Math/log (- 1 (rand)))) lambda))
(defn rand-geometric [p] (clojure.contrib.math/floor
                          (/ (Math/log (- 1 (rand)))
                             (Math/log (- 1 p)))))

(defn kmer-freqn [n inseqs]
  (->> (map #(freqn n %) inseqs)
       (apply merge-with +)
       (reduce-kv (fn [M k v]
                    (assoc M k (double (/ v (count inseqs)))))
                  {})))

(defn melt-temp
  "Finds the melting temp of a short seq of length len using the 2-4
  rule. Only returns those seq with temp >= 60. "

  [s]
  (let [m {\G 4 \C 4 \T 2 \A 2}
        scorefn (fn [s] (reduce + (map m s)))
        calc-temp (fn [len s]
                    (->> (str-partition len s)
                         (map-indexed (fn [i s] [(inc i) (scorefn s) len s]))
                         (filter #(>= (second %) 60))))]
    (->> (for [len (range 19 25)] (calc-temp len s))
         (remove empty?)
         first)))

(defn const-start
  "attempts to find the start of the constant region in a sequence"
  
  [s]
  (let [dists (map-indexed (fn [i x] [i (levenshtein prime3-const x)])
                           (str-partition 20 s))
        m (apply min (map second dists))]
    (-> (drop-while #(not= (second %) m) dists)
        ffirst inc)))

(defn rand-sequence
  ([n len]
     (rand-sequence n len (probs 1 "ACGT")))
  ([n len probs]
     (let [variable (fn [] (apply str (repeatedly len #(markov-step probs))))]
       (repeatedly n variable))))

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
       (reduce (fn [M {:keys [id hit_seq cs]}]
                 (let [sq (second
                           (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" hit_seq))
                       cval (get M sq [])]
                   (if sq (assoc M sq (conj cval id)) M)))
               {})))
