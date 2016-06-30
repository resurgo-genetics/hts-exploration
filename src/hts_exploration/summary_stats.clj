(ns hts-exploration.summary-stats
  (:require [clojure.contrib.string :as str]
            [clojure.core.reducers :as r]
            [iota])
  (:use edu.bc.bio.seq-utils
        edu.bc.utils
        edu.bc.utils.probs-stats
        hts-exploration.db-queries
        hts-exploration.globals
        hts-exploration.utils))


;;;get lengths of usable seqs in this file
(comment ;deprecated stuff
  (defn- all-usable-seqs
  "Gets all usable seqs in the file which only contain standard bases
  and the sequences of interest either in the forward or reverse
  direction"
  
  [f]
  (->> (iota/vec f)
       read-fasta-csv
       (r/map second)
       (r/remove #(re-find #"[^ACGTU]" %))
       (r/map #(or (re-find forward-re %)
                   (re-find reverse-re %))) 
       (r/remove nil?)))
  
  (defn -main []
    (let [infiles [S15-r1-qual-filter-csv S15-round10-qual-filter-csv S15-round11-qual-filter-csv]]
      (doseq [f infiles
              :let [f (all-usable-seqs f)
                    ldist (length-distribution f)
                    total-usable (reduce + (vals ldist))
                    total-parasites (count-parasites f)]]
        (prn ldist)
        (prn :total-usable total-usable)
        (prn :total-parasites total-parasites)
        (prn :fraction-parasite (double (/ total-parasites total-usable))))))

  )

(defn- count-parasites [inseqs]
  (let [rev-parasite (str (reverse-compliment cds) "T+"
                          (reverse-compliment prime5-const))
        parasite? (re-pattern (str parasite "|" rev-parasite))
        parasites (->> (r/filter #(re-find parasite? %) inseqs)
                       (into []))]
    (count parasites)))

(defn- length-distribution [inseqs]
  (frequencies (r/map count inseqs)))




(defn find-const-start [round]
  (letfn [(get-seqs [round]
            "Gets sequences from the database and returns the seq with the ID."
            (->> (round-all-usable-seqs round) sql-query
                 (reduce (fn [M {:keys [id hit_seq]}]
                           (let [sq (second
                                     (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" hit_seq))
                                 cval (get M sq [])]
                             (if sq (assoc M sq (conj cval id)) M)))
                         {})))
          (const-start [s]
            (let [dists (map-indexed (fn [i x] [i (levenshtein prime3-const x)])
                                     (str-partition 20 s))
                  m (apply min (map second dists))]
              (-> (drop-while #(not= (second %) m) dists)
                  ffirst inc)))]
    (let [sqs (get-seqs round)
          cnt (count sqs)]
      (->> (pxmap (fn  ([[s ids]]
                        (let [c (const-start s)]
                          (interleave ids (repeat c)))))
                  20 sqs)
           flatten
           (partition-all 2)
           doall))))

(defn var-base-rate [nmer round]
  (letfn [(get-seqs [round]
            (->> (round-all-usable-seqs round) sql-query
                 (reduce (fn [M {:keys [id hit_seq cs]}]
                           (let [sq (second
                                     (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" hit_seq))]
                             (if sq (assoc M sq cs) M)))
                         {})))]
    (let [sqs (get-seqs round)
          cnt (count sqs)]
      (->> (pxmap (fn [[s cs]] (->> s (str/take cs) (probs nmer))) 10 sqs)
           (apply merge-with +)
           (reduce (fn [M [k v]]
                     (assoc M k (double (/ v cnt))))
                   {})))))

(defn var-base-distribution [nmer round]
  (letfn [(get-seqs [round]
            (->> (round-all-usable-seqs round) sql-query
                 (reduce (fn [M {:keys [id hit_seq cs]}]
                           (let [sq (second
                                     (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" hit_seq))]
                             (if sq (assoc M sq cs) M)))
                         {})))]
    (let [sqs (get-seqs round)
          cnt (count sqs)]
      (->> (pxmap #(freqn nmer (str/take (-> % second) (first %))) 10 sqs)
           (apply merge-with +)
           probs))))

(defn var-motif-ratio
  "Finds the ratio of mean motif rate in the variable region for 2 separate rounds (ra,
  rb). Reports the top 10 hits"

  [n ra rb]
  (->> (merge-with / (var-base-rate n ra)
                   (var-base-rate n rb))
       (sort-by second > )
       (take 10)))

(defn find-base-pairs
  "Takes a sequence and structure and returns the base pairing positions"
  
  [s st]
  (loop [s s
         st st
         i (range (count s))
         stack []
         bps []]
    (if (seq s)
      (let [[stack bps] (case (first st)
                          \( [(conj stack (first i)) bps];open bp
                          \. [ stack bps] ;gap do nothing
                          \) [(pop stack) (conj bps [(peek stack) (first i)])])];close bp
        (recur (rest s) (rest st) (rest i) stack bps))
      bps)))

(defn xconj
  "Special purpose tree. Tree represents the secondary structure of an
  RNA sequence. This function is for adding base/base-pairs to the
  tree"
  
  [t b1 & [b2]]
  (cond
   (nil? t) {:b1 b1 :next nil :b2 b2 }
   :else {:b1 (t :b1)
          :next (xconj (t :next) b1 b2)
          :b2 (t :b2)}))

(defn xseq
  "traverse the tree in order and return the nodes"
  
  [t]
  (when t
    (->> (concat [(t :b1)] (xseq (t :next)) [(t :b2)])
         (remove nil?))))

(defn adjacent-GUGC?
  "Checks a structure for adjacent GUGC motif, which is known to be
  involved in binding to S15"
  
  [bp-set s]
  (some (fn [[i j]]
          (let [si (.charAt s i)
                sj (.charAt s j)
                sk (str/get s (inc i))
                sl (str/get s (dec j))]
            (and
             (contains? bp-set [(inc i) (dec j)]);k-l bp must exist
             (= \C si) (= \G sj)
             (= \T sk) (= \G sl))))
        bp-set))

(defn xwalk
  "lie xseq but stops the walk if adjacent GUGC is found"

  [t]
  (when t (if (adjacent-GUGC? t) :true (xwalk (t :next)))))

(let [s  "GGAAAACCACCAAAAGG"
      st "((....)).((....))"]
  (let [M (reduce (fn [M bp]
                    (let [k (or (->> (keys M)
                                     (filter (fn [[i j]] (->> [i bp j]
                                                             flatten
                                                             (apply <))))
                                     first)
                                bp)
                          v (get M k [])]
                      (assoc M k (conj v bp))))
                  (sorted-map)
                  (sort-by first (find-base-pairs s st)))
        bp-list (fn [M s]
                  (prn M)
                  (loop [k (range (count s))
                         foo []]
                    (if (seq k)
                      (let [i (first k)
                            j (M i)]
                        (cond
                          (nil? j) ;;left-bulge
                          (recur (rest k) (conj foo [(.charAt s i)]))
                        
                          (< j (last k)) ;right-bulge
                          (recur (butlast k)
                                 (conj foo [nil (.charAt s (last k))]))
                        
                          :else ;;base pairing
                          (recur (-> k rest butlast)
                                 (conj foo [(.charAt s i) (.charAt s j)]))))
                      foo)))
        build-tree (fn [L]
                     (reduce (fn [t nts]
                               (apply xconj t nts))
                             nil L))]
    (prn M)
    (map (fn [s p]
           (-> (bp-list (into {} p) s)
               build-tree
               ))
         (map (fn [[i j]] (let [s (subs s i (inc j))]
                           (prn :s s) s))
              (as-> (keys M) x
                (sort-by first x)
                (mapv first x)
                (conj x (->> (keys M) flatten (apply max)))
                (partition 2 1 x)))
         (vals M))
    ))

(let [s  "GGAGAAAACCCACCCAAAACAGG"
                                     st "((.(....))).(((....).))"]
  (let [M (reduce (fn [M bp]
                    (let [k (or (->> (keys M)
                                     (filter (fn [[i j]] (->> [i bp j]
                                                            flatten
                                                            (apply <))))
                                     first)
                                bp)
                          v (get M k [])]
                      (assoc M k (conj v bp))))
                  (sorted-map)
                  (sort-by first (find-base-pairs s st)))
        bp-list (fn [M s mini maxj]
                  (prn :M M)
                  (loop [k (range mini maxj)
                         foo []]
                    (if (seq k)
                      (let [i (first k)
                            j (M i)]
                        (cond
                         (nil? j) ;;left-bulge
                         (recur (rest k) (conj foo [(.charAt s i)]))
                         
                         (< j maxj) ;right-bulge
                         (recur (butlast k)
                                (conj foo [nil (.charAt s (last k))]))
                         
                         :else ;;base pairing
                         (recur (-> k rest butlast)
                                (conj foo [(.charAt s i) (.charAt s j)]))))
                      foo)))
        build-tree (fn [L]
                     (reduce (fn [t nts]
                               (apply xconj t nts))
                             nil L))]
    (prn M)
    (prn s)
    (map (fn [[lbound rbound] p]
           (prn :l lbound :r rbound :p p)
           (-> (bp-list (into {} p) s lbound rbound)
               ;build-tree
               prn))
         (as-> (keys M) x
               (sort-by first x)
               (mapv first x)
               (conj x (->> (keys M) flatten (apply max)))
               (if (zero? (first x)) x (apply vector 0 (rest x)))
               (partition 2 1 x))
         (vals M))
    ))

(let [s  "GGAGAAAACCCACTCAAAACAGG"
      st "((.(....))).(((....).))"]
  (let [M (set (find-base-pairs s st))
        adjacent-GUGC? (fn [bp-list s]
                         (some (fn [[i j]]
                                 (let [si (.charAt s i)
                                       sj (.charAt s j)
                                       sk (str/get s (inc i))
                                       sl (str/get s (dec j))]
                                   (and (= \C si) (= \G sj)
                                        (= \T sk) (= \G sl))))
                               bp-list))]
    (adjacent-GUGC? M s)))

