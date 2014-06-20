(ns hts-exploration.summary-stats
  (:require [clojure.core.reducers :as r]
            [iota])
  (:use edu.bc.bio.seq-utils
        edu.bc.utils.probs-stats
        hts-exploration.globals
        hts-exploration.utils))


(defn- str-take [n s] (apply str (take n s)))

(defn- str-take-last [n s] (apply str (take-last n s)))

;;;get lengths of usable seqs in this file
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
       (r/remove nil?)
       ))



(defn- count-parasites [inseqs]
  (let [rev-parasite (str (reverse-compliment cds) "T+"
                          (reverse-compliment prime5-const))
        parasite? (re-pattern (str parasite "|" rev-parasite))
        parasites (->> (r/filter #(re-find parasite? %) inseqs)
                       (into []))]
    (count parasites)))

(defn- length-distribution [inseqs]
  (rfrequencies (r/map count inseqs)))

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
