(ns hts-exploration.hts-utils.file
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [clojure.java.io :as io2]
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
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.sequtils.files
        [edu.bc.bio.seq-utils :only (markov-step)]
        edu.bc.utils.fold-ops))

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

(defn write-fasta
  "Write out a fasta given sequences. If inseq is not a pair then
  numerical names are given for each sequence."

  [outfile inseqs & {:keys [info append]
                     :or {info :both append false}}]
  (cond
    (true? append)
    (with-open [wrtr (io2/writer outfile :append true)]
      (doseq [[nm s] inseqs]
        (.write wrtr (str ">" nm "\n"))
        (.write wrtr s)
        (.write wrtr "\n")))
    (= info :both)
    (io/with-out-writer outfile
      (doseq [[nm s] inseqs]
        (println (str ">" nm))
        (println s)))
    (= info :data)
    (write-fasta outfile (map-indexed vector inseqs))
    )
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
