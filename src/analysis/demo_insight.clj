(ns analysis.demo-insight
  (:require [clojure.java.jdbc :as jdbc]
            [clojure.core.reducers :as r]
            [clojure.java.io :as io2]
            [clojure.contrib.string :as str]
            [edu.bc.fs :as fs]
            [clojure.pprint :as pp])
  (:use edu.bc.utils.fold-ops
        edu.bc.utils.probs-stats
        hts-exploration.db-queries
        hts-exploration.globals
        hts-exploration.utils))

;;;---------preface-------------------------------------------------------------------

(declare filter-clusters)
(defn cluster-vs-ident
  "Sample code to get cluster sizes and mean cluster identity from
  MySql database. Returns a vector of [cluster-id cluster-size
  cluster-ident]"

  []
  (->> (filter-clusters 100000 84 0);large clusters
       (map (juxt :cid :cluster_size :mean_ident))))

;;;----------------sql queries---------------------------------------------------

(defn- generic-get-seqs
  "Generic way to get sequences out of the database. It joins the
  sequences to cluster_ids on selex_ids so we are getting only usable
  sequences (filtered on quality, base usage, strand) which have been
  clustered."

  []
  (let [q [(str "select sr.selex_id as id, cluster_id as cid,
                       substring(sr.sequence, sr.usable_start+1, sr.length) as hitseq
                  FROM selex_reads as sr, cdhit
                 WHERE sr.selex_id=cdhit.selex_id
                   AND sr.usable=1
                   AND sr.strand=?") 1]];example
    (sql-query q)))

(defn- get-seq-by-id
  "Get a particular sequence using the selex_id."
  
  [sid]
  (let [q [(str "select substring(sr.sequence, sr.usable_start+1, sr.length) as hitseq
                   FROM selex_reads as sr
                  WHERE selex_id=?") sid]]
    (->> (sql-query q) (map :hitseq) first)))

(defn cluster-id->seqs
  "Gets distinct sequences with cluster-id cs ordered by frequency"

  [cs]
  (letfn [(q2 [cid]
            [(str "select sr.selex_id as sid,
                          substring(sr.sequence, sr.usable_start+1, sr.length) as hitseq
                     FROM cdhit, selex_reads as sr
                    WHERE sr.selex_id=cdhit.selex_id
                      AND cdhit.cluster_id=?") cid])]
    (->> (q2 cs) sql-query
         (map (juxt :sid :hitseq))
         freqs-id-seq
         (sort-by #(-> % second count) >);order by freq
         (map (juxt #(->> % second first) first)) ;[[id s] ...]
         )))

(defn filter-clusters
  "Gets all clusters (ordered by the number of distinct sequences)
  which have cluster_size > csize with the number of distinct seqs >
  dsize and mean_ident > ident"
  
  [csize ident dsize]
  (let [q [(str "select cdhit.cluster_id as cid,
                        cdhit.cluster_size,
                        count(distinct substring(sr.sequence, sr.usable_start+1, sr.length)) as hitseq,
                        avg (cdhit.clstr_iden) as mean_ident
                  FROM cdhit, selex_reads as sr
                 WHERE sr.selex_id=cdhit.selex_id
                   AND cluster_size > ?
              GROUP BY cdhit.cluster_id
                HAVING mean_ident > ?
                   AND hitseq > ?
              ORDER BY hitseq DESC") csize ident dsize]]
    (sql-query q)))

;;;-----------------------------------------------------------------------------

;;;----------------constants ---------------------------------------------------------
(def most-freq-seqs
  "The most frequent sequences in the database"
  
  (->> (generic-get-seqs)
       (mapv (juxt (juxt :id :cid) :hitseq))
       freqs-id-seq ;seq freq map
       (remove #(< (-> % second count) 100))
       (sort-by #(-> % second count) >);sort by frequencies
       (map (fn [[s ids]]
              (let [id (ffirst ids);selex_id
                    cnt (count ids);seq count
                    cids (distinct (map second ids))];unique clusters seq is found in
                [s cnt id cids])))
       (take 10)))

(def tested-seqs
  "Sequences which have been experimentally tested"

  {:5motif {:id 3196036 :cluster 62979 :csize 1
            :ident 100 :seq "ATCGAAAGGAGAATGGAATCGAGCAATCGA"} 
   :4motif {:id 455019 :cluster 51432 :csize 1
            :ident 100 :seq "GTCAAACAAAACAAAAGACGAAGACGCACT"}
   :2motif {:id 4077286 :cluster 1503 :csize 72
            :ident 85.19 :seq "AACGATTCGAAAGTGAAAGAAAGAGAAA"}
   :ci7-255 {:id 5192670 :cluster 2346 :csize 367
             :ident 94.83 :seq "CGATCACACGAGAACATCGGTGATTTGGTGTCAAT"}
   :ci7-73 {:id 5229242 :cluster 1887 :csize 562
            :ident 92.86 :seq "AGGCAAACCGATCCTAACGAATGCTTGGTG"}
   :ci17-156 {:id 5205749 :cluster 3580 :csize 728
              :ident 94.83 :seq "ACCCAAGACGGCTCTACAGTAAGATAGCCTA"}
   :bob {:id 3597107 :cluster 9131 :csize 61
         :ident 98.25 :seq "TCCTTCGCTTATTCGGAGTAGATCACGTGA"}
   :wka {:id 5267819 :cluster 9391 :csize 8825
         :ident 99.29 :seq "GCGGACAGCGAGACAGATCGAAGGTTTTGA"}})

;;;-------------------------------------------------------------------------------

;;;---------counting kmers--------------------------------------------------------

(defn sqs->kmer 
  "Turns a set of sequences into a kmer probs according to the
  probability that it will appear in the cluster"
  
  [k sqs]
  (->> sqs
       (vfold #(freqn k %));vector of freq maps
       (apply merge-with +)
       probs));normalize to num. kmers

(defn cluster->kmer
  "Takes a cluster number and turns the cluster into a set of
  sequences and then calculates the kmer probs"
  
  [k clstr]
  (->> (cluster-id->seqs clstr)
       (map second);keep seq, ignore id
       (sqs->kmer k)))

;;;--------------------------------------------------------------------------------

;;;----------examples--------------------------------------------------------------

(defn kmer-count-cluster [] (cluster->kmer 2 9131)) ;particular cluster example

(defn kmer-count-clusters
  "Calculating kmers for sequences within a set of clusters."

  []
  (let [ks (vec (range 2 3));kmer length=2
        cs [:tested (map #(-> % second :cluster) tested-seqs)]
        a (atom "")]
    (reduce (fn [V c]
              (if (keyword? c)
                (do (swap! a #(identity %2) c) V)
                (let [prob-map (->> (map #(cluster->kmer % c) ks) ;calc kmers
                             (apply merge))]
                  (conj V [@a c prob-map]))))
            [] (->> cs flatten distinct))))

;;;------code to do work----------------------------------------------------------
(let [kmers (vec (range 2 5));interesting kmer lengths
      interesting-clusters (concat [:tested (map #(-> % second :cluster) tested-seqs)]
                                   [:doubleton (map #(-> % last first) most-freq-seqs)]
                                   [:large (->> (filter-clusters 1000 84 1) (map :cid) (take 5))]
                                   [:ident (->> (filter-clusters 100 90 0) (map :cid))])
      ;;all possible kmer combinations
      combos (sort (mapcat #(add-letters % [\A \C \T \G]) kmers))]
  (with-open [wrtr (io2/writer "data/150211/large-cluster-kmers-demo.csv")]
    (.write wrtr (str "cluster,type," (str/join "," combos) "\n"))
    (reduce (fn [type c]
              (if (keyword? c)
                (str/as-str c)
                (let [prob-map (->> (map #(cluster->kmer % c) kmers)
                             (apply merge))]
                  (.write wrtr (str c "," type))
                  (doseq [i combos]
                    ;;make sure all kmer vectors print out in same order
                    (.write wrtr (str "," (get prob-map i 0))))
                  (.write wrtr "\n")
                  type)));keep type
            "" (->> interesting-clusters flatten distinct))))

