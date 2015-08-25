(ns analysis.cluster-cdhit-sql
 (:require [clojure.java.jdbc :as jdbc]
            [clojure.core.reducers :as r]
            [clojure.java.io :as io2]
            [clojure.contrib.string :as str]
            [edu.bc.fs :as fs]
            [clojure.java.shell :as shell]
            [clojure.pprint :as pp]
            [edu.bc.bio.sequtils.files :refer (read-seqs)]
            [hts-exploration.hts-utils.file :refer (write-fasta)])
  (:use edu.bc.utils.fold-ops
        edu.bc.utils.probs-stats
        hts-exploration.db-queries
        hts-exploration.globals
        hts-exploration.utils))



(defn cluster-ids-to-ensemble
  "Gets distinct sequences with cluster-id cs organized by frequency"

  [cs]
  (letfn [(q2 [cid]
            [(str "select sr.selex_id as sid,
                          substring(sr.sequence, sr.usable_start+1, sr.length) as hitseq
                     FROM cdhit, selex_reads as sr
                    WHERE sr.selex_id=cdhit.selex_id
                      AND cdhit.cluster_id=?") cid])]
    (->> (q2 cs) sql-query
         (map (juxt :hitseq :sid))
         (reduce (fn [M [s id]] (assoc M s (conj (get M s []) id))) {})
         (sort-by #(-> % second count) >);order by freq
         (map (juxt #(->> % second first) first)) ;[id s]
         )))

(defn q1
  "Gets all clusters (ordered by the number of distinct sequences)
  which have cluster_size > csize with the number of distinct seqs >
  dsize and mean_ident > ident"
  
  [csize ident dsize]
  [(str "select cdhit.cluster_id as cid,
                cdhit.cluster_size,
                count(distinct substring(sr.sequence, sr.usable_start+1, sr.length)) as hitseq,
                avg (cdhit.clstr_iden) as mean_ident
           FROM cdhit, selex_reads as sr
          WHERE sr.selex_id=cdhit.selex_id
            AND cluster_size > ?
       GROUP BY cdhit.cluster_id having mean_ident > ?
            AND hitseq > ?
       ORDER BY hitseq DESC") csize ident dsize])

(defn take-frac [p coll] (take (* p (count coll)) coll))

(def foo
  (future
    (let [wdir "data/150309/cluster"
          x (->> (sql-query (q1 100 90 100))
                 (map :cid)
                 (pmap (fn [clstr]
                         (let [dir (fs/join wdir clstr)
                               _ (when-not (fs/directory? dir)
                                   (fs/mkdir (fs/join wdir clstr)))
                               infasta (fs/join dir "cluster-seqs.fna")
                               seqs-to-ensemble (take 1000 (cluster-ids-to-ensemble clstr))
                               ps-files (->> (map first seqs-to-ensemble)
                                             (map #(fs/join dir (str % "_dp.ps"))))]
                           ;;make sure ps files exist for all seqs
                           (when-not (every? fs/exists? ps-files)
                             (write-fasta infasta seqs-to-ensemble :info :both) 
                             (shell/sh "RNAfold" "-p"
                                       "-P" "/usr/local/ViennaRNA/misc/rna_andronescu2007.par"
                                       "--noPS"
                                       :in ((shell/sh "cat" infasta) :out)
                                       :dir dir))
                           ;;calc pairwise ensemble distances
                           (let [probs (for [f ps-files] (last (parse-ps f)))]
                             (vector clstr (mean (pairwise ensemble-dist probs))))))))]
      (with-open [wrtr (io2/writer "data/150309/intracluster-dist.csv")]
        (.write wrtr "cluster,ensemble distance\n")
        (doseq [i x] (.write wrtr (str (str/join "," i) "\n")))))))

(let [wdir "data/150309/cluster"
      x (->> (sql-query (q1 100 90 100))
             (map :cid)
             (map (fn [clstr]
                    (let [dir (fs/join wdir clstr)
                          seqs-to-ensemble (take 5 (cluster-ids-to-ensemble clstr))
                          ps-files (->> (map first seqs-to-ensemble)
                                        (map #(fs/join dir (str % "_dp.ps"))))]         
                      (interleave (repeat clstr) ps-files))))
             flatten
             (partition-all 2))
      x (pmap (fn [[x1 x2]]
                (let [[id1 ps1] x1
                      [id2 ps2] x2
                      [s1 f1 bprob1] (parse-ps ps1)
                      [s2 f2 bprob2] (parse-ps ps2)]
                  [id1 id2 f1 f2 s1 s2 (levenshtein s1 s2) (ensemble-dist bprob1 bprob2)]))
              (pairwise vector x))]
  (with-open [wrtr (io2/writer "data/150309/intercluster-dist.csv")]
    (.write wrtr "cluster1,cluster2,ps file1,ps file2,seq1,seq2,levenshtein dist,ensemble dist\n")
    (doseq [i x] (.write wrtr (str (str/join "," i) "\n")))))
