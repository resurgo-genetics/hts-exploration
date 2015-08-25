(ns hts-exploration.db-queries
  (:require [clojure.java.jdbc :as jdbc]))

(def mysql-ds
  {:classname "com.mysql.jdbc.Driver"
   :subprotocol "mysql"
   :subname "//127.0.0.1:3306/s15selexdata"
   :user "root"
   :password "rna314rulz"})

(defn sql-query [stmt & {:keys [f p] :or {f identity p false}}]
  (let [q (partial  jdbc/query mysql-ds)]
    (when p (println stmt))
    (cond (string? stmt) (q [stmt])
          (vector? stmt) (q stmt)
          :else "Invalid query")))

(defn round-all-usable-seqs [round]
  [(str "SELECT sr.selex_id as id, SUBSTRING(sr.sequence, sr.usable_start+1, sr.length) as hit_seq,
                sr.const_start as cs
          FROM selex_reads as sr
         WHERE sr.usable=1
           AND sr.strand=1
           AND sr.round_number=?")
   round])

(defn round-all-usable-seqs-structs [round]
  [(str "SELECT sr.selex_id as id, SUBSTRING(sr.sequence, sr.usable_start+1, sr.length) as hit_seq,
                sr.const_start as cs, sr.struct_id as sid, ss.structure as hit_struct
          FROM selex_reads as sr, selex_structs as ss
         WHERE sr.usable=1
           AND sr.strand=1
           AND sr.struct_id=ss.struct_id
           AND sr.round_number=?")
   round])

(defn round-all-usable-seqs-centroid [round]
  [(str "SELECT sr.selex_id as id, SUBSTRING(sr.sequence, sr.usable_start+1, sr.length) as hit_seq,
                sr.const_start as cs, sr.centroid_id as sid, ss.structure as hit_struct
          FROM selex_reads as sr, selex_structs as ss
         WHERE sr.usable=1
           AND sr.strand=1
           AND sr.centroid_id=ss.struct_id
           AND sr.round_number=?")
   round])

"SELECT SQL_CALC_FOUND_ROWS foo.hit_seq FROM (SELECT SUBSTRING(sr.sequence, sr.usable_start+1, sr.length) as hit_seq FROM selex_reads as sr WHERE sr.usable=1 AND sr.strand=1) as foo WHERE foo.hit_seq REGEXP 'A[ACT][GC][ACG]A.{3,}T[CGT][GC][TGA]T' LIMIT 10; select found_rows();"


(defn add-const-start
  "Adds the start of the constant region to the mysql data. Requires
  the start position of the start and the id"
  
  [[cstart id] & more]
  (apply vector "UPDATE selex_reads SET const_start=? WHERE selex_id=?" [cstart id] more))

(defn update-const-start
  "Updates the selex database to add a constant start. Uses the data
  produced by find-const-start. Update operations are broken up into
  smaller so that there are no overflow errors."
  
  [data]
  (doseq [i (->> (map (fn [[id const-start]] [(dec const-start) id]) data)
                 (partition-all 5000))]
    (jdbc/execute! mysql-ds (apply add-const-start i) :multi? true) ))

(defn get-seqs-in-cluster
  "Takes a cluster ID and gets the distinct seqs in the cluster"

  [cid]
  [(str "select sr.selex_id as id, cdhit.cluster_size as csize,
            substring(sequence, usable_start+1,sr.length) as hitseq
            FROM selex_reads as sr, cdhit where
            sr.selex_id=cdhit.selex_id AND sr.usable=1 AND sr.strand=1
            AND cluster_id=?
           GROUP by hitseq;") cid])

(defn most-freqn-in-cluster [n cid]
  (letfn [(q2 [cid]
            [(str "SELECT substring(sr.sequence, usable_start+1, sr.length) as hitseq,
                          count(sr.selex_id) as cnt , sr.selex_id as id, cluster_id as cid
                     FROM selex_reads as sr, cdhit
                    WHERE sr.selex_id=cdhit.selex_id and cluster_id=?
                 group by hitseq
                 order by cnt DESC") cid])]
    (->> (q2 cid) sql-query (map (juxt :hitseq :id :cnt))
         (take n))))
