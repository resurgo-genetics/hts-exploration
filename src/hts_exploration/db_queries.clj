(ns hts-exploration.db-queries
  (:require [clojure.java.jdbc :as jdbc]
            [clojure.string :as str]))

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

(defn get-seqs-sid [selex-id-vec]
  (->> ["select sr.selex_id as sid, ss.seq_id as ssid, sk.cluster_id2 as cid, cluster_size as csize, round_number as r, ss.hitseq from selex_reads as sr JOIN selex_keys as sk ON sr.selex_id=sk.selex_id JOIN selex_seqs as ss ON sk.seq_id=ss.seq_id JOIN cdhit2 ON sk.selex_id=cdhit2.selex_id where sr.selex_id IN (" (str/join "," selex-id-vec) ");"]
       (apply str)
       sql-query
       (map (juxt :sid :ssid :cid :csize :r :hitseq))))

(defn get-seqs-ssid [seq-id-vec]
  (->> ["select sr.selex_id as sid, ss.seq_id as ssid, sk.cluster_id2 as cid, cluster_size as csize, round_number as r, ss.hitseq from selex_reads as sr JOIN selex_keys as sk ON sr.selex_id=sk.selex_id JOIN selex_seqs as ss ON sk.seq_id=ss.seq_id JOIN cdhit2 ON sk.selex_id=cdhit2.selex_id where sk.seq_id IN (" (str/join "," seq-id-vec) ");"]
       (apply str)
       sql-query
       (map (juxt :sid :ssid :cid :csize :r :hitseq))))

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
  (letfn [(q2 [cid]
            [(str "SELECT hitseq,
                          count(sr.selex_id) as cnt , ss.seq_id as ssid, sk.cluster_id2 as cid,
                          sr.round_number as rnd, sr.const_start as cs
                     FROM selex_reads as sr, selex_keys as sk, selex_seqs as ss, cdhit2
                    WHERE sr.selex_id=cdhit2.selex_id and sr.selex_id=sk.selex_id
                      and sk.seq_id=ss.seq_id and sk.cluster_id2=? and sr.usable=1 AND sr.strand=1
                 group by ssid
                 order by cnt DESC") cid])]
    (->> (q2 cid) sql-query (map (juxt :hitseq :ssid :cnt :rnd :cs)))))

(defn most-freqn-in-cluster [n cid]
  ((comp (partial take n) get-seqs-in-cluster) cid))

(defn filter-clusters
  "Gets all clusters (ordered by the number of distinct sequences)
  which have cluster_size > csize with the number of distinct seqs >
  dsize and mean_ident > ident"
  
  [csize ident]
  (letfn [(q2 [csize ident]
            [(str "select cdhit2.cluster_id as cid,
                          cdhit2.cluster_size as csize,
                          avg(cdhit2.clstr_iden) as mean_ident
                     FROM cdhit2
                    WHERE cluster_size > ?
                 GROUP BY cdhit2.cluster_id
                   HAVING mean_ident > ?
                 ORDER BY csize DESC") csize ident])]
    (->> (q2 csize ident) sql-query (map (juxt :cid :csize :mean_ident)))))

;;;create selex seqs table
"create table selex_seqs (seq_id int(11) NOT NULL AUTO_INCREMENT, hitseq VARCHAR(100) NOT NULL, primary key (seq_id));"

"INSERT INTO selex_seqs (hitseq) select distinct substr(sr.sequence, usable_start+1, sr.length) from selex_reads as sr,cdhit where sr.selex_id=cdhit.selex_id and usable=1 and strand=1"

"CREATE UNIQUE INDEX seqIndex ON selex_seqs (hitseq);"

;;;test query for joining tables based on sequence
"select sr.selex_id, ss.seq_id, ss.hitseq, substr(sr.sequence, usable_start+1, sr.length) as s from selex_reads as sr INNER JOIN selex_seqs as ss ON substr(sr.sequence, usable_start+1, sr.length)=ss.hitseq limit 10;"

;;;create selex keys table
"CREATE TABLE selex_keys select sr.selex_id, ss.seq_id, cluster_id from selex_reads as sr INNER JOIN selex_seqs as ss ON substr(sr.sequence, usable_start+1, sr.length)=ss.hitseq INNER JOIN cdhit on sr.selex_id=cdhit.selex_id;"

"CREATE INDEX selexIndex on selex_keys (selex_id);
"
"CREATE INDEX clusterIndex on selex_keys (cluster_id);"

;;;get distinct seqs which appear in multiple rounds. Count the seqs 
"select ss.seq_id, group_concat(distinct round_number), count(distinct round_number) as cnt, count(sk.selex_id) from selex_keys as sk, selex_seqs as ss, selex_reads as sr where sr.selex_id=sk.selex_id and sk.seq_id=ss.seq_id and round_number in (4, 11) group by ss.seq_id having cnt > 1 order by cnt DESC limit 10;"
