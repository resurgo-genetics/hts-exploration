(ns hts-exploration.preprocess
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.io :as io]
            [clojure.contrib.string :as str]
            [clojure.core.reducers :as r]
            [iota])
  (:use edu.bc.utils.fold-ops
        hts-exploration.globals
        hts-exploration.utils))


(defn fastq->csv [infile]
  (let [outfile (fs/replace-type infile ".csv")]
   (io/with-out-writer outfile
     (println "name,sequence,quality")
     (doseq [[nm inseq _ qual] (->> infile
                                    io/read-lines
                                    (partition-all 4))]
       (println (str/join "," [nm inseq qual]))))
   outfile))

(defn fastq->table
  ([round-number infile outfile] (fastq->table 0 round-number infile outfile))
  ([id-start round-number infile outfile]
     (let [r round-number
           c (atom id-start)]
       (io/with-out-writer outfile
         (println (str "selex_id name sequence quality quality_filter standard_bases"
                       " not_parasite parasite_type round_number usable strand"
                       " from to length"))
         (doseq [[nm inseq _ qual] (->> infile
                                        io/read-lines
                                        (partition-all 4)
                                        ;(take 100)
                                        )
                 :let [qfilter (if (qual-good? 20 qual) 1 0)
                       std-bases (if (standard-chars? inseq) 1 0)
                       strand (cond (re-find forward-re inseq) 1
                                    (re-find reverse-re inseq) -1
                                    :else 0)]]
           (let [[from to] (cond (= strand 1) (str-regex-index forward-re inseq)
                                 (= strand -1) (str-regex-index reverse-re inseq)
                                 (zero? strand) [-1 -1])
                 s (when-not (zero? strand)
                     (subs inseq from to))
                 len (- to from)
                 parasite-type (cond (nil? s) ""
                                     (= strand 1) (first (parasite? s))
                                     (= strand -1) (first (parasite-rev? s)))
                 not-parasite (if-not parasite-type 1 0)
                 usable? (if (= qfilter std-bases not-parasite 1) 1 0)]
             (swap! c inc)
             (->> [@c nm inseq qual qfilter std-bases
                   not-parasite (or parasite-type "") r usable? strand
                   from to len]
                  (apply prn-str ) print;used in s15-round11-mysql.table
                  ;;(str/join ",") prn;error checking
                  ))))
       outfile)))

(defn quality-filter [infile]
  (let [outfile (fs/replace-type infile ".qual-filter.csv")]
   (io/with-out-writer outfile
     (println "name,sequence,quality")
     (doseq [x (->> infile
                    iota/seq
                    (r/drop 1) ;drop header
                    (r/map #(vec (str/split #"," 3 %)))
                    qual-remove
                    (r/fold (fn ([] [])
                              ([l r] (conj l r)))
                            (fn ([] [])
                              ([V x] 
                                 (conj V (str/join "," x))))))];csv format
       (println x)))))

(comment
  ;;process raw SELEX data
  (->> (fastq->csv S15-round11-fastq-1)
       quality-filter)

  (->> (fastq->csv S15-round11-fastq-2)
       quality-filter)

  ;;number of parasites in each round
  (let [infiles [S15-round10-qual-filter-csv S15-round11-qual-filter-csv S15-r1-qual-filter-csv]
        get-usable (fn [f]
                     (->> (iota/vec f)
                          read-fasta-csv
                          (r/map second)
                          (r/map #(re-find forward-re %))
                          (r/remove nil?)
                          (into [])))
        parasite (str prime5-const "AA" cds)
        contain-parasite? (fn [s] (re-find (re-pattern parasite) s))]
    (doseq [f infiles
            :let [usable (get-usable f)
                  cnt (->> usable
                           (r/filter contain-parasite?)
                           (into []) 
                           count)]]
      (prn :file (fs/basename f))
      (prn :parasites cnt
           :total (count usable)
           :fract (double (/ cnt (count usable))))))

  "create table selex_reads (selex_id INT not NULL auto_increment, name varchar(70) not NULL, sequence VARCHAR(100) not NULL, base_quality varchar(100) not null, quality_filter INT(1), standard_bases INT(1), not_parasite INT(1), parasite_type varchar(100), round_number INT, usable INT(1), strand INT, primary key (selex_id, name));"
  
  "load data local INFILE '/home/peis/S15SELEXHTSdata/s15-round11-mysql.table' into table selex_test FIELDS TERMINATED by ' ' ENCLOSED BY '\"' LINES TERMINATED BY '\n' IGNORE 1 ROWS;"

  "alter table selex_test change id selex_id int;"

  ;;update some mistakes made in making the db
  (let [stmt (str "update selex_reads as sr SET not_parasite=1, parasite_type='', usable=1 where sr.name IN "
                  (->> (pr-str (map first foo))
                       (str/replace-re #"\"" "'")
                       (str/replace-re #"'\s'" "','" ))
                  " ;")]
    stmt) ;generate command

  ;;update mistakes made in reverse strand seqs in db
  (let [stmts (map (fn [part]
                     (str "UPDATE selex_reads as sr SET not_parasite=1, parasite_type='', usable=1 where sr.name IN "
                          (->> part
                               (map first)
                               pr-str
                               (str/replace-re #"\"\s\"" "\",\""))
                          " ;"))
                   (partition-into 100 foo2))]
    (doseq [stmt stmts]
      (clojure.java.jdbc/execute! mysql-test/mysql-ds [stmt])));update
                                   ;table
  ;;above code doens't work cuz it makes the same update to each
  ;;row. this corrects the discrepancies. 
  (let [mysql (clojure-csv.core/parse-csv (clojure.java.io/reader "/tmp/s15-round11-temp.csv"))
        tablefn-generator (map first (rest (clojure-csv.core/parse-csv (clojure.java.io/reader "/tmp/s15-round11-table.csv"))))
        to-change (->> (map #(if (= %1 %2) true [%1 %2]) 
                            tablefn-generator
                            (map #(str/join "," %)  mysql))
                       (remove true?)
                       (map first)
                       (map (fn [s] ((juxt #(take 2 %) #(take 2 (take-last 8 %))) (str/split #"," s))))
                       (map #(apply concat %)) 
                       (map (fn [[selex-id name not-parasite? parasite-type]]
                              [(Integer/parseInt selex-id) name 
                               (if (empty? parasite-type) 1 0)
                               (if (empty? parasite-type) "" parasite-type)]))
                                        ;(take 100)
                       (group-by (fn [[_ _ not-parasite? parasite-type]]
                                   [not-parasite? parasite-type])))]
    (clojure.java.jdbc/with-db-connection [conn mysql-ds]
      (doall
       (doseq [[[not-parasite? parasite-type] x] to-change
               y (partition-all 10000 x) 
               :let [in (str/join "," (map first y))
                     stmt (str "UPDATE selex_reads"
                               " SET not_parasite=" not-parasite?
                               ", parasite_type='" parasite-type "'"
                               " where selex_id IN (" in ");")]]
         (clojure.java.jdbc/execute! conn [stmt])))))
  
  ;;test integrity of data and database
  ;;forward seqs
  (let [M1 (->> (mysql-test/sql-query "select sr.name, sr.sequence, sr.usable_start, sr.usable_stop from selex_reads as sr where sr.usable=1 and sr.not_parasite=1 and sr.strand=1")
                (map (fn [m]
                       [(m :name) 
                        (subs (m :sequence) (m :usable_start) (m :usable_stop))]))
                (into #{}))
        M2 (->> S15-round11-qual-filter-csv
                get-usable 
                parasite-remove
                r/foldcat set)]
    [(count (clojure.set/difference M1 M2)) (count (clojure.set/difference M2 M1))])

  ;;reverse seqs
  (let [M1 (->> (mysql-test/sql-query "select sr.name, sr.sequence, sr.usable_start, sr.usable_stop from selex_reads as sr where sr.usable=1 and sr.not_parasite=1 and sr.strand=-1")
                (map (fn [m]
                       [(m :name) 
                        (subs (m :sequence) (m :usable_start) (m :usable_stop))]))
                (into #{}))
        M2 (->> (iota/vec S15-round11-qual-filter-csv)
                read-fasta-csv
                (r/reduce (fn [V [nm inseq]]
                            (let [sq (re-find reverse-re inseq)]
                              (if sq
                                (conj V [nm sq])
                                V)))
                          [])
                (remove #(parasite-rev? (second %)))
                set)]
    (clojure.set/difference M2 M1)))

(let [stmts (map (fn [part]
                                         (str "select * from selex_reads as sr where sr.name IN "
                                              (->> part
                                                   (map first)
                                                   pr-str
                                                   (str/replace-re #"\"\s\"" "\",\""))
                                              " ;"))
                                       (partition-into 100 foo2))]
                        (prn :count (count foo2)) (count
                         (mapcat (fn [stmt]
                                   (->> (mysql-test/sql-query stmt)
                                        (map (fn [m]
                                               [(m :name) 
                                                (subs (m :sequence) (m :usable_start) (m :usable_stop))]))
                                        #_(filter #(parasite-rev? (second %)))))
                                 (take 2 stmts))))
