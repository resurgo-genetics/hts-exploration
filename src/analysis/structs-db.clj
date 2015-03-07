(ns analysis.structs-db
  (:require [clojure.java.jdbc :as jdbc]
            [clojure.contrib.io :as io]
            [clojure.java.shell :as shell])
  (:use edu.bc.utils.fold-ops
        hts-exploration.db-queries
        hts-exploration.globals
        hts-exploration.utils))

(defn- get-seqs-specific
  "Specific to getting full length sequences in a round which have a
  variable region between 28 and 32 nt in length. This is specific to
  this ns unlike get-seqs which is used to simply store variable
  regions with ids."
  
  [round]
  (->> (round-all-usable-seqs round) sql-query
       (reduce (fn [M {:keys [id hit_seq cs]}]
                 (let [sq (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" hit_seq)
                       cur-val (get M (second sq) [])]
                   (if (and sq (<= 28 cs 32)) ;standard sequence
                     (assoc M (second sq) (conj cur-val [id (first sq)]))
                     M)))
               {})))

(defn get-doubletons
  "Specific to getting full length doubletons which have a variable
  region between 28 and 32 nt in length."
  
  [round]
  (filter #(> (->> % second count) 1) (get-seqs-specific round))) ;doubleton

(defn- rnafold-shell [fasta]
  ((shell/sh "RNAfold"
             "-P" "/usr/local/ViennaRNA/misc/rna_andronescu2007.par"
             "--noPS"
             :in ((shell/sh "cat" fasta) :out))
   :out))

(defn- centroid-shell [fasta]
  ((shell/sh "/home/peis/bin/centroid_fold-0.0.9-x86_64-linux/centroid_fold" fasta) :out))

(defn partition-seqs-to-fold [n sql-sqs]
  (->> sql-sqs vals (apply concat)
       (partition-into n)))

(defn distribute-folding
  "Distribute RNAfold to multiple cpus by breaking into multiple fasta
  files and then using RNAfold on each file. The results are combined
  into the same file again for parsing."

  [tmp-prefix sql-sqs]
  (doall
   (map-indexed (fn [c p]
                  (let [outfile (str tmp-prefix c ".fna")]
                    (io/with-out-writer outfile
                      (doseq [[id sq] p]
                        (println (str ">" id))
                        (println sq)))
                    outfile))
                sql-sqs)))

(defn execute-folding [shell-fold outfile & files-to-fold]
  (io/with-out-writer outfile
    (doseq [out (doall (pmap shell-fold files-to-fold))];parallel folding
      (println out))))

(defn read-structs 
  "Reads the file and attempts to parse it into a
  datastructure according to the lines: name, sequence,
  structure+energy"
  [f]
  (->> (io/read-lines f)
       (partition-by #(re-find #">" %) )
       (partition-all 2)
       (map #(apply  concat %))))

(defn get-struct-deltag [v] (vec (take-last 2 v)))

(defn split-rnafold-struct-line
  "After reading the file, this function splits the structure and
  deltaG from the last line. The input should have been read and put
  into the 3 element vector corresponding to the 3 lines. Returns a 4
  element vector."
  
  [[id s st]]
  (let [[st g] (str/split #" " 2 st)
        dg (Double/parseDouble (str/replace-re #"\(|\)" "" g))]
    [(str/drop 1 id) s st dg]))

(defn split-centroid-struct-line [[id s st]]
     (let [[st g] (str/split #" " 2 st)
           dg (Double/parseDouble (second (re-find #"e=(.*)\)" g)))]
       [(str/drop 1 id) s st dg]))

(defn check-db 
  "Checks the structues to be added to see if they exist in the
  database. Catagorizes them into ignore, update, and insert according
  to if the sequence exists and the MFE/centroid flags"

  [type structs-to-add]
  {:pre [(or (= type :centroid)
             (= type :MFE))]}
  (let [S (->> (jdbc/query mysql-ds "SELECT * FROM selex_structs;")
               (map (fn [x] [(x :structure) {:MFE (when (x :MFE) :MFE)
                                           :centroid (when (x :centroid) :centroid)}]))
               (into {}))]
    (group-by (fn [st]
                (if (contains? S st)
                  (if (= type (get-in S [st type]))
                    :ignore
                    :update)
                  :insert))
              structs-to-add)))

(defn insert-struct [data & {:keys [MFE centroid] :as x}]
  (jdbc/with-db-transaction [c mysql-ds]
    (doseq [i data
            :let [d (->> (interleave [:structure :deltaG] i)
                         (apply hash-map))
                  d (if centroid (assoc d :centroid 1) d)
                  d (if MFE (assoc d :MFE 1) d)]]
      (jdbc/insert! c :selex_structs d))))

(defn update-struct [data & {:keys [MFE centroid] :as x}]
  (jdbc/with-db-transaction [c mysql-ds]
    (doseq [i data
            :let [d (->> (interleave [:structure :deltaG] i)
                         (concat (filterv (complement (comp nil? second)) x));mfe/centroid flags
                         (apply hash-map ))]]
      )))

(defn selex-id->struct-id [struct-parser f]
      (let [data (->> f read-structs
                      (map struct-parser)
                      (map #(vector (% 2) (% 0)))
                      (reduce (fn [M [st id]]
                                (assoc M st (conj (get M st []) id)))
                              {}))
            ;;hashmap of structs to their struct_ids
            sql-structs (reduce (fn [M x]
                                  (assoc M (x :structure) (x :struct_id)))
                                {} (jdbc/query mysql-ds "SELECT * FROM selex_structs"))]
        (reduce (fn [M [st ids]]
                  ;;id maps to struct->struct_id
                  (assoc M ids (sql-structs st)))
                {} data)))

(defn- update-struct-id [mysql-column selex-ids struct-id]
  (str "UPDATE selex_reads SET " mysql-column "=" struct-id
       " WHERE selex_id IN (" (str/join "," selex-ids) ");"))
