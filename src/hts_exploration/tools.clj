(ns hts-exploration.tools
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [clojure.pprint :as pp])
  (:use edu.bc.utils
        hts-exploration.globals
        hts-exploration.utils))

(def cdhit-dir "/usr/local/cd-hit/")

(assert-tools-exist [(str cdhit-dir "cd-hit-est")
                     (str cdhit-dir "clstr2txt.pl")])

(defn run-cdhit [infile outfile & args]
  (runx "/usr/local/cd-hit/cd-hit-est" args))

(defn run-clstr2txt [infile outfile]
  (runx "/usr/local/cd-hit/clstr2txt.pl" infile :> outfile))

(comment
  (runx "/usr/local/cd-hit/cd-hit-est" "-i" "/home/peis/data/s15-round11-cdhit-cluster-7.fna" "-o" "/home/peis/data/s15-round11-cdhit-cluster-7.cdhit" "-mismatch" "-1" "-gap" "-1" "-gap-ext" "0" "-c" "0.95" "-n" "8")
  (runx "/usr/local/cd-hit/clstr2txt.pl" "/home/peis/data/s15-round11-cdhit-cluster-7.cdhit.clstr" :> "/home/peis/data/s15-round11-cdhit-cluster-7.cdhit.clusters")
  (hts-exploration.tools/make-cluster-fna 50 "/home/peis/data/s15-round11-cdhit-cluster-7-" "/home/peis/data/s15-round11-cdhit-cluster-7.cdhit.clusters"))

(defn- add-seq-to-cluster
  "Takes the cluster data from cdhit converted using clstr2txt and
  adds the sequence from a hashmap containing [name sequence]"

  [data clusters]
  (let [[header clusters] (split-at 1 clusters)
        header (map keyword (first header))]
    (prn :header header)
    (second 
     (reduce (fn [[clstr-rep M] clstr]
               (let [c (zipmap header clstr)
                     k (if (= "1" (c :clstr_rep)) (c :id) clstr-rep)
                     new-entry {:seq (data (c :id)) 
                                :len (Integer/parseInt (c :length))
                                :id (c :id) 
                                :clstr-rep (Integer/parseInt (c :clstr_rep))}
                     old-entry (get M [(c :clstr) k] [])]
                 [k (assoc M [(c :clstr) k] (conj old-entry new-entry))]))
             ["" {}] clusters))))

(defn- create-map
  "reads in a fasta-csv file and creates a hashmap of [name seq]
  pairs. must be used because cdhit only uses first 19chars from the
  name. To make seq names unique, numbers have been added on hence
  this function"
  
  [infile]
  (->> (get-usable infile) ;S15-round11-qual-filter-csv
       (map-indexed (fn [i [nm s]] 
                      [(str/take 19 (str i nm)) s]))
       (into {})))

(defn- read-clstr [cdhit-clstr]
  (->> (io/read-lines cdhit-clstr)
       (map #(str/split #"\t" %))))

(defn make-clstr-fna

  [infasta-csv cdhit-clstr]
  (let [data (create-map infasta-csv)
        clusters (read-clstr cdhit-clstr)] 
    (pp/pprint (take 3 (rest clusters)))
    (->> (add-seq-to-cluster data clusters)      
         (pxmap (fn [[rep-name vs]] 
                  [rep-name ;[cluster# rep seq name]
                   (distinct-hits (map (fn [m] [(m :id) (m :seq)]) vs))])
                15))))

(defn make-cluster-fna
  "makes fnas for each cluster found from cdhit by copying from the
  temp files. example out-prefix is
  \"/home/peis/data/s15-round11-cdhit-cluster-\" and clusters file is
  \"/home/peis/data/s15-round11-qual-filter.cdhit.clusters\" "
  
  [n out-prefix clusters]
  (let [cluster-fnas (->> (make-clstr-fna S15-round11-qual-filter-csv clusters)
                          (filter #(>= (-> % second count) n)))
        outfiles (map (fn [[k hits]]
                        (let [n (first k)] ;cluster number
                          (prn :k k :n n :count (count hits))
                          (str out-prefix n ".fna")))
                      cluster-fnas)]
    (dorun
     (map (fn [out [_ data]]
            (write-fasta out data))
          outfiles cluster-fnas))))
;;====================================================================

(defn run-clustalw [fna]
  (runx "clustalw" (str "-INFILE=" fna) "-TYPE=DNA" "-QUIET"))
