(ns analysis.combine-cdhit-doubleton
  (:require [clojure.java.jdbc :as jdbc]
            [clojure.core.reducers :as r]
            [clojure.contrib.io :as io]
            [clojure.contrib.string :as str]
            [edu.bc.fs :as fs]
            [clojure.java.shell :as shell]
            [clojure.pprint :as pp]
            [edu.bc.bio.sequtils.files :refer (read-seqs)])
  (:use edu.bc.utils.fold-ops
        edu.bc.utils.probs-stats
        hts-exploration.db-queries
        hts-exploration.globals
        hts-exploration.utils))

(let [get-seqs (fn [r] (->> (round-all-usable-seqs r) sql-query
                          (reduce (fn [V {:keys [id hit_seq]}]
                                    (let [sq (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" hit_seq)]
                                      (if sq  ;standard sequence
                                        (conj V [id (second sq)])
                                        V)))
                                  [])))]
  (io/with-out-writer "data/150122/s15-all-rounds-qual-filter.fna"
    (doseq [i [4 9 10 11]
            [nm s] (get-seqs i)]
      (println (str ">" nm))
      (println s))))

(def id->cluster
  (with-open [r (clojure.java.io/reader "data/150122/out.cluster")]
    (let [l (line-seq r)]
      (->> (rest l)
           (map #(str/split #"\t" %) )
           (reduce (fn [M x]
                     (let [clstr (read-string (second x)) ;cluster
                           size (read-string (nth x 2))
                           ident (read-string (nth x 5))
                           id (read-string (first x))]
                       (assoc M id [clstr size ident])))
                   {} )
           ))))

(def seq->id 
  (with-open [r (clojure.java.io/reader "data/150122/s15-all-rounds-qual-filter.fna")]
    (let [l (line-seq r)]
      (->> (partition-all 2 l)
           (reduce (fn [M v]
                     (assoc M (second v)
                            (read-string (str/drop 1 (first v)))))
                   {} )))))

(def most-freq-seqs
  (letfn [(get-seqs [r]
            (->> (round-all-usable-seqs r) sql-query
                 (reduce (fn [M {:keys [id hit_seq]}]
                           (let [sq (second 
                                     (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" hit_seq))]
                             (if sq  ;standard sequence
                               (assoc M sq (inc (get M sq 0)))
                               M)))
                         {})))]
    (->> [4 9 10 11]
         (map get-seqs)
         (apply merge-with +)
         (remove #(< (second %) 100))
         (sort-by second >)
         (map (fn [[s cnt]] (let [id (seq->id s)] [s cnt id (id->cluster id)])))
         (take 10))))

(def tested-seqs {:5motif {:id 3196036 :cluster 62979 :csize 1
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

;;largest clusters by doubletons
(map (fn [[s cnt]]
       [s (id->cluster (seq->id s)) cnt])
     most-freq-seqs)

;;doubleton cluster ids
(->> (map (fn [[s cnt]]
            [s (id->cluster (seq->id s)) cnt])
          most-freq-seqs)
     (map #(-> % second first))
     (str/join ","))

;;getting the percent identity of these doubleton clusters
(->> (map (fn [[s cnt]]
            [s (id->cluster (seq->id s)) cnt])
          most-freq-seqs)
     (map (fn [[s x cnt]] [(first x) (get (cluster->id) (first x))]))
     (map (fn [[clstr x]] [clstr (count x) (mean (map #(-> % second last) x))])))

;;largest clusters by cdhit
(->> (vals id->cluster) (map #(vec (take 2 %))) (into {}) (sort-by second >) (take 10))

(defn id->seq [ids]
  (if (> (count ids) 10000)
    (apply concat (map id->seq (partition-all 10000 ids)))
    (let [q (str "select sr.selex_id as id, 
                SUBSTRING(sr.sequence, sr.usable_start+1, sr.length) as hit_seq,
                  sr.const_start as cs from selex_reads as sr 
                WHERE selex_id IN (" (str/join "," ids) ");")]
     (jdbc/with-db-connection [c mysql-ds]
       (jdbc/query c [q])))))


(let [cluster->id
      (group-by (fn [[id [clstr size]]] clstr) id->cluster)
      other-clusters
      (->> (map (fn [[clstr x]] [clstr (mapv first x)]) cluster->id)
           (sort-by #(-> %  second count) >)
           (take 10))]
  (->> other-clusters
       (reduce (fn [M [clstr x]]
                 (let [sqs (take 10 x)]
                   (assoc M clstr (id->seq sqs))))
               {} )
       (reduce (fn [V [clstr x]]
                 (conj V [clstr
                          (->> (map (juxt :id :hit_seq) x)
                               (mapv (fn [[id s]]
                                       (second (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" s)))))]))
               [])
       pp/pprint))

;;identify large clusters with high %identity and >100 seqs
(let [interesting-clusters [2 0 7 8 3 11 15 6 34 17
                            9391 1920 464 2361 2371
                            38794 8747 20798 2917 2361]
      cluster->id
      (group-by (fn [[id [clstr size ident]]] clstr) id->cluster)]
  (->> cluster->id
       (map (fn [[clstr x]]
              [clstr (-> x first second second);size
               (mean (map #(-> % second last) x))]));mean %identity
       (filter #(and (> (second %) 100)
                     (> (last %) 90)))
       (sort-by (juxt second first) )))

(let [cluster->id
      (group-by (fn [[id [clstr size]]] clstr) id->cluster)
      largest-clusters
      (->> (map (fn [[clstr x]]
                  [clstr
                   (mapv (juxt first #(-> % second last)) x)])
                cluster->id)
           (sort-by #(-> % second count) >)
           #_(take-while #(> (second %) 1000));clusters larger than 1000
           #_(take-while #(not= (first %) 9391));most freq cluster id
           (map (fn [[clstr x]] [clstr (count x) (mean (map second x))])))]
  (io/with-out-writer "data/150122/largest-clusters.csv"
    (println "cluster,size,per identity")
    (doseq [i largest-clusters]
      (println (str/join "," i)))))

(defn cluster->id []
  (group-by (fn [[id [clstr size]]] clstr) id->cluster))

(defn clusters-ids->seqs 
  "Maps clusters-ids to a set of seqs and returns a vector
  of [cluster-id [seqs]]"
  
  [clusters]
  (->> clusters
       (map (juxt :clstr #(->> % :sids id->seq)))
       (reduce (fn [V [clstr x]]
                 (conj V [clstr
                          (->> (map (juxt :id :hit_seq) x)
                               (mapv (fn [[id s]]
                                       [id (second (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" s))])))]))
               [])))

(let [other-clusters (create-cluster-entries)]
  (->> other-clusters
       (filter #(and (> (:size %) 1000) (> (:ident %) 90)))
       clusters-ids->seqs
       (map (juxt first #(-> %
                             second distinct
                             count)))));[cluster unique-sqs]

(defrecord Cluster [clstr ident size sids])

(defn create-cluster-entries []
  (map (fn [[clstr x]]
         (let [sids (mapv first x)
               size (-> x first second second)
               mean-ident (mean (mapv #(-> % second last) x))]
           (Cluster. clstr mean-ident size sids))) 
       (cluster->id)))

(defn sqs->kmer 
  "Turns a set of sequences into a kmer probs according to the
  probability that it will appear in the cluster"
  
  [k sqs]
  (->> sqs
       (map #(freqn k %))
       (apply merge-with +)
       probs))

(defn cluster->kmer
  "Takes a cluster number and turns the cluster into a set of
  sequences and then calculates the kmer probs"
  
  [k cl]
  (->> (create-cluster-entries)
       (filter #(= cl (:clstr %)))
       clusters-ids->seqs;get seqs in cluster
       first;only 1 entry since its a specific cluster
       ((juxt first #(map (fn [[id s]] s) (second %))))
       second
       (sqs->kmer k)))

(let [other-clusters
      (->> (cluster->id)
           (map (fn [[clstr x]] [clstr (mapv first x)]) )
           (sort-by #(-> %  second count) >)
           (take 10)
           clusters-ids->seqs
           (map (juxt first #(map second %)))
           (map #(vector (first %) (->> % second (sqs->kmer 2))) ))
      interest (get-in tested-seqs [:ci7-73 :cluster])
      x (->> [[interest ((cluster->id) interest)]]
             (map (fn [[clstr x]] [clstr (mapv first x)]))
             clusters-ids->seqs
             (map (juxt first #(map second %)))
             first
             second
             (sqs->kmer 2)) ]
  (map #(vector (first %) (merge-with / x (second %))) other-clusters))

(let [x (cluster->kmer 2 (get-in tested-seqs [:wka :cluster]))
      y (cluster->kmer 2 (get-in tested-seqs [:bob :cluster]))]
  (pp/pprint
   (merge-with / y x)))

(let [other-clusters (create-cluster-entries)
      filtered-clusters
      (->> other-clusters
           (filter #(and (> (:size %) 1000) (> (:ident %) 90)))
           clusters-ids->seqs
           (map (juxt first (fn [[_ sqs]]
                              (->> (map (juxt second first) sqs)
                                   (into {}))))))]
  (doseq [f filtered-clusters
          :let [clstr (first f)
                cname (edu.bc.fs/join "data/150122/cluster" clstr)]]
    (when-not (edu.bc.fs/directory? cname)
      (edu.bc.fs/mkdirs cname))
    (io/with-out-writer (str cname "/cluster-seqs.fna")
      (doseq [[sq id] (second f)]
        (println (str ">" id))
        (println (str prime5-const sq cds))))))

(let [wdir "data/150122/cluster"
      dirs (edu.bc.fs/listdir wdir)]
  (doseq [dir dirs
          :let [dir (edu.bc.fs/join wdir dir)
                infasta (edu.bc.fs/join dir "cluster-seqs.fna")
                fasta ((shell/sh "cat" infasta) :out)
                rnafold-out (shell/sh "RNAfold" "-p"
                                      "-P" "/usr/local/ViennaRNA/misc/rna_andronescu2007.par"
                                      "--noPS"
                                      :in fasta
                                      :dir dir)]]
    (prn "working on" infasta)))

;;clojure.java.io alternative
(with-open [f (clojure.java.io/reader "data/150122/cluster/9391/4278649_dp.ps")]
  (let [ex-file (->> (line-seq f)
                     (drop-while (complement #(= % "%start of base pair probability data"))))
        bp (->> (filter #(re-find #"ubox" %) ex-file)
                (map #(vec (str/split #" " %)))
                (reduce (fn [M x]
                          (let [i (read-string (x 0))
                                j (read-string (x 1))
                                p (read-string (x 2))]
                            (assoc M [i j] (* p p))))
                        {}))]
    bp))

(defn parse-ps [f]
  (with-open [infile (clojure.java.io/reader f)]
    (let [ex-file (line-seq infile)
          sq (->> (drop-while (complement #(= % "/sequence { (\\")) ex-file)
                  second
                  (str/replace-re #"U" "T")
                  (str/butlast 1)
                  (str/drop 15)
                  (str/butlast 15))
          bprobs (drop-while (complement #(= % "%start of base pair probability data")) ex-file)
          ]
      [sq
       (fs/basename f)
       (->> (filter #(re-find #"ubox" %) bprobs)
            (map #(vec (str/split #" " %)))
            (reduce (fn [M x]
                      (let [i (read-string (x 0))
                            j (read-string (x 1))
                            p (read-string (x 2))]
                        (assoc M [i j] (* p p))))
                    {})
            )])))

(let [
            b (jna/make-cbuf (* 16 12))
            p (jna/pointer b)]
        (jna/invoke Float RNA/pf_fold "GGGGATTACCCC" b)
        (.getString p 0)
        )

;;proof this works and we can fold many seqs
(let [bs [\A \C \G \T]
      rand-seq (fn [n] (apply str (repeatedly n #(rand-nth bs))))
      ]
  (clojure.pprint/pprint
   (take-last 10
              (doall
               (map (fn [s]
                      (let [l (count s)
                            b (jna/make-cbuf (inc (* 16 l)))
                            p (jna/pointer b)
                            _ (jna/invoke Integer RNA/read_parameter_file "/usr/local/ViennaRNA/misc/rna_andronescu2007.par")
                            _ (jna/invoke Integer RNA/free_pf_arrays)
                            _ (jna/invoke Float RNA/pf_fold s b)
                            bp-array (.getDoubleArray
                                      (jna/invoke com.sun.jna.Pointer RNA/export_bppm) 0 (* l l))
                            iindx (.getIntArray (jna/invoke com.sun.jna.Pointer RNA/get_iindx l) 0 l)]
                        [s (.getString p 0)]))
                    (map (fn [_] (rand-seq 57)) (range 10000))))))
  )

;;trying to get bp probs from the matrix and storing in a hashmap
(let [bs [\A \C \G \T]
             rand-seq (fn [n] (apply str (repeatedly n #(rand-nth bs))))
             ]
         
         (->> (map (fn [_] (rand-seq 30)) (range 100))
              distinct 
              (map (fn [s] 
                     (let [l (count s)
                           b (jna/make-cbuf (inc (+ 32 (* 2 l))))
                           p (jna/pointer b)
                           _ (jna/invoke Integer RNA/read_parameter_file "/usr/local/ViennaRNA/misc/rna_andronescu2007.par")
                           _ (jna/invoke Long RNA/free_pf_arrays)
                           _ (jna/invoke Double RNA/pf_fold s b)
                           bp-array (.getDoubleArray
                                     (jna/invoke com.sun.jna.Pointer RNA/export_bppm) 0 (* l l))
                           iindx (fn [i len] (+ (bit-shift-right (* (- (inc len) i) (- len i)) 1) len 1))
                           foo (.getIntArray (jna/invoke com.sun.jna.Pointer RNA/get_iindx l) 0 l)]
                       (when (re-find #"[\(\{\}\)]" (.getString p 0))
                         [s (->> (for [i (range 1 l)
                                       j (range 3 (inc l))
                                       :let [bprob (aget bp-array (- (iindx i l) j))
                                             bprob2 (aget bp-array (- (aget foo i) j))]
                                       :when (and (> (Math/sqrt bprob) 1e-5)
                                                  (< i j))]
                                   [[i j] [bprob bprob2]])
                                 (into {}))]))))
              doall
              last)
         42
         )

(let [other-clusters
      (->> (cluster->id)
           (map (fn [[clstr x]]
                  (let [sids (mapv first x)
                        size (-> x first second second)
                        mean-ident (mean (mapv #(-> % second last) x))]
                    [clstr mean-ident size sids ])) )
           )
      filtered-clusters
      (->> other-clusters
           (filter #(and (> (nth % 2) 1000) (> (second %) 90)))
           (map (juxt first last))
           clusters-ids->seqs
           (map (juxt first (fn [[_ sqs]]
                              (->> (map second sqs)
                                   frequencies
                                   (sort-by second >)
                                   first)))))]
  filtered-clusters)

(let [other-clusters (create-cluster-entries)
      filtered-clusters
      (->> other-clusters
           (filter #(and (> (:size %) 1000) (> (:ident %) 90)))
           clusters-ids->seqs
           (map (juxt first (fn [[_ sqs]]
                              (->> (map second sqs)
                                   frequencies probs
                                   (sort-by second >)
                                   (take-while #(> (second %) 0.001)))))))]
  (io/with-out-writer "data/150122/cluster-seq-comp.csv"
    (println "cluster,seq,fraction of cluster")
    (doseq [[c xs] filtered-clusters
            x xs]
      (println (str/join "," (cons c x))))))

(let [other-clusters (create-cluster-entries)
      filtered-clusters
      (->> other-clusters
           (filter #(and (> (:size %) 1000) (> (:ident %) 90)))
           clusters-ids->seqs
           (map (juxt first (fn [[_ sqs]]
                              (->> (map second sqs)
                                   frequencies probs
                                   (sort-by second >)
                                   (take-while #(> (second %) 0.01)))))))
      xs (->> (map second filtered-clusters)
              (apply concat)
              (map first))
      ]
  (io/with-out-writer "data/150122/top-cluster-top-seq-dist.csv"
    (println (str/join "," (for [[c xs] filtered-clusters
                                 [x _] xs]
                             (str x "-" c))))
    (doseq [i (pairwise-distance #(double
                                   (/ (levenshtein %1 %2)
                                      (max (count %1)(count %2)))) xs)]
      (println (str/join "," i)))))

(defn ensemble-dist
  "Calculate the distance between 2 ensembles x and y. "
  
  [x y]
  (->> (set (concat (keys x) (keys y))) ;S
       (map (fn [k]
              (let [diff (- (get x k 0)
                            (get y k 0))] ;(PAij-PBij)
                (* diff diff))));sq
       (reduce +)));sum

(defn pairwise [f coll]
  (loop [x coll
         Vx []]
    (if (seq x)
      (recur (rest x)
             (concat Vx
                     (reduce (fn [Vy y]
                               (conj Vy (f (first x) y)))
                             [] (vec (rest x)))))
      Vx)))

(def intracluster-dist
  (let [wdir "data/150122/cluster"
        dirs (->> (edu.bc.fs/listdir wdir)
                  (map #(fs/join wdir %)))]
    (doall
     (pmap (fn [dir] 
             (let [probs (for [f (fs/directory-files dir ".ps")]
                           (last (parse-ps f)))]
               (vector dir (mean (pairwise ensemble-dist probs)))))
           dirs))))

;;intracluster data file
(with-open [out (clojure.java.io/writer "data/150122/intracluster-distance.csv")]
  (.write out "cluster,ensemble distance\n")
  (doseq [[clstr dist] intracluster-dist]
    (.write out (str (fs/basename clstr) "," dist "\n"))))

(let [other-clusters (create-cluster-entries)
      filtered-clusters
      (->> other-clusters
           (filter #(and (> (:size %) 1000) (> (:ident %) 90)))
           (map (juxt first last))
           clusters-ids->seqs
           (map (juxt first (fn [[_ sqs]]
                              (->> (map second sqs)
                                   frequencies
                                   (sort-by second >)
                                   first)))))
      foo (->> filtered-clusters
               (map (juxt first #(-> % second first)) );[cluster-id seq]
               (map (fn [[c sq]]
                      (let [ps (->> (fs/directory-files (fs/join "data/150122/cluster/" c) ".ps")
                                    (map parse-ps);[sq {probs}]
                                    (filter #(= (first %) sq))
                                    first)]
                        (conj ps c))))
               )]
  (pp/pprint
   (map vector (pairwise vector (map #(nth % 3) foo))
        (pairwise vector (map #(% 0) foo))
        (map #(apply levenshtein %) (pairwise vector (map #(% 0) foo)))
        (pairwise ensemble-dist (map #(nth % 2) foo)))))

;;printing the above results into a file
(with-open [out (clojure.java.io/writer "data/150122/intercluster-dist.csv")]
  (let [;;last part of the above section to get results
        xs (map #(->> % flatten (str/join ","))
                (map vector (pairwise vector (map #(nth % 3) foo))
                     (pairwise vector (map #(% 0) foo))
                     (map #(apply levenshtein %) (pairwise vector (map #(% 0) foo)))
                     (pairwise ensemble-dist (map #(nth % 2) foo))))]
    (.write out "cluster1,cluster2,seq1,seq2,levenshtein dist,ensemble dist\n")
    (doseq [i xs] (.write out (str i "\n")))))
