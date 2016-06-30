(ns edu.bc.utils.fold-ops
  (:require [clojure.contrib.io :as io]
            [clojure.string :as str]
            [clojure.java.shell :as shell]
            [me.raynes.conch :refer [programs let-programs with-programs] :as sh2]
            [edu.bc.fs :as fs])
  (:use refold
        [edu.bc.bio.sequtils.files :only [join-sto-fasta-file]]
        [slingshot.slingshot :only [throw+]]
        [edu.bc.utils :only [assert-tools-exist getenv]
         :exclude [get-tool-path]]))

(def ^{:private true} viennadir
  (if (fs/directory? "/usr/local/ViennaRNA/")
    "/usr/local/ViennaRNA/" 
    (fs/join (fs/homedir) "/bin/ViennaRNA/")))

(def ^{:private true} param-file
  (let [pfile (str viennadir "rna_andronescu2007.par")
        pfile (if (fs/exists? pfile)
                pfile
                (str viennadir "misc/rna_andronescu2007.par"))]
    (if (fs/exists? pfile)
      pfile
      (throw+ {:file pfile} "parameter file not exist" ))))

(defn- get-tool-path
  "Changes in the ViennaRNA file structure. This allows for trying to
  create possible custom locations."

  [toolset-type]
  (let [path-loc (if (fs/exists? (fs/join viennadir "Progs"))
                   (fs/join viennadir "Progs")
                   (fs/join viennadir "bin"))]
    (case toolset-type
      :RNAfold
      (or (getenv "RNAfold")
          path-loc)
      :RNAalifold
      (or (getenv "RNAalifold")
          path-loc)
      :RNAdistance
      (or (getenv "RNAdistance")
          path-loc)
      :RNAinverse
      (or (getenv "RNAinverse")
          path-loc))))

(def ^{:private true} tools
  ["RNAfold" "RNAalifold" "RNAdistance" "RNAinverse"])
(def ^{:private true} tools-exist?
  (assert-tools-exist
   (mapv #(fs/join (get-tool-path (keyword %)) %) tools)))

(defn inverse-fold
  "Given a target structure, it will use RNAinverse to find n
   sequences which fold into a similar structure. If :perfect? is
   true, only returns sequences which fold into identical structures
   else returns the first n sequences. Returns a list of sequences."
  
  [target n & {:keys [perfect? ncore]
               :or {perfect? false ncore 2}}]
  (let [inv-fold (fn [target n perfect?]
                   (->> (map (fn [[s ensemble]]
                               (if perfect?
                                 (when-not (re-find #"d=" s) (re-find #"\w+" s)) ;perfect match
                                 (re-find #"\w+" s))) ;take all output
                             ;;calls the RNAinverse to generate inverse-fold seqs
                             (->> ((shell/sh "RNAinverse"
                                             "-Fmp"
                                             (str "-R" n)
                                             "-P" param-file
                                             :in target)
                                   :out)
                                  str/split-lines
                                  (partition-all 2)))
                        flatten
                        (remove nil? ))) ;imperfect matches removed if
                                        ;they were nil
        ;;generate the proper number of distinct inverse-fold sequences
        inv-seq (loop [c 0
                       cand []]
                  (if (< c n)    
                    (recur (count cand) ;distinct candidate seqs
                           ;;add current list to newly generated ones
                           (->> (pmap (fn [_] (inv-fold target (min 10 (quot n ncore)) perfect?)) (range ncore))
                                (apply concat cand )
                                distinct))
                    (take n cand)))] 
    inv-seq))

(defn struct->matrix
  "creates array of bp locations. Array resembles a hash-map where the
  keys are base-pair locations and the value is 1 if
  present. Locations not present are not represented in the hash-map."
  
  [st]
  (reduce (fn [m kv] ;creates array of bp locations
            (assoc m kv 1))
          {} (make-pair-table st)))

(defn calc-bp-prob
  "Takes a collection of structures and finds the bp probs for each
  base-pair in the structure"

  [structures]
  (let [n (count structures)
        map-structures (map struct->matrix structures)]
    (reduce (fn [m [k v]]
              (assoc m k (/ v n)))
            {} (apply merge-with + map-structures))))

(defn structs->centroid
  "Converts a list of structures into a centroid structure"
  
  [structures]
  (let [len (count (first structures))]
    (->> (calc-bp-prob structures)
         (filter (fn [[[i j] p]] 
                   (and (< (+ i 3) j) ;proper base-pair
                        (>= p 0.5)))) ;appears >=50% of the time
         (reduce (fn [m [[i j] p]]
                   (assoc m i "(" j ")"))
                 (apply sorted-map (interleave (range len) (repeat "."))))
         vals
         (apply str ))))

(defn centroid-n-structs
  "Finds the centroid structure of suboptimal structures and a vector
  representation using 0's and 1's using a RNAmutants or RNAsubopt. s
   is the RNA sequence (case doesn't matter as it will be all
   upper-cased) and n is the number of suboptimal structures to
   consider."
  
  [s n]
  (let [;;s "AACGAUCCCGGG"
        ;;n 10
        s (.toUpperCase s)
        structures (do (declare fold)
                       (fold s {:foldmethod :RNAsubopt :n n}))]
    (let [centroid (structs->centroid structures)]
      [centroid
       (map #(if (= \. %) 0 1) (seq centroid))]))) ;vector representation

(ns-unmap 'edu.bc.utils.fold-ops 'fold)
(defmulti
  ^{:doc "folds an RNA(s) (inseqs) into a structure using the method
               specified by :foldmethod. possible :foldmethod
               are :RNAfold, :RNAfoldp, :RNAsubopt, :centroid."
    :arglists '([inseqs & args])}

  fold
  (fn [inseqs & args]
    ((or (first args) {}) :foldmethod)))

(defmethod fold :RNAfold [inseqs args]
  (let [parser (fn [out]
                 (->> out
                      rest
                      (take-nth 2)
                      (map #(rest (re-find #"([\.|\(]\S*[\.|\)]) \((.*)\)" %))) ;split struct
                                        ;and energy
                      (map (juxt first #(Double/parseDouble (second %))))))
        ]
    (sh2/with-programs [RNAfold cat]
      (parser (RNAfold "--noPS" "-P" param-file {:in (cat {:in inseqs}) :seq true})))))

(defmethod fold :RNAfoldp [s args]
  (let [ensemble-div  (->> ((shell/sh "RNAfold"
                                      "-p"
                                      "-P" param-file
                                      "--noPS"
                                      :in s )
                            :out)
                           str/split-lines
                           last
                           (re-find #"ensemble diversity (\d*.\d*)")
                           second
                           Double/parseDouble)]
    ;;(when (fs/exists? "dot.ps") (do (prn "del dot") (fs/rm "dot.ps")))
    ensemble-div))

(defmethod fold :RNAsubopt [s args] 
  (sh2/with-programs [RNAsubopt cat]
    (->> (RNAsubopt "-p" (:n args) "-P" param-file {:in (cat {:in s}) :seq true})
         (remove (partial re-find #"[^\(\)\.]")))))

(defmethod fold :centroid [s args]
  (first (centroid-n-structs s (args :n))))

(defmethod fold :default [s]
  (let [pfold (fn [n s] (->> ((juxt #(/ (count %) n) identity) s)
                            (apply partition-all  )
                            (pmap fold )
                            (apply concat )))]
    (cond
      (>= (count s) 100000) (pfold 30 s)
      (>= (count s) 10000) (pfold 10 s)
      (>= (count s) 1000) (pfold 2 s)
      :else (fold s {:foldmethod :RNAfold}))))

(defn fold-fasta
  "Uses RNAfold to fold all sequences in a given fasta file and
  returns the structures. This method should be used when the number
  of sequences to fold becomes larger (ie > 200)."

  [infasta]
  (let [fasta ((shell/sh "cat" infasta) :out)
        rnafold-out ((shell/sh "RNAfold"
                               "-P" param-file
                               "--noPS"
                               :in fasta)
                     :out)]
    (->> rnafold-out
         (str/split-lines)
         (partition-all 3)
         (map last)
         (map #(str/split % #" "))
         (map first))))

(defn align-fold
  "Takes a fasta file and outputs an alignment with structure in
  stockholm file format. The outfile is the input file renamed to have
  a .sto extension. Uses mxscarna to do the fold/align"

  ([fna]
     (align-fold fna (fs/replace-type fna ".sto")))
  
  ([fna outfile]
     (let [tmp (fs/tempfile)
           call (shell/sh "mxscarna" "-stockholm" fna)]
       (io/with-out-writer tmp
         (-> call :out println))
       (join-sto-fasta-file tmp outfile :origin "#=GF AU mxscarna")
       outfile)))

;;;(ns-unmap 'edu.bc.utils.fold-ops 'fold-aln)
(defmulti

  ^{:doc "multimethod for folding an alignment. can use
  either :RNAalifold or :centroid_alifold as the first argument or
  none in args. Defaults to :RNAalifold. Always requires a file
  name (aln)."
    :arglists '([aln] [foldprogram aln])}
  
  fold-aln (fn [& args]
             (first args)))

(defmethod fold-aln :RNAalifold [_ aln]
  (-> ((shell/sh "RNAalifold"
                 "-P" param-file
                 "-r" "--noPS" aln)
       :out)
      str/split-lines
      second
      (str/split #"\s")
      first))

(defmethod fold-aln :centroid_alifold [_ aln]
  (-> ((shell/sh "centroid_alifold" aln)
       :out)
      str/split-lines
      last
      (str/split #"\s")
      first))

(defmethod fold-aln :default [aln]
  {:pre [(fs/exists? aln)]}
  (fold-aln :RNAalifold aln))

(comment 
  (defn fold-aln [aln]
    (-> ((shell/sh "RNAalifold"
                   "-P" param-file
                   "-r" "--noPS" aln)
         :out)
        (str/split-lines)
        second
        (str/split #"\s")
        first)))

(defn bpdist
  "finds the distance between 2 structures. uses tree edit distance by
  default. when bpdist = \"-DP\", uses base pair distance. Using
  \"-Xm\" will do bpdist in a pairwise manner. Using 'Xf' will compare
  each structure with the first one"
  
  [structs & args]
  (let [parser (fn [out] (map #(Integer/parseInt (re-find #"\d+" %)) out))]
    (sh2/with-programs [RNAdistance cat]
      (parser ((apply partial RNAdistance args) {:in (cat {:in structs}) :seq true})))))

(defn bpdist-fasta
  "Uses RNAdistance on folded sequences in a given fasta file
  containing name and structure and returns the base-pair
  distance. This method should be used when the number of sequences
  becomes larger (ie > 200)."

  [infasta]
  (let [fasta ((shell/sh "cat" infasta) :out)
        rnafold-out ((shell/sh "RNAdistance" (if bpdist "-DP" "")
                               "-Xf"
                               :in fasta)
                     :out)]
    (->> rnafold-out
         (str/split-lines)
         rest
         (partition-all 2)
         (map (fn [[nm dist]]
                [(apply str (rest nm))
                 (Integer/parseInt (re-find #"\d+" dist))]))
         (into {}))))

(defn centroid-fold [fasta]
  (-> ((shell/sh "centroid_fold"
                   :in fasta)
       :out)
      (str/split-lines)
      second
      (str/split #"\s")
      first))
