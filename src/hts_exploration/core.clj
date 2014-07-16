(ns hts-exploration.core
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.io :as io]
            [clojure.contrib.string :as str]
            [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :as math]
            [clojure.pprint :as pp]
            [iota]
            [foldable-seq.core :as fseq]
            gibbs-sampler
            [smith-waterman :as sw]
            )
  (:use edu.bc.utils.fold-ops
        edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.utils.probs-stats
        edu.bc.utils
        hts-exploration.db-queries
        hts-exploration.globals
        hts-exploration.utils))

(set! *print-length* 10)
*print-length*

(def data-r1 (iota/vec S15-r1-qual-filter-csv))



(defn get-name-seq [data & {:keys [info]
                            :or {info :both}}]
  (r/reduce (fn [V [nm s qual]]
              (let [s (re-find forward-re s)]
                (case info
                  :both (conj V [nm s])
                  :data (conj V s)
                  :name (conj V nm))))
            [] data))

(defn- sample-data
  "Samples n sequences from the data. Returns a reducible"

  [n data]
  (as-> data d
        (r/filter (fn [[nm s qual]] (re-find forward-re s)) d)
        (r/remove #(re-find #"[^ACGTU]" (second %)) d);remove nonstandard bases
        (into [] d)
        (repeatedly n #(rand-nth d))
        distinct
        (vec d)))

(defn sample-n
  "Samples n sequences from data and processes into [name seq struct]
  type format"
  
  [n data]
  (as-> data d
        (sample-data n d)
        (qual-remove d)
        (get-name-seq d :info :both)
        (add-structure 30 d)))

(defn similarity-remove

  [sampled-data]
  (loop [i sampled-data
         V []]
    (if (> (count (seq i)) 1)
      (let [[nm s st] (first i)]
        (recur (rest i)
               (if (->> (for [[_ _ stj] (rest i)] 
                          (normalized-dist st stj))
                        (every? #(>= % 0.2)))
                 (conj V [nm s st])
                 V)))
      (conj V (vec (first i))))))

(defn get-regions
  "Gets the interesting part of the sequence and returns a vector of
  vector pairs [name region-of-interest]"

  [from to reducible-data]
  (r/fold (fn ([] [])
            ([l r] (concat l r)))
          (fn ([] [])
            ([V [nm s]]
               (conj V [nm (subs s from to)])))
          reducible-data))



(defn struct-sim [n ref-structure]
  (->> data-r1
       (sample-data n )
       get-name-seq
       (add-structure 10)
       (dist-filter 10 0.2 ref-structure)
       (get-regions (- (count prime5-const) 3)
                    (+ 30 (count prime5-const) ;var+const region
                       (count prime3-const)))))

(defn motif-found [sim-seqs]
  (future
    (as-> (range 10) r                 ;multiple random gibbs samples
          (pmap (fn [_]
                  (gibbs-sampler/-main (->> sim-seqs
                                            (map second)
                                            (remove #(re-find #"[^ACGTU]" %))
                                            (mapv #(->> (subs % 3)
                                                        (str/take 30))))))
                r)
          (remove nil? r)
          (last r)                     ;pick a random gibbs sample
          (second r) ;last state
          (r :baseprob))))

(def struct-sim-bob
  ^{:doc "Pulling sequences which have simlar structure to bob"}
  (time 
   (let [bob-st (fold bob)]
     (struct-sim 1000000 bob-st))))

(def motif-found-bob
  @(motif-found struct-sim-bob))

(def struct-sim-ref-seq-firmicute
  ^{:doc "Pulling sequences which have simlar structure to Gk ref-seq"}
  (time 
   (let [firm-st (fold ref-seq-firmicute)]
     (struct-sim firm-st))))

(def struct-sim-ref-seq-ec
  ^{:doc "Pulling sequences which have simlar structure to Ec ref-seq"}
  (time 
   (let [ec-st (fold ref-seq)]
     (struct-sim 1000000 ec-st))))

(def motif-found-ref-seq-ec @(motif-found struct-sim-ref-seq-ec))





(defn  most-common-struct
  "Sample sequences from the fasta file and fold them. Then use this
  sequence set to compare against a sampling of 1M reads. This should
  find the most common structures because the sampled sequences will
  be from common structures. The data comes from reading in a file
  using iota/vec"

  [outfile data]
  (future
    (let [samples1 (similarity-remove (sample-n 200 data))
          samples2 (sample-n 1e6 data)
          _ (prn :samples (count samples1) (count samples2))
          xxx (r/reduce
               (fn [V [nm sample sample-st]]
                 (let [len (count sample-st)
                       dist-map (reduce (fn [M [nm d]]
                                          (->> (relative-dist d len)
                                               (assoc M nm)))
                                        {}
                                        (calc-dist 30 sample-st samples2))
                       dist-distribution (frequencies (vals dist-map))
                       sim (vec (filter (fn [[nm s]]
                                          (< (dist-map nm) 0.2))
                                        samples2))]
                   (conj V [nm sample sample-st dist-distribution sim])))
               [] samples1) ]
      (io/with-out-writer outfile
        (prn :samples (count samples1) (count samples2))
        (doseq [[nm sampled-seq sampled-struct distribution similar-st] xxx]
          (prn nm)
          (prn sampled-seq) 
          (prn sampled-struct) 
          (prn distribution)
          (prn similar-st))))))

(comment (most-common-struct "/home/peis/data/foo.test" (read-fasta-csv data-r1))) ;typical case

(defn read-common-struct-hits [f]
  (let [f (rest (io/read-lines f))]
     (reduce (fn [M x]
               (let [[name sample-seq sample-struct distribution hits]
                     (map read-string x)
                     hits (distinct hits)
                     key-vec [:sample-seq :sample-struct :distribution :hits]]
                 (assoc M name
                        (->> [sample-seq sample-struct distribution hits]
                             (interleave key-vec )
                             (apply assoc {})))))
               {} (partition 5 f))))

(defn get-centroid

  [hits]
  (prn :centroid)
  (structs->centroid
   (r/fold 50
           (fn ([] [])
             ([l r] (concat l r)))
           (fn ([] [])
             ([V s]
                (conj V (fold s))))
           hits)))

(defn determine-clusters
  "data comes from reading output from common-struct-hits"
  
  [data]
  (r/reduce (fn [M [k sample-seq hits]]
              (prn :seq k :count (count hits) :sample (type hits))
              (let [structs (map last hits)
                    vector-structs (map struct->vector structs)]
                (assoc M k {:sample-seq sample-seq
                            :n (count hits)
                            :hit-names (map first hits)
                            :hit-seqs (map second hits)
                            :hit-structs structs
                            :centroid (structs->centroid structs)
                            :vector-structs vector-structs})))
            {}
            (->> data
                 (r/map #(vector (first %)
                                 (-> % second :sample-seq)
                                 (-> % second :hits vec)))
                 (r/remove #(zero? (-> % last count))))))

(defn positional-H
  "after finding the most common seqs according to structure
  simlarity, we look at the entropy per position to see if there are
  any strongly conserved residues. (ie clusters=foo-clusters)"

  [clusters]
    (doseq [k (keys clusters)
            :let [m (clusters k )]
            :when (> (m :n) 5000)] 
      (pp/pprint
       (let [{centroid :centroid 
              sample-seq :sample-seq
              hit-seqs :hit-seqs} (select-keys m [:centroid :sample-seq :hit-seqs])
              most-freq (reduce (fn [max-entry entry]
                                  (let [cur-max (second max-entry)
                                        cur-cnt (second entry)]
                                    (if (> cur-cnt cur-max)
                                      entry
                                      max-entry)))
                                [0 0] (frequencies (map count hit-seqs)))
              hit-seqs (filter #(= (first most-freq) (count %)) hit-seqs)
              hit-seqs-T (->> (edu.bc.utils/transpose hit-seqs)
                              (drop 13)
                              (take 40))
              H (map-indexed (fn [i coll]
                               [(+ 14 i)
                                (round 3 (entropy (probs 1 coll)))
                                (repeatedly 10 #(markov-step (probs 1 coll)));(probs 1 coll)
                                ])
                             hit-seqs-T)]
         [k centroid sample-seq most-freq (count hit-seqs) H]))))
(comment
  (def foo-test (read-common-struct-hits "/home/peis/data/foo.test"))
  (def foo-clusters (determine-clusters foo-test))

  ;;after finding the most common seqs according to structure
  ;;simlarity, we look at the entropy per position to see if there are
  ;;any strongly conserved residues
  
  
  ;;discovery of TAATACGAC as a common motif in sequences with
  ;;deletions. To identify how common this or similar motifs are, use
  ;;the following:
  (let [f (iota/vec S15-round10-fastq-qual-filter-csv)
        motif "TAATACGAC"
        mutant-motif (as-> motif m
                           (mutant-neighbor 3 m)
                           (conj m motif)
                           (map re-pattern m))
        contain-motif? (fn [s]
                         (->> (for [re mutant-motif]
                                (re-find re s))
                              (not-every? nil?)))
        usable (->> (read-fasta-csv f)
                    (r/map second)
                    (r/map #(re-find forward-re %))
                    (r/remove nil?)
                    (into []))
        data (map contain-motif? usable)]
    [(frequencies (map count (filter #(re-find #"TAATACGAC" %) usable)))
     (->> (map #(when %1 (count %2)) data usable)
          (remove nil?)
          frequencies)
     (count (filter true? data))
     (count usable)])

  ;;Counts of finding similar seqs to the motif according to the
  ;;length of the sequence
  [{96 558, 97 819, 98 4102, 99 5707, 79 52528, 80 20062, 81 13060,
    82 10528, 83 3984, 84 3139, 85 3266, 86 4847, 87 2614, 88 1403,
    89 1397, 90 934, 91 1492, 92 712, 93 1047, 94 2989, 95 541}
   {96 601, 97 852, 98 4218, 99 5815, 79 54906, 80 22877, 81 16961,
    82 15359, 83 10791, 84 12939, 85 15941, 86 22035, 87 34923,
    88 8889, 89 3388, 90 1600, 91 1838, 92 920, 93 1183, 94 3151, 95 596}
   239783
   631958]

  (defn- mutant-motif [n motif]
    (as-> motif m
          (mutant-neighbor n m)
          (conj m motif)
          (map re-pattern m)))
  
  (defn count-motifs

    [motif f]
    (let [;f S15-round10-qual-filter-csv
          motif (re-pattern motif) ;GAGC
          usable (map second (get-usable f))
          data (frequencies (map count (filter #(re-find motif %) usable)))
          total-nt (reduce #(+ %1 (- (count %2) 3)) 0 usable)
          cnt-hits (reduce + (vals data) )
          expectation (/ total-nt 256)
          zscore (/ (- cnt-hits expectation)
                    (Math/sqrt (* total-nt (/ 256) (/ 255 256))))
          cnt-usable (count usable)]
      [:data data
       :total-nt total-nt
       :cnt-hits cnt-hits
       :expectation (double expectation)
       :zscore (double zscore)
       :cnt-usable cnt-usable]))


  ;;pick favorite parasite seqs
  (let [f S15-round11-qual-filter-csv
        usable (->> (get-usable f)
                    (r/map second) 
                    (r/filter parasite?)
                    (r/filter #(= (count %) 79))
                    (into []) )
        p-table (->> usable 
                     transpose
                     (map (fn [coll]
                            (probs 1 coll))))
        make-seq (fn [] (apply str (map #(markov-step %) p-table)))
        seq-set (repeatedly 100 make-seq)
        p-seq-set (repeatedly 100 #(rand-nth usable))
        ]
    (pp/pprint
     [(->> seq-set
           (pxmap (fn [s]
                    [s (mean (map #(hamming-dist s %) usable))])
                  30)
           (sort-by second )
           (take 10))
      (->> p-seq-set
           (pxmap (fn [s]
                    [s (mean (map #(hamming-dist s %) usable))])
                  30)
           (sort-by second )
           (take 10))
      (->> usable
           rfrequencies
           (sort-by second > )
           (take 10))]))
  
  )

(->> foo-test
     determine-clusters
     (filter #(> (->> % (drop 2) first) 10000))
     (map butlast)
     pp/pprint)

(io/with-out-writer "/home/peis/data/foo.csv"
  (doseq [i foo-clusters
          j (last i)]
    (println (str/join "," (concat (butlast i) j)))))


(time
 (let [data (->> (get-usable S15-round11-qual-filter-csv)
                 (map-indexed (fn [i [nm s]] 
                                [(str/take 19 (str i nm)) s]))
                 (into {}))
       clusters (->> (io/read-lines "/home/peis/data/s15-round11-qual-filter.cdhit.clusters")
                     (map #(str/split #"\t" %)))
       [header clusters] (split-at 1 clusters)
       header (map keyword (first header))
       xxx (second 
            (reduce (fn [[clstr-rep M] clstr]
                      (let [c (zipmap header clstr)
                            k (if (= "1" (c :clstr_rep)) (c :id) clstr-rep)
                            new-entry {:seq (data (c :id)) 
                                       :len (Integer/parseInt (c :length))
                                       :id (c :id) 
                                       :clstr-rep (Integer/parseInt (c :clstr_rep))}
                            old-entry (get M [(c :clstr) k] [])]
                        [k (assoc M [(c :clstr) k] (conj old-entry new-entry))]))
                    ["" {}] clusters))
       yyy 0 #_(r/reduce (fn [M [k vs]]
                           (prn "working on entry:" k (count vs))
                           (let [s1 (data (second k))
                                 new-entries (r/fold 100 concat
                                                     (fn ([] [])
                                                       ([V m]
                                                          (let [s2 (m :seq)
                                                                s2-aligned (last (sw/sw s1 s2 :global true))]
                                                            (conj V (assoc m :seq s2-aligned)))))
                                                     vs)]
                             (assoc M k new-entries)))
                         {} (r/take 2 (r/filter #(>= 2000 (-> % second count) 1000) xxx)))] 
   (prn :header header) 
   (pp/pprint (take 3 clusters))
   (let [c (filter #(>= (-> % second count) 10000) xxx)
         outfiles (repeatedly (count c) fs/tempfile)
         out (pmap (fn [outfile [k vs]] 
                     [outfile k (distinct-hits (map (fn [m] [(m :id) (m :seq)]) vs))])
                   outfiles c)]
     (doseq [[outfile k vs] out]
       (write-fasta outfile vs))
     (prn outfiles)
     (map (fn [[_ k vs]] [k (count vs)]) out))))

(dorun
 (map (fn [cluster-num f]
        (fs/copy f (str "/home/peis/data/s15-round11-cdhit-cluster-" (ffirst cluster-num) ".fna")))
      '([["20" "346150@HWI-ST1129:5"] 10513] [["3" "60035@HWI-ST1129:52"] 100543] [["7" "109128@HWI-ST1129:5"] 23402]
          [["13" "161701@HWI-ST1129:5"] 12255] [["39" "29626@HWI-ST1129:52"] 11058] [["1" "4103@HWI-ST1129:529"] 104857]
            [["0" "7@HWI-ST1129:529:H9"] 19095] [["10" "123410@HWI-ST1129:5"] 21696] [["5" "67664@HWI-ST1129:52"] 8451]
              [["2" "12981@HWI-ST1129:52"] 7681] [["6" "91594@HWI-ST1129:52"] 7838] [["8" "113894@HWI-ST1129:5"] 48164]) 
      '("/tmp/-fs-6029566329720862683" "/tmp/-fs-4491974394548404381" "/tmp/-fs-2679861326937564582" "/tmp/-fs-4394718262706558267"
        "/tmp/-fs-8209872026325559064" "/tmp/-fs-5751431934686883617" "/tmp/-fs-5742671127339322398" "/tmp/-fs-1653243009906879314"
        "/tmp/-fs-6262556634941290504" "/tmp/-fs-6007613895371995308" "/tmp/-fs-6286049487709373225" "/tmp/-fs-1124450394497733987")))

;;;generate random sequences
(defn rand-sequence
  ([n len]
     (rand-sequence n len (probs 1 "ACGT")))
  ([n len probs]
     (let [variable (fn [] (apply str (repeatedly len #(markov-step probs))))]
       (repeatedly n variable))))

(defn rand-background [n]
  (map #(str prime5-const
             %
             prime3-const
             cds)
       (rand-sequence n 30)))

(defn chi-sq [obs expected]
  (let [S (clojure.set/union (set (keys obs)) (set (keys expected)))]
    (sum
     (map (fn [k] (let [obs (get obs k 0)
                       exp (get expected k 0)]
                   (/ (Math/pow (- obs exp) 2)
                      exp)))
          S))))

(defn kmer-freq [size inseqs]
  (->> inseqs
       (r/fold (fn ([] {})
                 ([l r] (merge-with + l r)))
               (fn ([] {})
                 ([M s]
                    (merge-with + M (freqn size s)))))))

(defn kmer-probs [m]
  (let [tot (sum (vals m))]
    (reduce (fn [M [k v]]
                (assoc M k (/ v tot)))
              {} m)))

(doseq [infile (fs/re-directory-files "/home/peis/data" #"cluster*.fna")]
                          (let [f (->> (iota/seq infile)
                                       (partition-all 2) 
                                       (r/map second) 
                                       (into []))
                                kmer4 (kmer-freq f)
                                background (->> (rand-sequence (count f)) 
                                                vec
                                                (kmer-freq 4))
                                pback (kmer-probs background)]
                            (prn infile)
                            (->> (clojure.set/union (set (keys kmer4)) 
                                                    (set (keys pback)))
                                 (reduce (fn [M x]
                                           (->> (- (kmer4 x) (background x))
                                                (assoc M x )))
                                         {}) 
                                 (sort-by val >)
                                 (take 5) pp/pprint)
                            (prn :chisq (chi-sq kmer4 background))
                            (prn :test (chi-sq background background))
                            (prn (jensen-shannon (kmer-probs kmer4) pback))))

(let [f (->> (iota/seq "/home/peis/data/s15-round11-cdhit-cluster-17-cmsearch.fna" )
             (partition-all 2) 
             (r/map second) 
             (into []))
      kmer4 (kmer-freq 4 f)
      tot (sum (vals kmer4))
      background (->> (get-usable S15-round11-qual-filter-csv)
                      parasite-remove
                      (r/map second)
                      (kmer-freq 4)
                      kmer-probs
                      (reduce (fn [M [kmer pr]]
                                (assoc M kmer (* pr tot)))
                              {} ))]
  (->> (clojure.set/union (set (keys kmer4)) 
                          (set (keys background)))
       (reduce (fn [M x]
                 (->> (- (kmer4 x) (background x))
                      double (assoc M x )))
               {}) 
       (sort-by val >)
       pp/pprint)
  (prn :chisq (chi-sq kmer4 background))
  (prn :test (chi-sq background background)))

(let [usable (->> "/home/peis/data/s15-round11-cdhit-cluster-17.fna"
                  io/read-lines
                  (partition-all 2)
                  (map vec)
                  (into {}))
      cmsearch (->> (io/read-lines "/home/peis/data/test.cmsearch.out")
                    (remove #(.startsWith % "#"))
                    (map #(vec (str/split #"\s+" %))))
      keep (reduce (fn [M entry]
                     (let [nm (first entry)
                           qstart (entry 5)
                           qend (entry 6)
                           hstart (Integer/parseInt (entry 7))
                           hend (Integer/parseInt (entry 8))
                           score (Double/parseDouble (entry 14))
                           eval (Double/parseDouble (entry 15))]
                       (if (and (> hend 40)
                                (< eval 1e-7))
                         (assoc M (str ">" nm) [(dec hstart) (dec hend)])
                         M)))
                   {} cmsearch)
      out (->> (keys keep)
               (select-keys usable )
               (mapv (fn [[nm s]]
                       [(str ">a" (re-find #"\d+" nm))
                        (apply subs s (keep nm))]))
               distinct-hits
               (sort-by val))]
  (io/with-out-writer "/home/peis/data/s15-round11-cdhit-cluster-17-cmsearch-2.fna"
    (doseq [[nm s] out]
      (println nm)
      (println s))))

(let [data (->> "/home/peis/data/s15-round11-cdhit-cluster-17-cmsearch-2-052914.sto"
                read-sto 
                :seqs
                (map #(str/replace-re #"\." "" %))
                (map str/upper-case))
      get-energy (fn [s]
                   (->> ((clojure.java.shell/sh "RNAfold"
                                                "--noPS"
                                                :in s )
                         :out)
                        str/split-lines
                        second
                        (str/split #" ")
                        last
                        (re-find #"\-*\d*.\d+")
                        (Double/parseDouble)))
      zscore (fn [s] 
               (let [p (probs 1 s)
                     energies (map get-energy (rand-sequence 1000 (count s) p))]
                 (/ (- (get-energy s) (mean energies))
                    (sd energies))))
      z (pxmap zscore 10 data)]
  [(mean z) (/ -11.22 (mean (map get-energy data)))])
;;[-1.474822732828705 0.8164149910022694]

cmalign --noprob -o s15-round11-cdhit-cluster-17-cmsearch-2.sto test.cm s15-round11-cdhit-cluster-17-cmsearch-2.fna

;;;identification of high indels
(time 
 (let [stmt "SELECT sr.sequence, sr.usable_start, sr.usable_stop 
                                              FROM selex_reads as sr 
                                             WHERE sr.usable=1 
                                               AND sr.strand=1 
                                             LIMIT 100000;"
       st (str prime3-const cds)
       aligned-seqs (->> (get-db-seq stmt)
                         align-to-const)
       _ (prn :count (count aligned-seqs))
       xxx (->> (map #(apply count-indel %) aligned-seqs)
                (map vector aligned-seqs)
                (filter (fn [[_ [_ ins-lens _ del-lens]]]
                          (or (some #(> % 7) ins-lens)
                              (some #(> % 5) del-lens)))))]
   (prn :cntxxx (count xxx))
   (doseq [[i j] xxx
           :let [[s1 s2] i]]
     (prn s1)
     (prn s2)
     (println j))))

(def foo 
  (future 
    (time 
     (let [st prime3-const
           stmt "SELECT sr.sequence, sr.usable_start, sr.usable_stop 
                                              FROM selex_reads as sr 
                                             WHERE sr.usable=1 
                                               AND sr.strand=1;"
           aligned-seqs (-> stmt get-db-seq align-to-const)
           _ (prn :count (count aligned-seqs))
           [mrates lrates overall totals] (->> aligned-seqs count-mutations get-mutation-rates )
           mrates (mapv #(/ % (count st)) mrates)
           overall (/ overall (count st))
           t-matrix (mutation-matrix aligned-seqs)]
       {:overall overall
        :totals totals
        :mrates mrates
        :lrates lrates
        :tmatrix t-matrix}))))

(def mutation-rates-primer-region
  {:overall 0.007994311914897497,
   :totals {:sub 0.4827822120866591, :del 0.3372862029646522, :ins 0.1799315849486887},
   :mrates [1.2865038634616492E-4 4.794764044738977E-5 8.98790370363618E-5],
   :lrates [{3 0.0038022813688212928, 2 0.03802281368821293, 1 0.9581749049429658} {16 6.76132521974307E-4, 9 6.76132521974307E-4, 4 0.0074374577417173765, 3 0.016903313049357674, 2 0.08384043272481406, 1 0.8904665314401623}],
   :tmatrix {\A {\T 0.2690812995491994, \G 0.6841287113321933, \C 0.04678998911860718}, \C {\T 0.40845990527519194, \G 0.3512983831455169, \A 0.2402417115792912}, \G {\T 0.6201950659781984, \C 0.05220883534136546, \A 0.32759609868043604}, \T {\G 0.14688362632573004, \C 0.5754758099665843, \A 0.2776405637076856}, \- {\T 0.12369760116307243, \G 0.32117761085534285, \C 0.20002423067603584, \A 0.3551005573055488}}})

(def mutation-rates-const-region
  {:overall 0.08580696500078834,
   :totals {:sub 0.3287608603446824, :del 0.5896483534297089, :ins 0.08159078622560872},
   :mrates [0.028209971637225224 0.0070010577380476115 0.05059593562551548],
   :lrates [{9 1.9289193229493178E-5, 8 1.1573515937695906E-4, 7 2.1218112552442495E-4,
             6 0.0014466894922119883, 5 0.0074359839899696195, 4 0.009519216858754882,
             3 0.025558181029078458, 2 0.09354294256642716, 1 0.862149780585427}
            {1 0.758103669855776, 2 0.1681388884069711, 3 0.0528878811617451,
             4 0.01526715519783911, 5 0.0035205205779632494, 6 0.0011543784306058417,
             7 6.752780183659606E-4, 8 2.3354476919771366E-4, 9 4.003624614817948E-6,
             10 4.003624614817948E-6, 11 1.3345415382726495E-6, 13 2.669083076545299E-6,
             14 2.669083076545299E-6, 15 2.669083076545299E-6, 17 1.3345415382726495E-6}],
   :tmatrix {\A {\T 0.33469268407802627, \G 0.4642048803068977, \C 0.20110243561507593},
             \C {\T 0.3532442874834455, \G 0.3930922401972083, \A 0.25366347231934616},
             \G {\T 0.39803274512067444, \C 0.18781273262132134, \A 0.41415452225800414},
             \T {\G 0.22932979622636543, \C 0.5287716493922594, \A 0.24189855438137509},
             \- {\T 0.33998587240298017, \G 0.2370278272296705, \C 0.2566282771490891, \A 0.16635802321826024}}})

(defn simulate-seq []
  (let [primer-mrate mutation-rates-primer-region
        const-mrate mutation-rates-const-region
        bprobs (apply hash-map (interleave [\A \C \G \T] (repeat 0.25)))]
    (str/upper-case
     (str (simulate2 prime5-const primer-mrate)
          (generate-rand-seq 30 bprobs)
          (simulate2 prime3-const const-mrate)
          (simulate2 cds primer-mrate)))))

(def background-seqs (future (doall (repeatedly 1e6 simulate-seq))))

(def background-stats
  (let [cnt (count @background-seqs)
        kmerfn (fn [n inseqs]
                 (->> (map #(freqn n %) inseqs)
                      (apply merge-with +)
                      (reduce-kv (fn [M k v]
                                   (assoc M k (double (/ v (count inseqs)))))
                                 {})))
        [kmer1 kmer2 kmer3 kmer4 kmer5] (pmap #(kmerfn % @background-seqs) (range 1 6))]
    {:n cnt
     :kmer1 kmer1
     :kmer2 kmer2
     :kmer3 kmer3
     :kmer4 kmer4
     :kmer5 kmer5
    :len-dist (probs 1 (map count @background-seqs))}))

(let [sqs (->> (sql-query "SELECT SUBSTRING(sr.sequence, sr.usable_start+1, sr.length) as hit_seq FROM selex_reads as sr WHERE sr.usable=1 AND sr.strand=1 ;")
                                     (map :hit_seq))
                            re #"T[CGT][GC][TGA]T" ;#"A[ACT][CG][ACG]A"
                            fun (fn [sqs] (->> sqs
                                               ;(filter #(re-find re %) ) 
                                               (map (fn [s]
                                                      (->> (str-re-pos re s)
                                                           (remove #(<= 15 (first %) 45)))))))
                            normalize (fn [x] (mean (map count x)))]
                        [(normalize (fun sqs)) (normalize (fun @background-seqs))
                         (frequencies (mapcat #(map first %)(fun sqs)))
                         (reduce (fn [M [v k]]
                                   (assoc-in M [k v] (inc (get-in M [k v] 0))))
                                 {} (apply concat (fun sqs)))])

(let [sqs (->> (round-all-usable-seqs 11) sql-query (mapv :hit_seq ))
      kmer (pmap #(kmer-freqn % sqs) [2 3 4 5])
      chisq (fn [Mobs Mexp]
              (as-> (merge-with - Mobs Mexp) x
                    (merge-with (fn [num denom] (/ (sqr num) denom)) x Mexp)))]
  (map #(->> (chisq %1 %2) (sort-by val > ) (take 10))
       kmer
       (map background-stats [:kmer2 :kmer3 :kmer4 :kmer5])))

(def mean-sd ((juxt mean sd) (map #(-> % fold second) (take 10000 @background-seqs))))
(let [sqs (->> (round-all-usable-seqs 11) sql-query (mapv :hit_seq )
                                     (filter #(re-find #"A[ACT][CG][ACG]A.{3,}T[CGT][GC][TGA]T" %)))
                            foo (->> (map #(vector % (count (re-seq #"A[ACT][CG][ACG]A" %))) sqs)
                                     (group-by second)
                                     (sort-by key >)
                                     (map #(->> % second (sort-by count >)));sort on length
                                     )
                            trivial-zscore (fn [s] (/ (- (-> s fold second) (first mean-sd)) (second mean-sd)))
                            chisq (fn [Mobs Mexp]
                                    (as-> (merge-with - Mobs Mexp) x
                                          (merge-with (fn [num denom] (/ (sqr num) denom)) x Mexp)))]
                        (for [xxx foo]
                          (reduce (fn [V data]
                                    ;;data is formated [seq motif-cnt [chisq-kmer2 chisq-kmer3 ... kmer5]]
                                    (if (and (<= (-> data first count) 87)
                                             (>= (-> data last first) 25)
                                             (>= (-> data last second) 83)
                                             (>= (-> data last third) 293))
                                      (let [dist (levenshtein (first data) bob)
                                            z (trivial-zscore (first data))]
                                        (if (<= dist 27) 
                                          (->> (conj V (conj data dist z))
                                               (sort-by last )
                                               (take 3))
                                          V))
                                      V)) []
                                  (map (fn [[s cnt]]
                                         [s cnt 
                                          (map #(apply + (map second (chisq %1 %2)))
                                               (map #(kmer-freqn % [s]) [2 3 4 5])
                                               (map background-stats [:kmer2 :kmer3 :kmer4 :kmer5]))])
                                       xxx))))

["TGCGTAACGTACACTATCGAAAGGAGAATGGAATCGAGCAATCGATCATTCTATCTTAGGATTTAAAAATGTCTCTAAGTACT" 5 ;pick this one
  (20.58710780740946;most deviation from expected dimer count...
   100.35018189484066
   389.83382282396497
   1514.1072691229554) 22]
["TGCGTAACGTACACTGTCAAACAAAACAAAAGACGAAGACGCACTTCATTCAATACTTGGACTCTTAAAATGTCTCTAAGTACT" 4 (25.450923249798738 89.96136254358372 403.93468076920476 1608.6484871597056) 27 -0.5247170019887462]
["TGCGTAACGTACACTCACGAAGAGGACGGAAGACAGATGAAGAGCTCGTTCTATACTTTGGAGATTTAAAATGTCTCTAAGTACT" 3 (29.843538435035306 93.35265995533078 394.51988629953615 1588.9677584109838) 22 -0.8349109319076299]
["TGCGTAACGTACACTAACGATTCGAAAGTGAAAGAAAGAGAAATCATTCTTATACTTTGGAGAGTTAAAATGTCTCTAAGTACT" 2 (27.700575628411325 86.44315532626388 414.62469863857905 1431.9806485636368) 22 -1.004107620954294]
["TGCGTAACGTACACTGGGAGCCGCCCCACCCAGGCGCCCTCGGTGTCATTCTATACGCTTTGGGGTTTTAAAATGTCTCTAAGTACT" 1 (34.66052698680395 96.3740191847706 334.6083922999456 1343.6005201916003) 23 -3.1754651303864803]
["TGCGTAACGTACACTGCGGACAGCGAGACAGATCGAAGGTTTTGATCATTCTATACTTTGGATTTTAAAATGTCTCTAAGTACT" 3955];most frequent in round11
["TGCGTAACGTACACTGCGGGAACAGACCCAACCTACCCTGCGGTGTCTTCTATTACTTTGAGTTTTAAAATGTCTCTAAGTACT" 1085];second
["TGCGTAACGTACACTGTGACGAAGACAAAGACTAGGTTACTGACTTCATTCTATAACTTGGTTTTAAAATGTCTCTAAGTACT" 696];third
"TGCGTAACGTACACTCGATCACACGAGAACATCGGTGATTTGGTGTCAATTCTATATACTTTGGGAGTTTTTAAAATGTCTCTAAGTACT" ;cluster17-108
"TGCGTAACGTACACTACCCAAGACGGCTCTACAGTAAGATAGCCTATCATTCTATATGCTTTGGAGTTTTTAAAATGTCTCTAAGTACT" ;cluster17-157
"TGCGTAACGTACACTGGGCAGATCGCACACACGTCTTGCTCGGTGTCATTCTATATACTTTGGGAGTTTTAAAATGTCTCTAAGTACT" ;cluster7-2147
"TGCGTAACGTACACTCGATCACACGAGAACATCGGTGATTTGGTGTCAATCCTATATACTTTGGGAGTTTTAAAATGTCTCTAAGTACT" ;cluster7-255
"TGCGTAACGTACACTCGGCAAATCCACTAACGGACTACTGGGTGATCATTCAATATACTTTTGGAGTTTTAAAATGTCTCTAAGTACT" ;cluster7-637
"TGCGTAACGTACACTAGGCAAACCGATCCTAACGAATGCTTGGTGTCATTCTATATACCTTGGAGTGTTTTAAAATGTCTCTAAGTACT" ;cluster7-73



