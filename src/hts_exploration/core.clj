(ns hts-exploration.core
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.io :as io]
            [clojure.contrib.string :as str]
            [clojure.core.reducers :as r]
            
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
        hts-exploration.globals
        hts-exploration.utils))

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



(time
 (let [st prime3-const
       get-subs (fn [c] (markov-step (get t-matrix c)))
       {:keys [mrates lrates tmatrix]} mutation-rates-const-region]
   (pp/pprint 
    (doall
     (for [i (range 10)
           :let [get-locs (fn [s rate] (-> (mut-locations rate (count s)) reverse vec))
                 [srate irate drate] mrates]]
       (let [add-subs (fn [s loc] (str-replace-at loc (get-subs (str/get s loc)) s))
             add-ins (fn [s loc] 
                       (let [size (-> (first lrates) markov-step)]
                         (prn :ins loc)
                         (prn :before @del-locs) 
                         (as-> size n
                               (repeatedly n #(get-subs \-))
                               (apply str n)
                               (str/lower-case n)
                               (str-insert-at loc n s))))
             add-del (fn [s loc] (prn :del loc) 
                       (-> (second lrates) markov-step
                           (str-remove-at loc s)))]
         (as-> prime3-const st
               (reduce add-del st (get-locs st drate))
               (reduce add-subs st (get-locs st srate))
               (reduce add-ins st (get-locs st irate)))))))))

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
    (time (let [st prime3-const
                stmt "SELECT sr.sequence, sr.usable_start, sr.usable_stop 
                                              FROM selex_reads as sr 
                                             WHERE sr.usable=1 
                                               AND sr.strand=1;"
                aligned-seqs (-> stmt get-db-seq align-to-const)
                _ (prn :count (count aligned-seqs))
                [mrates lrates overall] (->> aligned-seqs count-mutations get-mutation-rates )
                mrates (mapv #(/ % (count st)) mrates)
                t-matrix (mutation-matrix aligned-seqs)]
            {:overall overall
             :mrates mrates
             :lrates lrates
             :tmatrix t-matrix}))))

(def mutation-rates-primer-region
  {:mrates [1.984575348257437E-4 8.07315775278278E-5 8.127829249979953E-5] ;sub,ins,del
   :lrates [{3 7.524454477050414E-4, 2 0.012791572610985704, 1 0.9864559819413092};ins lengths
            {4 7.473841554559044E-4, 3 7.473841554559044E-4, 2 0.014200298953662182, 1 0.984304932735426}],;del lens
   :tmatrix {\A {\T 0.27785434114279656, \G 0.5916463479274886, \C 0.1304993109297148},
             \C {\T 0.3725585345141332, \G 0.39196894334586924, \A 0.23547252213999756},
             \G {\T 0.3318593853636966, \C 0.05328321910236568, \A 0.6148573955339376},
             \T {\G 0.18420299532576553, \C 0.4722884670418774, \A 0.34350853763235717},
             \- {\T 0.07349451201423911, \G 0.1832542272322753, \C 0.1261495105309997, \A 0.617101750222486}}})

(def mutation-rates-const-region
  {:mrates [0.03223329251605598 0.009255832498587937 0.055500963713136585]
   :lrates [{1 0.8612541764542815, 2 0.09524504296823706, 3 0.02335896350982652,
             4 0.010891609156830418, 5 0.004501086972380689, 6 0.0024073884941420215,
             7 0.0013277112301025692, 8 6.127897985088781E-4, 9 2.6262419936094777E-4,
             10 1.0213163308481302E-4, 11 2.1885349946745648E-5, 12 1.4590233297830432E-5}
            {9 1.2165982938425528E-6, 8 5.718011981059998E-5, 7 1.7275695772564248E-4,
             6 3.69845881328136E-4, 5 0.0010474911309984378, 4 0.0031911373247490156,
             3 0.019055579076455904, 2 0.12037023519278217, 1 0.8557345577178562}],
   :tmatrix {\A {\T 0.33324249405940765, \G 0.44204819655407135, \C 0.2247093093865209},
             \C {\T 0.3450192114665972, \G 0.40237560831433644, \A 0.2526051802190663},
             \G {\T 0.46401157839942464, \C 0.17467195036453304, \A 0.36131647123604227},
             \T {\G 0.22885529437082425, \C 0.4907117501616078, \A 0.28043295546756797},
             \- {\T 0.25210274423815193, \G 0.16376058545660183, \C 0.1908705929051287, \A 0.3932660774001174}}})
