(ns hts-exploration.background-set
  (:require
   [clojure.contrib.string :as str]
            [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :as math]
            [smith-waterman :as sw]
            [edu.bc.bio.seq-utils :refer (markov-step generate-rand-seq)])
  (:use [clojure.contrib.core :only [dissoc-in]]
        hts-exploration.globals
        edu.bc.utils
        edu.bc.utils.probs-stats
        hts-exploration.utils))

;;;--------derived from counting ----------------------------------
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
;;;------------------------------------------------------

(defn starting-gaps [s]
  (if-let [gaps (re-find #"^\-+" s)] (count gaps) 0))

(defn- -align
  [alignfn inseqs]
  (->> inseqs
       (r/map alignfn)
       (r/fold concat
               (fn ([] [])
                 ([V [_ s1 s2]]
                    (let [d (starting-gaps s1)]
                      (conj V [(str/drop d s1) (str/drop d s2)])))))))

(defn count-subs
  "Count total substitutions between seq1 and seq2. Does not count the
  number of indels."
  
  [s1 s2]
  (->> (map (fn [x1 x2]
              (if (or (= x1 \-)
                      (= x2 \-))
                0
                (if (= x1 x2) 0 1)))
            s1 s2)
       (reduce +)))

(defn count-indel [s1 s2]
  (let [helper (fn [x] (map count (re-seq #"\-+" x)))
        ins-len (helper s1)
        del-len (helper s2)]
    [(count ins-len) ;insertion
     ins-len
     (count del-len) ;deletion
     del-len]))

(defn align-to-const [inseqs]
  (let [const prime3-const
        get-str (fn [s]
                  (let [s (->> s (str-take-last 52) (str/butlast 15))
                        rtstr (str-take-last 15 s)
                        re #"AATGTC"]
                    (if-let [start-codon (re-find re rtstr)]
                      (str/butlast  (- 15 (.indexOf rtstr start-codon)) s)
                      s)))
        to-const (fn [x] (first (sw/sw const x :gap-open-weight -2
                                      :gap-ext-weight -1 :anchor-right true)))]
    (->> inseqs
         (r/map get-str)
         (-align to-const))))

(defn align-to-primer [inseqs]
  (let [primers (str prime5-const cds)
        to-primer (fn [x] (first (sw/sw primers x :gap-open-weight -2)))]
    (->> inseqs
         (r/map #(str (str/take 15 %) (str-take-last 15 %)))
         (-align to-primer))))

(defn count-mutations
  "Counts the number of subs,ins, and dels when comparing each set of
  aligned sequences. Returns a collection of lists containing
  substition rate, insertion rate, insertion lengths, deletion rate,
  deletion lengths for each aligned seq pair."
  
  [aligned-seqs]
  (map cons (map #(apply count-subs %) aligned-seqs)
       (map #(apply count-indel %) aligned-seqs)))

(defn get-mutation-rates
  "Combines the mutation rates for each individual pair together to
  get a mean mutation rate. Returns a vector of mutation rates and
  indel lengths as [[sub insertion deletion] [insertion deletion]]."

  [counts]
  (let [[subs ins ins-lens dels del-lens] (transpose counts)
        ins-lens (probs 1 (flatten ins-lens))
        del-lens (probs 1 (flatten del-lens))
        overall (->> (transpose [subs ins dels]) (map #(apply + %)) mean)
        totals (map #(apply + %) [subs ins dels])
        ]
    [(mapv mean [subs ins dels])
     [ins-lens del-lens]
     overall
     (->> (map #(double (/ % (apply + totals))) totals)
          (interleave [:sub :ins :del])
          (apply hash-map))]))

(defn mutation-matrix [aligned-seqs]
  (let [mutations (fn [s1 s2]
                    (->> (map (fn [x1 x2]
                                (when (and (not= x1 x2)
                                           (not= x2 \-))
                                  [x1 x2])) 
                              s1 s2)
                         (remove nil?)))
        init-transitions (as-> (repeat {\A 0.1 \C 0.1 \G 0.1 \T 0.1}) x
                               (interleave [\A \C \G \T \-] x)
                               (apply assoc {} x)
                               (reduce #(dissoc-in %1 %2);remove self transitions
                                       x (map #(repeat 2 %) [\A \C \G \T])))]
    (->> (mapcat #(apply mutations %) aligned-seqs)
         (reduce (fn [M [wt mut]] 
                   (assoc-in M [wt mut] (inc (get-in M [wt mut]))))
                 init-transitions)
         (reduce (fn [M [k m]] (assoc M k (probs m))) {}))))

(defn mut-locations
  "Generates locations for any mutations based on the assumption that
  mutations follow a Poisson distribution. Here the positions are
  determined by using the exponential distribution to predict the
  number of nucleotides until the next mutation. Returns a collection
  of positions which will not exceed the length."

  [rate len]
  (->> (repeatedly #(rand-exponential rate))
       (map math/floor)
       (reductions +)
       (take-while #(< % len) )))

(defn rand-mutations [s rate] (-> (mut-locations rate (count s)) reverse vec))

(defmulti ^{:doc "Makes mutations of a type to the a sequence"
            :arglists '([s loc mutation-rates mutant-type])}

  add-mutation (fn [s loc mutation-rates mutant-type] mutant-type))

(defmethod add-mutation :sub [s loc Mrates mutant-type]
  (let [{tmatrix :tmatrix} Mrates
        replacement (->> (str/get s loc) tmatrix markov-step)]
    (str-replace-at loc replacement s)))

(defmethod add-mutation :ins [s loc Mrates mutant-type]
  (let [{:keys [lrates tmatrix]} Mrates
        ins-len (-> lrates first markov-step)]
    (as-> ins-len n
          (repeatedly n #(markov-step (tmatrix \-)))
          (apply str n)
          (str/lower-case n)
          (str-insert-at loc n s))))

(defmethod add-mutation :del [s loc Mrates mutant-type]
  (-> (Mrates :lrates) second markov-step ;del length
      (str-remove-at loc s)))

(defn calc-mutation-rates
  "Calculates mutation rates for either the primers or the constant
  region in round11 of the SELEX data. Returns a nested-map with the
  mutation rates and indel lengths."
  
  [st alignfun] 
  (let [stmt "SELECT sr.sequence, sr.usable_start, sr.usable_stop 
                                              FROM selex_reads as sr 
                                             WHERE sr.usable=1 
                                               AND sr.strand=1
                                               AND sr.round_number=11;"
        aligned-seqs (-> stmt get-db-seq alignfun)
        _ (prn :count (count aligned-seqs))
        [mrates lrates overall totals] (->> aligned-seqs count-mutations get-mutation-rates )
        mrates (mapv #(/ % (count st)) mrates)
        overall (/ overall (count st))
        t-matrix (mutation-matrix aligned-seqs)]
    {:overall overall
     :totals totals
     :mrates mrates
     :lrates lrates
     :tmatrix t-matrix}))

(defn simulate [start-seq mutation-rates]
  (let [[srate irate drate] (mutation-rates :mrates)
        mutate (fn [type s loc] 
                 (add-mutation s loc mutation-rates type))
        add-del (partial mutate :del)
        add-subs (partial mutate :sub)
        add-ins (partial mutate :ins)]
    (as-> start-seq st
          (reduce add-del st (rand-mutations st drate))
          (reduce add-subs st (rand-mutations st srate))
          (reduce add-ins st (rand-mutations st irate)))))

(defn simulate2
  "Takes a starting sequence and a map of mutation rates. Applies the
  mutation rates to the starting seqence and returns a new sequence
  which has been mutated."

  [start-seq mutation-rates]
  (let [mrate (mutation-rates :overall)
        mtype (mutation-rates :totals)
        mutate (fn [[s last-pos] [loc type]] 
                 (let [loc (if (>= loc (count s))
                             (dec (count s)) loc)
                       ;(prn :s s :loc loc :type type)
                       newstr (add-mutation s loc mutation-rates type)]
                   [newstr loc]))
        ;;get a list of mutation location and type
        locs (->> (rand-mutations start-seq mrate)
                  (map (fn [loc] [loc (markov-step mtype)]))
                  #_(group-by second))]
    (first (reduce #(mutate %1 %2) [start-seq -1] locs))))

(defn simulate-seq
  "Simulates a SELEX sequence based on the observed mutation rates
  seen in the HTS sequence data. The variable region is completely
  random."
  
  []
  (let [primer-mrate mutation-rates-primer-region
        const-mrate mutation-rates-const-region
        bprobs (apply hash-map (interleave [\A \C \G \T] (repeat 0.25)))]
    (str/upper-case
     (str (simulate2 prime5-const primer-mrate)
          (generate-rand-seq 30 bprobs)
          (simulate2 prime3-const const-mrate)
          (simulate2 cds primer-mrate)))))
