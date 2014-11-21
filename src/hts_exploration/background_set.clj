(ns hts-exploration.background-set
  (:require
   [clojure.contrib.string :as str]
            [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :as math]
            [smith-waterman :as sw])
  (:use [clojure.contrib.core :only [dissoc-in]]
        hts-exploration.globals
        edu.bc.utils
        edu.bc.utils.probs-stats
        [edu.bc.bio.seq-utils :only (markov-step)]
        hts-exploration.utils))

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
