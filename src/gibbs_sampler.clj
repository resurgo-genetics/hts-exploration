(ns gibbs-sampler
  (:require [clojure.core.reducers :as r]
            [clojure.contrib.io :as io]))

(def counter (atom 0))
(def inspectable-state (atom {}))

(def wordlength 6)
(def pseudocount 0.5)
(def percent-converge 0.0001)

(def ^{:dynamic true} *sequence-vec* nil)
(def state {:baseprob {} :samples {}})

(defn- transpose [coll] (apply map vector coll))

(defn- sample
  "samples a kmer from each string in a collection. Returns a map of
  kmer strings where the key is the sequence number the kmer is
  sampled from"
  
  [coll]
  (->> coll
       (map-indexed (fn [k s]
                      (let [maxi (- (count s) wordlength)
                            i (rand-int maxi)
                            j (+ i wordlength)]
                        [k (subs s i j)])))
       (into {})))

(defn- remove-seq
  "Removes a random sequence from a map where k=sequence number and
  v=sampled kmer. Returns a vector of which k sequence was removed and
  m without the removed sequence"

  [m]
  (let [k-to-remove (rand-int (count m))]
    [k-to-remove (dissoc m k-to-remove)]))

(defn rfrequencies [coll]
  (r/fold (fn ([] {})
            ([l r] (merge-with + l r)))
          (fn ([] {})
            ([M v]
               (assoc M v (inc (get M v 0)))))
          coll))

(defn- pseudo-prob
  "Determines the probability of a single position of the kmer and
  adds in a pseudocount so that all bases are represented"

  [coll]
  (let [len (+ (* 4 pseudocount)
               (count coll)
               )]
    (->> coll
         frequencies
         (merge-with + {\A 0.5 \C 0.5 \G 0.5 \T 0.5})
         (map (fn [[k v]] [k (/ v len)]))
         (into {}))))

(defn- probabilities
  "Calculates the base probabilies for each position of the
  kmer. Returns a vector of maps where each map is the base
  probabilies for that position."

  [coll]
  (reduce (fn [V col]
            (conj V (pseudo-prob col)))
          [] (transpose coll)))

(defn score
  "Determines the highest scoring kmer in a sequence based on the
  probability table"

  [s prob-table]
  (let [kmers (partition wordlength 1 s)] ;all possible kmers in sequence
    (reduce (fn [[top-score top-seq :as top] cur-kmer]
              (let [cur-score (reduce * (map #(get %1 %2) prob-table cur-kmer))]
                (if (> cur-score top-score)
                  [cur-score (apply str cur-kmer)]
                  top)))
            [-1 ""] kmers)))

(defn- converge?
  "Checks for convergence of the gibbs sampling when the delta base
  probabilies is less than a threshold. By default the threshold is
  0.001"

  [[pstate cstate] & {:keys [thr]
                      :or {thr 0.001}}]
  ;;finds delta across every position in the kmer
  (->> (map (fn [ps cs] 
              (merge-with (fn [l r] (Math/abs (- l r))) ps cs));delta
            (pstate :baseprob) (cstate :baseprob))
       (apply merge-with max);max delta for each base is stored
       vals
       (every? #(< % thr))));all base probs must have delta < thr

(defn- single-round
  "Single round of the gibbs sampler. Returns a vector of the previous
  state and the new current state."

  [states]
  (let [prev-state (second states)
        sampled-kmer (prev-state :samples)
        [kseq temp-samples] (remove-seq sampled-kmer)
        baseprob (probabilities (vals temp-samples))
        [score seq-to-add] (score (*sequence-vec* kseq) baseprob)]
    (swap! counter inc)
    (when (zero? (mod @counter 1000))
      (swap! inspectable-state assoc :baseprob baseprob))
    [prev-state
     (assoc state ;create new state
       :baseprob baseprob 
       :samples (assoc temp-samples kseq seq-to-add))]))

(defn -main [inseqs & {:keys [limit]
                       :or {limit 100000}}]
  (reset! counter 0)
  (reset! inspectable-state {})
  (binding [*sequence-vec* inseqs]
    (->> [state (assoc state :samples (sample *sequence-vec*))] ;initial state
         (iterate single-round)    ;repeatedly executing the sampling
         (drop 10)                 ;Burn in period
         (take-while #(complement (converge? % :thr 0.0001))) ;while not converged
         (take limit)            ;hard limit on number of iterations
         last)))

;;;test cases
(comment
 (def test-case
   (mapv #(.toUpperCase %)
         ["gatgacgtggtttacgaccccaTTTAGTagtcaaccgcagtgagtgagtc"
          "ttgaaaccagacgtttcgccccTATTACagactcacaaccacatgatgac"
          "ctggcggcgtagcgatgcgctgGTTACTctgaaaacggtctatgcaaatt"]))

 (def example-file
   (->> (io/read-lines "/home/kitia/Downloads/gibbs-test.txt")
        (remove empty?)
        rest
        (mapv #(.toUpperCase %))))

 (-main test-case)
 (-main example-file))
