(ns analysis.context
  (:require [clojure.java.jdbc :as jdbc]
            [clojure.core.reducers :as r]
            [clojure.core.match :refer [match]]
            [clojure.java.io :as io2]
            [clojure.data.csv :as csv]
            [semantic-csv.core :as sc]
            [clojure.contrib.string :as str]
            [edu.bc.fs :as fs]
            [clojure.java.shell :as shell]
            [clojure.pprint :as pp]
            [edu.bc.bio.sequtils.files :refer (read-seqs)]
            [hts-exploration.hts-utils.file :refer (write-fasta)]
            [edu.bc.utils :refer (sum transpose log pxmap reversev count-if)]
            [edu.bc.bio.seq-utils :refer (markov-step markov-step-order-n
                                                      generate-rand-seq dint-seq-shuffle)]
            [instaparse.core :as insta]
            [clojure.zip :as zip]
            [edu.bc.utils.fold-ops :refer (fold bpdist)]
            [hts-exploration.background-set :as bg]
            [refold :as refold])
  (:use edu.bc.utils.probs-stats
        hts-exploration.db-queries
        hts-exploration.globals
        hts-exploration.utils)
  (:import [java.util Random]))

;;;---------structure stuff------------------------------------------------------------
(defn fold-subopt-sequence [s] (fold s {:foldmethod :RNAsubopt :n 1000}))

(def parse-struct
  (insta/parser "S = S S | U | H
                 H = '(' (M|U|H) ')' | epsilon
                 M =  H H
                 U = '.' (S|U|H|epsilon)"))

(def parse-struct2
  (insta/parser "<S> =  S S | unpaired | stem | epsilon
                 <H> = stem | loop | M | H H
                 unpaired = <'.'>
                 loop = <'.'> 
                 stem = <'('> H <')'>
                 M = loop* stem loop* (stem loop*)+"))

(defn struct->context
  "Takes a structure and turns it into a context. The alphabet used
  here is u=unpaired, l=loop, p=paired. It is hard coded in but could
  be modified later to include more."
  
  [st]
  (let [stack (atom [])]
    (map (fn [c]
           (cond (and (not (empty? @stack)) (= c \.)) \l
                 (= c \() (do (swap! stack conj c) \p)
                 (= c \)) (do (swap! stack pop) \p)
                 :else \u))
         (seq st))))

(defn calc-struct-context
  "Take a set of structure contexts, presumably the subopt structures
  and returns probabilities for each structure context at each
  position. Returns a map."
  
  [incontexts & {:keys [alphabet]
                 :or {\p 1 \u 1 \l 1}}]
  (->> (transpose incontexts)
       (map (partial freqn 1) )
       (map #(merge-with + % alphabet))
       (map probs)
       (map-indexed vector)))

(defn check-multiloop
  "Looks for a multiloop in a parsed structure"

  [loc]
  (if (zip/end? loc)
    :end
    (if (= (zip/node loc) :M)
      :M
      (recur (zip/next loc)))))

(defn get-helicies
  "Returns helicies. ignores the closing base pairs of the multiloops"
  
  [loc]
  (let [get-remaining (comp vec                               ; fn to get subtree and remaining helicies
                         (partial apply concat)
                         (juxt #(vector (zip/node %)) zip/rights))
        floc (fn [mloc] (->> mloc                                ; multiloop location
                            (iterate zip/next)
                            (drop-while #(not= (zip/node %) :M)); goto multiloop
                            first zip/next                      ; subtree in multiloop
                            get-remaining))]
    (->> loc zip/next get-remaining
         (mapcat #(if (= (check-multiloop (zip/vector-zip %)) :M)
                    (floc (zip/vector-zip %))
                    [%]))
         (filterv #(= (first %) :stem)))))

(defn check-3helix
  "Looks for 3 separate hairpins in a parsed structure"
  
  [loc]
  (count (get-helicies loc)))

(defn bp-counter
  "Takes a seq of nodes from the tree and counts up the continous
  stacks to get the max helix size. Doesn't include bulges on both
  sides of the helix, just the 5' side of it. Should also integrate
  with the continuous function but this is going to have to do for
  now."
  [helix]
  (->> (flatten helix)
       (partition-by #{:loop})
       (reduce (fn [[next-fn cur-max cnt] x]
                 (let [c (count x)
                       next-cnt (next-fn cnt c)]
                   (cond (and (< c 4) (= :loop (first x))) [+ cur-max cnt]; small loop continue counting
                         (= :stem (first x)) [max (max cur-max next-cnt) next-cnt]; cur-max or new stem
                         :else [max cur-max 0]))); reset stem length counter
               [max 0 0])
       second))

(defn bp-counter2
  "only counts base pairs in a stem. A stem is defined as a stack of
  base pairs with bulges/internal loops <= 3"
  [tree]
  (letfn [(hwalk [cur-max cnt loc]
            (if (zip/end? loc)
              cur-max
              (if (and (zip/branch? loc)
                       (= :stem (-> loc zip/down zip/node)))
                (let [child-seq (zip/children loc)
                      next-cnt (inc cnt)
                      [lbulge rbulge :as loop-cnt]
                      (->> ((juxt rest reverse) child-seq) ; check if we are in a bulge
                           (map #(take-while #{[:loop]} %)); get bulges
                           (map count))]
                  (if (and (some #{:stem} (flatten (rest child-seq))); more bp in helix?
                           (every? #(<= % 3) loop-cnt))              ; bulge <= 3
                    (recur cur-max next-cnt
                           (->> child-seq rest    ; continue progressing on helix
                                (drop lbulge)     ; remove any 5' bulge
                                (drop-last rbulge); remove any 3' bulge
                                vec zip/vector-zip zip/next))
                    (recur (max cur-max next-cnt) 0 (zip/next loc))))
                (recur cur-max cnt (zip/next loc)))))]
    (->> tree get-helicies (map #(hwalk 0 0 (zip/vector-zip %))))))

(defn get-gugc-pos
  "Checks a structure for the gugc base pair and reports locations where GG is paired to (CT)|(TC)"
  ([s] (get-gugc-pos s (-> s fold ffirst)))
  ([s st]
   (let [pairs (refold/make-pair-table st)]
     (->> (str-re-pos #"(?=(GG))" s)
          (reduce (fn [[hit cnt :as x] y]
                    (let [i (first y)
                          j (get pairs i 0)
                          k (inc i)
                          l (get pairs k 0) ]
                      (if (and (not (or (zero? l) (zero? j)))
                               (= (dec j) l))
                        (if (or (= [\C \T] [(str/get s l) (str/get s j)])
                                (= [\T \C] [(str/get s l) (str/get s j)]))
                          [(conj hit [i j k l]) (inc cnt)]
                          x)
                        x)))
                  [[] 0])))))

(defn get-gcgc-pos
  "Checks a structure for the gugc base pair and reports locations where GG is paired to CC"
  ([s] (get-gcgc-pos s (-> s fold ffirst)))
  ([s st]
   (let [pairs (refold/make-pair-table st)]
     (->> (str-re-pos #"(?=(GG))" s)
          (reduce (fn [V y]
                    (let [i (first y)
                          j (get pairs i 0)
                          k (inc i)
                          l (get pairs k 0) ]
                      (if (and (not (or (zero? l) (zero? j)))
                               (= (dec j) l))
                        (conj V [i j k l])
                        V)))
                  [])
          (reduce (fn [[hit cnt :as x] [i j k l :as y]]
                    (if (= [\C \C] [(str/get s l) (str/get s j)])
                      [(conj hit y) (inc cnt)]
                      x))
                  [[] 0])))))

(defn get-guac-pos
  "Checks a structure for the gugc base pair and reports locations
  where GA and AG is paired to UC and CU, respectively."
  ([s] (get-guac-pos s (-> s fold ffirst)))
  ([s st]
   (let [pairs (refold/make-pair-table st)]
     (->> ((juxt #(str-re-pos #"GA" %) #(str-re-pos #"AG" %)) s);check ga, and ag
          (apply merge )
          (reduce (fn [[hit cnt :as x] y]
                    (let [i (first y)
                          j (get pairs i 0)
                          k (inc i)
                          l (get pairs k 0) ]
                      (if (and (not (or (zero? l) (zero? j)))
                               (= (dec j) l))
                        (if (= (if (= (second y) "GA")
                                 [\T \C]
                                 [\C \T])
                               [(str/get s l) (str/get s j)])
                          [(conj hit [i j k l]) (inc cnt)]
                          x)
                        x)))
                  [[] 0])))))

(defn valid-bp-stacks
  "Returns base pair stacks of size n"
  [n]
  (let [valid-pairs (mapcat (juxt identity #(reversev %)) [[\G \C] [\A \T] [\G \T]])
        helper (fn this [n i cur]
                 (if (< i n)
                   (let [new (for [c cur x valid-pairs] (concatv c x))]
                     (this n (inc i) new))
                   cur))]
    (mapv #(mapv vec (partition 2 %)) (helper (dec n) 0 valid-pairs))))

(defn valid-bp-stacks-str [n sep]
  (->> (valid-bp-stacks n)
       (mapv #(->> % (map (partial apply str)) (str/join sep)))))

(defn continuous?
  "Checks a collection of numbers to see if they are continuous with a
  max gap of x (<= -(+ x 1) % (+ x 1)). x is the number of gaps to allow"
  [x coll] ; continuous seq of base pairs?
  (let [x (inc x)]
    (->> (sort coll)
         (partition 2 1)
         (map #(apply - %))
         (every? #(<= (- x) % x)); should be (<= -4 x 4) for bulges of 3
         )))

(defn potential-stacks
  "Checks the stack location to identify the helix length and where in
  the helix the stack is found."
  [max-gap m stack-locs]
  (let [i (apply min stack-locs); stack start position
        j (apply max stack-locs); stack end position
        valid-stack? (fn [l r]
                       (let [m (->> m
                                    ;;keep bp that are nested in [i j] or that [i j] is nested within
                                    (filter (fn [[k l]] (or (<= i k l j) (<= k i j l))))
                                    (into (sorted-map)))
                             m2 (reduce (fn [x [k l]] ; max stem size containing [i j]
                                          (let [check-continuous? (fn [m]
                                                                    (and (continuous? max-gap (keys m))
                                                                         (continuous? max-gap (vals m))))]
                                            (if (check-continuous? (assoc x k l))
                                              (assoc x k l)
                                              x)))
                                        (select-keys m [i])
                                        (sort-by #(Math/abs (- (first %) i)) m))];abs dist from i
                         [(count (filter #(< % i) (keys m2)))
                          (count (filter #(> % i) (keys m2)))
                          (count m2)]))
        stacks (valid-stack? i j)]
    stacks))

(defn get-motif-pos
  "Checks a structure for the base pairs (stack1(i,j) stack2(k,l)
  where (< i k l j)) and reports locations where stacks are adjacent
  and proper pairs. "
  ([max-gap s] (get-motif-pos 2 max-gap s (-> s fold ffirst)))
  ([max-gap s st] (get-motif-pos 2 max-gap s st))
  ([stack-size max-gap s st]
   (let [continuous-bp-pos? (fn [& ps]
                              (let [ps (transpose (mapv second ps))]
                                (every? (partial continuous? 0) ps)))
         pairs (->> (refold/make-pair-table st)
                    (remove (fn [[i j]] (> i j)))
                    (into (sorted-map) ))
         stack-locs-vec (->> pairs
                             (map (fn [x] [(mapv #(str/get s %) x) x]))
                             (partition stack-size 1)
                             (filter #(apply continuous-bp-pos? %)))]
     (reduce (fn [M x]
               (let [locs-vec (apply concat (map second x))
                     stack-str (apply str (apply concat (map first x)))
                     p-stack (potential-stacks max-gap pairs locs-vec) ; flanked by bp-stack
                     existing-val (get M stack-str [])]
                 (assoc M stack-str (conj existing-val p-stack))))
             {} stack-locs-vec))))

(defn motif-stack?
  "takes a the 2 bp stack (motif-bp1, motif-bp2), seq s and returns if
  the seq has the base stack in question."
  ([motif-bp1 motif-bp2 s]
   (motif-stack? motif-bp1 motif-bp2 s (-> s fold ffirst)))
  ([motif-bp1 motif-bp2 s st]
   (->> (get-motif-pos motif-bp1 motif-bp2 s st)
        ((complement empty?)))))

(defn calc-motif-probs
  "Takes a vector of motifs (valid-bp-stacks) in the format [bp1 bp2]
  where bp1 and bp2 are the first and second base-pairs as a vector of
  chars (ie [\\G \\C]). The inseqs are just sequences to calculate if
  they have the motifs. The return is a map of the motifs and how
  prevalent they are in the inseqs"

  [motifs-vec inseqs]
  (let [instructs (map first (fold inseqs))
        data (mapv (fn [[a b]]
                     [(apply str (concat a b))
                      (double
                       (/ (->> (vfold (partial motif-stack? a b) 30 inseqs instructs)
                               (count-if true?))
                          (count inseqs)))])
                   motifs-vec)]
    (into (sorted-map) data)))

(defn multi-or-3helix
  "checks for :M or 3 stems, really needs to be rewritten. Entry is a
  parsed tree containing the propper flags already"
  [entry]
  (or (>= (last entry) 3)
      (and (= (last entry) 1)
           (= (second entry) :M))))

(defn st-context-tree
  "takes a set of structures, presumably subopt structures of a
  sequences s, and checks for multiloops or 3 helix. returns a vector
  of structure, multiloop, 3helix."

  [struct-set]
  (r/fold 50
          (fn ([] []) ([l r] (concat l r)))
          (fn ([] [])
            ([V st]
             (let [tree (->> (parse-struct2 st)
                             vec zip/vector-zip)]
               (->> ((juxt check-multiloop check-3helix) tree)
                    (apply vector st)
                    (conj V)))))
          (vec struct-set)))

(defn count-multiloops [st-context]
  (->> (filter #(= (second %) :M) st-context)
       count))

(defn count-3helix [st-context]
  (->> (filter #(= (last %) 3) st-context)
       count))

(defn calc-struct-features
  "Takes seqs and corresponding structs to get the mfe,
  multiloop/helix counts, motif positions and returns a vector."
  
  [s [st nrg]]
  [nrg
   (if (re-find #"[^\.]" st)
     (let [loc (->> (parse-struct2 st)
                    vec zip/vector-zip )
           helix-lengths (bp-counter2 loc)]
       (vector (if (= :M (check-multiloop loc)) 1 0)
               (check-3helix loc)
               (apply max helix-lengths)
               (mean helix-lengths)
               (apply min helix-lengths)
               (str/join "," helix-lengths)))
     [0 0 0 0])
   (mapv (fn [[a b]]
           (let [x (get-motif-pos a b s st)
                 y (->> (remove empty? x)
                        (map first))
                 z (map (fn [foo]
                          (match [foo]
                                 [[_ :5prime]] :1
                                 [[_ :3prime]] :2))
                        y)]
             [(count (get-motif-pos a b s st))
              (count-if #(= % :1) z); stem before
              (count-if #(= % :2) z)]));stem after
         valid-bp-stacks)])

(defn calc-struct-features2
  "Takes seqs and corresponding structs to get the mfe,
  multiloop/helix counts, motif positions and returns a vector."
  
  [stack-size max-gap s [st nrg]]
  [nrg
   (if (re-find #"[^\.]" st)
     (let [loc (->> (parse-struct2 st)
                    vec zip/vector-zip )
           helix-lengths (bp-counter2 loc)]
       (vector (if (= :M (check-multiloop loc)) 1 0)
               (check-3helix loc)
               (apply max helix-lengths)
               (mean helix-lengths)
               (apply min helix-lengths)
               (str/join "," helix-lengths)))
     [0 0 0 0])
   (let [motif-pos (get-motif-pos stack-size max-gap s st)]
     (mapv (fn [stack-str]
             (let [x (get motif-pos stack-str [["." "." "."]])
                   all-abs-helix-starts (str/join "," (map #(nth % 0) x))
                   all-abs-helix-ends (str/join "," (map #(nth % 1) x))
                   all-abs-helix-lengths (str/join "," (map #(nth % 2) x))]
               [(if (= x [["." "." "."]]) 0 (count x))
                all-abs-helix-starts
                all-abs-helix-ends
                all-abs-helix-lengths]))
           (valid-bp-stacks-str stack-size "")))
   (let [dinucleotides (add-letters 2 [\A \C \G \T])
         m (probs 2 s)]
     (mapv #(get m % 0) dinucleotides))])

(defn write-header [wrtr stack-size]
  (.write wrtr "sid,ssid,cid,csize,round,hitseq,mfe,multiloop,helix,max helix size,mean helix size,min helix size,helix lengths,")
  (.write wrtr (str/join ","
                         (for [vbs (valid-bp-stacks-str stack-size "")
                               h ["" ".all.abs.helix.start"
                                  ".all.abs.helix.end"
                                  ".all.abs.helix.length"]]
                           (str vbs h))))
  (.write wrtr ",")
  (.write wrtr (str/join "," (add-letters 2 [\A \C \G \T])))
  (.write wrtr "\n"))

(defn explode-calc-features
  "If there are multple NCMs/stack motifs that appear in the same seq,
  then it will create multiple entries for the seq, each with one
  start;stop;length entry."
  
  [stack-size infile outfile]
  (with-open [rdr (io2/reader infile)]
    (let [const-header ["ssid" "round" "mfe" "multiloop" "helix"]
          data (sc/mappify {:keyify false} (csv/read-csv rdr))
          out (fn [wrtr line-map]
                (doseq [sk (valid-bp-stacks-str stack-size "")]; sk = ncm stack key
                  (let [const (cons sk (map line-map const-header))
                        explode (->> [".all.abs.helix.start" ".all.abs.helix.end"
                                      ".all.abs.helix.length"]
                                     (mapv (partial str sk)); make ncm specific keys
                                     (map (fn [sstk] ; function returns vecs of starts, stops, lengths
                                            (str/split #"," (get line-map sstk))))
                                     transpose ; return tuples in form of start;stop;length
                                     (remove #(= ["." "." "."] %))); NAs
                        line (mapv (partial concatv const) explode)]
                    (csv/write-csv wrtr line))))]
      (with-open [wrtr (io2/writer outfile)]
        (.write wrtr (str/join "," ["stack.type" "ssid" "round" "mfe" "multiloop"
                                    "helix" "all.abs.helix.start" "all.abs.helix.end"
                                    "all.abs.helix.length"]))
        (.write wrtr "\n")
        (doseq [m data] (out wrtr m))))))

(defn comparpare-saves [f1 f2]
  (let [headers (->> (valid-bp-stacks-str 2 "")
                     (mapcat (fn [x]
                               (for [h ["" ".all.abs.helix.start"
                                        ".all.abs.helix.end"
                                        ".all.abs.helix.length"]]
                                 (keyword (str x h))))))
        read (fn [f]
               (with-open [rdr (io2/reader f)]
                 (doall (set
                         (map (fn [m]
                                (mapv #(get m %) (cons :sid headers)))
                              (sc/mappify (csv/read-csv rdr)))))))
        [f1 f2] (map read [f1 f2])]
    [(count (clojure.set/difference f2 f1))
     (count (clojure.set/difference f1 f2))]))

;;;-----------------------------------------------------------------------------------

;;;---------test utilities -----------------------------------------------------------

(clojure.test/deftest test-get-helicies
  (clojure.test/is (= (get-helicies (-> "((((.((..)).((...))..))))..(((...)))"
                                        parse-struct2 vec zip/vector-zip))
                      '([:stem [:stem [:loop] [:loop]]] [:stem [:stem [:loop] [:loop] [:loop]]] [:stem [:stem [:stem [:loop] [:loop] [:loop]]]])))
  (clojure.test/is (= (get-helicies (-> "((((.((..)).((...))..))))." parse-struct2 vec zip/vector-zip))
                      [[:stem [:stem [:loop] [:loop]]] [:stem [:stem [:loop] [:loop] [:loop]]]]))
  (clojure.test/is (= (get-helicies (-> "((((.((..))..))))." parse-struct2 vec zip/vector-zip))
                      [[:stem [:stem [:stem [:stem [:loop] [:stem [:stem [:loop] [:loop]]] [:loop] [:loop]]]]]]))
  (clojure.test/is (= (get-helicies (-> "((((..)))).((...))..(((...)))" parse-struct2 vec zip/vector-zip))
                      [[:stem [:stem [:stem [:stem [:loop] [:loop]]]]] [:stem [:stem [:loop] [:loop] [:loop]]] [:stem [:stem [:stem [:loop] [:loop] [:loop]]]]])))

(clojure.test/deftest test-get-motif-pos-2
  (clojure.test/is (= {"ATGC" [[5 1 7]], "CGGC" [[1 2 4] [1 8 10]], "CGTG" [[6 3 10]], "GCTA" [[2 1 4] [2 7 10]], "TGTA" [[0 6 7]], "TATA" [[1 5 7]], "GTCG" [[0 3 4]], "ATCG" [[0 9 10] [5 4 10]], "TAGT" [[2 4 7]], "GTGC" [[3 3 7]], "CGAT" [[4 5 10] [8 1 10]], "GCAT" [[4 2 7]]}
                      (let [s "TGCGTAACGTACACTGACCATACATGCAGGTGAACGTGAACAACGTCATTCTATATTTTGGAGTTTTAAAATGTCTCTAAGTACT"
                            [st nrg] (first (fold s))]
                        (prn s)
                        (prn st)
                        (get-motif-pos 2 3 s st)))))

;;;-----------------------------------------------------------------------------------

;;;---------background seq generation-------------------------------------------------


(defn my-simulate-seq
  "Simulates a selex sequence taking into account the const-region
  mutation rate. primer region mutation rate is ignored because we
  expect perfect primers in our data currently (20150805)."
  
  ([] (my-simulate-seq (zipmap [\A \C \G \T] (repeat 0.25))))
  ([bprobs] (my-simulate-seq 0 bprobs))
  ([order bprobs] (let [;primer-mrate bg/mutation-rates-primer-region
                  const-mrate bg/mutation-rates-const-region]
              (str/upper-case
               (str prime5-const
                    (generate-rand-seq 30 order bprobs)
                    (bg/simulate2 prime3-const const-mrate)
                    cds))))
  ([order bprobs init-bprobs]
   (let [;primer-mrate bg/mutation-rates-primer-region
         const-mrate bg/mutation-rates-const-region]
     (str/upper-case
      (str prime5-const
           (generate-rand-seq 30 order bprobs init-bprobs)
           (bg/simulate2 prime3-const const-mrate)
           cds)))))
;;;------------------------------------------------------------------------------

;;;------Calculations------------------------------------------------------------
(defn rand-normal [mu sigma]
  (let [r (Random.)]
    (+ mu (* sigma (.nextGaussian r )))))

(defn pssm
  "Scoring matrix for bases."

  [inseqs]
  (->> (transpose inseqs)
       (map (partial freqn 1) )
       (map #(merge-with + % {\A 1 \C 1 \G 1 \T 1}))
       (map probs)
       (map-indexed vector)))

(defn calc-kmer-score [k inseq]
  (sum ))

(defn ktuples [k s]
  (let [st-context (->> (fold-subopt-sequence s)
                        (map struct->context)
                        calc-struct-context
                        (map #(->> % second (apply max-key val) key)))]
    (freqn 1 (map vector
                  (partition k 1 s)
                  (partition k 1 st-context)))))

(defn get-seqs-context [cluster-ids]
  (mapcat (fn [c] (->> (get-seqs-in-cluster c) sql-query )) cluster-ids))



M (apply hash-map
         (for [i (map seq (hts-exploration.utils/addletters k ["A" "C" "G" "T"]))
               j (map seq (hts-exploration.utils/addletters k [\l \u \p]))]
           [[i j] 0]))

(def ktuple4
  (let [q (fn [n len]
           [(str "SELECT sr.selex_id as id, cdhit.cluster_id as cid, cluster_size as csize, sr.const_start as cs, substring(sequence, usable_start+1, sr.length) as hitseq FROM selex_reads as sr, cdhit WHERE sr.selex_id=cdhit.selex_id and cluster_size=1 and sr.usable=1 and sr.strand=1 and sr.length=? LIMIT ?") len n])
       q2 (fn [cid]
            [(str "select cdhit.cluster_size as csize,
                          substring(sequence, usable_start+1,sr.length) as hitseq
                     FROM selex_reads as sr, cdhit
                    where sr.selex_id=cdhit.selex_id 
                      AND sr.usable=1 AND sr.strand=1
                      AND cluster_id=?;") cid])
       k 4
       singletons (fn [k len]
                    (->> (q 1000 len) sql-query (mapv (juxt :cs :hitseq))
                         (r/fold 50
                                 (fn ([] {}) ([& m]
                                              (do (prn :l (count l) :r (count r))
                                                  (apply merge-with + m))))
                                 (fn ([] {}) ([m [cs s]]
                                              (let [s (str (str/take 15 s)
                                                           (->> (str/drop 15 s)
                                                                (str/take cs)
                                                                (dint-seq-shuffle 1)
                                                                keys first)
                                                           (str/drop (+ 15 cs) s))])
                                              (merge-with + m (ktuples k s)))))
                         ((fn [x] (prn :afterfold) x))
                         probs))        ;~57k singletons
       doubletons  (->> (q1 100 88 10) sql-query
                        (map :cid)
                        (map (fn [c]
                               [c (->> (q2 c) sql-query (map :hitseq)
                                       frequencies (sort-by second >) ffirst)])))
       ston-map (into {} (map (fn [[cid s]]
                                (let [l (count s)]
                                  [l (singletons k l)]))
                              doubletons))]
   (->> doubletons
        (map (fn [[cid s]]
               (let [L (count s)
                     stons (get ston-map L)]
                 [cid (merge-with (comp edu.bc.utils/log /) (probs (ktuples k s)) stons)])))
        doall)))


(def ston-map
  (let [q (fn [n len]
            [(str "SELECT sr.selex_id as id, cdhit.cluster_id as cid, cluster_size as csize, sr.const_start as cs, substring(sequence, usable_start+1, sr.length) as hitseq FROM selex_reads as sr, cdhit WHERE sr.selex_id=cdhit.selex_id and cluster_size=1 and sr.usable=1 and sr.strand=1 and sr.length=? LIMIT ?") len n])
        q2 (fn [cid]
                               [(str "select cdhit.cluster_size as csize,
                          substring(sequence, usable_start+1,sr.length) as hitseq
                     FROM selex_reads as sr, cdhit
                    where sr.selex_id=cdhit.selex_id 
                      AND sr.usable=1 AND sr.strand=1
                      AND cluster_id=?;") cid])
                          k 4
                          singletons (fn [k len]
                                       (->> (q 1000 len) sql-query (mapv (juxt :cs :hitseq))
                                            (r/fold 100
                                                    (fn ([] {}) ([& m]
                                                                 (do (prn :l (count l) :r (count r))
                                                                     (apply merge-with + m))))
                                                    (fn ([] {}) ([m [cs s]]
                                                                 (let [s (if (<= 28 cs 32)
                                                                           (str (str/take 15 s)
                                                                                (->> (str/drop 15 s)
                                                                                     (str/take cs)
                                                                                     (dint-seq-shuffle 1)
                                                                                     keys first)
                                                                                (str/drop (+ 15 cs) s))
                                                                           s)]
                                                                   (merge-with + m (ktuples k s))))))
                                            ((fn [x] (prn :afterfold) x))
                                            probs))        ;~57k singletons
                          doubletons  (->> (q1 100 88 10) sql-query
                                           (map :cid)
                                           (map (fn [c]
                                                  [c (->> (q2 c) sql-query (map :hitseq)
                                                          frequencies (sort-by second >) ffirst)])))]
                      (->> doubletons
                           (map second);get seqs
                           (map count) ;doubleton lengths
                           distinct
                           (pmap (fn [l] [l (singletons k l)]))
                           (into {}))))

(let [doubletons  (->> (q1 100 88 10) sql-query
                                         (map :cid)
                                         (map (fn [c] [c (ffirst (most-freqn-in-cluster 1 c))])))
                        foo (fn [sts] (->> sts
                                           (map struct->context)
                                           calc-struct-context
                                           (map #(->> % second (apply max-key val) key))))
                        out (map (fn [[cid s]]
                                   (let [struct-set (fold-subopt-sequence s)
                                         st-context (foo struct-set)
                                         other-context (->> (st-context-tree (vec struct-set))
                                                            (filter multi-or-3helix)
                                                            count)]
                                     [cid other-context s (apply str st-context)]))
                                 doubletons)]
                    (with-open [wrtr (io2/writer "data/150325/ktuple-contexts.txt")]
                      (doseq [[c o s st] out]
                        (.write wrtr (str c " " o "\n"))
                        (.write wrtr (str s "\n"))
                        (.write wrtr (str st "\n")))))

(defn build-fastas-multiple-clusters [& cids]
  (let [get-seqs (comp #(map (juxt :id :hitseq) %)
                    sql-query
                    get-seqs-in-cluster )]
    (mapcat (fn [c]
              (->> (get-seqs c)
                   (map (fn [[id s]] [(str id "-" c) s]))))
            cids)))


(comment "this stuff needs to be moved to its own ns or something. This is just a quick piece of code to find the zscore of the kmer counts to normalize the counts. Also, I have been having trouble getting a signal above background likely because we don't know what the background looks like. This time, we are only analyzing 2 clusters (9391, 2) because it will simplify the number of sequences we are looking at. I picked cluster 2 because it is the largest cluster but has low percent idienty"
         
         (defn zscore [x mu sigma] (if (zero? sigma) 0 (/ (- x mu) sigma)))
         (defn mean-kmer [n totals]
           (reduce (fn [M k]
                     (assoc M k (double (/ (totals k) n))))
                   {} (keys totals)))
         (defn sd-kmer [S totals]
           (into {}
                 (map (fn [k]
                        [k (sd (map (fn [m] (get m k 0)) S))])
                      (keys totals))))
         (defn create-set [k coll] (mapv (fn [[cid id s]] [cid id (freqn k s)]) coll))
         (defn attach-label [lab coll] (map #(vec (cons lab %)) coll))
         (let [q (fn [cid] (->> (get-seqs-in-cluster cid) sql-query (map (juxt :id :hitseq)) ))
               mostfreq (attach-label 9391 (q 9391))
               background (attach-label 2 (q 2))
               S (create-set 2 (concat mostfreq background))
               totals (->> (map last S) (apply merge-with +))
               mu (mean-kmer (count S) totals)
               sdev (sd-kmer (map last S) totals)]
           (take 3
                 (map (fn [[cid id m]]
                        [cid id
                         (into (sorted-map)
                               (mapv (fn [[k v]]
                                       [k (zscore v (mu k) (sdev k))])
                                     m))])
                      S)))

         (let [q (fn [cid] (->> (get-seqs-in-cluster cid) sql-query (map (juxt :id :hitseq)) ))
               mostfreq (attach-label 9391 (q 9391))
               background (attach-label 2 (q 2))
               S (create-set 2 (concat mostfreq background))
               totals (->> (map last S) (apply merge-with +))
               mu (mean-kmer (count S) totals)
               sdev (sd-kmer (map last S) totals)
               data (r/fold (fn ([] [])
                              ([& vs] (apply concat vs)))
                            (fn ([] [])
                              ([V [cid id m]]
                               (conj V
                                     [cid id
                                      (into (sorted-map)
                                            (map (fn [[k v]]
                                                   [k (zscore v (mu k) (sdev k))])
                                                 m))])))
                            S)]
           (with-open [wrtr (io2/writer "data/150703/zscore-2clusters.csv")]
             (let [kmers (add-letters 2 ["A" "C" "G" "T"])]
               (.write wrtr (str  "cluster id,seq id," (str/join "," kmers) "\n"))
               (doseq [[cid id d] data]
                 (.write wrtr (str cid "," id ","
                                   (str/join "," (map #(get d % ".") kmers)) "\n")))))))
