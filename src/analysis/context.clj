(ns analysis.context
  (:require [clojure.java.jdbc :as jdbc]
            [clojure.core.reducers :as r]
            [clojure.core.match :refer [match]]
            [clojure.java.io :as io2]
            [clojure.contrib.string :as str]
            [edu.bc.fs :as fs]
            [clojure.java.shell :as shell]
            [clojure.pprint :as pp]
            [edu.bc.bio.sequtils.files :refer (read-seqs)]
            [hts-exploration.hts-utils.file :refer (write-fasta)]
            [edu.bc.utils :refer (sum transpose log pxmap reversev)]
            [edu.bc.bio.seq-utils :refer (markov-step generate-rand-seq dint-seq-shuffle)]
            [instaparse.core :as insta]
            [clojure.zip :as zip]
            [edu.bc.utils.fold-ops :refer (fold bpdist)]
            [hts-exploration.background-set :as bg]
            [refold :as refold])
  (:use edu.bc.utils.probs-stats
        hts-exploration.db-queries
        hts-exploration.globals
        hts-exploration.utils))

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

(defn check-3helix
  "Looks for 3 separate hairpins in a parsed structure"
  
  [loc]
  (let [rights (->> loc zip/down zip/rights)
        down (->> loc zip/down zip/node)]
    (->> (conj rights down)
         (map first )
         (filter #(= % :stem))
         count)))

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

(defn valid-stack? [l r st]
  (let [valid? (fn [m]
                 (and (>= (reduce + (vals m)) 8) ;min num of bases to look at
                      (>= (get (probs m) \( 0) 0.75)))]
    (->> (subs st l r) (freqn 1) valid?)))

(defn potential-stacks
  [st stack-locs]
  (->> stack-locs
       (map #(apply min %));stack start position 
       (mapcat (fn [n]
                 (let [before (max 0 (- n 12))
                       after (min (+ n 12 2) (count st))]
                   ((juxt (partial valid-stack? before n);bp stack prior
                          (partial valid-stack? (+ 2 n) after));;+2 for 2bp stack and 12 after
                    st))))))

(defn get-motif-pos
  "Checks a structure for the base pairs (stack1 stack2) and reports locations
  where stacks are adjacent and proper pairs. "
  ([stack1 stack2 s] (get-motif-pos stack1 stack2 s (-> s fold ffirst)))
  ([stack1 stack2 s st]
   (let [valid-bp-pos? (fn [[[bp1 p1] [bp2 p2]]]
                         (let [[i j] p1
                               [k l] p2]
                           (and (= (inc i) k)
                                (= (dec j) l)
                                (match [bp1 bp2]
                                       [stack1 stack2] true
                                       [(stack1 :<< reversev) (stack2 :<< reversev)] true
                                       [stack2 stack1] true
                                       [(stack2 :<< reversev) (stack1 :<< reversev)] true
                                       :else false))))]
     (->> (refold/make-pair-table st)
          (remove (fn [[i j]] (> i j)))
          (sort-by first)
          (map (fn [x] [(mapv #(.charAt s %) x) x])); get bases at i,j locations
          (partition-all 2 1)
          (filter valid-bp-pos?)
          (map (fn [[[_ p1] [_ p2]]] (vec (concat p1 p2))))))))

(defn gugc-stack? [s]
  (let [st (-> s fold ffirst)]
    (->> (get-gugc-pos s st)
         first ;ignore counts
         (potential-stacks st)
         (some true?))))

(defn motif-stack?
  "takes a seq s and a set of motif locations i.e. (get-gugc-pos s st)"
  ([motif-bp1 motif-bp2 s]
   (motif-stack? motif-bp1 motif-bp2 s (-> s fold ffirst)))
  ([motif-bp1 motif-bp2 s st]
   (->> (get-motif-pos motif-bp1 motif-bp2 s st)
        (potential-stacks st)
        (some true?))))

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

;;;-----------------------------------------------------------------------------------

;;;---------sql queries---------------------------------------------------------------

(defn get-seqs-in-cluster
  "Takes a cluster ID and gets the distinct seqs in the cluster"

  [cid]
  [(str "select sr.selex_id as id, cdhit.cluster_size as csize,
            substring(sequence, usable_start+1,sr.length) as hitseq
            FROM selex_reads as sr, cdhit where
            sr.selex_id=cdhit.selex_id AND sr.usable=1 AND sr.strand=1
            AND cluster_id=?
           GROUP by hitseq;") cid])

(defn q1
  "Gets all clusters (ordered by the number of distinct sequences)
  which have cluster_size > csize with the number of distinct seqs >
  dsize and mean_ident > ident"
  
  [csize ident dsize]
  [(str "select cdhit.cluster_id as cid,
                cdhit.cluster_size,
                count(distinct substring(sr.sequence, sr.usable_start+1, sr.length)) as cnt,
                avg (cdhit.clstr_iden) as mean_ident
           FROM cdhit, selex_reads as sr
          WHERE sr.selex_id=cdhit.selex_id
            AND cluster_size > ?
            AND sr.strand=1
       GROUP BY cdhit.cluster_id having mean_ident > ?
            AND cnt > ?
       ORDER BY cnt DESC") csize ident dsize])

(defn most-freqn-in-cluster [n cid]
  (letfn [(q2 [cid]
            [(str "SELECT substring(sr.sequence, usable_start+1, sr.length) as hitseq,
                          count(sr.selex_id) as cnt , sr.selex_id as id, cluster_id as cid
                     FROM selex_reads as sr, cdhit
                    WHERE sr.selex_id=cdhit.selex_id and cluster_id=?
                 group by hitseq
                 order by cnt DESC") cid])]
    (->> (q2 cid) sql-query (map (juxt :hitseq :id :cnt))
         (take n))))

;;;------------------------------------------------------------------------------

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
