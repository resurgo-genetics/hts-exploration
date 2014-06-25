(ns smith-waterman)

(def ^{:dynamic true :private true} global)
(def ^{:dynamic true :private true} match)
(def ^{:dynamic true :private true} mis)
(def ^{:dynamic true :private true} gap-open)
(def ^{:dynamic true :private true} gap-ext)

(defn- array-keys
  "positions of the array to work on"

  [s1 s2]
  (for [i (range (count s1)) ;initialize scoring array. similar to a sparse matrix
        j (range (count s2))]
    [i j]))

(defn- init-array
  "initial conditions of the array to fill in for the gapped row/col"
  
  [s1 s2 locations gap global]
  (reduce (fn [m [i j]]
            (assoc m [i j]
                   (cond
                    (= i j 0) [0 [0 0] "-" "-"]
                    (= i 0) [(if global (* gap j) 0) :l "-" (.charAt s2 j)]
                    (= j 0) [(if global (* gap i) 0) :u (.charAt s1 i) "-"])))
          {} (filter #(or (zero? (first %))
                          (zero? (second %))) locations)))

(defn- dir->coord [dir i j]
  (case dir
    :d [(dec i) (dec j)]
    :u [(dec i) j]
    :l [i (dec j)]))

(defn- get-score [m dir i j] (->> (dir->coord dir i j) m))

(defn- maxa
  "Determine max score and direction for position [i j]"
  [coll]
  (->> coll
       (filter #(= (apply max (map second coll)) (second %)))
       (sort-by first)
       first))


(defn- gapfn [recur-from]
  (case recur-from
    :d gap-open
    :u gap-ext
    :l gap-ext))

(defn- fill-array
  
  [locations s1 s2]
  (reduce (fn [m [i j]];;score array format
            (let [[d dfrom] (get-score m :d i j) ;;score match/mismatch (diagonal)
                  [u ufrom] (get-score m :u i j) ;;score deletion (above)
                  [l lfrom] (get-score m :l i j) ;;score insertion (left)
                  aa1 (.charAt s1 i) ;;current char in s1
                  aa2 (.charAt s2 j) ;;current char in s2
                  [from score] ;;chooses from d, u, l and scores associated with it.
                  (maxa [(if (= aa1 aa2 ) 
                           [:d (+ d match)]
                           [:d (+ d mis)])
                         [:u (+ u (gapfn ufrom))]
                         [:l (+ l (gapfn lfrom))]])]
              (assoc m [i j] ;;insertion of the best score into the matrix
                     (case from
                       :d [score :d aa1 aa2]
                       :u [score :u aa1 "-"]
                       :l [score :l "-" aa2]))))
          (init-array s1 s2 locations gap-open global)
          (remove #(or (zero? (first %))
                       (zero? (second %))) locations)))

(defn- trace-back

  [score start-loc H]
  (loop [loc start-loc
         aln_s1 ""
         aln_s2 ""]
    (let [[_ dir a1 a2] (get H loc) ;;stores the next location [score[i j] to go to in H]
          next-coord (apply dir->coord dir loc)] 
      (if (not= [0 0] next-coord)
        (recur next-coord
               (str a1 aln_s1) ;;builds strings up from the right to left
               (str a2 aln_s2))
        (if (= "-" a1 a2)
          [score aln_s1 aln_s2]
          [score (str a1 aln_s1) (str a2 aln_s2)])))))

(defn sw [seq1 seq2 & {:keys [global
                              anchor-right
                              match-weight
                              mismatch-weight
                              gap-open-weight
                              gap-ext-weight]
                       :or {global false
                            anchor-right false
                            match-weight 2
                            mismatch-weight -1
                            gap-open-weight -1
                            gap-ext-weight -1}}]
  (binding [global global
            match match-weight  ;match
            mis mismatch-weight ;mismatch
            gap-open gap-open-weight
            gap-ext gap-ext-weight]
    (let [s1 (str "-" seq1)
          s2 (str "-" seq2)
          H (fill-array (array-keys s1 s2) s1 s2) ;creates score matrix
          start (cond global
                      (get H (mapv dec [(count s1) (count s2)]))
                      anchor-right
                      (->> (filter #(= (ffirst %) (dec (count s1))) H)
                           (sort-by #(-> % second first) >)
                           first
                           second)
                      :else
                      (-> (sort-by #(-> % second first) H) ;finds highest value in matrix
                          last
                          second))
          start-locs (map first (filter #(= start (val %)) H))] ;;starts traceback from this highest value
      (mapv #(trace-back (first start) % H) start-locs))))

(def test-case ([12 "a-cacacta" "agcacac-a"]))
