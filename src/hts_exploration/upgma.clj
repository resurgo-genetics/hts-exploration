(ns hts-exploration.upgma
  (:require clojure.set))

(def ^:dynamic *M*)
(defn min-entry [m] (first (sort-by second m)))

(defn distfn [c1 c2]
  (let [c1 (flatten c1)
        c2 (flatten c2)]
    (/ (reduce + (for [i c1 j c2] (*M* (mapv vector (sort [i j])))))
      (* (count c1) (count c2)))))

(defn remove-keys-contain [x coll]
  (->> coll
       (remove (fn [[k v]]
                 (some #(contains? (set (flatten x)) %) (flatten k))))
       (into {})))

(defn re-calc [min m distfn]
  (let [rkc (remove-keys-contain min m)]
    ;(prn :remove-keys rkc)
    (reduce (fn [m x] (let [x (if (vector? x) x [x])
                           d (distfn min x)]
                       (prn :min min :x x :dist d)
                       (assoc m [min x] d)))
            rkc
            (->> rkc keys (apply concat) distinct))))

(defn bfs-eager [tree]
  (loop [ret []
         queue (conj clojure.lang.PersistentQueue/EMPTY tree)]
    (if (seq queue)
      (let [[node & children] (peek queue)]
        (recur (conj ret node)
               (into (pop queue) children)))
      ret)))

(defn tree-traversal [tree]
  (loop [t tree
         stack []]
    (let [x (first (pop t))]
      (prn :pop x :stack stack :count (count (flatten x)))
      (if (>= (count (flatten x)) 2)
        (recur x
               (conj stack (peek t)))
        [t stack]))))

(defn add-weights [weights tree]
  (let [S (tree-traversal tree)
        add-edge (fn [node new-edge] (conj node (/ (weights new-edge) 2)))]
    (second
     (reduce (fn [[k V] x]
               (let [nextk [k x]]
                 [nextk [V (add-edge x nextk)]]))
             [(first S) (add-edge (first S) (first S))]
             (reverse (second S))))))

(defn -main [distfn dist-matrix]
  (binding [*M* dist-matrix]
    (loop [c 10
           next-min (min-entry *M*)
           m *M*
           tree (apply hash-map next-min)]
      (prn :>>> next-min)
      (let [cluster (re-calc (first next-min) m distfn)]
        (prn :cluster cluster)
        (if (or (empty? cluster) (neg? c)) 
          (let [foo (clojure.set/difference
                     (set (flatten (keys dist-matrix)))
                     (set (flatten (first next-min))))
                end-tree [(first next-min) (vec foo)]]
            
            [(assoc tree end-tree (apply distfn end-tree))
             end-tree])
          (recur (dec c) 
                 (min-entry cluster)
                 cluster
                 (apply assoc tree (min-entry cluster))))))))

(let [s "ABCDEFG"
      M (apply hash-map
               (interleave 
                (distinct
                 (for [i s j s :when (not= i j)] (mapv vector (sort [i j]))))
                [19 27 8 33 18 13 31 18 36 1 13 
                 26 41 32 29 31 17 14 35 28 12]))]
  (apply add-weights (-main distfn M)))
