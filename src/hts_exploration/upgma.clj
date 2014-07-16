(ns hts-exploration.upgma)

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

(defn tree-traversal [weights tree]
  (loop [t tree
         q []]
    (if (seq tree)
      (recur (first t) (conj q tree))
      (weights tree))))

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
                     (set (flatten (first next-min))))]
            (distfn (first next-min) (vec foo))
            [tree [(first next-min) (vec foo)]])
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
  (-main distfn M))
