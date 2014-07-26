(ns hts-exploration.hmm
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :as math]
            [clojure.math.combinatorics :as combo])
  #_(:use hts-exploration.globals
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.sequtils.files
        [edu.bc.bio.seq-utils :only (markov-step)]))

(defn sum [coll] (reduce + coll))

(def transitions {:s {:s 0.3 :t 0.7} :t {:s 0.1 :t 0.9}})

(def emissions {:s {:A 0.4 :B 0.6} :t {:A 0.5 :B 0.5}})

(def initial {:s 0.85 :t 0.15})

(def path [:A :B :B :A])

(defn alpha [path pos state]
  (if (zero? pos)
    (* (initial state)
       (get-in emissions [state (path pos)]))
    (sum
     (map (fn [s] (* (alpha path (dec pos) s)
                    (get-in transitions [s state])
                    (get-in emissions [state (path pos)])))
          (keys initial)))))

(defn beta [path pos state]
  (if (= (-> path count dec) pos)
    1
    (sum
     (map (fn [u] (* (get-in transitions [state u])
                    (get-in emissions [u (path (inc pos))])
                    (beta path (inc pos) u)))
          (keys initial)))))

(defn gamma [path pos state]
  (* (alpha path pos state)
     (beta path pos state)
     (/ (sum
         (map (fn [statej]
                (* (alpha path pos statej)
                   (beta path pos statej)))
              (keys initial))))))

(defn psi [path pos state1 state2]
  (let [fun (fn [path pos state1 state2]
              (* (alpha path pos state1)
                 (get-in transitions [state1 state2])
                 (beta path (inc pos) state2)
                 (get-in emissions [state2 (path (inc pos))])
                 ))]
    (/ (fun path pos state1 state2)
       (sum
        (for [k (keys initial)
              l (keys initial)]
          (fun path pos k l))))))

(reduce (fn [astar [state1 state2]]
          (update-in astar [state1 state2] (fn [_] (/ (sum 
                                                      (for [i (range (dec (count path)))]
                                                        (psi path i state1 state2)))
                                                     (sum
                                                      (for [i (range (dec (count path)))]
                                                        (gamma path i state1)))))))
        transitions
        (for [k (keys initial)
              l (keys initial)]
          [k l]))

(reduce (fn [bstar [state emit]]
          (update-in bstar [state emit] (fn [_] (/ (sum 
                                                   (for [i (range (count path))]
                                                     (* (if (= (path i) emit) (gamma path i state) 0))))
                                                  (sum
                                                   (for [i (range (count path))]
                                                     (gamma path i state)))))))
        emissions
        (for [s (keys initial)
              e [:A :B]]
          [s e]))
