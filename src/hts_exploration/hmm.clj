(ns hts-exploration.hmm
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :as math])
  (:use hts-exploration.globals
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.sequtils.files
        [edu.bc.bio.seq-utils :only (markov-step)]))

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
         (map (fn [j]
                (* (alpha path j state)
                   (beta path j state)))
              (range (count path)))))))
