(ns learn-net-eval.sandbox
  (:require [clojure.contrib.io :as io]
            [clojure.contrib.string :as str]
            [net-eval.core :as ne]
            [clojure.core.reducers :as r]
            [clojure.pprint :as pp]
            [iota]
            [foldable-seq.core :as fseq]
            gibbs-sampler)
  (:use [edu.bc.utils.fold-ops]
        hts-explortation.globals))

(def babs "136.167.49.162")
(def roz "136.167.54.218")

(ne/deftask foo [st data]
  (->> data
       ;;(sample-data 100 )
       (dist-filter 0.2 st)
       (get-regions (- (count prime5-const) 3)
                    (+ 30 (count prime5-const);var+const region
                       (count prime3-const)))
       ))




(ne/deftask sum-and-print-task [x y]
 (let [s (+ x y)]
  (do
   (println s)
   (r/fold + (r/map inc (range 10))))))

(ne/deftask toy-fold-ex [s]
  (prn (type s))
  (r/fold + (r/map inc s)))

(comment
  (def response (ne/net-eval [[babs 4013 #'sum-and-print-task 3 4]
                              [babs 4013 #'toy-fold-ex (vec (range 10))]]))

  (println (map deref response)))
