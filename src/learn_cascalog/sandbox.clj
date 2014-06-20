(ns learn-cascalog.sandbox
  (:require [clojure.contrib.io :as io]
            [clojure.contrib.string :as str]
            [cascalog.logic.def :as d]
            [cascalog.more-taps :as mt]
            [net-eval.core :as ne])
  (:use [cascalog.api]
        [cascalog.playground]
        [edu.bc.utils.fold-ops]
        globals))

(bootstrap-emacs)

(defn- dist [ref-st st]
  (/ (bpdist ref-st st :bpdist true)
     (count ref-st)))

(d/defmapfn find-leader [s] (re-find #"TGCG.{70,90}GTACT" s))
(d/deffilterfn similar? [ref-st st] (< (dist ref-st st) 0.1))

(defn -main
  "Some sort of example program to run"

  [f]
  (let [data (cascalog.more-taps/hfs-delimited f :delimiter "\t" :skip-header? true)
        ref-st (fold ref-seq)]
    (?<- (stdout)
         [?seq ?leader ?st]
         (data ?name ?seq ?qual)
         (find-leader ?seq :> ?leader)
         (fold ?leader :> ?st)
         (similar? ref-st ?st))))



