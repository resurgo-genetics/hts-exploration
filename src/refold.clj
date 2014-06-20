(ns refold
  (:require [clojure.contrib.string :as str]))

;;creates a list of where the base pairings are located
;;according to index number
(defn make_pair_table [structure]
  (loop [c (rest (str/split #"" structure))
         s [] ;;stack structure which stores the location of the
         ;;opening ( so that closing parens pop location
         table {}
         i 0]
    ;;(prn "c=" c s table)
    (if(seq c)
      (do 
        (cond
         ;;opening bracket for structure stores location in stack
         (or (= (first c) "(") (= (first c) "<"))
         ;;(recur (rest c) (cons i s) (cons i table) (inc i))
         (recur (rest c) (cons i s) table (inc i))

         ;;closing paren pops location off stack and creates a key
         ;;value pair where the 2 base paired bases are put into the map
         (or (= (first c) ")") (= (first c) ">"))
         (recur (rest c)
                (rest s)
                (assoc table (first s) i
                       i (first s))
                (inc i))

         ;;gap does nothing
         :else (recur (rest c) s table (inc i))))
      table ;;returns a map
      )))


(defn remove-gaps
  "attempts to remove gaps from a sequence and corresponding structure
   and tries to ensure equal number of ( and ). Takes 2 strings, must
   be upper-case "

  ([inseq] (remove-gaps inseq (apply str (take (count inseq) (cycle ".")))))
  
  ([inseq struct]
     (let [table (make_pair_table struct)
           pairs {"AU" 1 "GC" 1 "CG" 1 "UA" 1 "AT" 1
                  "TA" 1 "GU" 1 "GT" 1 "TG" 1 "UG" 1}] ;;possible base pairings
       (loop [i 0
              st struct]
         (if (< i (count inseq))
           (let [c (.substring inseq i (inc i)) ;;char at position i
                 compi (get table i -1)] ;;where i is base paired to
             ;; (prn "i=" i "compi=" compi "c=" c "length=" (count st))
             ;; (prn inseq)
             ;; (prn st)
             (cond
              ;;encounters gap in sequence = removal of that from
              ;;structure if that location is base paired then its
              ;;complement is turned into a "." in the structure
              (= c ".")
              (recur (inc i)
                     (cond
                      (and (< compi i)
                           (pos? compi)
                           (not= "x" (subs st compi (inc compi))))
                      (str (subs st 0 compi) "." (subs st (inc compi) i) "x" (subs st (inc i)))
                      
                      (and (pos? compi)
                           (not= "x" (subs st compi (inc compi))))
                      (str (subs st 0 i) "x" (subs st (inc i) compi) "." (subs st (inc compi)))
                      
                      :else
                      (str (subs st 0 i) "x" (subs st (inc i)))))
             
              ;;encounters a location where base pair was thought to
              ;;occur but there is no proper base pairing due because
              ;;sequence doesn't contain necessary bases
              (and (> compi i) (not (contains? pairs (str c (.substring inseq compi (inc compi))))))
              (recur (inc i)
                     (str (subs st 0 i) "." (subs st (inc i) compi) "." (subs st (inc compi))))
              :else
              (recur (inc i)
                     st)))
         [(str/replace-re  #"\." "" inseq) (str/replace-re #"x" "" st)])))))


