(ns analysis.motif-ratio-rounds
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.io :as io]
            [clojure.contrib.string :as str]
            [clojure.core.reducers :as r]
            [clojure.java.jdbc :as jdbc]
            [clojure.math.numeric-tower :as math]
            [clojure.pprint :as pp]
            )
  (:use edu.bc.utils.fold-ops
        edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.utils.probs-stats
        edu.bc.utils
        hts-exploration.db-queries
        hts-exploration.globals
        hts-exploration.utils))

(defn kmerfn
  "Takes a length k and a set of sequences and computes the mean
  probability of seeing the kmer in a sequence."
  
  [k inseqs]
  (->> (map #(probs k %) inseqs)
       (apply merge-with +)
       (reduce (fn [M [kmer v]]
                 (assoc M kmer (double (/ v (count inseqs)))))
               {})))

(defn motif-ratio
  "Finds the ratio of mean motif rate in 2 separate rounds (ra,
  rb). Reports the top 10 hits"

  [n ra rb]
  (time
   (let [get-seqs (fn [round]
                    (->> (round-all-usable-seqs round) sql-query
                         (map :hit_seq )
                         (map #(second
                                (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" %)))
                         (remove nil?)
                         distinct))
         sqsa (get-seqs ra)
         sqsb (get-seqs rb)]
     (->> (merge-with / (kmerfn n sqsa)
                      (kmerfn n sqsb))
          (sort-by val > )
          (take 10)))))

;;round 11 vs round 10
(motif-ratio 6 11 10)
"Elapsed time: 980305.595233 msecs"
(["GTGTCT" 2.06174109652782] ["GCGGTG" 2.0616517848451674] ["CCTGCG" 2.042251885207106] 
 ["CTGCGG" 1.9855611802607638] ["TGCGGT" 1.9823356229670237] ["GCGGAC" 1.966339850816003] 
 ["TTGATC" 1.9568510477633874] ["GTTTTG" 1.9135521648495453] ["TTTGAT" 1.890558758524509] 
 ["TGTCTT" 1.874874631767203])

;;after forcing all seqs to be distinct
(["GGTGTC" 1.6428656241397235] ["CGGGTA" 1.6195919644356367] ["CCGGGT" 1.5897515542640006] 
 ["GCGGTG" 1.5739453570424213] ["TGGGTG" 1.5727188947041364] ["GGGTAT" 1.5714839436422858]
 ["CGGTGT" 1.5652489716303113] ["TGGTGT" 1.5467607633840224] ["TTGTGG" 1.535166940421112]
 ["ACTCTT" 1.4982984993353305])

;;round 11 vs round 4
"Elapsed time: 2516346.633465 msecs"
(["TTGGTT" 25.987668137754557] ["TGGTTT" 25.89882809444897] ["TTTGGT" 20.17223022804577]
 ["CTTGGT" 19.287400487754503] ["CTTCTT" 17.46616736742761] ["GGTTTA" 17.435844990074866]
 ["TTCTTT" 16.394382056124638] ["TTGGTG" 13.530774021665612] ["TGGTGT" 12.088138188711829]
 ["TCTTTG" 11.876481436902408])

;;round 11 vs round 9
(["TTGGTT" 25.91787330674479] ["TGGTTT" 25.845936514710676] ["TTTGGT" 19.877741563560754]
 ["CTTGGT" 19.24089275573208] ["CTTCTT" 17.549219268199415] ["GGTTTA" 17.50697341968236]
 ["TTCTTT" 16.634552819390308] ["TTGGTG" 13.529265555173582] ["TCTTTG" 11.941030978936078]
 ["TGGTGT" 11.848737115678052])

(defn locate-motifs
  "attempts to locate the position of motifs in round 11 based on the
  most prevalent motifs found when comparing different rounds"

  [res]
  (time
   (let [get-seqs (fn [round]
                    (->> (round-all-usable-seqs round) sql-query
                         (map :hit_seq )
                         (map #(second
                                (re-find #"TGCGTAACGTACACT(.*)ATGTCTCTAAGTACT" %)))
                         (remove nil?)
                         distinct))
         sqs11 (get-seqs 11)
         cnt (count sqs11)]
     (reduce (fn [M re]
               (let [sqs (->> (map #(vector % (str-re-pos (re-pattern re) %)) sqs11)
                              (remove #(-> % second empty?)))
                     motifs-per-seq (->> (map second sqs) (map count) frequencies)
                     cutoff (fn [cnt thr] (<= (/ cnt (reduce + (vals motifs-per-seq)))
                                            thr))]
                 (assoc M re {"starts" (->>
                                        (pxmap (fn [[s m]]
                                                 (let [ks (keys m)
                                                       c (const-start s)]
                                                   (map #(/ (inc %) c) ks)))
                                               20 sqs)
                                        (apply concat) frequencies
                                        (remove #(cutoff (second %) 0.05))
                                        (into (sorted-map)))
                              "motifs per seq" motifs-per-seq
                              ;; length of sequence fragment so that we understand the starts as a relative position in teh variable region
                              "lengths" (->> (map first sqs) (map count) frequencies (into (sorted-map)))
                              "percent containing motif" (double (/ (count sqs) cnt))})))
             {}
             res))))

(locate-motifs (mapv first [["GTGTCT" 2.06174109652782] ["GCGGTG" 2.0616517848451674]
                          ["CCTGCG" 2.042251885207106] ["CTGCGG" 1.9855611802607638]
                          ["TGCGGT" 1.9823356229670237] ["GCGGAC" 1.966339850816003]
                          ["TTGATC" 1.9568510477633874] ["GTTTTG" 1.9135521648495453]
                          ["TTTGAT" 1.890558758524509] ["TGTCTT" 1.874874631767203]]))

(locate-motifs (mapv first  '(["GGTGTC" 1.6428656241397235] ["CGGGTA" 1.6195919644356367]
                              ["CCGGGT" 1.5897515542640006] ["GCGGTG" 1.5739453570424213]
                              ["TGGGTG" 1.5727188947041364] ["GGGTAT" 1.5714839436422858]
                              ["CGGTGT" 1.5652489716303113] ["TGGTGT" 1.5467607633840224]
                              ["TTGTGG" 1.535166940421112] ["ACTCTT" 1.4982984993353305] )))
"Elapsed time: 671047.354179 msecs"
{"GGGTAT" {"starts" {2/31 519, 3/31 355, 4/31 353}, 
           "motifs per seq" {1 5819, 2 7}, 
           "lengths" {49 43, 50 52, 51 75, 52 134, 53 230, 54 404, 
                      55 660, 56 1194, 57 2424, 58 422, 59 109, 60 41, 61 13,
                      62 9, 63 6, 64 3, 65 3, 66 2, 68 2}, 
           "percent containing motif" 0.0127320050613326}, 
 "GGTGTC" {"starts" {27/31 8459}, 
           "motifs per seq" {1 12982, 2 41},
           "lengths" {49 85, 50 70, 51 155, 52 419, 53 751, 54 1189, 
                      55 1564, 56 2564, 57 4702, 58 1119, 59 250, 60 68,
                      61 35, 62 21, 63 12, 64 9, 65 4, 66 4, 67 1, 68 1}, 
           "percent containing motif" 0.02846016167417344}, 
 "TGGGTG" {"starts" {12/31 288, 24/31 635, 25/31 1180}, 
           "motifs per seq" {1 5069, 2 9}, 
           "lengths" {49 57, 50 45, 51 95, 52 190, 53 260, 54 399, 55 638, 56 982, 
                      57 1734, 58 490, 59 113, 60 41, 61 16, 62 10, 63 5, 66 1, 67 2}, 
           "percent containing motif" 0.01109734323746085},
 "TTGTGG" {"starts" {8/31 128, 9/31 307, 10/31 117, 23/31 149},
           "motifs per seq" {1 2310, 2 1}, 
           "lengths" {49 19, 50 22, 51 111, 52 104, 53 129, 54 218, 55 300,
                      56 508, 57 659, 58 156, 59 57, 60 17, 61 5, 62 4, 63 2}, 
           "percent containing motif" 0.005050405715197329},
 "CCGGGT" {"starts" {1/31 799, 2/31 1103}, 
           "motifs per seq" {1 7476, 2 2}, 
           "lengths" {49 35, 50 35, 51 85, 52 162, 53 285, 54 473, 55 929, 
                      56 1478, 57 3041, 58 644, 59 188, 60 55, 61 27, 62 17, 
                      63 13, 64 8, 66 3}, 
           "percent containing motif" 0.01634224748517768},
 "TGGTGT" {"starts" {26/31 2705, 45/31 808, 46/31 3182}, 
           "motifs per seq" {1 12308, 2 78}, 
           "lengths" {49 52, 50 68, 51 152, 52 316, 53 379, 54 615, 55 1346, 
                      56 2535, 57 5420, 58 1091, 59 267, 60 75, 61 34, 62 7, 
                      63 10, 64 4, 65 5, 66 4, 67 1, 68 4, 69 1}, 
           "percent containing motif" 0.02706807667175859},
 "GCGGTG" {"starts" {23/31 211, 24/31 319, 25/31 1113}, 
           "motifs per seq" {1 3762, 2 7}, 
           "lengths" {49 98, 50 47, 51 83, 52 156, 53 347, 54 581, 55 437, 56 598, 
                      57 1074, 58 260, 59 48, 60 18, 61 10, 62 5, 63 4, 64 2, 66 1}, 
           "percent containing motif" 0.008236685045685302},
 "ACTCTT" {"starts" {27/29 1007, 14/15 2600, 29/31 8072},
           "motifs per seq" {1 16761, 2 17}, 
           "lengths" {49 515, 50 883, 51 1468, 52 2673, 53 3158, 54 3447, 
                      55 2435, 56 1061, 57 768, 58 197, 59 89, 60 38, 61 24, 
                      62 9, 63 5, 64 6, 65 1, 69 1}, 
           "percent containing motif" 0.03666625144508039},
 "CGGGTA" {"starts" {1/31 1089, 2/31 828, 3/31 998}, 
           "motifs per seq" {1 10256, 2 13}, 
           "lengths" {49 60, 50 86, 51 105, 52 202, 53 368, 54 690, 55 1210,
                      56 2091, 57 4270, 58 813, 59 208, 60 79, 61 39, 62 21, 
                      63 13, 64 6, 65 4, 66 3, 68 1}, 
           "percent containing motif" 0.02244163404991838},
 "CGGTGT" {"starts" {24/31 431, 25/31 486, 26/31 3499}, 
           "motifs per seq" {1 7954, 2 11}, 
           "lengths" {49 46, 50 60, 51 90, 52 248, 53 520, 54 884, 55 866, 
                      56 1573, 57 2807, 58 609, 59 164, 60 45, 61 17, 62 16,
                      63 7, 64 7, 65 3, 66 2, 67 1}, 
           "percent containing motif" 0.01740652597211022}}

(let [sql-data (fn [r] (->> ["select hitseq, const_start as cs from selex_keys as sk join selex_reads as sr on sk.selex_id=sr.selex_id join selex_seqs as ss on sk.seq_id=ss.seq_id where round_number=?" r]
                           sql-query
                           (map (juxt :hitseq :cs))
                           (mapv (fn [[s cs]] (subs s 15 (+ 15 cs))))))
      r4-10 (future (kmerfn 6 (sql-data 10)))
      r4-11 (future (kmerfn 6 (sql-data 11)))]
  (->> (merge-with / @r4-11 @r4-10)
       (sort-by second >)
       (take 10)
       vec))

(comp #(str/take cs %) #(str/drop 15 %))
