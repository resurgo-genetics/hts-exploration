(ns hts-exploration.globals
  (require [edu.bc.fs :as fs])
  (use [edu.bc.bio.seq-utils :only [reverse-compliment]]))

(def forward-re #"TGCG.{70,90}GTACT")
(def reverse-re #"AGTAC.{70,90}CGCA")

(def ex-file "/home/kitia/data/index46_TCCCGA_combined_R1_001.filtered.fastq")
(def S15-r1-csv "/home/peis/S15SELEXHTSdata/index46_TCCCGA_combined_R1_001.csv"); round 4
(def S15-r1-qual-filter-csv (fs/replace-type S15-r1-csv ".qual-filter.csv"))

(def S15-round4-fastq "/home/peis/S15SELEXHTSdata/index46_TCCCGA_combined_R1_001.fastq")

(def S15-round9-fastq "/home/peis/S15SELEXHTSdata/oto9952/Index47_Ot9952_Michelle.Meyer_400-595b_PCR_amplicon_1M_08142013/index47_TCGAAG_combined_R1_001.fastq")

(def S15-round10-fastq "/home/peis/S15SELEXHTSdata/Sequencing Data/400-595c/Index8_OtA2903_Michelle.Meyer_400-595c_PCR_1M_05022014/index8_ACTTGA_L001-L002_R1_001.fastq")
(def S15-round10-qual-filter-csv (fs/replace-type S15-round10-fastq ".qual-filter.csv"))

(def S15-round11-fastq
  "/home/peis/S15SELEXHTSdata/Sequencing Data/400-595c/Index9_OtA2904_Michelle.Meyer_400-595c_PCR_1M_05022014/index9_GATCAG_L001-L002_R1_001.fastq")
(def S15-round11-qual-filter-csv
  (fs/replace-type S15-round11-fastq ".qual-filter.csv"))

(def t7 "TAATACGACTCACTATAgg")
(def prime5-const "TGCGTAACGTACACT")
(def variable-region "GGGATCGCTGAATTAGAGATCGGCGTCCTT");stem
(def prime3-const "TCATTCTATATACTTTGGAGTTTTAAA");after stem
(def cds "ATGTCTCTAAGTACT") ;pseduonoted partner

(def parasite [(str prime5-const "A+" cds)
               "TGCGT.*TGCGT"
               "TGCGCCTTCGTATGT" ;random seq doesn't fit
               "TAATACGACTCACTATA"  ;very similar to t7
               "AAGTACTGAAG"]) ;random seq doesn't fit

(def parasite-rev (flatten
                   ["AGTACTTAGAGACATT+AGTGTACGTTACGCA"
                    "ACGCA.*ACGCA"
                    (->> parasite (drop 2)
                         (map reverse-compliment))]))

(def ref-seq
  (str prime5-const variable-region prime3-const cds))

(def ref-seq-firmicute
  "UCAGUGUAUGCGAACCGUUGCUUGGCCAGGCGACUCACCGACGCCCGCUCGGCAACCGGGGAUCGAAGACUAGGGAGGUGAACCAUGAUGGCAUUGA")

(def bob
  (str ;"TAATACGACTCACTATAGG"
       "TGCGTAACGTACACT"
       "TCCTTCGCTTATTCGGAGTAGATCACGTGA"
       "TCATTCTGTATGCTTTGGAGTTTTAAAATGTCTCTAAGTACT"))

(def rna-8-11-P
  (str "TGCGTAACGTACACT"
       "CCCTACTCGTGGATTGGACTCTATAATAGAT"
       "CATTCTATATACTGTGGAGCTTAAA"
       "ATGTCTCTAAGTACT"))

(def rna-8-12-35
  (str "TGCGTAACGTACACT"
       "AAAGACGGAAGGCCAAAACCTATGCTTACCTTTACTTTGGAGTTTTAAA"
       "ATGTCTCTAAGTACT"))

(def round4-count 2060481)
(def round9-count 2089600)
(def round10-count 481766)
(def round11-count 545144)

(def tested-seqs [{:name :5motif :id 3196036 :cluster 62979 :csize 1
                   :ident 100 :seq "ATCGAAAGGAGAATGGAATCGAGCAATCGA"} 
                  {:name :4motif :id 455019 :cluster 51432 :csize 1
                   :ident 100 :seq "GTCAAACAAAACAAAAGACGAAGACGCACT"}
                  {:name :2motif :id 4077286 :cluster 1503 :csize 72
                   :ident 85.19 :seq "AACGATTCGAAAGTGAAAGAAAGAGAAA"}
                  {:name :ci7-255 :id 1895496 :cluster 2346 :csize 367
                   :ident 94.83 :seq "CGATCACACGAGAACATCGGTGATTTGGTGTCAAT"}
                  {:name :ci7-73 :id 2774645 :cluster 1887 :csize 562
                   :ident 89.83 :seq "AGGCAAACCGATCCTAACGAATGCTTGGTG"}
                  {:name :ci17-156 :id 2152795 :cluster 3580 :csize 728
                   :ident 94.92 :seq "ACCCAAGACGGCTCTACAGTAAGATAGCCTA"}
                  {:name :bob :id 3597107 :cluster 9131 :csize 61
                   :ident 98.25 :seq "TCCTTCGCTTATTCGGAGTAGATCACGTGA"}
                  {:name :wka :id 5267819 :cluster 9391 :csize 8825
                   :ident 99.29 :seq "GCGGACAGCGAGACAGATCGAAGGTTTTGA"}
                  {:name :A :id 11025540 :cluster 1079 :csize 1022
                   :ident "not calculated" :seq "none" :reason "cdhit cluster"}
                  {:name :B :id 10951764 :cluster 1920 :csize 3718
                   :ident "not calculated" :seq "none" :reason "cdhit cluster"}
                  {:name :C :id 10999049 :cluster 2917 :csize 1648 :reason "multiton"}
                  {:name :D :id 11029916 :cluster 464 :csize 4845 :reason "multiton"}
                  {:name :E :id 10594516 :cluster 62378 :csize 62378 :reason "singleton large cluster"}
                  {:name :F :id 1707334 :cluster 77409 :csize 1 :reason "singleton singleton cluster"}
                  {:name :G :id 10632040 :cluster 697
                   :csize 490 :reason "cluster contains structure similar to wka"}
                  {:name :H :id 39485 :cluster 303 :csize 1912 :reason "ratio less than 1"}])
