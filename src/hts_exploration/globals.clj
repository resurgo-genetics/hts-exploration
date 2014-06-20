(ns hts-exploration.globals
  (require [edu.bc.fs :as fs])
  (use [edu.bc.bio.seq-utils :only [reverse-compliment]]))

(def forward-re #"TGCG.{70,90}GTACT")
(def reverse-re #"AGTAC.{70,90}CGCA")

(def ex-file "/home/kitia/data/index46_TCCCGA_combined_R1_001.filtered.fastq")
(def S15-r1-csv "/home/peis/S15SELEXHTSdata/index46_TCCCGA_combined_R1_001.csv")
(def S15-r1-qual-filter-csv (fs/replace-type S15-r1-csv ".qual-filter.csv"))

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

