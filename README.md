# hts-exploration

A Clojure library containing methods for analyzing HTS-SELEX data found in the SRA --- SUB1661164

## Usage

The primary methods for analyzing the HTS reads for 2\_2 nucleotide cyclic motifs (NCM) is found in the analysis.context namespace. The function is calc-struct-features2 which takes the arguments stack-size, max gaps in current helix, sequence and structure tuple. 

```clojure
;;; for the 2\_2 NCM, the number of base-pairs stacked is 2. The max bulge size is 3.
;;; The format of the data must be a sequence (s) and a structure tuple containing
;;; [structure free-energy] where the structure is represented in dot-bracket notation.
;;; This tuple is the structure line output from secondary structure prediction programs such
;;; as RNAfold

(calc-struct-features2 2 3 s [st nrg]); 2\_2 NCM, 3 base bulge/interior loop

(calc-struct-features2 3 3 s [st nrg]); 3\_3 NCM, 3 base bulge/interior loop

```
## License

Copyright © 2014 FIXME

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
