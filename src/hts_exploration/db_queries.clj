(ns hts-exploration.db-queries)

(defn round-all-usable-seqs [round]
  [(str "SELECT sr.selex_id as id, SUBSTRING(sr.sequence, sr.usable_start+1, sr.length) as hit_seq
          FROM selex_reads as sr
         WHERE sr.usable=1
           AND sr.strand=1
           AND sr.round_number=?")
   round])

"SELECT SQL_CALC_FOUND_ROWS foo.hit_seq FROM (SELECT SUBSTRING(sr.sequence, sr.usable_start+1, sr.length) as hit_seq FROM selex_reads as sr WHERE sr.usable=1 AND sr.strand=1) as foo WHERE foo.hit_seq REGEXP 'A[ACT][GC][ACG]A.{3,}T[CGT][GC][TGA]T' LIMIT 10; select found_rows();"
