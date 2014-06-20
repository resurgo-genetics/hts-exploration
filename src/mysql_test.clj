(ns mysql-test
  (:require [clojure.java.jdbc :as jdbc])
  (:use [edu.bc.bio.sequtils.files]))




(def mysql-ds
  {:classname "com.mysql.jdbc.Driver"
   :subprotocol "mysql"
   :subname "//127.0.0.1:3306/biosql"
   :user "root"
   :password "rna314rulz"})

(defn sql-query [stmt & {:keys [f p] :or {f identity p false}}]
  (let [q (partial jdbc/query mysql-ds)]
    (when p (println stmt))
    (cond (string? stmt) (q [stmt])
          (vector? stmt) (q stmt)
          :else "Invalid query")))

(def genome-counts
     "select count(*) from
      (select distinct tx.taxon_id
              from bioentry as be,
                   taxon as tx,
                   ancestor as an
              where be.taxon_id=tx.taxon_id
              and tx.ncbi_taxon_id=an.ncbi_taxon_id
              and be.name regexp \"^NC\"
              and be.description not regexp \"plasmid\"
              and an.ancestors regexp \"/taxon/\") as foo")

(def query-operon
  "SELECT be.bioentry_id, be.name, tx.taxon_id, tx.ncbi_taxon_id, opl.*
   FROM taxon as tx, bioentry as be, ancestor as an, operon as op, operon_loc as opl
   WHERE tx.taxon_id=be.taxon_id
     AND tx.ncbi_taxon_id=an.ncbi_taxon_id
     and an.ancestors REGEXP 'Gamma*'
     and op.bioentry_id=be.bioentry_id
     and be.description not REGEXP \"plasmid\"
     and op.operon_id=opl.operon_id
   LIMIT 3;")

(def query-operon2
  "SELECT be.accession, tx.taxon_id, tx.ncbi_taxon_id, opl.*,be.description
   FROM taxon as tx, bioentry as be, ancestor as an, operon as op, operon_loc as opl
   WHERE tx.taxon_id=be.taxon_id
     AND tx.ncbi_taxon_id=an.ncbi_taxon_id
     and be.description not REGEXP \"plasmid\"
     and an.ancestors REGEXP 'enterobacteria'
     and op.bioentry_id=be.bioentry_id and op.operon_id=opl.operon_id LIMIT 3")

(def query-operon3
  "SELECT *
   FROM bioentry as be, seqfeature as sf, seqfeature_qualifier_value as sq,
        location as loc
   WHERE be.bioentry_id=sf.bioentry_id
     AND sf.seqfeature_id=sq.seqfeature_id
     AND be.description not REGEXP 'plasmid'
     AND loc.seqfeature_id=sf.seqfeature_id
     AND sq.value REGEXP '^rp[sml]'
     AND be.bioentry_id=738
   LIMIT 3;")

(def query-operon4
  ["SELECT *
   FROM bioentry as be, seqfeature as sf, seqfeature_qualifier_value as sq,
        location as loc
   WHERE be.bioentry_id=sf.bioentry_id
     AND sf.seqfeature_id=sq.seqfeature_id
     AND be.description not REGEXP 'plasmid'
     AND loc.seqfeature_id=sf.seqfeature_id
     AND sq.value REGEXP ?
     AND be.bioentry_id = ?
   LIMIT 3;"
   "rpl" 738]) ;args

(defn map->entry [sqlout]
  (let [{:keys [accession start_pos end_pos strand]} sqlout]
    (make-entry [accession [start_pos end_pos] strand])))

(->> (sql-query query-operon2)
     (map #(select-keys % [:accession :start_pos :end_pos :strand]))
     (map map->entry)
     (map gen-name-seq))
