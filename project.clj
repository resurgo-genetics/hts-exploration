(defproject hts-exploration "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  
  :repl-options {:port 4011}
  :profiles {:babs {:jvm-opts ["-Xmx10g"]}
             :kitia {:jvm-opts ["-Xmx1g"]}}

  :plugins [[cider/cider-nrepl "0.13.0"]
            ]
  
  :global-vars {*print-length* 100}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [org.clojure/tools.nrepl "0.2.12"]
                 [org.clojure/core.async "0.2.385"]
                 [org.clojure/core.match "0.3.0-alpha4"]
                 [mlabs.jars/clojure-contrib "1.2.0-mlab"]
                 [medley "0.7.0"]
                 [me.raynes/fs "1.4.4"]
                 [slingshot "0.12.1"]
                 
                 ;; math stuff
                 [org.clojure/math.numeric-tower "0.0.4"]
                 [net.mikera/core.matrix "0.42.0"]

                 ;; shell stuff
                 [me.raynes/conch "0.8.0"]
                 [clj-shell "0.1.0"]
                 
                 ;;; dealing with nested stuff
                 [com.rpl/specter "0.7.1"]
                 [instar "1.0.10"]
                 
                 ;; exporting data should figure out one that works
                 [org.clojure/data.csv "0.1.3"]
                 [clojure-csv/clojure-csv "2.0.1"]
                 [semantic-csv "0.1.0"]; handle csvs in a better way
                 [cheshire "5.5.0"]; json
                 
                 [instaparse "1.4.1"]
                 
                 [org.clojure/java.jdbc "0.3.3"]
                 [mysql/mysql-connector-java "5.1.18"]
                                  
                 [cascalog "2.0.0"]
                 [org.apache.hadoop/hadoop-core "1.1.2"]

                 ;;file reading in parallel
                 [iota "1.1.1"]
                 [foldable-seq "0.2"]

                 ;;pdf manipulation
                 [clj-pdf "1.11.20"]
                 [camelot "0.2.0"]

                 ;;graph lib
                 [aysylu/loom "0.5.0"]
                 
                 ;;graph visualization
                 [rhizome "0.2.1"]
                 ])
