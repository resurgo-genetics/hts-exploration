(defproject hts-exploration "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  
  :repl-options {:port 4011}
  :profiles {:babs {:jvm-opts ["-Xmx10g"]}
             :user {:jvm-opts ["-Xmx4g"]}}

  :plugins [[cider/cider-nrepl "0.8.2"]]
  
  :global-vars {*print-length* 100}
  :dependencies [[org.clojure/clojure "1.7.0"]
                 [org.clojure/tools.nrepl "0.2.10"]
                 [org.clojure/core.match "0.3.0-alpha4"]
                 [mlabs.jars/clojure-contrib "1.2.0-mlab"]
                 [clojure-csv/clojure-csv "2.0.1"]
                 [org.clojure/math.numeric-tower "0.0.4"]
                 [me.raynes/fs "1.4.4"]
                 [me.raynes/conch "0.8.0"]
                 
                 [instaparse "1.4.1"]
                 [net.mikera/core.matrix "0.28.0"]
                 
                 [org.clojure/java.jdbc "0.3.3"]
                 [mysql/mysql-connector-java "5.1.18"]
                                  
                 [cascalog "2.0.0"]
                 [clj-shell "0.1.0"]
                 [org.apache.hadoop/hadoop-core "1.1.2"]

                 ;;for net-eval before jaring it
                 [slingshot "0.10.3"]

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
