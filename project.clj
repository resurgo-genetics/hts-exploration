(defproject hts-exploration "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}

  :jvm-opts ["-Xmx4g"]
  
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [mlabs.jars/clojure-contrib "1.2.0-mlab"]
                 [clojure-csv/clojure-csv "2.0.1"]
                 [org.clojure/math.numeric-tower "0.0.4"]
                 [me.raynes/fs "1.4.4"]

                 [org.clojure/java.jdbc "0.3.3"]
                 [mysql/mysql-connector-java "5.1.18"]
                                  
                 [cascalog "2.0.0"]
                 [clj-shell "0.1.0"]
                 [org.apache.hadoop/hadoop-core "1.1.2"]

                 ;;for net-eval before jaring it
                 [org.clojure/tools.nrepl "0.2.2"]
                 [slingshot "0.10.3"]

                 ;;file reading in parallel
                 [iota "1.1.1"]
                 [foldable-seq "0.2"]])
