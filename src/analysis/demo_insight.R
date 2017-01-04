library(ggplot2)

setwd("~/bin/hts-exploration/data/150211/")

d <- read.csv("large-cluster-kmers-demo.csv")
d2 <- d[, which(!apply(d,2,function(x) { all(x==0) }))] #remove cols of all 0
d2 <- subset(d2, type != "ident")
p <- prcomp(d2[,c(-1,-2)], center=T, scale.=T)

color_values <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c") #maintain same colors for categories
names(color_values) <- c("doubleton", "tested", "ident", "large")


##Are tested sequences different or similar to other untested sequences?

##pc1 vs pc2
ggplot(data.frame(name=d2[,1], type=d2[,2],x=p$x[,1], y=p$x[,2]), aes(x=x, y=y, color=type))+
    geom_point()+theme_bw()+
    geom_text(aes(label=name))+
    xlab("PC1")+ylab("PC2")+
    scale_color_manual(name="Filter", values=color_values, 
                       breaks=c("doubleton","tested","ident","large"), 
                       labels=c("Doubleton","Tested", "Identity+Size", "Large"))+
    theme(legend.position="right",
          legend.text = element_text(size=18))

##are any kmer vectors similar to each other?
##if any vectors are similar to known binders then they could be candidates for further testing
rownames(d2) <- apply(d2[,c(1,2)],1,function(x) { paste(x, collapse="_") })
plot(hclust(dist(scale(d2[,c(-1,-2)], center=T, scale=T))))

