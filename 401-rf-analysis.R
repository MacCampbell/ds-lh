#!  /usr/local/bin/RScript

# This is the first bit of random forest analysis. We'll 
# 1, get genotypes and replace nas.
# 2, we will choose the mtry value that minimizes the OOB-ER.
library(tidyverse)
library(randomForest)
genos<-read_tsv("outputs/300/202.geno.gz", col_names = FALSE)

#Removing last column and transposing
genos<-genos %>% select(-X205)
trans<-as_tibble(t(genos))
trans<-trans[-1:-2,]
names<-read_tsv("bamlists/202.bamlist",col_names=FALSE) %>% mutate(Names=gsub("bams/","",X1))
names <- names %>% mutate(Sample=gsub(".sort-n.fixmate-m.sort.markdup-r.bam","",Names)) %>% select(Sample)

phenos <- read_tsv("phenos/202.phenos", col_names=FALSE) %>% rename(Pop=X1) %>% mutate(Phenotype=Pop)

meta<-as_tibble(cbind(names,phenos))

input<-as_tibble(cbind(meta,trans))

#Using df for generalization
df<-input
df[df < 0] <- NA

length(which(input$Phenotype==0)) # 61 here today
length(which(input$Phenotype==1)) # 141 here today
sample_size <- c(61,61)

results_optimization <- matrix(data=NA , nrow = 0, ncol = 3)
for (i in seq(from = 100, to = 1000 , by = 100)){  # values of ntree
  print(i)
  for (j in c(71,143,513,1016,1712,5135)){    #values of mtry based on on all our loci
    ## We will run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, and p, where p is the number of loci
    
    rf_ij <- randomForest(x = na.roughfix(df[,3:5138]), y = as.factor(df$Phenotype), importance=TRUE ,proximity=TRUE, ntree=i, mtry=j, strata=as.factor(df$Phenotype), sampsize=sample_size)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$err.rate,1)[1]))
  }
}

##Now let's make some outputs

pdf("outputs/400/400-mtry.pdf")

results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","OOB_ER")

# Now plot results to see if there's a plateau
# I adjusted the y axis to be smaller
plot(results_optimization$ntree[results_optimization$mtry == 71],results_optimization$OOB_ER[results_optimization$mtry == 71], type="l", col="black", xlab="ntree",ylab="OOB-ER",ylim=c(0,.25))
lines(results_optimization$ntree[results_optimization$mtry == 143],results_optimization$OOB_ER[results_optimization$mtry == 143], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 513],results_optimization$OOB_ER[results_optimization$mtry == 513], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 1027],results_optimization$OOB_ER[results_optimization$mtry == 1016], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 1712],results_optimization$OOB_ER[results_optimization$mtry == 1712], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 5135],results_optimization$OOB_ER[results_optimization$mtry == 5135], col="red")

dev.off()

#Look for plateaus in the graphic...

#With the whole dataset, mtry of 1027 (.2p) produces a low OOB
write.csv(results_optimization, "outputs/400/results_optimization.csv")