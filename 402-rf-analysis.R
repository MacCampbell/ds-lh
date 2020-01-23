#! /usr/local/bin/RScript
#This is the second bit of a Random Forest analysis, basically using the tutorial form.

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

#renaming my objects to match
class_data_corrected <- df
class_data_corrected<-class_data_corrected %>% rename(resistance=Phenotype)
###########################################################################################################################################
###########################################################################################################################################
# Now begin the full Random Forest analyses

# Recall that even though we optimized mtry, we must now run a larger number of trees in order to achieve convergence of importance values between forests.
# As a starting point, we will grow 25,000 trees and increase if necessary. We do not need to worry about this increase in ntree affecting our mtry optimization,
# since the OOB-ER reached a plateau for a given mtry value after about 400 trees.

rf_all_1 = randomForest(x = na.roughfix(class_data_corrected[,3:5135]), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=1027, ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_all_1,file="outputs/400/rf_all_1.Rdata")

rf_all_2 = randomForest(x = na.roughfix(class_data_corrected[,3:5135]), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=1027, ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_all_2,file="outputs/400/rf_all_2.Rdata")

#Check correlation of locus importance values between forests 
importance_rf_all_1<-data.frame(importance(rf_all_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_all_1)<-c("importance")
importance_rf_all_2<-data.frame(importance(rf_all_2,type=1))
colnames(importance_rf_all_2)<-c("importance")

cor(importance_rf_all_1,importance_rf_all_2) # A correlation of 0.98 for locus importance values between forests is extremely good, so we'll use 25,000 trees for the remaining forests

rf_all_3 = randomForest(x = na.roughfix(class_data_corrected[,3:5135]), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=1027, ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_all_3,file="outputs/400/rf_all_3.Rdata")
importance_rf_all_3<-data.frame(importance(rf_all_3,type=1))
colnames(importance_rf_all_3)<-c("importance")

############################################################################################################################################
############################################################################################################################################

# The predictive ability of classification trees is measured by the out-of-bag error rate.  An error rate is calculated for each tree within a forest. 
# We will use the error rate from the last tree in the forest, which takes all previous trees into account and thus represents the error rate after the model stabilizes/converges

rf_all_1_err.rate <- rf_all_1$err.rate[25000]
rf_all_2_err.rate <- rf_all_2$err.rate[25000]
rf_all_3_err.rate <- rf_all_3$err.rate[25000]

#Combine importance (mean decrease in accuracy) values of each locus across the three forests
importance_rf_all <-cbind(rownames(importance_rf_all_1),importance_rf_all_1,importance_rf_all_2, importance_rf_all_3)
colnames(importance_rf_all)<-c("Variable","Importance1","Importance2", "Importance3")

# Export importance values for future reference
write.csv(importance_rf_all,file="rf_importance_values_all_loci_classification_tutorial.csv",row.names=FALSE)

############################################################################################################################################
############################################################################################################################################

# Now conduct RF on subsets of the data to identify a group of loci that may be predictive of disease resistance. 
# For each subset, we will use mtry=p since that is the optimal setting that we previously found.

##### Best 2% 

names_best_2perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.98))]
names_best_2perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.98))]
names_best_2perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.98))]
names_best_2perc_unique<-unique(c(names_best_2perc_1,names_best_2perc_2,names_best_2perc_3))

# Extract genotypes 
genotypes_2perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_2perc_unique]

# Now conduct RF on this subset
rf_2perc_1 = randomForest(x = genotypes_2perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_2perc_1,file="rf_2perc_1.Rdata")

rf_2perc_2 = randomForest(x = genotypes_2perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_2perc_2,file="rf_2perc_2.Rdata")

rf_2perc_3 = randomForest(x = genotypes_2perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_2perc_3,file="rf_2perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_2perc_1_err.rate <- rf_2perc_1$err.rate[25000]
rf_2perc_2_err.rate <- rf_2perc_2$err.rate[25000]
rf_2perc_3_err.rate <- rf_2perc_3$err.rate[25000]

rm(rf_2perc_1,rf_2perc_2,rf_2perc_3) # remove the objects to save memory in R 


##### Best 3% 

names_best_3perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.97))]
names_best_3perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.97))]
names_best_3perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.97))]
names_best_3perc_unique<-unique(c(names_best_3perc_1,names_best_3perc_2,names_best_3perc_3))

# Extract genotypes
genotypes_3perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_3perc_unique]

# Now conduct RF on this subset
rf_3perc_1 = randomForest(x = genotypes_3perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_3perc_1,file="rf_3perc_1.Rdata")

rf_3perc_2 = randomForest(x = genotypes_3perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_3perc_2,file="rf_3perc_2.Rdata")

rf_3perc_3 = randomForest(x = genotypes_3perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_3perc_3,file="rf_3perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_3perc_1_err.rate <- rf_3perc_1$err.rate[25000]
rf_3perc_2_err.rate <- rf_3perc_2$err.rate[25000]
rf_3perc_3_err.rate <- rf_3perc_3$err.rate[25000]

rm(rf_3perc_1,rf_3perc_2,rf_3perc_3) # remove the objects to save memory in R 

##### Best 4% 

names_best_4perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.96))]
names_best_4perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.96))]
names_best_4perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.96))]
names_best_4perc_unique<-unique(c(names_best_4perc_1,names_best_4perc_2,names_best_4perc_3))

# Extract genotypes
genotypes_4perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_4perc_unique]

# Now conduct RF on this subset
rf_4perc_1 = randomForest(x = genotypes_4perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_4perc_1,file="outputs/400/rf_4perc_1.Rdata")

rf_4perc_2 = randomForest(x = genotypes_4perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_4perc_2,file="outputs/400/rf_4perc_2.Rdata")

rf_4perc_3 = randomForest(x = genotypes_4perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_4perc_3,file="ouputs/400/rf_4perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_4perc_1_err.rate <- rf_4perc_1$err.rate[25000]
rf_4perc_2_err.rate <- rf_4perc_2$err.rate[25000]
rf_4perc_3_err.rate <- rf_4perc_3$err.rate[25000]

rm(rf_4perc_1,rf_4perc_2,rf_4perc_3) # remove the objects to save memory in R 

##### Best 5% 

names_best_5perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.95))]
names_best_5perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.95))]
names_best_5perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.95))]
names_best_5perc_unique<-unique(c(names_best_5perc_1,names_best_5perc_2,names_best_5perc_3))

# Extract genotypes
genotypes_5perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_5perc_unique]

# Now conduct RF on this subset
rf_5perc_1 = randomForest(x = genotypes_5perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_5perc_1,file="outputs/400/rf_5perc_1.Rdata")

rf_5perc_2 = randomForest(x = genotypes_5perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_5perc_2,file="outputs/400/rf_5perc_2.Rdata")

rf_5perc_3 = randomForest(x = genotypes_5perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_5perc_3,file="outputs/400/rf_5perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_5perc_1_err.rate <- rf_5perc_1$err.rate[25000]
rf_5perc_2_err.rate <- rf_5perc_2$err.rate[25000]
rf_5perc_3_err.rate <- rf_5perc_3$err.rate[25000]

rm(rf_5perc_1,rf_5perc_2,rf_5perc_3) # remove the objects to save memory in R 

##### Best 10% 

names_best_10perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.90))]
names_best_10perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.90))]
names_best_10perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.90))]
names_best_10perc_unique<-unique(c(names_best_10perc_1,names_best_10perc_2,names_best_10perc_3))

# Extract genotypes
genotypes_10perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_10perc_unique]

# Now conduct RF on this subset
rf_10perc_1 = randomForest(x = genotypes_10perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_10perc_1,file="rf_10perc_1.Rdata")

rf_10perc_2 = randomForest(x = genotypes_10perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_10perc_2,file="rf_10perc_2.Rdata")

rf_10perc_3 = randomForest(x = genotypes_10perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_10perc_3,file="rf_10perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_10perc_1_err.rate <- rf_10perc_1$err.rate[25000]
rf_10perc_2_err.rate <- rf_10perc_2$err.rate[25000]
rf_10perc_3_err.rate <- rf_10perc_3$err.rate[25000]

rm(rf_10perc_1,rf_10perc_2,rf_10perc_3) # remove the objects to save memory in R 

##### Best 20% 

names_best_20perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.80))]
names_best_20perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.80))]
names_best_20perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.80))]
names_best_20perc_unique<-unique(c(names_best_20perc_1,names_best_20perc_2,names_best_20perc_3))

# Extract genotypes
genotypes_20perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_20perc_unique]

# Now conduct RF on this subset
rf_20perc_1 = randomForest(x = genotypes_20perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_20perc_1,file="rf_20perc_1.Rdata")

rf_20perc_2 = randomForest(x = genotypes_20perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_20perc_2,file="rf_20perc_2.Rdata")

rf_20perc_3 = randomForest(x = genotypes_20perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_20perc_3,file="rf_20perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_20perc_1_err.rate <- rf_20perc_1$err.rate[25000]
rf_20perc_2_err.rate <- rf_20perc_2$err.rate[25000]
rf_20perc_3_err.rate <- rf_20perc_3$err.rate[25000]

rm(rf_20perc_1,rf_20perc_2,rf_20perc_3) # remove the objects to save memory in R 

##### Best 30% 

names_best_30perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.70))]
names_best_30perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.70))]
names_best_30perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.70))]
names_best_30perc_unique<-unique(c(names_best_30perc_1,names_best_30perc_2,names_best_30perc_3))

# Extract genotypes
genotypes_30perc<-class_data_corrected[,colnames(class_data_corrected) %in% names_best_30perc_unique]

# Now conduct RF on this subset
rf_30perc_1 = randomForest(x = genotypes_30perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_30perc_1,file="rf_30perc_1.Rdata")

rf_30perc_2 = randomForest(x = genotypes_30perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_30perc_2,file="rf_30perc_2.Rdata")

rf_30perc_3 = randomForest(x = genotypes_30perc, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_30perc_3,file="rf_30perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_30perc_1_err.rate <- rf_30perc_1$err.rate[25000]
rf_30perc_2_err.rate <- rf_30perc_2$err.rate[25000]
rf_30perc_3_err.rate <- rf_30perc_3$err.rate[25000]

rm(rf_30perc_1,rf_30perc_2,rf_30perc_3) # remove the objects to save memory in R 


# Now combine all of the error rates from the subsets and for all loci to identify a group for the backward purging approach

All_initial_err.rate <- rbind(cbind(rf_all_1_err.rate,rf_all_2_err.rate,rf_all_3_err.rate),
                              cbind(rf_2perc_1_err.rate,rf_2perc_2_err.rate,rf_2perc_3_err.rate),
                              cbind(rf_3perc_1_err.rate,rf_3perc_2_err.rate,rf_3perc_3_err.rate),
                              cbind(rf_4perc_1_err.rate,rf_4perc_2_err.rate,rf_4perc_3_err.rate),
                              cbind(rf_5perc_1_err.rate,rf_5perc_2_err.rate,rf_5perc_3_err.rate),
                              cbind(rf_10perc_1_err.rate,rf_10perc_2_err.rate,rf_10perc_3_err.rate),
                              cbind(rf_20perc_1_err.rate,rf_20perc_2_err.rate,rf_20perc_3_err.rate),
                              cbind(rf_30perc_1_err.rate,rf_30perc_2_err.rate,rf_30perc_3_err.rate))

# Plot error rates for the various subsets
All_initial_err.rate<-data.frame(All_initial_err.rate)
All_initial_err.rate$Number_loci<-c(1000,length(names_best_2perc_unique),length(names_best_3perc_unique),length(names_best_4perc_unique),length(names_best_5perc_unique),length(names_best_10perc_unique),length(names_best_20perc_unique),length(names_best_30perc_unique))
rownames(All_initial_err.rate)<-c("All","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
All_initial_err.rate$Average<-apply(All_initial_err.rate[,1:3],1,mean)

# Write error rates to file for future reference
write.csv(All_initial_err.rate,file="All_initial_err_rate_classification_tutorial.csv")

pdf("outputs/400/error-rates.pdf")
# Plot error rates as well
par(mar=c(5,6,3,3))
plot(All_initial_err.rate$Number_loci,All_initial_err.rate$Average,log="x", pch=19,xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5)
dev.off()
# Based on this table and plot, the best 5% of loci have the lowest error rate
# As a conservative measure, I'll run backward purging RF with the best 10% loci

#################### Backward purging approach
names_purging <- names_best_10perc_unique

genotypes_purging<-class_data_corrected[,colnames(class_data_corrected) %in% names_purging]

rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_purging_1,file="rf_purging_1.Rdata")
rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_purging_2,file="rf_purging_2.Rdata")
rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_purging_3,file="rf_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[length(names_purging)]]<-names_purging
error_rate_best<-data.frame(V1=1:length(names_purging),V2=1:length(names_purging),V3=1:length(names_purging))
rownames(error_rate_best)<-1:length(names_purging)
error_rate_best[length(names_purging),] <- c(rf_purging_1$err.rate[25000],rf_purging_2$err.rate[25000],rf_purging_3$err.rate[25000])


for (i in 1:(length(names_purging)-2)){  # RF cannot be conducted with 1 locus, which is why the loop is from 1:length(names_purging)-2
  print(i)
  imp_purging_1<-data.frame(importance(rf_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging)
  colnames(imp_purging_1)<-"Mean_Decrease_Accuracy1"
  rownames(imp_purging_2)<-colnames(genotypes_purging)
  colnames(imp_purging_2)<-"Mean_Decrease_Accuracy2"
  rownames(imp_purging_3)<-colnames(genotypes_purging)
  colnames(imp_purging_3)<-"Mean_Decrease_Accuracy3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[length(names_purging)-i]]<-names_keep
  genotypes_purging<-class_data_corrected[,colnames(class_data_corrected) %in% names_keep]
  rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
  rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
  rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
  error_rate_best[length(names_purging)-i,] <- c(rf_purging_1$err.rate[25000],rf_purging_2$err.rate[25000],rf_purging_3$err.rate[25000])
}


error_rate_best$Average<-apply(error_rate_best,1,mean)
write.csv(error_rate_best, file="Backward_purging_OOB-ER_classification_tutorial.csv") # Save the error rates

# Now plot the backward purging results. Omit error rates from one locus since RF cannot be conducted with just one locus
plot(seq(2,nrow(error_rate_best),1),error_rate_best$Average[-c(1)],xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5,pch=16)

# Which group of loci yields the lowest error rate?
which(error_rate_best$Average==min(error_rate_best$Average[-c(1)])) #34 loci have the lowest OOB-ER

# Export the names of the predictor loci
write.csv(names_all_iterations[[34]],file="Predictor_loci_classification_tutorial.csv")
