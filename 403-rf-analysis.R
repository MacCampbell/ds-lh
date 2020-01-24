#! /usr/local/bin/RScript

#Third part of random forest analysis
library(tidyverse)
library(randomForest)

# Retrive outputs from previous script
importance_rf_all <-read.csv(file = "outputs/400/rf_importance_values_all_loci_classification_tutorial.csv")
load("outputs/400/class_data_corrected.rda")
load("outputs/400/sample_size.rda")

#And recomputing importance
load("outputs/400/rf_all_1.Rdata")
load("outputs/400/rf_all_2.Rdata")
load("outputs/400/rf_all_3.Rdata")

rf_all_1_err.rate <- rf_all_1$err.rate[25000]
rf_all_2_err.rate <- rf_all_2$err.rate[25000]
rf_all_3_err.rate <- rf_all_3$err.rate[25000]

importance_rf_all_1<-data.frame(importance(rf_all_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_all_1)<-c("importance")

importance_rf_all_2<-data.frame(importance(rf_all_2,type=1))
colnames(importance_rf_all_2)<-c("importance")

cor(importance_rf_all_1,importance_rf_all_2) # A correlation of 0.98 for locus importance values between forests is extremely good, so we'll use 25,000 trees for the remaining forests

# Mac 01232020 -> my importance was 0.93325

importance_rf_all_3<-data.frame(importance(rf_all_3,type=1))
colnames(importance_rf_all_3)<-c("importance")

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
rf_2perc_1 = randomForest(x = na.roughfix(genotypes_2perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_2perc_1,file="outputs/400/rf_2perc_1.Rdata")

rf_2perc_2 = randomForest(x = na.roughfix(genotypes_2perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_2perc_2,file="outputs/400/rf_2perc_2.Rdata")

rf_2perc_3 = randomForest(x = na.roughfix(genotypes_2perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_2perc_3,file="outputs/400/rf_2perc_3.Rdata")

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
rf_3perc_1 = randomForest(x = na.roughfix(genotypes_3perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_3perc_1,file="outputs/400/rf_3perc_1.Rdata")

rf_3perc_2 = randomForest(x = na.roughfix(genotypes_3perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_3perc_2,file="outputs/400/rf_3perc_2.Rdata")

rf_3perc_3 = randomForest(x = na.roughfix(genotypes_3perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_3perc_3,file="outputs/400/rf_3perc_3.Rdata")

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
rf_4perc_1 = randomForest(x = na.roughfix(genotypes_4perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_4perc_1,file="outputs/400/rf_4perc_1.Rdata")

rf_4perc_2 = randomForest(x = na.roughfix(genotypes_4perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_4perc_2,file="outputs/400/rf_4perc_2.Rdata")

rf_4perc_3 = randomForest(x = na.roughfix(genotypes_4perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_4perc_3,file="outputs/400/rf_4perc_3.Rdata")

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
rf_5perc_1 = randomForest(x = na.roughfix(genotypes_5perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_5perc_1,file="outputs/400/rf_5perc_1.Rdata")

rf_5perc_2 = randomForest(x = na.roughfix(genotypes_5perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_5perc_2,file="outputs/400/rf_5perc_2.Rdata")

rf_5perc_3 = randomForest(x = na.roughfix(genotypes_5perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
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
rf_10perc_1 = randomForest(x = na.roughfix(genotypes_10perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_10perc_1,file="outputs/400/rf_10perc_1.Rdata")

rf_10perc_2 = randomForest(x = na.roughfix(genotypes_10perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_10perc_2,file="outputs/400/rf_10perc_2.Rdata")

rf_10perc_3 = randomForest(x = na.roughfix(genotypes_10perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_10perc_3,file="outputs/400/rf_10perc_3.Rdata")

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
rf_20perc_1 = randomForest(x = na.roughfix(genotypes_20perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_20perc_1,file="outputs/400/rf_20perc_1.Rdata")

rf_20perc_2 = randomForest(x = na.roughfix(genotypes_20perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_20perc_2,file="outputs/400/rf_20perc_2.Rdata")

rf_20perc_3 = randomForest(x = na.roughfix(genotypes_20perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_20perc_3,file="outputs/400/rf_20perc_3.Rdata")

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
rf_30perc_1 = randomForest(x = na.roughfix(genotypes_30perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_30perc_1,file="outputs/400/rf_30perc_1.Rdata")

rf_30perc_2 = randomForest(x = na.roughfix(genotypes_30perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_30perc_2,file="outputs/400/rf_30perc_2.Rdata")

rf_30perc_3 = randomForest(x = na.roughfix(genotypes_30perc), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
save(rf_30perc_3,file="outputs/400/rf_30perc_3.Rdata")

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
write.csv(All_initial_err.rate,file="outputs/400/All_initial_err_rate_classification_tutorial.csv")

pdf("outputs/400/error-rates.pdf")
# Plot error rates as well
par(mar=c(5,6,3,3))
plot(All_initial_err.rate$Number_loci,All_initial_err.rate$Average,log="x", pch=19,xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5)
dev.off()
