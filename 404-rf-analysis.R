#! /usr/local/bin/RScript
# Backward purging and selection of variants

library(randomForest)
# Based on this table and plot, the best 5% of loci have the lowest error rate
# As a conservative measure, I'll run backward purging RF with the best 10% loci

# Mac notes that all his loci have about the same error rate (0.0)

#Loading and creating necessary objects

importance_rf_all <-read.csv(file = "outputs/400/rf_importance_values_all_loci_classification_tutorial.csv")
load("outputs/400/class_data_corrected.rda")
load("outputs/400/sample_size.rda")

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

names_best_10perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.90))]
names_best_10perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.90))]
names_best_10perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.90))]
names_best_10perc_unique<-unique(c(names_best_10perc_1,names_best_10perc_2,names_best_10perc_3))

#################### Backward purging approach
names_purging <- names_best_10perc_unique

genotypes_purging<-class_data_corrected[,colnames(class_data_corrected) %in% names_purging]

#rf_purging_1 = randomForest(x=na.roughfix(genotypes_purging), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_purging_1,file="outputs/400/rf_purging_1.Rdata")

#rf_purging_2 = randomForest(x=na.roughfix(genotypes_purging), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_purging_2,file="outputs/400/rf_purging_2.Rdata")

#rf_purging_3 = randomForest(x=na.roughfix(genotypes_purging), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
#save(rf_purging_3,file="outputs/400/rf_purging_3.Rdata")

load("outputs/400/rf_purging_1.Rdata")
load("outputs/400/rf_purging_2.Rdata")
load("outputs/400/rf_purging_3.Rdata")

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
    #table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
   # Here there is an error: Error in sample(x = dont_keep, n = 1) : unused argument (n = 1)
    # replace size = 1, 01/23/2020
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep, size=1)]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[length(names_purging)-i]]<-names_keep
  genotypes_purging<-class_data_corrected[,colnames(class_data_corrected) %in% names_keep]
  rf_purging_1 = randomForest(x=na.roughfix(genotypes_purging), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
  rf_purging_2 = randomForest(x=na.roughfix(genotypes_purging), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
  rf_purging_3 = randomForest(x=na.roughfix(genotypes_purging), y = as.factor(class_data_corrected$resistance), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000, strata=as.factor(class_data_corrected$resistance), sampsize=sample_size)
  error_rate_best[length(names_purging)-i,] <- c(rf_purging_1$err.rate[25000],rf_purging_2$err.rate[25000],rf_purging_3$err.rate[25000])
}


error_rate_best$Average<-apply(error_rate_best,1,mean)
write.csv(error_rate_best, file="outputs/400/Backward_purging_OOB-ER_classification_tutorial.csv") # Save the error rates

pdf("outputs/400/error-rates-line-307.pdf")
# Now plot the backward purging results. Omit error rates from one locus since RF cannot be conducted with just one locus
plot(seq(2,nrow(error_rate_best),1),error_rate_best$Average[-c(1)],xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5,pch=16)

# Which group of loci yields the lowest error rate?
which(error_rate_best$Average==min(error_rate_best$Average[-c(1)])) #34 loci have the lowest OOB-ER

dev.off()

num<-length(which(error_rate_best$Average==min(error_rate_best$Average[-c(1)])))

# Export the names of the predictor loci
# I made this 100 so that it can be trimmed as necessary based on the output of "which(error_rate_best$Average==min(error_rate_best$Average[-c(1)]))"
write.csv(names_all_iterations[[100]],file="outputs/400/Predictor_loci_classification_tutorial.csv")

write.csv(names_all_iterations[[num]],file="outputs/400/Predictor_loci_classification_tutorial_num.csv")