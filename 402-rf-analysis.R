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
save(sample_size, file="outputs/400/sample_size.rda")
#renaming my objects to match
class_data_corrected <- df
class_data_corrected<-class_data_corrected %>% rename(resistance=Phenotype)

save(class_data_corrected, file="outputs/400/class_data_corrected.rda")
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

# Mac 01232020 -> my importance was 0.93325

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
write.csv(importance_rf_all,file="outputs/400/rf_importance_values_all_loci_classification_tutorial.csv",row.names=FALSE)


