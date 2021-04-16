#Random forest analysis R code to obtain best SNPs based on the phenotypic traits.
#Cisel Kemahli Aytekin, PhD [mckemahli@gmail.com]

library(RRF)
#library(assigner) #https://github.com/eriqande/assigner
library(randomForest)
library(ggplot2)
library(ggthemes)
library("WGCNA")
#-----------

# Read in data files. Geno file consists of the genotypes of each SNP. The data also includes trait of each indvidual that you want to identify for random forest. 
# You can mangae your data based on your geno file.

data_btr <- read.table("btr_rf.geno", header=T)
rownames(data_btr)=data_btr[,1]
data_btr$Site <- NULL
data_btr$Major <- NULL
data_btr$Minor <- NULL
data_btr2=transposeBigData(data_btr)

#--------Testing/training sets. Randomly select 2/3 of samples for the training set, and put the remaining in the test set

samp <- sample(nrow(data_btr2), 0.66 * nrow(data_btr2))
train <- data_btr2[samp, ]
test <- data_btr2[-samp, ]

#--------------------

# ran RF using 250, 500, 1,000, 2,000, 4,000, 8,000 and 10,000 trees, 10 times each
# minimum node size of 5. Other parameters default

# ntree = number of trees
# mtry - Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3)
# nodesize = Minimum size of terminal nodes. 

n_btr=ncol(train)-1
ntrees <- rep(c(250,500,1000,2000,4000,8000,10000),each=10)
oob <- rep(0,times=70)
ntree_oob <- data.frame(ntrees,oob)

# Run 7 different tree numbers 10 times each, and record OOB error for each run
for (i in 1:10) {
  print(paste0("ntree = ", ntree_oob$ntrees[i], ", replicate ", i))
  rf_125tree=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=10000, nodesize=5, importance=TRUE)
  ntree_oob$oob[i] = as.numeric(rf_125tree$err.rate[125,1])
  
  print(paste0("ntree = ", ntree_oob$ntrees[10+i], ", replicate ", i))
  rf_250tree=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=250, nodesize=5, importance=TRUE)
  ntree_oob$oob[10+i] = as.numeric(rf_250tree$err.rate[250,1])
  
  print(paste0("ntree = ", ntree_oob$ntrees[20+i], ", replicate ", i))
  rf_500tree=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=500, nodesize=5, importance=TRUE)
  ntree_oob$oob[20+i] = as.numeric(rf_500tree$err.rate[500,1])
  
  print(paste0("ntree = ", ntree_oob$ntrees[30+i], ", replicate ", i))
  rf_1000tree=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=1000, nodesize=5, importance=TRUE)
  ntree_oob$oob[30+i] = as.numeric(rf_1000tree$err.rate[1000,1])
  
  print(paste0("ntree = ", ntree_oob$ntrees[40+i], ", replicate ", i))
  rf_2000tree=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=2000, nodesize=5, importance=TRUE)
  ntree_oob$oob[40+i] = as.numeric(rf_2000tree$err.rate[2000,1])
  
  print(paste0("ntree = ", ntree_oob$ntrees[50+i], ", replicate ", i))
  rf_4000tree=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=4000, nodesize=5, importance=TRUE)
  ntree_oob$oob[50+i] = as.numeric(rf_4000tree$err.rate[4000,1])
  
  print(paste0("ntree = ", ntree_oob$ntrees[60+i], ", replicate ", i))
  rf_8000tree=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=8000, nodesize=5, importance=TRUE)
  ntree_oob$oob[60+i] = as.numeric(rf_8000tree$err.rate[8000,1])
  
}

# go with number of trees with stabilized out of bag error

ggplot(ntree_oob,aes(group=ntrees, x=as.factor(ntrees), y=oob)) + geom_boxplot() + theme_few() #plot OOB error over all runs

#Choose number of trees to proceed with
chosenNTree <- 2000

# mtry tested at default, half default, twice default. go with lowest error. Try mtry starting from 2.
# Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3)

mtry <- rep(c(2,4,8,16),each=10)
oob <- rep(0,times=40)
mtry_oob <- data.frame(mtry,oob)


for (i in 1:10) {
  print(paste0("mry = ", mtry_oob$mtry[i], ", replicate ", i))
  rf_2mtry=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=chosenNTree, nodesize=5, importance=TRUE, mtry=2)
  mtry_oob$oob[i] = as.numeric(rf_2mtry$err.rate[chosenNTree,1])
  
  print(paste0("mry = ", mtry_oob$mtry[10+i], ", replicate ", i))
  rf_4mtry=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=chosenNTree, nodesize=5, importance=TRUE, mtry=4)
  mtry_oob$oob[10+i] = as.numeric(rf_4mtry$err.rate[chosenNTree,1])
  
  print(paste0("mry = ", mtry_oob$mtry[20+i], ", replicate ", i))
  rf_8mtry=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=chosenNTree, nodesize=5, importance=TRUE, mtry=8)
  mtry_oob$oob[20+i] = as.numeric(rf_8mtry$err.rate[chosenNTree,1])
  
  print(paste0("mry = ", mtry_oob$mtry[30+i], ", replicate ", i))
  rf_16mtry=randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=chosenNTree, nodesize=5, importance=TRUE, mtry=16)
  mtry_oob$oob[30+i] = as.numeric(rf_16mtry$err.rate[chosenNTree,1])
  
}

ggplot(mtry_oob,aes(group=mtry, x=as.factor(mtry), y=oob)) + geom_boxplot() + theme_few() #plot OOB error over all runs


#run the constructed model on the testing set to classify individuals
pred <- predict(rf_8mtry, newdata = test)
table(pred, test$Behavior)

# 5 separate runs to make 5 lists, ranked by MDA. Panels created from lists by choosing 10 different MDA thresholds (all positive), to make 40-700 loci/panel
# For example, SNPs consistently ranked within the top 800 loci in all five lists were aggregated to form a consensus panel of 67 SNPs

rf_mtry_r1 <- randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=chosenNTree, nodesize=5, importance=TRUE, mtry=chosenmtry)
rf_mtry_r2 <- randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=chosenNTree, nodesize=5, importance=TRUE, mtry=chosenmtry)
rf_mtry_r3 <- randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=chosenNTree, nodesize=5, importance=TRUE, mtry=chosenmtry)
rf_mtry_r4 <- randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=chosenNTree, nodesize=5, importance=TRUE, mtry=chosenmtry)
rf_mtry_r5 <- randomForest(train[,1:n_btr], train$Behavior, data=train, proximity=TRUE, ntree=chosenNTree, nodesize=5, importance=TRUE, mtry=chosenmtry)

print(rf_mtry_r1)
print(rf_mtry_r2)
print(rf_mtry_r3)
print(rf_mtry_r4)
print(rf_mtry_r5)

par(mfcol=c(3,2))
importance(rf_mtry_r1, type=1)
importance(rf_mtry_r2, type=1)
importance(rf_mtry_r3, type=1)
importance(rf_mtry_r4, type=1)
importance(rf_mtry_r5, type=1)

#write.table(importance(rf_mtry_r1, type=1), "BTR_rf1_1.txt", quote=F, row.names=T)
#write.table(importance(rf_18mtry_r2, type=1), "BTR_rf2_1.txt", quote=F, row.names=T)
#write.table(importance(rf_18mtry_r3, type=1), "BTR_rf3_1.txt", quote=F, row.names=T)
#write.table(importance(rf_18mtry_r4, type=1), "BTR_rf4_1.txt", quote=F, row.names=T)
#write.table(importance(rf_18mtry_r5, type=1), "BTR_rf5_1.txt", quote=F, row.names=T)

#After determining best ntree and mtry values for the dataset, run 5 separate randomForest. After this, based on the OOB error rates, the SNPs that have all 0>= values were eliminated. 
#Until oob error rate is stabilized and minimized for all runs, this step can be repeated. At the end, the best SNPs that identifes the searched pattern will be obtained. 
