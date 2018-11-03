setwd("D:/NB_deeplearning/autoencoder")
load("TARGET.NBL.expr.Rdata")
load("cnv.mat.Rdata")

cnv.matrix <- read.table("TARGET_CNV.all_data_by_genes.txt",header = TRUE,sep = "\t",row.names = 1)
cnv.matrix <- cnv.matrix[,c(-1,-2)]
colnames(cnv.matrix) <- gsub("\\.","-",colnames(cnv.matrix))
# load("discovery.clinical.Rdata")

library(ANN2)
# train.target.matrix <- TARGET.NBL.expr.mat[row.names(TARGET.NBL.expr.mat) %in% overlap.genes,grep("^High",TARGET.clinical[,"COG.Risk.Group"])]
gene.matrix <- TARGET.NBL.expr.mat[,grep("^High",TARGET.clinical[,"COG.Risk.Group"])]
train.gene.matrix <- gene.matrix[,match(intersect(colnames(cnv.mat),colnames(gene.matrix)),colnames(gene.matrix))]
# train.gene.matrix <- train.gene.matrix[apply(train.gene.matrix>0,1,sum)>38,]

# train.cnv.matrix <- cnv.mat[,match(intersect(colnames(cnv.mat),colnames(gene.matrix)),colnames(cnv.mat))]
# train.cnv.matrix.dup <- train.cnv.matrix[apply(train.cnv.matrix>0,1,sum)>19,]
# train.cnv.matrix.del <- train.cnv.matrix[apply(train.cnv.matrix<0,1,sum)>19,]
train.cnv.matrix <- cnv.matrix[,match(intersect(colnames(cnv.matrix),colnames(gene.matrix)),colnames(cnv.matrix))]
train.cnv.mat <- cnv.mat[,match(intersect(colnames(cnv.mat),colnames(gene.matrix)),colnames(cnv.mat))]
train.cnv.matrix.dup <- train.cnv.matrix[apply(train.cnv.mat>0,1,sum)>19,]
train.cnv.matrix.del <- train.cnv.matrix[apply(train.cnv.mat<0,1,sum)>19,]
train.target.matrix.merge <- rbind(train.gene.matrix,train.cnv.matrix.del,train.cnv.matrix.dup)

# aeNN <- autoencoder(t(train.target.matrix.merge),hiddenLayers = c(1000,500,1000), lossFunction = "log",standardize = TRUE,
#                     learnRate = 1e-06,maxEpochs = 10, batchSize = 32, momentum = 0.2,L1 = 1e-4, L2 = 1e-3)
# rX <- reconstruct(aeNN, t(train.target.matrix.merge),mahalanobis = F)

# autoencoder.object <- autoencode(X.train=train.target.matrix,nl = 5,N.hidden = c(500,100,500), unit.type = "tanh",
#                                  optim.method = c("BFGS", "L-BFGS-B", "CG"), epsilon = 0.001,lambda = 0.0002,
#                                  max.iterations = 10, beta = 6,rho = 0.01,
#                                  rescale.flag = c(F, T), rescaling.offset = 0.001)

# selected.feature <- t(encode(aeNN,t(train.target.matrix.merge)))

library(survival)
efs_time <- as.numeric(TARGET.clinical[match(colnames(train.target.matrix.merge),TARGET.clinical[,"TARGET.USI"]),"Event.Free.Survival.Time.in.Days"])
efs_event <- as.character(TARGET.clinical[match(colnames(train.target.matrix.merge),TARGET.clinical[,"TARGET.USI"]),"First.Event"])
efs_event <- as.numeric(efs_event!="Censored")
# cox_efs_pvalue <- apply(selected.feature,1,function(x){y=coxph(Surv(efs_time,efs_event)~x);y=summary(y)$sctest[3]})

os_time <- as.numeric(TARGET.clinical[match(colnames(train.target.matrix.merge),TARGET.clinical[,"TARGET.USI"]),"Overall.Survival.Time.in.Days"])
os_event <- as.character(TARGET.clinical[match(colnames(train.target.matrix.merge),TARGET.clinical[,"TARGET.USI"]),"Vital.Status"])
os_event <- as.numeric(os_event=="Dead")
# cox_pvalue <- apply(selected.feature,1,function(x){y=coxph(Surv(os_time,os_event)~x);y=summary(y)$sctest[3]})
# choose_cox_matrix <- selected.feature[cox_pvalue<0.05,]

os_all_pvalue <- apply(train.target.matrix.merge,1,function(x){y=coxph(Surv(os_time,os_event)~x);y=summary(y)$sctest[3]})
efs_all_pvalue <- apply(train.target.matrix.merge,1,function(x){y=coxph(Surv(efs_time,efs_event)~x);y=summary(y)$sctest[3]})
train.selected.matrix <- train.target.matrix.merge[efs_all_pvalue<0.05|os_all_pvalue<0.05,]

aeNN1 <- autoencoder(t(train.selected.matrix),hiddenLayers = c(500,100,500), lossFunction = "log",standardize = TRUE, 
                    learnRate = 1e-06,maxEpochs = 10, batchSize = 32, momentum = 0.2,L1 = 1e-4, L2 = 1e-3)
selected.feature <- t(encode(aeNN1,t(train.selected.matrix)))
cox_pvalue <- apply(selected.feature,1,function(x){y=coxph(Surv(os_time,os_event)~x);y=summary(y)$sctest[3]})
cox_efs_pvalue <- apply(selected.feature,1,function(x){y=coxph(Surv(efs_time,efs_event)~x);y=summary(y)$sctest[3]})

selected.surv.feature <- selected.feature[cox_pvalue<0.05|cox_efs_pvalue<0.05,]
dim(selected.surv.feature)

##################
hc <- hclust(dist(t(selected.surv.feature)))
plot(hc)

surv.kmeans <- kmeans(t(selected.surv.feature),centers=2)
center = surv.kmeans$centers
surv.cluster <- surv.kmeans$cluster
plot(t(selected.surv.feature),col=surv.cluster,pch=20)

############
library(survcomp)
fit <- coxph(Surv(os_time, os_event) ~ surv.cluster)
cindex <- concordance.index(predict(fit),surv.time = os_time, surv.event = os_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper

fit <- coxph(Surv(efs_time, efs_event) ~ surv.cluster)
cindex <- concordance.index(predict(fit),surv.time = efs_time, surv.event = efs_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper
############

############
library(cluster)
library(fpc)
sil <- silhouette(surv.cluster,dist(t(selected.surv.feature)))

surv.kmeans <- kmeans(t(selected.surv.feature),centers=2)
center = surv.kmeans$centers
surv.cluster <- surv.kmeans$cluster
summary(silhouette(surv.cluster,dist(t(selected.surv.feature))))
round(calinhara(t(selected.surv.feature),surv.cluster),digits=2)
############


##########DGE
train.overlap.gene.matrix <- train.gene.matrix[row.names(train.gene.matrix) %in% intersect(row.names(SEQC.gene.expr.high),row.names(all.target.expr.filter)),]
library(limma)
design.cluster <- sub("1","One",surv.cluster)
design.cluster <- sub("2","Two",design.cluster)
cluster.mat <- factor(design.cluster, levels = c("One","Two"))
cluster.matrix <- model.matrix(~0 + cluster.mat)
colnames(cluster.matrix) <- levels(cluster.mat)
cluster.fit <- lmFit(train.selected.overlap.matrix,cluster.matrix)
cont.cluster <- makeContrasts(One-Two,levels = cluster.matrix)
cluster.fit2 <- contrasts.fit(cluster.fit, cont.cluster)
cluster.fit2 <- eBayes(cluster.fit2)

##decidetest
diff.res <- decideTests(cluster.fit2,method = "global",adjust.method = "BH",
                        p.value = 0.01,lfc = log2(1))
summary(diff.res)

One.diff <- topTable(cluster.fit2,coef=1,n=nrow(cluster.fit2),lfc=log2(1))
One.diff.final <- One.diff[One.diff[,"adj.P.Val"]<0.05,]
nrow(One.diff.final)
up.sample.One <- One.diff.final[One.diff.final[,"logFC"] > 0.2,]
nrow(up.sample.One)
down.sample.One <- One.diff.final[One.diff.final[,"logFC"] < -0.2,]
nrow(down.sample.One)
class.label <- c(row.names(up.sample.One),
                 row.names(down.sample.One))
class.label <- c(row.names(up.sample.One)[order(up.sample.One[,1],decreasing = TRUE)][1:1000],
                 row.names(down.sample.One)[order(down.sample.One[,1])][1:1000])

# train.gene.fvalue <- apply(train.selected.matrix,1,function(x){summary(aov(surv.cluster~x))[[1]]$`F value`[1]})
# all.feature.fvalue <- apply(train.target.matrix.merge,1,function(x){summary(aov(surv.cluster~x))[[1]]$`F value`[1]})
# svm.train.matrix <- train.selected.matrix[order(train.gene.fvalue,decreasing = TRUE)[1:50],]
class.label <- c(intersect(row.names(up.sample.One),row.names(svm.train.matrix)),
                 intersect(row.names(down.sample.One),row.names(svm.train.matrix)))
length(class.label)
length(intersect(row.names(up.sample.One),row.names(svm.train.matrix)))
length(intersect(row.names(down.sample.One),row.names(svm.train.matrix)))

##########################xgboost
library(xgboost)
dtrain <- xgb.DMatrix(t(train.overlap.gene.matrix), label = as.numeric(surv.cluster==1))
param <- list(max_depth = 8, eta = 0.1, silent = 1, nthread = 2,
              objective = "binary:logistic", eval_metric = "auc",
              lambda = 1e-4,alpha = 1e-9)

bst <- xgb.train(param,dtrain,nround = 30)
# bst <- xgboost(data = t(train.overlap.gene.matrix),label = as.numeric(surv.cluster==1),max_depth = 15, eta = 0.1, silent = 0, 
#                objective = "binary:logistic",nround = 100)
pred <- as.numeric(predict(bst,t(train.overlap.gene.matrix))<0.5) + 1
table(pred,surv.cluster)
# names <- row.names(train.overlap.gene.matrix)
# importance.matrix <- xgb.importance(names,model = bst)
# dim(importance.matrix)
bst <- xgb.cv(param,dtrain,nround = 20,nfold = 10)

##############################
#random Forest
library(randomForest)
tumor.rf <- data.frame(cbind(t(train.gene.matrix),surv.cluster))
tumor.rf$surv.cluster[tumor.rf$surv.cluster==1] <- "One"
tumor.rf$surv.cluster[tumor.rf$surv.cluster==2] <- "Two"
tumor.rf$surv.cluster <- factor(tumor.rf$surv.cluster)
n<-length(names(tumor.rf))
Error<-NULL
set.seed(123)
for(i in 1:(n-1)){
  rf<-randomForest(surv.cluster ~ .,data = tumor.rf,mtry=i,importance=TRUE,proximity=TRUE)
  err<-mean(rf$err.rate)
  Error[i]<-err
}
m=which.min(Error)
rf.res <- randomForest(surv.cluster ~ .,data = tumor.rf,mtry=10,importance=TRUE,proximity=TRUE)
plot(rf.res)
rf.res <- randomForest(surv.cluster ~ .,data = tumor.rf,mtry=10,ntree=100,importance=TRUE,proximity=TRUE)

##############################

MYCN <- TARGET.clinical[match(colnames(train.target.matrix.merge),TARGET.clinical[,"TARGET.USI"]),"MYCN.status"]
Gender <- TARGET.clinical[match(colnames(train.target.matrix.merge),TARGET.clinical[,"TARGET.USI"]),"Gender"]
Race <- TARGET.clinical[match(colnames(train.target.matrix.merge),TARGET.clinical[,"TARGET.USI"]),"Race"]
Stage <- TARGET.clinical[match(colnames(train.target.matrix.merge),TARGET.clinical[,"TARGET.USI"]),"INSS.Stage"]
Histology <- TARGET.clinical[match(colnames(train.target.matrix.merge),TARGET.clinical[,"TARGET.USI"]),"Histology"]

fisher.test(MYCN,surv.cluster)
fisher.test(Gender,surv.cluster)
fisher.test(Race,surv.cluster)
fisher.test(Stage,surv.cluster)
fisher.test(Histology,surv.cluster)

#####
library(GGally)
library(ggplot2)
library(gtable)
library(survival)
sdf <- survdiff(Surv(efs_time, efs_event) ~ surv.cluster)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
sdf <- survdiff(Surv(os_time, os_event) ~ surv.cluster)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
kmsurvival <- survfit(Surv(efs_time,efs_event) ~ surv.cluster)
ggsurv(kmsurvival,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="EFS p-value 2.20e-7")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))
summary(coxph(Surv(efs_time, efs_event) ~ surv.cluster))$concordance

kmsurvival <- survfit(Surv(os_time,os_event) ~ surv.cluster)
ggsurv(kmsurvival,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="OS p-value 2.80e-8")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
summary(coxph(Surv(os_time, os_event) ~ surv.cluster))$concordance

############
library(survcomp)
fit <- coxph(Surv(os_time, os_event) ~ surv.cluster)
cindex <- concordance.index(predict(fit),surv.time = os_time, surv.event = os_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper

fit <- coxph(Surv(efs_time, efs_event) ~ surv.cluster)
cindex <- concordance.index(predict(fit),surv.time = efs_time, surv.event = efs_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper
############

####################################PCA
p <- prcomp(as.matrix(train.selected.matrix))
pca.select.feature <- t(p$rotation[,1:100])
pca.cox_pvalue <- apply(pca.select.feature,1,function(x){y=coxph(Surv(os_time,os_event)~x);y=summary(y)$sctest[3]})
pca.cox_efs_pvalue <- apply(pca.select.feature,1,function(x){y=coxph(Surv(efs_time,efs_event)~x);y=summary(y)$sctest[3]})

pca.selected.surv.feature <- pca.select.feature[pca.cox_pvalue<0.05|pca.cox_efs_pvalue<0.05,]

hc <- hclust(dist(t(pca.selected.surv.feature)))
plot(hc)

surv.kmeans.pca <- kmeans(t(pca.selected.surv.feature),centers=2)
center.pca = surv.kmeans.pca$centers
surv.cluster.pca <- surv.kmeans.pca$cluster
table(surv.cluster.pca)
plot(t(pca.selected.surv.feature),col=surv.cluster.pca,pch=20)

sdf.pca <- survdiff(Surv(efs_time, efs_event) ~ surv.cluster.pca)
pvalue <- 1 - pchisq(sdf.pca$chisq, length(sdf.pca$n) - 1)
pvalue
sdf.pca <- survdiff(Surv(os_time, os_event) ~ surv.cluster.pca)
pvalue <- 1 - pchisq(sdf.pca$chisq, length(sdf.pca$n) - 1)
pvalue
kmsurvival.pca <- survfit(Surv(efs_time,efs_event) ~ surv.cluster.pca)
ggsurv(kmsurvival.pca,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="EFS p-value 0.068")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))
summary(coxph(Surv(efs_time, efs_event) ~ surv.cluster.pca))$concordance

kmsurvival.pca <- survfit(Surv(os_time,os_event) ~ surv.cluster.pca)
ggsurv(kmsurvival.pca,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="OS p-value 0.012")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))
summary(coxph(Surv(os_time, os_event) ~ surv.cluster.pca))$concordance

############
library(survcomp)
fit <- coxph(Surv(os_time, os_event) ~ surv.cluster.pca)
cindex <- concordance.index(predict(fit),surv.time = os_time, surv.event = os_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper

fit <- coxph(Surv(efs_time, efs_event) ~ surv.cluster.pca)
cindex <- concordance.index(predict(fit),surv.time = efs_time, surv.event = efs_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper
############

############
library(iClusterPlus)
gene.efs.pvalue <- apply(train.gene.matrix,1,function(x){y=coxph(Surv(efs_time,efs_event)~x);y=summary(y)$sctest[3]})
gene.os.pvalue <- apply(train.gene.matrix,1,function(x){y=coxph(Surv(os_time,os_event)~x);y=summary(y)$sctest[3]})
cnv.efs.pvalue <- apply(rbind(train.cnv.matrix.del,train.cnv.matrix.dup),1,function(x){y=coxph(Surv(efs_time,efs_event)~x);y=summary(y)$sctest[3]})
cnv.os.pvalue <- apply(rbind(train.cnv.matrix.del,train.cnv.matrix.dup),1,function(x){y=coxph(Surv(os_time,efs_event)~x);y=summary(y)$sctest[3]})
icluster.train.gene.matrix <- train.gene.matrix[gene.efs.pvalue<0.05|gene.os.pvalue<0.05,]
icluster.train.cnv.matrix <- rbind(train.cnv.matrix.del,train.cnv.matrix.dup)[cnv.efs.pvalue<0.05|cnv.os.pvalue<0.05,]
iCluster.NB <- iClusterPlus(t(icluster.train.gene.matrix),t(icluster.train.cnv.matrix),
                              type = c("gaussian","gaussian"),K=1,alpha=c(1,1),
                              lambda=c(0.03,0.03),
                              n.burnin=100,n.draw=200,maxiter=20,
                              sdev=0.05,eps=1.0e-4)
surv.cluster.iCluster <- iCluster.NB$clusters

sdf.iCluster <- survdiff(Surv(efs_time, efs_event) ~ surv.cluster.iCluster)
pvalue <- 1 - pchisq(sdf.iCluster$chisq, length(sdf.iCluster$n) - 1)
pvalue
sdf.iCluster <- survdiff(Surv(os_time, os_event) ~ surv.cluster.iCluster)
pvalue <- 1 - pchisq(sdf.iCluster$chisq, length(sdf.iCluster$n) - 1)
pvalue
kmsurvival.iCluster <- survfit(Surv(efs_time,efs_event) ~ surv.cluster.iCluster)
ggsurv(kmsurvival.iCluster,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="EFS p-value 1.22e-4")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))
summary(coxph(Surv(efs_time, efs_event) ~ surv.cluster.iCluster))$concordance

kmsurvival.iCluster <- survfit(Surv(os_time,os_event) ~ surv.cluster.iCluster)
ggsurv(kmsurvival.iCluster,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="OS p-value 3.76e-5")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))
summary(coxph(Surv(os_time, os_event) ~ surv.cluster.iCluster))$concordance

############
library(survcomp)
fit <- coxph(Surv(os_time, os_event) ~ surv.cluster.iCluster)
cindex <- concordance.index(predict(fit),surv.time = os_time, surv.event = os_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper

fit <- coxph(Surv(efs_time, efs_event) ~ surv.cluster.iCluster)
cindex <- concordance.index(predict(fit),surv.time = efs_time, surv.event = efs_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper
############

