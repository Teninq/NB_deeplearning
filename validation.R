setwd("D:/NB_deeplearning/autoencoder")
load("EnsID2Symbol.Rdata")
# load("GSE73517.expr.Rdata")
load("SEQC.NB.expr.Rdata")
source("D:/sourceRScript/NTP.R")
library(GGally)
library(ggplot2)
library(gtable)
library(survival)
library(survcomp)

row.names(SEQC.gene.expr) <- EnsID2Symbol[match(row.names(SEQC.gene.expr),EnsID2Symbol[,1]),2]
SEQC.NB.clinical.high <- SEQC.NB.clinical[grep("HR",SEQC.NB.clinical[,"HighRisk"]),]
SEQC.gene.expr.high <- SEQC.gene.expr[,match(SEQC.NB.clinical.high[,1],substr(colnames(SEQC.gene.expr),6,nchar(colnames(SEQC.gene.expr))))]
SEQC.gene.expr.high <- SEQC.gene.expr.high[apply(SEQC.gene.expr.high>0,1,sum)>50,]
# SEQC.gene.expr.high <- SEQC.gene.expr.high[apply(SEQC.gene.expr.high>0,1,sum)>18,]

##########NTP
SEQC.gene.expr.high <- SEQC.gene.expr.high[match(row.names(train.selected.overlap.matrix),row.names(SEQC.gene.expr.high)),]
gene.class <- c(rep("One",113),rep("Two",295))
gene.class.mat <- factor(gene.class, levels = c("One","Two"))
gene.class.matrix <- model.matrix(~0 + gene.class.mat)
colnames(gene.class.matrix) <- levels(gene.class.mat)
row.names(gene.class.matrix) <- class.label
SEQC.expr.filter <- log(SEQC.gene.expr.high[match(class.label,row.names(SEQC.gene.expr.high)),]+1)
SEQC.pre <- NTP_Fuc(SEQC.expr.filter,gene.class.matrix)
SEQC.label <- SEQC.pre[SEQC.pre[,6]<0.05,2]
table(SEQC.label)
SEQC.expr.filter.NTP <- SEQC.expr.filter[,SEQC.pre[,6]<0.05]
SEQC.NB.clinical.high.NTP <- SEQC.NB.clinical.high[SEQC.pre[,6]<0.05,]

SEQC.efs_time <- as.numeric(SEQC.NB.clinical.high.NTP[match(substr(colnames(SEQC.expr.filter.NTP),6,nchar(colnames(SEQC.expr.filter.NTP))),
                                                            SEQC.NB.clinical.high.NTP[,1]), "EFS_d"])
SEQC.efs_event <- as.numeric(SEQC.NB.clinical.high.NTP[match(substr(colnames(SEQC.expr.filter.NTP),6,nchar(colnames(SEQC.expr.filter.NTP))),
                                                             SEQC.NB.clinical.high.NTP[,1]),"E_EFS_HR"])


SEQC.os_time <- as.numeric(SEQC.NB.clinical.high.NTP[match(substr(colnames(SEQC.expr.filter.NTP),6,nchar(colnames(SEQC.expr.filter.NTP))),
                                                           SEQC.NB.clinical.high.NTP[,1]),"OS_d"])
SEQC.os_event <- as.numeric(SEQC.NB.clinical.high.NTP[match(substr(colnames(SEQC.expr.filter.NTP),6,nchar(colnames(SEQC.expr.filter.NTP))),
                                                            SEQC.NB.clinical.high.NTP[,1]),"F_OS_HR"])

sdf <- survdiff(Surv(SEQC.efs_time, SEQC.efs_event) ~ SEQC.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
sdf <- survdiff(Surv(SEQC.os_time, SEQC.os_event) ~ SEQC.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue


train.selected.overlap.matrix <- train.selected.matrix[row.names(train.selected.matrix) %in% intersect(row.names(SEQC.gene.expr.high),
                                                                                                     row.names(all.target.expr.filter)),]
# train.selected.overlap.matrix <- train.selected.matrix[row.names(train.selected.matrix) %in% row.names(SEQC.gene.expr.high),]
#########################SVM
###############ANOVA
train.gene.fvalue <- apply(train.selected.overlap.matrix,1,function(x){summary(aov(surv.cluster~x))[[1]]$`F value`[1]})
# train.gene.pvalue <- apply(train.selected.overlap.matrix,1,function(x){summary(aov(surv.cluster~x))[[1]]$`Pr(>F)`[1]})
svm.train.matrix <- train.selected.overlap.matrix[order(train.gene.fvalue,decreasing = TRUE)[1:56],]
SEQC.expr.filter <- scale(log(SEQC.gene.expr.high[match(row.names(svm.train.matrix),row.names(SEQC.gene.expr.high)),]+1))
svm.train.matrix.filter <- scale(svm.train.matrix)


SEQC.gene.expr.high.filter <- SEQC.gene.expr.high[match(row.names(svm.train.matrix),row.names(SEQC.gene.expr.high)),]
####SVM
library(e1071)
accuracy <- rep(0,10)
for(i in 1:10){
  a=summary(svm(t(svm.train.matrix.filter),surv.cluster,type="C-classification",
                cost=500,cross = 5,kernel="linear",probability=TRUE,scale=FALSE))
  accuracy[i] <- a$tot.accuracy
}
mean(accuracy)

rocall <- rep(0,10)
for(i in 1:10){
  samplelabel <- sample(190,19,replace = FALSE)
  a=svm(t(svm.train.matrix.filter[,-samplelabel]),surv.cluster[-samplelabel],type="C-classification",
        cost=500,cross = 5,kernel="linear",probability=TRUE,scale=FALSE)
  pre.label=predict(a,t(svm.train.matrix.filter[,samplelabel]),decision.values=TRUE)
  auc <- roc(surv.cluster[samplelabel],as.numeric(attr(pre.label,"decision.values")))
  rocall[i] <- auc$auc
}
mean(unlist(rocall))

library(pROC)
model.svm<-svm(t(svm.train.matrix.filter),surv.cluster,type="C-classification",
               cost=500,cross = 5,kernel="linear",probability=TRUE,scale=FALSE)
pre.svm=predict(model.svm,t(svm.train.matrix.filter),decision.values=TRUE)
roca <- roc(surv.cluster,as.numeric(attr(pre.svm,"decision.values")))
table(surv.cluster,pre.svm)
pre.svm=predict(model.svm,t(SEQC.expr.filter),decision.values=TRUE)
table(pre.svm)
SEQC.label <- pre.svm

svm.train.matrix <- train.selected.overlap.matrix[order(train.gene.fvalue,decreasing = TRUE)[1:57],]
svm.train.matrix.filter <- scale(svm.train.matrix)
SEQC.expr.filter.new <- scale(log(SEQC.gene.expr.high[match(row.names(svm.train.matrix),row.names(SEQC.gene.expr.high)),]+0.01))
model.svm<-svm(t(svm.train.matrix.filter),surv.cluster,type="C-classification",
               cost=500,cross = 5,kernel="linear",probability=TRUE,scale=FALSE)
pre.svm=predict(model.svm,t(svm.train.matrix.filter),decision.values=TRUE)
roc(surv.cluster,as.numeric(attr(pre.svm,"decision.values")))
table(surv.cluster,pre.svm)
pre.svm.1=predict(model.svm,t(SEQC.expr.filter.new),decision.values=TRUE)
table(pre.svm.1)
SEQC.label.1 <- pre.svm.1
table(SEQC.label,SEQC.label.1)

sdf <- survdiff(Surv(SEQC.efs_time, SEQC.efs_event) ~ SEQC.label.1)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
sdf <- survdiff(Surv(SEQC.os_time, SEQC.os_event) ~ SEQC.label.1)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue

####naive byase
svm.train.matrix <- train.selected.overlap.matrix[order(train.gene.fvalue,decreasing = TRUE)[1:46],]
svm.train.matrix.filter <- scale(svm.train.matrix)

library(caret)
for(i in 1:10){
  a=train(t(svm.train.matrix.filter),as.factor(surv.cluster), method = "nb",
          trControl = trainControl(method = "cv", number = 10))
  accuracy[i] <- a$results[2,5]
}
mean(accuracy)

library(naivebayes)
bayes_model <- naiveBayes(t(svm.train.matrix.filter),as.factor(surv.cluster))
roc(surv.cluster,as.numeric(predict(bayes_model,t(svm.train.matrix.filter),type = "raw")[,1]))
pre.bys=predict(bayes_model,t(svm.train.matrix.filter),decision.values=TRUE)
# predict(bayes_model,t(svm.train.matrix.filter),type = "raw")
table(surv.cluster,pre.bys)
rocb <- roc(surv.cluster,as.numeric(predict(bayes_model,t(svm.train.matrix.filter),type = "raw")[,1]))
pre.bys=predict(bayes_model,t(SEQC.expr.filter),decision.values=TRUE)
table(pre.bys)
# SEQC.label <- pre.bys

for(i in 1:10){
  samplelabel <- sample(190,95,replace = FALSE)
  a=naiveBayes(t(svm.train.matrix.filter[,-samplelabel]),as.factor(surv.cluster[-samplelabel]))
  pre.label=predict(a,t(svm.train.matrix.filter[,samplelabel]),decision.values=TRUE)
  auc <- roc(surv.cluster[samplelabel],as.numeric(pre.label))
  rocall[i] <- auc$auc
}
mean(unlist(rocall))

####logistic
svm.train.matrix <- train.selected.overlap.matrix[order(train.gene.fvalue,decreasing = TRUE)[1:44],]
svm.train.matrix.filter <- scale(svm.train.matrix)
for(i in 1:10){
  a=train(t(svm.train.matrix.filter),as.factor(surv.cluster), method = "glm",
          trControl = trainControl(method = "cv", number = 10))
  accuracy[i] <- a$results[1,2]
}
mean(accuracy)

svm.train.matrix.filter.bys <- data.frame(t(train.selected.overlap.matrix[order(train.gene.fvalue,decreasing = TRUE)[1:44],]))
svm.train.matrix.filter <- scale(svm.train.matrix.filter.bys)
svm.train.matrix.filter.bys$surv <- surv.cluster-1
logit <- glm(surv ~ .,svm.train.matrix.filter.bys,family = binomial(link = "logit"),control=list(maxit=100))
pre.logit <- predict(logit,type='response')
table(surv.cluster,as.numeric(pre.logit>0.5)+1)
rocc <- roc(surv.cluster,pre.logit)
pre.logit <- predict(pre.logit,t(SEQC.expr.filter),decision.values=TRUE)
table(pre.logit)
# SEQC.label <- pre.logit

for(i in 1:10){
  samplelabel <- sample(190,19,replace = FALSE)
  a <- glm(surv ~ .,svm.train.matrix.filter.bys[-samplelabel,],family = binomial(link = "logit"),control=list(maxit=100))
  pre.label=predict(a,type='response')
  auc <- roc(surv.cluster[samplelabel],pre.label[samplelabel])
  rocall[i] <- auc$auc
}
mean(unlist(rocall))


sdf <- survdiff(Surv(SEQC.efs_time, SEQC.efs_event) ~ SEQC.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
sdf <- survdiff(Surv(SEQC.os_time, SEQC.os_event) ~ SEQC.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue


############xgboost
SEQC.expr.filter <- scale(log(SEQC.gene.expr.high[match(row.names(train.selected.overlap.matrix),row.names(SEQC.gene.expr.high)),]+1))
dtrain <- xgb.DMatrix(t(scale(train.selected.overlap.matrix)), label = as.numeric(surv.cluster==1))
param <- list(max_depth = 5, eta = 0.2, silent = 1, nthread = 2,
              objective = "binary:logistic", eval_metric = "auc",
              lambda = 3,alpha = 2,min_child_weight = 8)
bst <- xgb.train(param,dtrain,nround = 3)
# param <- list(max_depth = 8, eta = 1, silent = 1, nthread = 2,
#               objective = "binary:logistic", eval_metric = "auc",
#               lambda = 3,min_child_weight = 8)
# my_etas <- list(eta = seq(0.001,0.5,by = 0.001))
# bst <- xgb.cv(param,dtrain,nround = 500,nfold = 10,callbacks = list(cb.reset.parameters(my_etas)))
bst <- xgb.cv(param,dtrain,nround = 20,nfold = 10)
  # bst <- xgb.train(param,dtrain,nround = 500,callbacks = list(cb.reset.parameters(my_etas)))
names <- row.names(train.selected.overlap.matrix)
importance.matrix <- xgb.importance(names,model = bst)
dim(importance.matrix)
pred <- as.numeric(predict(bst,t(scale(train.selected.overlap.matrix)))<0.5) + 1
roc(surv.cluster,predict(bst,t(scale(train.selected.overlap.matrix))))
rocd <- roc(surv.cluster,predict(bst,t(scale(train.selected.overlap.matrix))))
table(pred,surv.cluster)

ggroc(list(svm=roca,naive_bayes=rocb,logistic=rocc,xgboost=rocd),size=1.5,linetype = 1)

# SEQC.label <- as.numeric(predict(bst,t(SEQC.expr.filter)) < 0.5) + 1
table(as.numeric(predict(bst,t(SEQC.expr.filter)) < 0.5) + 1)

SEQC.efs_time <- as.numeric(SEQC.NB.clinical.high[match(substr(colnames(SEQC.expr.filter),6,nchar(colnames(SEQC.expr.filter))),
                                                        SEQC.NB.clinical.high[,1]), "EFS_d"])
SEQC.efs_event <- as.numeric(SEQC.NB.clinical.high[match(substr(colnames(SEQC.expr.filter),6,nchar(colnames(SEQC.expr.filter))),
                                                         SEQC.NB.clinical.high[,1]),"E_EFS_HR"])


SEQC.os_time <- as.numeric(SEQC.NB.clinical.high[match(substr(colnames(SEQC.expr.filter),6,nchar(colnames(SEQC.expr.filter))),
                                                       SEQC.NB.clinical.high[,1]),"OS_d"])
SEQC.os_event <- as.numeric(SEQC.NB.clinical.high[match(substr(colnames(SEQC.expr.filter),6,nchar(colnames(SEQC.expr.filter))),
                                                        SEQC.NB.clinical.high[,1]),"F_OS_HR"])


sdf <- survdiff(Surv(SEQC.efs_time, SEQC.efs_event) ~ SEQC.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
sdf <- survdiff(Surv(SEQC.os_time, SEQC.os_event) ~ SEQC.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue


kmsurvival <- survfit(Surv(SEQC.efs_time,SEQC.efs_event) ~ SEQC.label)
ggsurv(kmsurvival,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="EFS p-value 7.93e-5")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))
summary(coxph(Surv(SEQC.efs_time, SEQC.efs_event) ~ SEQC.label))$concordance

kmsurvival <- survfit(Surv(SEQC.os_time,SEQC.os_event) ~ SEQC.label)
ggsurv(kmsurvival,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="OS p-value 4.62e-6")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
summary(coxph(Surv(SEQC.os_time, SEQC.os_event) ~ SEQC.label))$concordance

############
fit <- coxph(Surv(SEQC.efs_time, SEQC.efs_event) ~ SEQC.label)
cindex <- concordance.index(predict(fit),surv.time = SEQC.efs_time, surv.event = SEQC.efs_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper

fit <- coxph(Surv(SEQC.os_time, SEQC.os_event) ~ SEQC.label)
cindex <- concordance.index(predict(fit),surv.time = SEQC.os_time, surv.event = SEQC.os_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper
############


sdf <- survdiff(Surv(SEQC.efs_time, SEQC.efs_event) ~ SEQC.label.1)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
sdf <- survdiff(Surv(SEQC.os_time, SEQC.os_event) ~ SEQC.label.1)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue

kmsurvival <- survfit(Surv(SEQC.efs_time,SEQC.efs_event) ~ SEQC.label.1)
ggsurv(kmsurvival,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="EFS p-value 5.66e-6")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))

kmsurvival <- survfit(Surv(SEQC.os_time,SEQC.os_event) ~ SEQC.label.1)
ggsurv(kmsurvival,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="OS p-value 1.28e-5")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


###################################
##############xgboost + svm
train.overlap.gene.matrix.1 <- train.gene.matrix[row.names(train.gene.matrix) %in% row.names(SEQC.gene.expr.high),]
dtrain <- xgb.DMatrix(t(train.overlap.gene.matrix.1), label = as.numeric(surv.cluster==1))
param <- list(max_depth = 8, eta = 0.8055, silent = 1, nthread = 2,
              objective = "binary:logistic", eval_metric = "auc",
              lambda = 1e-9,alpha = 1e-15)
# bst <- xgb.cv(param,dtrain,nround = 50,nfold = 10)
bst <- xgb.train(param,dtrain,nround = 50)
names <- row.names(train.overlap.gene.matrix.1)
importance.matrix <- xgb.importance(names,model = bst)
dim(importance.matrix)
svm.train.matrix.filter <- scale(train.overlap.gene.matrix.1[match(importance.matrix$Feature,row.names(train.overlap.gene.matrix.1)),])

####
dtrain <- xgb.DMatrix(t(train.selected.overlap.matrix), label = as.numeric(surv.cluster==1))
# param <- list(max_depth = 8, eta = 0.7944, silent = 1, nthread = 2,
#               objective = "binary:logistic", eval_metric = "auc",
#               lambda = 1e-9,alpha = 1e-7)
param <- list(max_depth = 8, eta = 0.02, silent = 1, nthread = 2,
              objective = "binary:logistic", eval_metric = "auc",
              lambda = 12,min_child_weight = 9)
# bst <- xgb.cv(param,dtrain,nround = 250,nfold = 10)
bst <- xgb.train(param,dtrain,nround = 500)
names <- row.names(train.selected.overlap.matrix)
importance.matrix <- xgb.importance(names,model = bst)
dim(importance.matrix)


# SEQC.expr.filter <- scale(log(SEQC.gene.expr.high[match(row.names(target.test.expr),row.names(SEQC.gene.expr.high)),]+1))
# svm.train.matrix.filter <- scale(train.overlap.gene.matrix.1[match(row.names(target.test.expr),row.names(train.overlap.gene.matrix.1)),])
####
SEQC.expr.filter <- scale(log(SEQC.gene.expr.high[match(importance.matrix$Feature,row.names(SEQC.gene.expr.high)),]+0.01))
svm.train.matrix.filter <- scale(train.selected.overlap.matrix[match(importance.matrix$Feature,row.names(train.selected.overlap.matrix)),])
####
model.svm<-svm(t(svm.train.matrix.filter),surv.cluster,type="C-classification",
               cost=5,cross=10,kernel="linear",probability=TRUE,scale=FALSE)
pre.svm=predict(model.svm,t(svm.train.matrix.filter),decision.values=TRUE)
table(surv.cluster,pre.svm)
pre.svm=predict(model.svm,t(SEQC.expr.filter),decision.values=TRUE)
table(pre.svm)
SEQC.label <- pre.svm
