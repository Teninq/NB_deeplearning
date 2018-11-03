setwd("D:/NB_deeplearning/autoencoder")
cnv.matrix <- read.table("TARGET_CNV.all_data_by_genes.txt",header = TRUE,sep = "\t",row.names = 1)
cnv.matrix <- cnv.matrix[,c(-1,-2)]
colnames(cnv.matrix) <- gsub("\\.","-",colnames(cnv.matrix))
test.cnv.matrix <- cnv.matrix[,!(colnames(cnv.matrix) %in% colnames(train.cnv.matrix))]

ALL.TARGET.samples <- read.table("TARGET_NBL_Clinical.txt",header = T,sep = "\t")
TARGET.CNV.clinic <- ALL.TARGET.samples[match(colnames(test.cnv.matrix),ALL.TARGET.samples$TARGET.USI),]
cnv_efs_time <- as.numeric(TARGET.CNV.clinic[,"Event.Free.Survival.Time.in.Days"])
cnv_efs_event <- as.character(TARGET.CNV.clinic[,"First.Event"])
cnv_efs_event <- as.numeric(cnv_efs_event!="Censored")

cnv_os_time <- as.numeric(TARGET.CNV.clinic[,"Overall.Survival.Time.in.Days"])
cnv_os_event <- as.character(TARGET.CNV.clinic[,"Vital.Status"])
cnv_os_event <- as.numeric(cnv_os_event=="Dead")

train.cnv.matrix.new <- rbind(train.cnv.matrix.dup,train.cnv.matrix.del)
os_cnv_pvalue <- apply(train.cnv.matrix.new,1,function(x){y=coxph(Surv(cnv_os_time,cnv_os_event)~x);y=summary(y)$sctest[3]})
efs_cnv_pvalue <- apply(train.cnv.matrix.new,1,function(x){y=coxph(Surv(cnv_efs_time,cnv_efs_event)~x);y=summary(y)$sctest[3]})
train.cnv.selected.matrix <- train.cnv.matrix.new[os_cnv_pvalue<0.05|efs_cnv_pvalue<0.05,]
test.cnv.selected.matrix <- test.cnv.matrix[match(row.names(train.cnv.selected.matrix),row.names(test.cnv.matrix)),]

##################
# dtrain <- xgb.DMatrix(t(train.cnv.matrix), label = as.numeric(surv.cluster==1))
dtrain <- xgb.DMatrix(t(train.cnv.selected.matrix), label = as.numeric(surv.cluster==1))
param <- list(max_depth = 10, eta = 0.2, silent = 1, nthread = 2,
              objective = "binary:logistic", eval_metric = "auc",
              lambda = 5,alpha = 1,min_child_weight = 6)
bst <- xgb.cv(param,dtrain,nround = 15,nfold = 10)

bst <- xgb.train(param,dtrain,nround = 350)
names <- row.names(train.cnv.selected.matrix)
importance.matrix <- xgb.importance(names,model = bst)
dim(importance.matrix)
pred <- as.numeric(predict(bst,t(train.cnv.selected.matrix))<0.5) + 1
roc1 <- roc(surv.cluster,predict(bst,t(train.cnv.selected.matrix)))
table(pred,surv.cluster)
cnv.label <- as.numeric(predict(bst,t(test.cnv.selected.matrix)) < 0.5) + 1
table(cnv.label)

rocall <- rep(0,10)
for(i in 1:10){
  samplelabel <- sample(190,19,replace = FALSE)
  dtrain <- xgb.DMatrix(t(train.cnv.selected.matrix[,-samplelabel]), label = as.numeric(surv.cluster==1)[-samplelabel])
  param <- list(max_depth = 10, eta = 0.2, silent = 1, nthread = 2,
                objective = "binary:logistic", eval_metric = "auc",
                lambda = 5,alpha = 1,min_child_weight = 6)
  bst <- xgb.train(param,dtrain,nround = 350)
  auc <- roc(surv.cluster[samplelabel],predict(bst,t(train.cnv.selected.matrix[,samplelabel])))
  rocall[i] <- auc$auc
}
mean(unlist(rocall))

sdf <- survdiff(Surv(cnv_efs_time, cnv_efs_event) ~ cnv.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
sdf <- survdiff(Surv(cnv_os_time, cnv_os_event) ~ cnv.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue

##################
train.cnv.fvalue <- apply(train.cnv.selected.matrix,1,function(x){summary(aov(surv.cluster~x))[[1]]$`F value`[1]})
svm.train.cnv.matrix <- train.cnv.selected.matrix[order(train.cnv.fvalue,decreasing = TRUE)[1:30],]
svm.test.cnv.matrix <- test.cnv.selected.matrix[match(row.names(svm.train.cnv.matrix),row.names(test.cnv.selected.matrix)),]

for(i in 1:10){
  a=summary(svm(t(svm.train.cnv.matrix),surv.cluster,type="C-classification",gamma = 2,
                cost=20,cross = 10,kernel="polynomial",probability=TRUE,scale=FALSE))
  accuracy[i] <- a$tot.accuracy
}
mean(accuracy)

model.svm<-svm(t(svm.train.cnv.matrix),surv.cluster,type="C-classification",
               cost=50,cross = 10,kernel="linear",probability=TRUE,scale=FALSE)
pre.svm=predict(model.svm,t(svm.train.cnv.matrix),decision.values=TRUE)
table(surv.cluster,pre.svm)
roc(surv.cluster,as.numeric(attr(pre.svm,"decision.values")))
roc2 <- roc(surv.cluster,as.numeric(attr(pre.svm,"decision.values")))
pre.svm=predict(model.svm,t(svm.test.cnv.matrix),decision.values=TRUE)
table(pre.svm)
cnv.label <- pre.svm

rocall <- rep(0,10)
for(i in 1:10){
  samplelabel <- sample(190,19,replace = FALSE)
  a=svm(t(svm.train.cnv.matrix[,-samplelabel]),surv.cluster[-samplelabel],type="C-classification",
        cost=20,cross = 10,kernel="linear",probability=TRUE,scale=FALSE)
  pre.label=predict(a,t(svm.train.cnv.matrix[,samplelabel]),decision.values=TRUE)
  auc <- roc(surv.cluster[samplelabel],as.numeric(attr(pre.label,"decision.values")))
  rocall[i] <- auc$auc
}
mean(unlist(rocall))

####naive byase
library(caret)
svm.train.cnv.matrix <- train.cnv.selected.matrix[order(train.cnv.fvalue,decreasing = TRUE)[1:24],]
for(i in 1:10){
  a=train(t(svm.train.cnv.matrix),as.factor(surv.cluster), method = "nb",
          trControl = trainControl(method = "cv", number = 10))
  accuracy[i] <- a$results[2,4]
}
mean(accuracy)

library(naivebayes)
svm.train.cnv.matrix <- train.cnv.selected.matrix[order(train.cnv.fvalue,decreasing = TRUE)[1:24],]
bayes_model <- naiveBayes(t(svm.train.cnv.matrix),as.factor(surv.cluster))
pre.bys=predict(bayes_model,as.data.frame(t(svm.train.cnv.matrix)))
table(surv.cluster,pre.bys)
roc(surv.cluster,as.numeric(predict(bayes_model,t(svm.train.cnv.matrix),type = "raw")[,1]))
roc3 <- roc(surv.cluster,as.numeric(predict(bayes_model,t(svm.train.cnv.matrix),type = "raw")[,1]))
pre.bys=predict(bayes_model,t(svm.test.cnv.matrix),decision.values=TRUE)
table(pre.bys)
SEQC.label <- pre.bys

for(i in 1:10){
  samplelabel <- sample(190,95,replace = FALSE)
  a=naiveBayes(t(svm.train.cnv.matrix[,-samplelabel]),as.factor(surv.cluster[-samplelabel]))
  pre.label=predict(a,t(svm.train.cnv.matrix[,samplelabel]),decision.values=TRUE)
  auc <- roc(surv.cluster[samplelabel],as.numeric(pre.label))
  rocall[i] <- auc$auc
}
mean(unlist(rocall))

####logistic
svm.train.cnv.matrix <- train.cnv.selected.matrix[order(train.cnv.fvalue,decreasing = TRUE)[1:15],]
for(i in 1:10){
  a=train(t(svm.train.cnv.matrix),as.factor(surv.cluster), method = "glm",
          trControl = trainControl(method = "cv", number = 10))
  accuracy[i] <- a$results[1,2]
}
mean(accuracy)

svm.train.cnv.matrix.bys <- data.frame(t(train.cnv.selected.matrix[order(train.cnv.fvalue,decreasing = TRUE)[1:15],]))
svm.train.cnv.matrix.bys$surv <- surv.cluster-1
logit <- glm(surv ~ .,svm.train.cnv.matrix.bys,family = binomial(link = "logit"))
pre.logit <- predict(logit,type='response')
table(surv.cluster,as.numeric(pre.logit>0.5)+1)
roc(surv.cluster,pre.logit)
roc4 <- roc(surv.cluster,pre.logit)
pre.logit <- predict(logit,data.frame(t(svm.test.cnv.matrix)),type='response')
table(pre.logit)
SEQC.label <- pre.logit

for(i in 1:10){
  samplelabel <- sample(190,19,replace = FALSE)
  a <- glm(surv ~ .,svm.train.cnv.matrix.bys[-samplelabel,],family = binomial(link = "logit"),control=list(maxit=100))
  pre.label=predict(a,type='response')
  auc <- roc(surv.cluster[samplelabel],pre.label[samplelabel])
  rocall[i] <- auc$auc
}
mean(unlist(rocall))

ggroc(list(svm=roc2,naive_bayes=roc3,logistic=roc4,xgboost=roc1),size=1.5,linetype = 1)

kmsurvival <- survfit(Surv(cnv_efs_time, cnv_efs_event) ~ cnv.label)
ggsurv(kmsurvival,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="EFS p-value 0.011")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))
summary(coxph(Surv(cnv_efs_time, cnv_efs_event) ~ cnv.label))$concordance

kmsurvival <- survfit(Surv(cnv_os_time, cnv_os_event) ~ cnv.label)
ggsurv(kmsurvival,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="OS p-value 0.036")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
summary(coxph(Surv(cnv_os_time, cnv_os_event) ~ cnv.label))$concordance

############
fit <- coxph(Surv(cnv_efs_time[!is.na(cnv_efs_time)],cnv_efs_event[!is.na(cnv_efs_time)]) ~ cnv.label[!is.na(cnv_efs_time)])
cindex <- concordance.index(predict(fit),surv.time = cnv_efs_time[!is.na(cnv_efs_time)],
                            surv.event = cnv_efs_event[!is.na(cnv_efs_time)],method = "noether")
cindex$c.index; cindex$lower; cindex$upper

fit <- coxph(Surv(cnv_os_time[!is.na(cnv_os_time)], cnv_os_event[!is.na(cnv_os_time)]) ~ cnv.label[!is.na(cnv_os_time)])
cindex <- concordance.index(predict(fit),surv.time = cnv_os_time[!is.na(cnv_os_time)],
                            surv.event = cnv_os_event[!is.na(cnv_os_time)],method = "noether")
cindex$c.index; cindex$lower; cindex$upper
############

