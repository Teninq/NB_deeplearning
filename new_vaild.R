setwd("D:/NB_deeplearning/autoencoder")
all.target <- read.table("TARGET_NBL_gene_expr.txt",header = T,row.names = 1,sep = "\t")
ens2symb <- read.table("GRCh37.59.ens2symb.txt",header = F)
ALL.TARGET.samples <- read.table("TARGET_NBL_Clinical.txt",header = T,sep = "\t")
# load("EnsID2Symbol.Rdata")
all.target.expr <- all.target[row.names(all.target) %in% ens2symb$V1,]
all.target.expr <- all.target.expr[!duplicated(ens2symb[match(row.names(all.target.expr),ens2symb$V1),2]),]
row.names(all.target.expr) <- ens2symb[match(row.names(all.target.expr),ens2symb$V1),2]
all.target.expr.filter <- all.target.expr[apply(all.target.expr>0,1,sum)>16,]
table(importance.matrix$Feature %in% row.names(all.target.expr.filter))

all.target.expr.filter <- all.target.expr.filter[,grep("01A$",colnames(all.target.expr.filter))]
colnames(all.target.expr.filter) <- gsub("\\.","-",substr(colnames(all.target.expr.filter),1,nchar(colnames(all.target.expr.filter))-4))
all.target.expr.filter <- all.target.expr.filter[,!(colnames(all.target.expr.filter) %in% colnames(train.overlap.gene.matrix))]
new.target.samples <- ALL.TARGET.samples[match(colnames(all.target.expr.filter),ALL.TARGET.samples[,"TARGET.USI"]),]
new.target.samples <- new.target.samples[grep("^High",new.target.samples$COG.Risk.Group),]

new_efs_time <- as.numeric(new.target.samples[,"Event.Free.Survival.Time.in.Days"])
new_efs_event <- as.character(new.target.samples[,"First.Event"])
new_efs_event <- as.numeric(new_efs_event!="Censored")
# cox_efs_pvalue <- apply(selected.feature,1,function(x){y=coxph(Surv(efs_time,efs_event)~x);y=summary(y)$sctest[3]})

new_os_time <- as.numeric(new.target.samples[,"Overall.Survival.Time.in.Days"])
new_os_event <- as.character(new.target.samples[,"Vital.Status"])
new_os_event <- as.numeric(new_os_event=="Dead")

#################ANOVA
all.target.expr.aov <- scale(log(all.target.expr.filter[match(row.names(svm.train.matrix),row.names(all.target.expr.filter)),
                                                        colnames(target.test.expr) %in% new.target.samples$TARGET.USI]+0.1))
pre.svm=predict(model.svm,t(all.target.expr.aov),decision.values=TRUE)
table(pre.svm)
TARGET.label <- pre.svm

sdf <- survdiff(Surv(new_efs_time,new_efs_event) ~ TARGET.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
sdf <- survdiff(Surv(new_os_time,new_os_event) ~ TARGET.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue


target.test.expr <- scale(all.target.expr.filter[na.omit(match(importance.matrix$Feature,row.names(all.target.expr.filter))),])
target.test.expr.filter <- target.test.expr[,colnames(target.test.expr) %in% new.target.samples$TARGET.USI]
#################SVM
svm.train.matrix.filter <- scale(train.overlap.gene.matrix.1[match(row.names(target.test.expr),row.names(train.overlap.gene.matrix.1)),])
# svm.train.matrix.filter <- scale(train.overlap.gene.matrix[match(row.names(target.test.expr),row.names(train.overlap.gene.matrix)),])
model.svm<-svm(t(svm.train.matrix.filter),surv.cluster,type="C-classification",
               cost=50,cross=10,kernel="linear",probability=TRUE,scale=FALSE)
pre.svm=predict(model.svm,t(svm.train.matrix.filter),decision.values=TRUE)
table(surv.cluster,pre.svm)
pre.svm=predict(model.svm,t(target.test.expr.filter),decision.values=TRUE)
table(pre.svm)
TARGET.label <- pre.svm

##################xgboost
all.target.train <- all.target.expr.filter[match(row.names(train.selected.overlap.matrix),row.names(all.target.expr.filter)),
                                           colnames(target.test.expr) %in% new.target.samples$TARGET.USI]
TARGET.label <- as.numeric(predict(bst,t(all.target.train)) < 0.5) + 1
table(TARGET.label)
##################

library(GGally)
library(ggplot2)
library(gtable)
library(survival)
sdf <- survdiff(Surv(new_efs_time,new_efs_event) ~ TARGET.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
sdf <- survdiff(Surv(new_os_time,new_os_event) ~ TARGET.label)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue


kmsurvival <- survfit(Surv(new_efs_time,new_efs_event) ~ TARGET.label)
ggsurv(kmsurvival,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="EFS p-value 0.03")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))
summary(coxph(Surv(new_efs_time, new_efs_event) ~ TARGET.label))$concordance

kmsurvival <- survfit(Surv(new_os_time,new_os_event) ~ TARGET.label)
ggsurv(kmsurvival,size.est = 1.2,
       cens.size = 3,
       cens.shape = 3,
       surv.col=c("#E41A1C","#377EB8"),
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="OS p-value 0.05")+
  ylim(0,1)+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
summary(coxph(Surv(new_os_time, new_os_event) ~ TARGET.label))$concordance

############
library(survcomp)
fit <- coxph(Surv(new_os_time, new_os_event) ~ TARGET.label)
cindex <- concordance.index(predict(fit),surv.time = new_os_time, surv.event = new_os_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper

fit <- coxph(Surv(new_efs_time, new_efs_event) ~ TARGET.label)
cindex <- concordance.index(predict(fit),surv.time = new_efs_time, surv.event = new_efs_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper
############


