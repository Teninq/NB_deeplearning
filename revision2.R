setwd("D:/NB_deeplearning/autoencoder")
load("EnsID2Symbol.Rdata")
load("SEQC.NB.expr.Rdata")
load("train.selected.RData")
load("surv.cluster.RData")

row.names(SEQC.gene.expr) <- EnsID2Symbol[match(row.names(SEQC.gene.expr),EnsID2Symbol[,1]),2]
SEQC.NB.clinical.high <- SEQC.NB.clinical[grep("HR",SEQC.NB.clinical[,"HighRisk"]),]
SEQC.gene.expr.high <- SEQC.gene.expr[,match(SEQC.NB.clinical.high[,1],substr(colnames(SEQC.gene.expr),6,nchar(colnames(SEQC.gene.expr))))]

res <- cbind(apply(SEQC.gene.expr.high,2,function(x){y=sum(x>0);sum(x>0.01)/y}),
             apply(SEQC.gene.expr.high,2,function(x){y=sum(x>0);sum(x>0.05)/y}),
             apply(SEQC.gene.expr.high,2,function(x){y=sum(x>0);sum(x>0.1)/y}),
             apply(SEQC.gene.expr.high,2,function(x){y=sum(x>0);sum(x>0.5)/y}),
             apply(SEQC.gene.expr.high,2,function(x){y=sum(x>0);sum(x>1)/y}))
colnames(res) <- c("0.01","0.05","0.1","0.5","1")

SEQC.gene.expr.high.change <- SEQC.gene.expr.high + as.numeric(SEQC.gene.expr.high==0)
pseudo.value <- max(apply(SEQC.gene.expr.high.change,2,min))
pseudo.value1 <- min(apply(SEQC.gene.expr.high.change,2,min))
SEQC.gene.expr.high[SEQC.gene.expr.high<pseudo.value] <- pseudo.value1

all.target <- read.table("TARGET_NBL_gene_expr.txt",header = T,row.names = 1,sep = "\t")
ens2symb <- read.table("GRCh37.59.ens2symb.txt",header = F)
ALL.TARGET.samples <- read.table("TARGET_NBL_Clinical.txt",header = T,sep = "\t")
# load("EnsID2Symbol.Rdata")
all.target.expr <- all.target[row.names(all.target) %in% ens2symb$V1,]
all.target.expr <- all.target.expr[!duplicated(ens2symb[match(row.names(all.target.expr),ens2symb$V1),2]),]
row.names(all.target.expr) <- ens2symb[match(row.names(all.target.expr),ens2symb$V1),2]
all.target.expr.filter <- all.target.expr[apply(all.target.expr>0,1,sum)>16,]

train.selected.overlap.matrix <- train.selected.matrix[row.names(train.selected.matrix) %in% intersect(row.names(SEQC.gene.expr.high),
                                                                                                       row.names(all.target.expr.filter)),]

library(e1071)
library(pROC)
train.gene.fvalue <- apply(train.selected.overlap.matrix,1,function(x){summary(aov(surv.cluster~x))[[1]]$`F value`[1]})
svm.train.matrix <- train.selected.overlap.matrix[order(train.gene.fvalue,decreasing = TRUE)[1:57],]
svm.train.matrix.filter <- scale(svm.train.matrix)
SEQC.expr.filter.new <- scale(log(SEQC.gene.expr.high[match(row.names(svm.train.matrix),row.names(SEQC.gene.expr.high)),]))
model.svm<-svm(t(svm.train.matrix.filter),surv.cluster,type="C-classification",
               cost=500,cross = 5,kernel="linear",probability=TRUE,scale=FALSE)
pre.svm=predict(model.svm,t(svm.train.matrix.filter),decision.values=TRUE)
roc(surv.cluster,as.numeric(attr(pre.svm,"decision.values")))
table(surv.cluster,pre.svm)
pre.svm.1=predict(model.svm,t(SEQC.expr.filter.new),decision.values=TRUE)
table(pre.svm.1)
SEQC.label.1 <- pre.svm.1


SEQC.efs_time <- as.numeric(SEQC.NB.clinical.high[match(substr(colnames(SEQC.expr.filter.new),6,nchar(colnames(SEQC.expr.filter.new))),
                                                        SEQC.NB.clinical.high[,1]), "EFS_d"])
SEQC.efs_event <- as.numeric(SEQC.NB.clinical.high[match(substr(colnames(SEQC.expr.filter.new),6,nchar(colnames(SEQC.expr.filter.new))),
                                                         SEQC.NB.clinical.high[,1]),"E_EFS_HR"])


SEQC.os_time <- as.numeric(SEQC.NB.clinical.high[match(substr(colnames(SEQC.expr.filter.new),6,nchar(colnames(SEQC.expr.filter.new))),
                                                       SEQC.NB.clinical.high[,1]),"OS_d"])
SEQC.os_event <- as.numeric(SEQC.NB.clinical.high[match(substr(colnames(SEQC.expr.filter.new),6,nchar(colnames(SEQC.expr.filter.new))),
                                                        SEQC.NB.clinical.high[,1]),"F_OS_HR"])

library(survcomp)
sdf <- survdiff(Surv(SEQC.efs_time, SEQC.efs_event) ~ SEQC.label.1)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue
sdf <- survdiff(Surv(SEQC.os_time, SEQC.os_event) ~ SEQC.label.1)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
pvalue


library(GGally)
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

############
fit <- coxph(Surv(SEQC.efs_time, SEQC.efs_event) ~ SEQC.label.1)
cindex <- concordance.index(predict(fit),surv.time = SEQC.efs_time, surv.event = SEQC.efs_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper

fit <- coxph(Surv(SEQC.os_time, SEQC.os_event) ~ SEQC.label.1)
cindex <- concordance.index(predict(fit),surv.time = SEQC.os_time, surv.event = SEQC.os_event,method = "noether")
cindex$c.index; cindex$lower; cindex$upper
############
