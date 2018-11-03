
SEQC.gene.expr.overlap <- SEQC.gene.expr[match(intersect(row.names(gene.matrix),row.names(SEQC.gene.expr)),row.names(SEQC.gene.expr)),]
gene.matrix.overlap <- gene.matrix[match(intersect(row.names(gene.matrix),row.names(SEQC.gene.expr)),row.names(gene.matrix)),]

SEQC.gene.expr.overlap.norm <- scale(log(SEQC.gene.expr.overlap+0.1))
gene.matrix.overlap.norm <- scale(gene.matrix.overlap)

library(sva)
batch = rep(c("SEQC","TARGET"),c(ncol(SEQC.gene.expr.overlap.norm),ncol(gene.matrix.overlap.norm)))
expr.all <- cbind(SEQC.gene.expr.overlap.norm,gene.matrix.overlap.norm)
expr.combat = ComBat(expr.all, batch)

SEQC.gene.expr.overlap.combat <- expr.combat[,1:ncol(SEQC.gene.expr.overlap.norm)]
gene.matrix.overlap.combat <- expr.combat[,(ncol(SEQC.gene.expr.overlap.norm)+1):ncol(expr.combat)]

SEQC.gene.expr.overlap.combat.high <- SEQC.gene.expr.overlap.combat[,match(SEQC.NB.clinical.high[,1],
                                                                             substr(colnames(SEQC.gene.expr),6,nchar(colnames(SEQC.gene.expr))))]
SEQC.gene.expr.overlap.combat.high.filter <- SEQC.gene.expr.overlap.combat.high[apply(SEQC.gene.expr.overlap.combat.high>0,1,sum)>50,]

train.gene.matrix.overlap.combat <- gene.matrix.overlap.combat[,match(intersect(colnames(cnv.matrix),colnames(gene.matrix.overlap.combat)),
                                                                      colnames(gene.matrix.overlap.combat))]

train.gene.matrix.overlap.combat.surv <- train.gene.matrix.overlap.combat[row.names(train.gene.matrix.overlap.combat) %in% 
                                                                              row.names(train.selected.overlap.matrix),]
train.gene.combat.fvalue <- apply(train.gene.matrix.overlap.combat.surv,1,function(x){summary(aov(surv.cluster~x))[[1]]$`F value`[1]})
train.gene.matrix.overlap.combat.select <- train.gene.matrix.overlap.combat.surv[order(train.gene.combat.fvalue,decreasing = TRUE)[1:56],]
SEQC.expr.filter.new <- SEQC.gene.expr.overlap.combat.high.filter[match(row.names(train.gene.matrix.overlap.combat.select),
                                                                        row.names(SEQC.gene.expr.overlap.combat.high.filter)),]
model.svm<-svm(t(train.gene.matrix.overlap.combat.select),surv.cluster,type="C-classification",
               cost=500,cross = 5,kernel="linear",probability=TRUE,scale=FALSE)
pre.svm=predict(model.svm,t(svm.train.matrix.filter),decision.values=TRUE)
roc(surv.cluster,as.numeric(attr(pre.svm,"decision.values")))
table(surv.cluster,pre.svm)
pre.svm.1=predict(model.svm,t(SEQC.expr.filter.new),decision.values=TRUE)
table(pre.svm.1)
SEQC.label.1 <- pre.svm.1
table(SEQC.label,SEQC.label.1)