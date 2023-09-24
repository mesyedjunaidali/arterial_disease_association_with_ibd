setwd(".")

list.of.packages <- c("easypackages", "PRROC", "e1071", "pacman","readr","randomForest","caret","e1071","mltools","boot","umap","kernlab")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library("easypackages")
libraries(list.of.packages)

mean.fun <- function(d, i){    m <- mean(d$data[i])
n <- length(i)
v <- (n-1)*var(d$data[i])/n^2
c(m, v)
}
ci <- function(data){
  k <- as.data.frame( cbind(seq(length(data)),data)  )
  if(all(k[,2]==k[,2][1])){
    out  <- c(k[,2][1],k[,2][1],k[,2][1])
  }else{
    names(k) <- c("id","data")
    boot <- boot(k, mean.fun, R=999 );
    boot.ci <- boot.ci(boot, type = c("stud"))
    out <- c(boot.ci$student[4],boot.ci$t0,boot.ci$student[5])
  }
  names(out) <- c("lower","mean","upper")
  return(out)
}
m_std <- function(confmats){
mymcc <- myf1 <- myacc  <- mytpr <- mytnr <- myppv <- mynpv <- rep(0,length(confmats))

metrics <- data.frame(measure=c("MCC","F1","ACC","TPR","TNR","PPV","NPV"),mean=rep(0,7),sd=rep(0,7))

for(i in 1:length(confmats)){
  mymcc[i] <- mcc(confusionM = confmats[[i]])
  TP <- confmats[[i]][1,1]
  TN <- confmats[[i]][2,2]
  FP <- confmats[[i]][1,2]
  FN <- confmats[[i]][2,1]
  myf1[i] <- ifelse(2*TP+FP+FN>0,2*TP/(2*TP+FP+FN),NA)
  myacc[i] <- ifelse(TP+TN+FP+FN>0,(TP+TN)/(TP+TN+FP+FN),NA)
  #myprauc[i] <- 
  #myrocauc[i] <- 
  mytpr[i] <- ifelse(TP+FN>0,TP/(TP+FN),NA)
  mytnr[i] <- ifelse(TN+FP>0,TN/(TN+FP),NA)
  myppv[i] <- ifelse(TP+FP>0,TP/(TP+FP),NA)
  mynpv[i] <- ifelse(TN+FN>0,TN/(TN+FN),NA)
}  
metrics[metrics$measure=="MCC","mean"] <- round(mean(mymcc,na.rm=TRUE),3)
metrics[metrics$measure=="MCC","sd"] <- round(sd(mymcc,na.rm=TRUE),3)
metrics[metrics$measure=="F1","mean"] <- round(mean(myf1,na.rm=TRUE),3)
metrics[metrics$measure=="F1","sd"] <- round(sd(myf1,na.rm=TRUE),3)
metrics[metrics$measure=="ACC","mean"] <- round(mean(myacc,na.rm=TRUE),3)
metrics[metrics$measure=="ACC","sd"] <- round(sd(myacc,na.rm=TRUE),3)
metrics[metrics$measure=="TPR","mean"] <- round(mean(mytpr,na.rm=TRUE),3)
metrics[metrics$measure=="TPR","sd"] <- round(sd(mytpr,na.rm=TRUE),3)
metrics[metrics$measure=="TNR","mean"] <- round(mean(mytnr,na.rm=TRUE),3)
metrics[metrics$measure=="TNR","sd"] <- round(sd(mytnr,na.rm=TRUE),3)
metrics[metrics$measure=="PPV","mean"] <- round(mean(myppv,na.rm=TRUE),3)
metrics[metrics$measure=="PPV","sd"] <- round(sd(myppv,na.rm=TRUE),3)
metrics[metrics$measure=="NPV","mean"] <- round(mean(mynpv,na.rm=TRUE),3)
metrics[metrics$measure=="NPV","sd"] <- round(sd(mynpv,na.rm=TRUE),3)


return(metrics)
}

## Datasets

arterial_event <- read_csv("../data/journal.pone.0201991.s001_EDITED_binary_event.csv")
## 30 arterial event (1) + 60 controls (0)
## 11 boolean features
## class: TARGET_arterial_event

event_type <- read_csv("../data/journal.pone.0201991.s001_EDITED_event_type_binary.csv")

## 22 ACS (0) + 8 Stroke (1)
## 11 boolean features
## class: TARGET_type_0ACS_1stroke

n_preds <- 11

# Tasks
#  A1) dimostrare che: il machine learning puo' predire sani e malati
#  A2) analizzare le caratteristiche piu' importanti nella distinzione tra malati e sani generate dal machine learning
#  B1) dimostrare che: il machine learning puo' predire infarto e sindrome coronarica acuta.
#  B2) analizzare le caratteristiche piu' importanti nella distinzione tra infarto e sindrome coronarica acuta generate dal machine learning
# Dovresti occuparti delle fasi A1 e B1, facendo classificazioni binarie su entrambi i datasets.
# Io ho gia' provato a fare le due classificazioni ma mi vengono risultati non eccellenti 
# MCC = +0.227 per il primo dataset e MCC =+0.121 per il secondo dataset con Random Forests).

################## UMAP  Analysis ###################################

UMAP_ANALYSIS <- FALSE

if(UMAP_ANALYSIS == TRUE){

    custom.config <- umap.defaults
    custom.config$n_neighbors=10
    custom.config$min_dist <- 0.01
    tmp <- umap(d=arterial_event[,1:n_preds],config = custom.config)
    png("Arterial_Event_UMAP.png")
    plot(tmp$layout,
	col=arterial_event$TARGET_arterial_event+3,
	pch=19, 
	asp=1,
	main="Arterial Event", 
	xlab="UMAP 1st feature",
	ylab="UMAP 2nd feature")
    legend(-4, 2, legend=c("Control", "Arterial Event"), col=c(3,4), pch=19, cex=0.8)
    dev.off()


    custom.config <- umap.defaults
    custom.config$n_neighbors=20
    custom.config$min_dist <- .01
    tmp <- umap(d=event_type[,1:n_preds],config = custom.config)
    png("Event_Type_UMAP.png")
    plot(tmp$layout,
	col=event_type$TARGET_type_0ACS_1stroke+3,
	pch=19, 
	asp=1,
	main="Event Type", 
	xlab="UMAP 1st feature",
	ylab="UMAP 2nd feature")
    legend(-2, 1.5, legend=c("ACS", "Stroke"), col=c(3,4), pch=19, cex=0.8)
    dev.off()

}

################# Machine Learning #############################

####### setting up #################
N <- 100      # num of exps
prop_tr_ts <- 80  #  training / test
event_type$TARGET_type_0ACS_1stroke <- factor(event_type$TARGET_type_0ACS_1stroke)
arterial_event$TARGET_arterial_event <- factor(arterial_event$TARGET_arterial_event)


ae_cl_one <- which(arterial_event$TARGET_arterial_event==1,arr.ind = TRUE)
ae_cl_zero <- which(arterial_event$TARGET_arterial_event==0,arr.ind = TRUE)
n_ae_ts_one <- floor((100-prop_tr_ts)*length(ae_cl_one)/100)
n_ae_ts_zero <- floor((100-prop_tr_ts)*length(ae_cl_zero)/100)

et_cl_one <- which(event_type$TARGET_type_0ACS_1stroke==1,arr.ind = TRUE)
et_cl_zero <- which(event_type$TARGET_type_0ACS_1stroke==0,arr.ind = TRUE)
n_et_ts_one <- floor((100-prop_tr_ts)*length(et_cl_one)/100)
n_et_ts_zero <- floor((100-prop_tr_ts)*length(et_cl_zero)/100)




## RF on the whole Arterial_Event - 100 runs to consider different seed
confmats <- list()
for(j in 1:N) confmats[[j]] <- randomForest(TARGET_arterial_event ~ ., data=arterial_event)$confusion[1:2,1:2]

cat("100 runs, arterial event, Random Forests: 
")
m_std(confmats)
#measure  mean    sd
#     MCC 0.240 0.031
#      F1 0.789 0.009
#     ACC 0.690 0.012
#     TPR 0.723 0.008
#     TNR 0.559 0.031
#     PPV 0.867 0.015
#     NPV 0.337 0.025

## RF on the whole Event_Type - 100 runs to consider different seed
confmats <- list()
for(j in 1:N) confmats[[j]] <- randomForest(TARGET_type_0ACS_1stroke ~ ., 
                                            data=event_type)$confusion[1:2,1:2]

cat("100 runs, event type, Random Forests: 
")
m_std(confmats)
#measure  mean    sd
#     MCC 0.046 0.073
#      F1 0.798 0.009
#     ACC 0.679 0.017
#     TPR 0.743 0.014
#     TNR 0.310 0.073
#     PPV 0.862 0.009
#     NPV 0.178 0.062

## Gaussian SVM on the whole Arterial_Event - 100 runs to consider different seed
confmats <- list()
for(j in 1:N){
  out <- ksvm(TARGET_arterial_event ~ ., data=arterial_event, type="C-svc",C=100)
  confmats[[j]] <- matrix(confusionMatrix(data=predict(out,arterial_event[,1:n_preds]),
                                   reference=arterial_event$TARGET_arterial_event)$table,2,2)
}
cat("100 runs, arterial event, Gaussian SVM: 
")
m_std(confmats)
#measure  mean    sd
#     MCC 0.900 0.000
#      F1 0.967 0.000
#     ACC 0.956 0.000
#     TPR 0.984 0.003
#     TNR 0.899 0.007
#     PPV 0.951 0.003
#     NPV 0.966 0.007

## Linear SVM on the whole Event_Type - 100 runs to consider different seed
confmats <- list()
for(j in 1:N){
  out <- ksvm(TARGET_type_0ACS_1stroke ~ ., data=event_type, type="C-svc",C=1000, kernel="vanilladot")
  confmats[[j]] <- matrix(confusionMatrix(data=predict(out,event_type[,1:n_preds]),
                                          reference=event_type$TARGET_type_0ACS_1stroke)$table,2,2)
}

cat("100 runs, event type, linear SVM: 
")
m_std(confmats)
#measure  mean sd
#     MCC 0.830  0
#      F1 0.955  0
#     ACC 0.933  0
#     TPR 0.955  0
#     TNR 0.875  0
#     PPV 0.955  0
#     NPV 0.875  0



## RF on the Arterial_Event - 100 runs x  MonteCarlo stratified splits
confmats <- list()
for(j in 1:N){
TS_O <- sample(ae_cl_one)[1:n_ae_ts_one]
TS_Z <- sample(ae_cl_zero)[1:n_ae_ts_zero]
TR_O <- setdiff(ae_cl_one,TS_O)
TR_Z <- setdiff(ae_cl_zero,TS_Z)

TR <- sample(c(TR_O,TR_Z))
TS <- sample(c(TS_O,TS_Z))


out <- randomForest(TARGET_arterial_event ~ ., 
             data=arterial_event[TR,], 
             xtest=arterial_event[TS,1:n_preds],
             ytest=arterial_event$TARGET_arterial_event[TS]
             )
confmats[[j]] <- out$test$confusion[1:2,1:2]
}

cat("MonteCarlo, arterial event, Random Forests: 
")
m_std(confmats)
#measure  mean    sd
#     MCC 0.257 0.224
#      F1 0.787 0.061
#     ACC 0.692 0.083
#     TPR 0.732 0.065
#     TNR 0.593 0.244
#     PPV 0.861 0.103
#     NPV 0.355 0.199


## RF on the Event_Type - 100 runs x  MonteCarlo stratified splits
confmats <- list()
for(j in 1:N){
  TS_O <- sample(et_cl_one)[1:n_et_ts_one]
  TS_Z <- sample(et_cl_zero)[1:n_et_ts_zero]
  TR_O <- setdiff(et_cl_one,TS_O)
  TR_Z <- setdiff(et_cl_zero,TS_Z)
  
  TR <- sample(c(TR_O,TR_Z))
  TS <- sample(c(TS_O,TS_Z))
  
  
    out <- randomForest(TARGET_type_0ACS_1stroke ~ ., 
                      data=event_type[TR,], 
                      xtest=event_type[TS,1:n_preds],
                      ytest=event_type$TARGET_type_0ACS_1stroke[TS]
  )
  confmats[[j]] <- out$test$confusion[1:2,1:2]
}

cat("MonteCarlo, event type, Random Forests: 
")
m_std(confmats)
#measure  mean    sd
#     MCC 0.244 0.512
#      F1 0.854 0.128
#     ACC 0.774 0.177
#     TPR 0.857 0.117
#     TNR 0.473 0.457
#     PPV 0.872 0.176
#     NPV 0.380 0.488


## Gaussian SVM on the Arterial_Event - 100 runs x  MonteCarlo stratified splits
confmats <- list()
for(j in 1:N){
  TS_O <- sample(ae_cl_one)[1:n_ae_ts_one]
  TS_Z <- sample(ae_cl_zero)[1:n_ae_ts_zero]
  TR_O <- setdiff(ae_cl_one,TS_O)
  TR_Z <- setdiff(ae_cl_zero,TS_Z)
  
  TR <- sample(c(TR_O,TR_Z))
  TS <- sample(c(TS_O,TS_Z))
  
  out <- ksvm(TARGET_arterial_event ~ ., data=arterial_event[TR,], type="C-svc", C=100)
  confmats[[j]] <- matrix(confusionMatrix(data=predict(out,arterial_event[TS,1:n_preds]),
                                          reference=arterial_event$TARGET_arterial_event[TS])$table,2,2)
}
cat("MonteCarlo, arterial event, Gaussian SVM: 
")
m_std(confmats)
#measure  mean    sd
#     MCC 0.286 0.232
#      F1 0.771 0.080
#     ACC 0.689 0.094
#     TPR 0.801 0.129
#     TNR 0.465 0.207
#     PPV 0.755 0.072
#     NPV 0.560 0.225


## Linear SVM on the Event_Type - 100 runs x  MonteCarlo stratified splits
confmats <- list()
for(j in 1:N){
  TS_O <- sample(et_cl_one)[1:n_et_ts_one]
  TS_Z <- sample(et_cl_zero)[1:n_et_ts_zero]
  TR_O <- setdiff(et_cl_one,TS_O)
  TR_Z <- setdiff(et_cl_zero,TS_Z)
  
  TR <- sample(c(TR_O,TR_Z))
  TS <- sample(c(TS_O,TS_Z))
  
  out <- ksvm(TARGET_type_0ACS_1stroke ~ ., data=event_type[TR,], type="C-svc", 
              C=1000, kernel="vanilladot")
  confmats[[j]] <- matrix(confusionMatrix(data=predict(out,event_type[TS,1:n_preds]),
                                                reference=event_type$TARGET_type_0ACS_1stroke[TS])$table,2,2)
}
cat("MonteCarlo, event type, linear SVM: 
")
m_std(confmats)
#measure  mean    sd
#     MCC 0.070 0.420
#      F1 0.764 0.153
#     ACC 0.664 0.175
#     TPR 0.745 0.219
#     TNR 0.340 0.476
#     PPV 0.834 0.130
#     NPV 0.231 0.314

