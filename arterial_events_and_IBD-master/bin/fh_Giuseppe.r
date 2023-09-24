setwd(".")


list.of.packages <- c("easypackages", "PRROC", "e1071", "randomForest","class", "gmodels", "formula.tools", "dplyr", "pastecs", "ROSE", "readr", "randomForest", "caret", "e1071", "mltools", "boot", "umap", "kernlab")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library("easypackages")
libraries(list.of.packages)

mean.fun <- function(d, i)
{    m <- mean(d$data[i])
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
  names(out) <- c("lower","Mean","upper")
  return(out)
}


#A) file: "journal.pone.0201991.s001_EDITED_binary_event.csv"
#target: "TARGET_arterial_event"
#numero caratteristiche: 12
#numero soggetti: 90
#classi: 30 pazienti (classe 1) e 60 sani (classe 0)

a_data <- read_csv("../data/journal.pone.0201991.s001_EDITED_binary_event.csv")
head(a_data)
colnames(a_data)

#B) file: "journal.pone.0201991.s001_EDITED_event_type_binary.csv"
#target: "TARGET_type_0ACS_1stroke"
#numero soggetti: 30
#numero caratteristiche: 12
#classi: 22 pazienti con sindrome coronarica acuta (ACS, classe 0) e 8
#pazienti con infarto (stroke, classe 1)

b_data <- read_csv("../data/journal.pone.0201991.s001_EDITED_event_type_binary.csv")
head(b_data)
colnames(b_data)

n_preds <- 11

#In questo progetto vorrei:
#  A1) dimostrare che: il machine learning puo' predire sani e malati
#  A2) analizzare le caratteristiche piu' importanti nella distinzione tra malati e sani generate dal machine learning
#  B1) dimostrare che: il machine learning puo' predire infarto e sindrome coronarica acuta.
#  B2) analizzare le caratteristiche piu' importanti nella distinzione tra infarto e sindrome coronarica acuta generate dal machine learning

#Dovresti occuparti delle fasi A1 e B1, facendo classificazioni binarie su entrambi i datasets.
I# o ho gia' provato a fare le due classificazioni ma mi vengono risultati non eccellenti 
# MCC = +0.227 per il primo dataset e MCC =+0.121 per il secondo dataset con Random Forests).

################## UMAP ###################################

z <- umap(a_data[,1:n_preds])
plot(z$layout,col=a_data$TARGET_arterial_event+2,pch=19, main="File A")

w <- umap(b_data[,1:n_preds])
plot(w$layout,col=b_data$TARGET_type_0ACS_1stroke+2,pch=19, main="File B")

################# RANDOM FOREST #############################

## on the whole dataset - repeat 100 times to allow for different seeds
N <- 100
MCC <- rep(0,N)
for(j in 1:N){
a_data$TARGET_arterial_event <- factor(a_data$TARGET_arterial_event)
out <- randomForest(TARGET_arterial_event ~ ., data=a_data, )
MCC[j] <- mcc(confusionM = out$confusion[1:2,1:2])
}
#summary(MCC)

cat("[Event] Repeat 100 times Random Forest mean MCC = ")
cat(summary(MCC)[c("Mean")][[1]])
cat("\n")

##   Min.  1st Qu.  Median Mean    3rd Qu.    Max. 
## 0.1768  0.2118  0.2357  0.2349  0.2609  0.2946 
ci(MCC)



#   lower      mean     upper 
#0.2289506 0.2349335 0.2410874 

MCC <- rep(0,N)
for(j in 1:N){
b_data$TARGET_type_0ACS_1stroke<- factor(b_data$TARGET_type_0ACS_1stroke)
out <- randomForest(TARGET_type_0ACS_1stroke ~ ., data=b_data)
MCC[j] <- mcc(confusionM = out$confusion[1:2,1:2])
}
(cat("[Type]  Repeat 100 times Random Forest mean MCC = "))
cat(summary(MCC)[c("Mean")][[1]], "
")

(cat("\n"))

##    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##-0.06742 -0.01478 -0.01478  0.05340  0.13484  0.26382 
ci(MCC)
##  lower       mean      upper 
##0.03918573 0.05339849 0.06892570 


### try stratified MonteCarlo Split
a_data$TARGET_arterial_event <- factor(a_data$TARGET_arterial_event)
cl_one <- which(a_data$TARGET_arterial_event==1,arr.ind = TRUE)
cl_zero <- which(a_data$TARGET_arterial_event==0,arr.ind = TRUE)
prop_tr_ts <- 80
n_ts_one <- floor((100-prop_tr_ts)*length(cl_one)/100)
n_ts_zero <- floor((100-prop_tr_ts)*length(cl_zero)/100)

N <- 100
MCC <- rep(0,N)
for(j in 1:N){
TS_O <- sample(cl_one)[1:n_ts_one]
TS_Z <- sample(cl_zero)[1:n_ts_zero]
TR_O <- setdiff(cl_one,TS_O)
TR_Z <- setdiff(cl_zero,TS_Z)

TR <- sample(c(TR_O,TR_Z))
TS <- sample(c(TS_O,TS_Z))


out <- randomForest(TARGET_arterial_event ~ ., 
             data=a_data[TR,], 
             xtest=a_data[TS,1:n_preds],
             ytest=a_data$TARGET_arterial_event[TS]
             )
MCC[j] <- mcc(confusionM = out$test$confusion[1:2,1:2])
}

(cat("[Event] MonteCarlo Random Forest mean MCC = "))
cat(summary(MCC)[c("Mean")][[1]], "
")

(cat("\n"))



##  Min.   1st Qu.  Median  Mean   3rd Qu.  Max. 
##-0.3162  0.1250  0.3162  0.2740  0.4725  0.8771 

ci(MCC)
##  lower      mean     upper 
##0.2241028 0.2740252 0.3267426



b_data$TARGET_type_0ACS_1stroke <- factor(b_data$TARGET_type_0ACS_1stroke)
cl_one <- which(b_data$TARGET_type_0ACS_1stroke==1,arr.ind = TRUE)
cl_zero <- which(b_data$TARGET_type_0ACS_1stroke==0,arr.ind = TRUE)
prop_tr_ts <- 80
n_ts_one <- floor((100-prop_tr_ts)*length(cl_one)/100)
n_ts_zero <- floor((100-prop_tr_ts)*length(cl_zero)/100)

N <- 100
MCC <- rep(0,N)
for(j in 1:N){
  TS_O <- sample(cl_one)[1:n_ts_one]
  TS_Z <- sample(cl_zero)[1:n_ts_zero]
  TR_O <- setdiff(cl_one,TS_O)
  TR_Z <- setdiff(cl_zero,TS_Z)
  
  TR <- sample(c(TR_O,TR_Z))
  TS <- sample(c(TS_O,TS_Z))
  
  
  out <- randomForest(TARGET_type_0ACS_1stroke ~ ., 
                      data=b_data[TR,], 
                      xtest=b_data[TS,1:n_preds],
                      ytest=b_data$TARGET_type_0ACS_1stroke[TS]
  )
  MCC[j] <- mcc(confusionM = out$test$confusion[1:2,1:2])
}

(cat("[Type] MonteCarlo Random Forest mean MCC = "))
cat(summary(MCC)[c("Mean")][[1]], "
")

(cat("\n"))

##  Min.   1st Qu.  Median  Mean   3rd Qu.  Max. 
##-0.4082 -0.2500  0.0000  0.1617  0.6124  1.0000 
ci(MCC)
##  lower      mean     upper 
##0.06809962 0.16169600 0.26522236



############# SVM ###########################

## on the whole dataset - repeat 100 times to allow for different seeds
N <- 100
MCC <- rep(0,N)
for(j in 1:N){
  a_data$TARGET_arterial_event <- factor(a_data$TARGET_arterial_event)
  out <- ksvm(TARGET_arterial_event ~ ., data=a_data, type="C-svc",C=100)
  MCC[j] <- mcc(pred=predict(out,a_data[,1:11]),actual=a_data$TARGET_arterial_event)
}

(cat("[Event] Repeat 100 times kernel SVM mean MCC = "))
cat(summary(MCC)[c("Mean")][[1]], "
")

(cat("\n"))

##   Min.  1st Qu.  Median Mean    3rd Qu.    Max. 
##  0.8995  0.8995  0.8995  0.8996  0.8995  0.9014 
ci(MCC)
#   lower      mean     upper 
# 0.8995109 0.8995838 0.8997140

MCC <- rep(0,N)
for(j in 1:N){
  b_data$TARGET_type_0ACS_1stroke<- factor(b_data$TARGET_type_0ACS_1stroke)
  out <- ksvm(TARGET_type_0ACS_1stroke ~ ., data=b_data, type="C-svc",C=1000, kernel="vanilladot")
  MCC[j] <- mcc(pred=predict(out,a_data[,1:11]),actual=b_data$TARGET_type_0ACS_1stroke)
}

(cat("[Type] Repeat 100 times kernel SVM mean MCC = "))
cat(summary(MCC)[c("Mean")][[1]], "
")

(cat("\n"))

##    Min.    1st Qu.   Median     Mean  3rd Qu.     Max. 
##   0.05913 0.05913 0.05913 0.05913 0.05913 0.05913
ci(MCC)
##  lower       mean      upper 
##0.05913124 0.05913124 0.05913124


### try stratified MonteCarlo Split
a_data$TARGET_arterial_event <- factor(a_data$TARGET_arterial_event)
cl_one <- which(a_data$TARGET_arterial_event==1,arr.ind = TRUE)
cl_zero <- which(a_data$TARGET_arterial_event==0,arr.ind = TRUE)
prop_tr_ts <- 80
n_ts_one <- floor((100-prop_tr_ts)*length(cl_one)/100)
n_ts_zero <- floor((100-prop_tr_ts)*length(cl_zero)/100)

N <- 100
MCC <- rep(0,N)
for(j in 1:N){
  TS_O <- sample(cl_one)[1:n_ts_one]
  TS_Z <- sample(cl_zero)[1:n_ts_zero]
  TR_O <- setdiff(cl_one,TS_O)
  TR_Z <- setdiff(cl_zero,TS_Z)
  
  TR <- sample(c(TR_O,TR_Z))
  TS <- sample(c(TS_O,TS_Z))
  
  
  out <- ksvm(TARGET_arterial_event ~ ., data=a_data[TR,], type="C-svc", C=100)
  MCC[j] <- mcc(preds = predict(out,a_data[TS,1:n_preds]),actuals = a_data$TARGET_arterial_event[TS])      
}

(cat("[Event] MonteCarlo kernel SVM mean MCC = "))
cat(summary(MCC)[c("Mean")][[1]], "
")

(cat("\n"))

##  Min.   1st Qu.  Median  Mean   3rd Qu.  Max. 
##-0.1581  0.1157  0.3162  0.2976  0.5000  0.6447 

ci(MCC)
##  lower      mean     upper 
##0.2522547 0.2976431 0.3456693 



b_data$TARGET_type_0ACS_1stroke <- factor(b_data$TARGET_type_0ACS_1stroke)
cl_one <- which(b_data$TARGET_type_0ACS_1stroke==1,arr.ind = TRUE)
cl_zero <- which(b_data$TARGET_type_0ACS_1stroke==0,arr.ind = TRUE)
prop_tr_ts <- 80
n_ts_one <- floor((100-prop_tr_ts)*length(cl_one)/100)
n_ts_zero <- floor((100-prop_tr_ts)*length(cl_zero)/100)

N <- 100
MCC <- rep(0,N)
for(j in 1:N){
  TS_O <- sample(cl_one)[1:n_ts_one]
  TS_Z <- sample(cl_zero)[1:n_ts_zero]
  TR_O <- setdiff(cl_one,TS_O)
  TR_Z <- setdiff(cl_zero,TS_Z)
  
  TR <- sample(c(TR_O,TR_Z))
  TS <- sample(c(TS_O,TS_Z))
  
  
  out <- ksvm(TARGET_type_0ACS_1stroke ~ ., data=b_data[TR,], type="C-svc", C=1000, kernel="vanilladot")
  MCC[j] <- mcc(preds = predict(out,b_data[TS,1:n_preds]),actuals = b_data$TARGET_type_0ACS_1stroke[TS])      
  
}

(cat("[Type] MonteCarlo kernel SVM mean MCC = "))
cat(summary(MCC)[c("Mean")][[1]], "
")

(cat("\n"))

##  Min.   1st Qu.  Median  Mean   3rd Qu.  Max. 
##-0.6124 -0.2500  0.4082  0.2456  0.6124  1.0000
ci(MCC)
##  lower      mean     upper 
##0.1554203 0.2455931 0.3392529



