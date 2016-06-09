library(limma)
library(marray)
library(arrayQuality)     #opcja wgrAnia pakietu: source("http://www.bioconductor.org/biocLite.R"), a pozniej  biocLite("arrayQuality")
library(e1071)
library(vsn)      # stworzony dla wersji nizszych od wersji 2.10.0
library(dyebias)
library(convert)
library(sda)
library(samr)
library(st)
require(party)
require(rpart)
require(ipred)
library(pROC)
library(ROCR)
library(MLInterfaces)
library(randomForest)
setwd('/Users/joannazyprych/Documents/obliczenia do prac naukowych/luiza praca o normalizacjach RNA-seq/wyniki do publikacji')
load('dane_potrzebne_do_wykresow.RData' )
bledy_klasyfikacji<-list()

for (zestaw in c("Cheung","Bodymap","Leukemia"))
{
  
  if (zestaw=='Cheung') {
    f=30
    data<-DANE$Cheung$counts
    L<-DANE$Cheung$L
    lista_auc<-list(TMM=wyniki$Cheung$TMM$kolejnosc[1:f],UQ=wyniki$Cheung$UQ$kolejnosc[1:f],DES=wyniki$Cheung$DES$kolejnosc[1:f],EBS=wyniki$Cheung$EBS$kolejnosc[1:f],PS=wyniki$Cheung$PS$kolejnosc[1:f],RD=wyniki$Cheung$RD$kolejnosc[1:f])
    index<-lista_auc
  }
  if (zestaw=='Bodymap') {
    f=12
    data<-DANE$Bodymap$counts
    L<-DANE$Bodymap$L
    lista_auc<-list(TMM=wyniki$Bodymap$TMM$kolejnosc[1:f],UQ=wyniki$Bodymap$UQ$kolejnosc[1:f],DES=wyniki$Bodymap$DES$kolejnosc[1:f],EBS=wyniki$Bodymap$EBS$kolejnosc[1:f],PS=wyniki$Bodymap$PS$kolejnosc[1:f],RD=wyniki$Bodymap$RD$kolejnosc[1:f])
    index<-lista_auc
  }
  if (zestaw=='Leukemia') {
    f=20
    data<-DANE$Leukemia$counts
    L<-DANE$Leukemia$L
    lista_auc<-list(TMM=wyniki$Leukemia$TMM$kolejnosc[1:f],UQ=wyniki$Leukemia$UQ$kolejnosc[1:f],DES=wyniki$Leukemia$DES$kolejnosc[1:f],EBS=wyniki$Leukemia$EBS$kolejnosc[1:f],PS=wyniki$Leukemia$PS$kolejnosc[1:f],RD=wyniki$Leukemia$RD$kolejnosc[1:f])
    index<-lista_auc
  }
mis_classified1<-numeric()
mis_classified2<-numeric()
mis_classified3<-numeric()
mis_classified4<-numeric()
mis_classified5<-numeric()
for(r in 1: length(index)){
  cat(paste(r,'... \n'))
  dane<-data
  dane_k<-dane[index[[r]],]                 #indeksy genow roznicujacych to index 
  dane_k<-t(dane_k) 

  dane<-data.frame(dane_k,L)
  # Naive Bayes
  

  True<-as.numeric(L)-1
  poprawne<-0
  for (i in 1:(dim(dane)[1])){
    training_set = c(1:(dim(dane)[1]))[-i]
    klasyfikator = MLearn(L~ ., data = dane, naiveBayesI, training_set )
    if(testPredictions(klasyfikator)==True[i])
      poprawne<-poprawne+1;
  }
  mis_classified1[r]<-(dim(dane)[1])-poprawne
  
  # knn
  
  poprawne<-0
  for (i in 1:(dim(dane)[1])){
    training_set= c(1:(dim(dane)[1]))[-i]
    klasyfikator = MLearn(L~ ., data=dane, knnI(k=3,l=0), training_set)
    if(testPredictions(klasyfikator)==True[i])
      poprawne<-poprawne+1;
  }
  mis_classified2[r]<-(dim(dane)[1])-poprawne
  
  # SVM
  
  poprawne<-0
  for (i in 1:(dim(dane)[1])){
    training_set= c(1:(dim(dane)[1]))[-i]
    klasyfikator = MLearn(L~ ., data=dane, svmI, training_set)
    if(testPredictions(klasyfikator)==True[i])
      poprawne<-poprawne+1;
  }
  mis_classified3[r]<-(dim(dane)[1])-poprawne
  
  # nnet - Fit single-hidden-layer neural network, possibly with skip-layer connections
  
  mis_classified4a<-numeric()
  for(q in 1:6){
    poprawne<-0
    for (i in 1:(dim(dane)[1])){
      training_set= c(1:(dim(dane)[1]))[-i]
      klasyfikator = MLearn(L~ ., data=dane, nnetI, training_set, size=3, decay=.01 )  #, size = 5, decay = 0.01, MaxNWts = 2000)   #
      if(testPredictions(klasyfikator)==True[i])
        poprawne<-poprawne+1;
    }
    mis_classified4a[q]<-(dim(dane)[1])-poprawne
  }
  mis_classified4[r]<-mean(mis_classified4a)
  
  #random forest
  
  mis_classified5a<-numeric()
  for(q in 1:6){
    poprawne<-0
    for (i in 1:(dim(dane)[1])){
      training_set= c(1:(dim(dane)[1]))[-i]
      klasyfikator = MLearn(L~ ., data=dane, randomForestI,training_set)
      if(testPredictions(klasyfikator)==True[i])
        poprawne<-poprawne+1;
    }
    mis_classified5a[q]<-((dim(dane)[1]))-poprawne
  }
  mis_classified5[r]<-mean(mis_classified5a)
}
bledy<-matrix(c(mis_classified1,mis_classified2,mis_classified3,mis_classified4,mis_classified5),5,length(index),byrow=T)
  colnames(bledy)<-names(index)

rownames(bledy)<-c('NB','knn','svm','nnet','rf')
bledy_klasyfikacji[[zestaw]]<-bledy
write.csv(bledy_klasyfikacji[[zestaw]],file=paste('bledy_klasyfikacji',zestaw,'.csv'))
}
bledy_klasyfikacji$Cheung<-bledy_klasyfikacji$Cheung/41
bledy_klasyfikacji$Bodymap<-bledy_klasyfikacji$Bodymap/16
bledy_klasyfikacji$Leukemia<-bledy_klasyfikacji$Leukemia/27
lapply(bledy_klasyfikacji,function(xcolMeans(x)*100)
