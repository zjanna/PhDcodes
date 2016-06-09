 # Kod wykorzystany do obliczen prezentowanych w pracy doktorskiej Joanny Zyprych-Walczak
 ###########################################################################################################
# wczytanie pakietów potrzebnych do programu
# jesli jakis pakiet wczesniej nie zostal zainstalowany to nalezy przed jego wczytaniem uzyc polecenia -- np. install.packages('limma')   

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

# 1. Funkcja wczytania danych  z plikow gpr

# setwd('F:\\program do doktoratu\\dane')         okreslenie sciezki miejsca pracy - tutaj beda zapisywane wszystkie wykresy

setwd('F:\\program_do_doktoratu')

wczytanie_danych<-function(targets,spottypes){        #targets -'nazwa pliku', spottypes -'nazwa pliku'
    targets<- readTargets(targets)
    RG<-read.maimages(targets,source="genepix",columns=list(R="F647 Mean",G="F555 Mean",Rb="B647 Median",Gb="B555 Median"),annotation=c("Block","Row","Column","ID","Name"),wt.fun=wtflags(weight=0, cutoff=-50))
    spottypes<-readSpotTypes(spottypes)
    wyniki<-list(targets=targets,RG=RG,spottypes=spottypes)
    attach(wyniki)
                                             }
#wywolanie funkcji wczytania danych  z plikow gpr dla przykladowych plików

wczytanie_danych('przykladowy_targets.txt','Bialaczki2_ test_SpotTypes.txt')
RG$genes$Status=controlStatus(spottypes, RG)
MA<-normalizeWithinArrays(RG, method="none")          #zamiana danych RGlist na MAlist

###################################################################################################################

#2. Analiza obrazu, korekta tla, normalizacja, uzupelnianie brakujacych danych
#wczytanie danych jako obiekt marray konieczne do dalszych wykresow
#ngr -	liczba bloków w wierszach na calej mikromacierzy
#ngc -  liczba bloków w kolumnach na calej mikromacierzy
#nsr -	liczba wierszy w blokach
#nsc - 	liczba kolumn w blokach
#skip 	liczba wierszy w pliku z danymi, które nie zawieraja bezposrednio danych i ktore nalezy przy wczytywaniu ominac

obiekt_marray<-function(nazwa_targets,ngr,ngc,nsr,nsc,skip){
  sample.file <-  file.path(nazwa_targets)
  targets <- read.marrayInfo(sample.file)
  layout <- read.marrayLayout(fname=NULL, ngr=ngr, ngc=ngc, nsr=nsr, nsc=nsc)
  n.spots <- maNspots(layout)
  first.file <- file.path(maInfo(targets)$FileName[1])
  gnames <- read.marrayInfo(first.file, info.id=c(4,5), labels=4, skip=skip)
  names(maInfo(gnames)) <- c("Name", "ID")
  maInfo(gnames)$ID <- as.character(maInfo(gnames)$ID)
  controls <- targets@maInfo[,3]
  maControls(layout) <- as.factor(controls)
  data.raw <- read.GenePix(fnames = maInfo(targets)$FileName,
                          name.Gf ="F555 Mean",
                          name.Gb ="B555 Median",
                          name.Rf ="F647 Mean",
                          name.Rb ="B647 Median",
                          name.W ="Flags",
                          layout = layout,
                          gnames = gnames,
                          targets = targets,
                          skip=skip, sep = "\t", quote = '"', DEBUG=TRUE)
  return(data.raw=data.raw)
                          }
#wywolanie funkcji

surowe.dane<-obiekt_marray('przykladowy_targets.txt',10,4,9,9,31)

#wykresy intensywnosci - typu: plot densities - Rysunek 4.1
nazwa<-paste("wykresy_intensywnosci",".jpeg",sep="")
jpeg(nazwa)
par(mfrow=c(3,1))
plotDensities(RG,log=TRUE)   #wykres dla obu kanalów jednoczesnie
plotDensities(RG,log=TRUE,singlechannels=((ncol(RG)+1):(ncol(RG)+ncol(RG)) ))  #wykres dla kanalu zielonego
plotDensities(RG,log=TRUE,singlechannels= (1:ncol(RG) ))                #wykres dla kanalu zielonego
dev.off()

#zestaw wykresów diagnostycznych dla kazdej macierzy osobno - Rysunek 4.2

for(i in 1:dim(RG)[2]){
maQualityPlots(RG[,i])
}

#boxploty dla danych przed normalizacj¹    - Rysunek 4.4

nazwa<-paste("boxploty_przed_normalizacja",".jpeg",sep="")
jpeg(nazwa)
boxplot(MA$M~col(MA$M),names=colnames(MA$M), col=rainbow(12))
dev.off()

# ten sam rysunek tylko bez zapisu do pliku

boxplot(MA$M~col(MA$M),names=colnames(MA$M), col=rainbow(12))

# MA plots - kady plik zawiera 6 MA wykresów w 3 wierszach i 2 kolumnach  - Rysunek 4.5

plotMA3by2(RG, prefix="MAplots", device="jpeg")

# obrazy t³a dla koloru czerwonego i zielonego   - Rysunek 4.6

for(i in 1:dim(surowe.dane)[2]){
nazwa<-paste("obraz_t³a_macierz_",i,".jpeg",sep="")
jpeg(nazwa)
par(mfrow=c(1,2))
imageplot(log2(RG$Rb[,i]), RG$printer, low="white", legend=TRUE,mar=c(2,2,4,4),ncolors = 123,high="red",main=c("obraz_t³a_macierz",i,"\n"))
imageplot(log2(RG$Gb[,i]), RG$printer, low="white", legend=TRUE,mar=c(2,2,4,4),ncolors = 123,high="green",main=c("obraz_t³a_macierz",i,"\n"))
dev.off()
}

#wykresy rozrzutu - Rysunek 4.7

dane<-MA$M

# uzupelnianie brakujacych danych dla surowych danych w celu wyznaczenia macierzy korelacji, metoda median, mean

uzupelnianie_danych<-function(metoda,dane){
  if(sum(is.na(dane))!=0){
    switch(metoda,
    'Median'=(dane<-impute(dane,'median')),
    'Mean'  =(dane<-impute(dane,'mean')))
                      }
                      }

#przykladowe wywolanie funkcji 

dane<- uzupelnianie_danych('Mean',dane)
M<-dane
Mb<-log2(RG$Rb/RG$Gb)

#do korekty tla   - wykresy rozrzutu M vs Mb

for(i in 1:dim(RG)[2]){
      nazwa<-paste("wykres_rozrzutu_macierz_",i,".jpeg",sep="")
      jpeg(nazwa)
      plot(Mb[,i],M[,i], xlab='Mb', ylab='M', cex.lab = 1.5,cex.main = 1.5,main=c("macierz",i))
      dev.off()
			                 }
cor(Mb,M,method="spearman")

#macierz korelacji pomiedzy probami - heatmap

nazwa<-paste("heatmap_dla_wszystkich_mikromacierzy.jpeg",sep="")
jpeg(nazwa)
heatmap.2(cor(as.matrix(dane)),col=redgreen(75),  scale="none",
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1)
dev.off()

#ten sam rysunek tylko bez zapisu do pliku

heatmap.2(cor(as.matrix(dane)),col=redgreen(75),  scale="none",
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1)
          
###############################################################NORMALIZACJA########################

# normalizacja danych

normalizacja=TRUE      #jesli zachodzi potrzeba normalizacji - ustawic na TRUE

# normalizacja marray loess

if(normalizacja){
    data.normPLS=maNorm(surowe.dane, norm="loess")
    MAlsm<-data.normPLS@maM          #nazwa po loess z marray
    MAloessM_M<-data.normPLS@maM
    MAloessM_A<-data.normPLS@maA
    
    # normalizacja loess i printtiploess limma
    
    MAls<-normalizeWithinArrays(RG, method="loess")
    MAlsp<- normalizeWithinArrays(RG, method="printtiploess")
    
    # vsn  - opcja dla wersji ponizej R 2.10.0  przy zainstalowanym pakiecie 'vsn'
    
    data.norm=justvsn(RG, backgroundsubtract=FALSE)
    r=assayData(data.norm)$R
    g=assayData(data.norm)$G
    asd=assayData(data.norm)
    Avsn=(asd$R+asd$G)/2
    Mvsn=(asd$R-asd$G)
                          }
                          
# boxploty dla normalizacji - Rysunek 4.8

boxplot=TRUE    # jesli TRUE to wyrysowane zostana boxploty

if(boxplot){
       nazwa<-paste("boxploty dla poszczegolnych normalizacji.jpeg",sep="")
       jpeg(nazwa)
       par(mfrow=c(2,2))
       boxplot(MAls$M~col(MAls$M), col=rainbow(12),main="po normalizacji loess",xlab='mikromacierze', ylab='M', las=2)
  		 boxplot(MAlsm~col(MAlsm), col=rainbow(12),main="po normalizacji loess marray", xlab='mikromacierze', ylab='M',las=2)
		   boxplot(MAlsp$M~col(MAlsp$M), col=rainbow(12),main="po normalizacji printtiploess", xlab='mikromacierze', ylab='M',las=2)
		   boxplot(Mvsn~col(Mvsn), col=rainbow(12),main="po normalizacji vsn", xlab='mikromacierze', ylab='M',las=2)
       dev.off()
          }
          
#ten sam rysunek tylko bez zapisu do pliku

if(boxplot){
       par(mfrow=c(2,2))
       boxplot(MAls$M~col(MAls$M), col=rainbow(12),main="po normalizacji loess",xlab='mikromacierze', ylab='M', las=2)
  		 boxplot(MAlsm~col(MAlsm), col=rainbow(12),main="po normalizacji loess marray", xlab='mikromacierze', ylab='M',las=2)
		   boxplot(MAlsp$M~col(MAlsp$M), col=rainbow(12),main="po normalizacji printtiploess", xlab='mikromacierze', ylab='M',las=2)
		   boxplot(Mvsn~col(Mvsn), col=rainbow(12),main="po normalizacji vsn", xlab='mikromacierze', ylab='M',las=2)
          }
          
####################################PORÓWNANIA NORMALIZACJI################################
# porównanie normalizacji loess, printtiploess, loess marray, vsn -opcjonalnie

BiasLS<-numeric()         # blad dla metody loess
BiasLSP<-numeric()        # blad dla printtiploess
BiasloessM<-numeric()     # blad dla loess(marray)
BiasVSN<-numeric()        # blad dla vsn
VarianceLS<-numeric()      #wariancja
VarianceLSP<-numeric()
VarianceloessM<-numeric()
VarianceVSN<-numeric()
nazwy_kontroli<-unique(RG$genes$Status)
spajki<-list()
for(i in 1:(length(nazwy_kontroli))){
  spajki[[i]]<-which(RG$genes$Status==nazwy_kontroli[i])
  BiasLS[i]<-sqrt(((sum(MAls$M[spajki[[i]],]^2,na.rm=TRUE))/(length(spajki[[i]])*dim(RG)[2])))
  VarianceLS[i]<-(1/(length(spajki[[i]])*dim(RG)[2]-1))*sum((MAls$M[spajki[[i]],]-mean(MAls$M[spajki[[i]],],na.rm=TRUE))^2,na.rm=TRUE)
  BiasLSP[i]<-sqrt(((sum(MAlsp$M[spajki[[i]],]^2,na.rm=TRUE))/(length(spajki[[i]])*dim(RG)[2])))
  VarianceLSP[i]<-(1/(length(spajki[[i]])*dim(RG)[2]-1))*sum((MAlsp$M[spajki[[i]],]-mean(MAlsp$M[spajki[[i]],],na.rm=TRUE))^2,na.rm=TRUE)
  BiasloessM[i]<-sqrt(((sum(MAlsm[spajki[[i]],]^2,na.rm=TRUE))/(length(spajki[[i]])*dim(RG)[2])))
  VarianceloessM[i]<-(1/(length(spajki[[i]])*dim(RG)[2]-1))*sum((MAlsm[spajki[[i]],]-mean(MAlsm[spajki[[i]],],na.rm=TRUE))^2,na.rm=TRUE)
  BiasVSN[i]<-sqrt(((sum(Mvsn[spajki[[i]],]^2,na.rm=TRUE))/(length(spajki[[i]])*dim(RG)[2])))
  VarianceVSN[i]<-(1/(length(spajki[[i]])*dim(RG)[2]-1))*sum((Mvsn[spajki[[i]],]-mean(Mvsn[spajki[[i]],],na.rm=TRUE))^2,na.rm=TRUE)
                                    }
                                    
#opcja z vsn

bledy<-matrix(c(BiasLS,BiasLSP,BiasloessM,BiasVSN, VarianceLS,VarianceLSP,VarianceloessM,VarianceVSN),8,length(nazwy_kontroli),byrow=T)
rownames(bledy)<- c('BiasLOESS','BiasPRINTTIPLOESS','BiasLOESSmarray','BiasVSN', 'VarianceLOESS','VariancePRINTTIPLOESS','VarianceLOESSmarray','VarianceVSN')
colnames(bledy)<- nazwy_kontroli
bledy<-cbind(c('BiasLOESS','BiasPRINTTIPLOESS','BiasLOESSmarray','BiasVSN', 'VarianceLOESS','VariancePRINTTIPLOESS','VarianceLOESSmarray','VarianceVSN'),bledy)
write.xls(bledy,file='bledy_normalizacji_wedlug_kontroli.xls')

#opcja bez vsn

bledy<-matrix(c(BiasLS,BiasLSP,BiasloessM, VarianceLS,VarianceLSP,VarianceloessM),6,length(nazwy_kontroli),byrow=T)
rownames(bledy)<- c('BiasLOESS','BiasPRINTTIPLOESS','BiasLOESSmarray', 'VarianceLOESS','VariancePRINTTIPLOESS','VarianceLOESSmarray')
colnames(bledy)<- nazwy_kontroli
write.xls(bledy,file='bledy_normalizacji_wedlug_kontroli_bez_vsn.xls')

##################ZGODNOSC NORMALIZACJI#############################################
#zgodnosc normalizacji - Rysunek 4.9    - wybrac nalezy ktore porownania nas interesuja

loessMvsVSN=TRUE
loessMvsPrinttiploess=TRUE
loessMvsLoess=TRUE
PrinttiploessvsLoess=TRUE
vsnvsPrinttiploess=TRUE
vsnvsLoess=TRUE

#loessM vs VSN

#a - numer macierzy dla której chcemy sprawdziæ zgodnoœæ

for(a in 1: (dim(RG)[2])){
    if(loessMvsVSN){
        m<-mean(MAloessM_M[,a]-Mvsn[,a],na.rm=T)
        s<-sd(MAloessM_M[,a]-Mvsn[,a],na.rm=T)
        nazwa<-paste("porownania_metod_loessMvsVSN_macierz_",a,".jpeg",sep="")
        jpeg(nazwa)
        plot((MAloessM_A[,a]+Avsn[,a])/2,MAloessM_M[,a]-Mvsn[,a], xlab='œrednia',ylab='ró¿nica',main=c('loess(marray) - vsn', '\n', 'macierz',a))
        abline(h=m,col='green')
        abline(h=m-2*s,lty=2,col='red')
        abline(h=2*s+m,lty=2,col='red')
        dev.off()  
                            }
                            
        #loessM vs Printtiploess
     if(loessMvsPrinttiploess){
        m<-mean(MAloessM_M[,a]-MAlsp$M[,a],na.rm=T)
        s<-sd(MAloessM_M[,a]-MAlsp$M[,a],na.rm=T)
        nazwa<-paste("porownania_metod_loessMvsPrinttiploess_macierz_",a,".jpeg",sep="")
        jpeg(nazwa)
        plot((MAloessM_M[,a]+MAlsp$M[,a])/2,MAloessM_M[,a]-MAlsp$M[,a], main=c('loess(marray) - printtiploess','\n', 'macierz',a),xlab='œrednia',ylab='ró¿nica',ylim=c(-2,2))
        abline(h=m,col='green')
        abline(h=m-2*s,lty=2,col='red')
        abline(h=2*s+m,lty=2,col='red')
        dev.off()
                              }
                              
        #loessM vs loess
        
    if(loessMvsLoess){
        m<-mean(MAloessM_M[,a]-MAls$M[,a],na.rm=T)
        s<-sd(MAloessM_M[,a]-MAls$M[,a],na.rm=T)
        nazwa<-paste("porownania_metod_loessMvsLoess_macierz_",a,".jpeg",sep="")
        jpeg(nazwa)
        plot((MAloessM_M[,a]+MAls$M[,a])/2,MAloessM_M[,a]-MAlsp$M[,a], main=c('loess(marray) - loess (limma)', '\n', 'macierz',a),xlab='œrednia',ylab='ró¿nica',ylim=c(-2,2))
        abline(h=m,col='green')
        abline(h=m-2*s,lty=2,col='red')
        abline(h=2*s+m,lty=2,col='red')
        dev.off()
                      }
                      
        #Printtiploess vs Loess
        
    if(PrinttiploessvsLoess){
        m<-mean(MAlsp$M[,a]-MAls$M[,a],na.rm=T)
        s<-sd(MAlsp$M[,a]-MAls$M[,a],na.rm=T)
        nazwa<-paste("porownania_metod_PrinttiploessvsLoess_macierz_",a,".jpeg",sep="")
        jpeg(nazwa)
        plot((MAlsp$M[,a]+MAls$M[,a])/2,MAlsp$M[,a]-MAls$M[,a], main=c('loess(limma) - printtiploess', '\n', 'macierz',a),xlab='œrednia',ylab='ró¿nica',ylim=c(-2,2))
        abline(h=m,col='green')
        abline(h=m-2*s,lty=2,col='red')
        abline(h=2*s+m,lty=2,col='red')
        dev.off()
                            }
                            
        #vsn vs Printtiploess
        
    if(vsnvsPrinttiploess){
        m<-mean(Mvsn[,a]-MAlsp$M[,a],na.rm=T)
        s<-sd(Mvsn[,a]-MAlsp$M[,a],na.rm=T)
        nazwa<-paste("porownania_metod_vsnvsPrinttiploess_macierz_",a,".jpeg",sep="")
        jpeg(nazwa)
        plot((Mvsn[,a]+MAlsp$M[,a])/2,Mvsn[,a]-MAlsp$M[,a], main=c('vsn - printtiploess', '\n', 'macierz',a),xlab='œrednia',ylab='ró¿nica')
        abline(h=m,col='green')
        abline(h=m-2*s,lty=2,col='red')
        abline(h=2*s+m,lty=2,col='red')
        dev.off()
                          }
                          
        #vsn vs Loess(limma)
        
    if(vsnvsLoess){
        m<-mean(Mvsn[,a]-MAls$M[,a],na.rm=T)
        s<-sd(Mvsn[,a]-MAls$M[,a],na.rm=T)
        nazwa<-paste("porownania_metod_vsnvsLoess_macierz_",a,".jpeg",sep="")
        jpeg(nazwa)
        plot((Mvsn[,a]+MAls$M[,a])/2,Mvsn[,a]-MAls$M[,a], main=c('vsn - loess(limma)', '\n', 'macierz',a),xlab='œrednia',ylab='ró¿nica')
        abline(h=m,col='green')
        abline(h=m-2*s,lty=2,col='red')
        abline(h=2*s+m,lty=2,col='red')
        dev.off()
                  }
                              }
                              
##################UZUPELNIANIE BRAKUJACYCH DANYCH############################################

# uzupelnianie brakujacych danych po normalizacji - nalezy wybrac TRUE dla wybranej  do dalszych analiz metody normalizacji

metoda='Mean' # nalezy wybrac interesujaca nas metode uzupelniania danych, do wyboru: Mean, Median
loess=FALSE
loessM=FALSE
printtiploess=TRUE
vsn=FALSE
  if(loess){
      dane<-MAls$M
      dane<- uzupelnianie_danych(metoda,dane)
      
      #usuwanie wierszy empty i buffer
      
      empty_buffer<-c(which(RG$genes$Name=='empty'),which(RG$genes$Name=='buffer'))
      dane<-dane[-empty_buffer,]
            }
  if(loessM){
      dane<-MAlsm
      dane<- uzupelnianie_danych(metoda,dane)
      
      #usuwanie wierszy empty i buffer
      
      empty_buffer<-c(which(RG$genes$Name=='empty'),which(RG$genes$Name=='buffer'))
      dane<-dane[-empty_buffer,]
            }
  if(printtiploess){
      dane<-MAlsp$M
      dane<- uzupelnianie_danych(metoda,dane)
      
      #usuwanie wierszy empty i buffer
      
      empty_buffer<-c(which(RG$genes$Name=='empty'),which(RG$genes$Name=='buffer'))
      dane<-dane[-empty_buffer,]
                  }
  if(vsn){
      dane<-Mvsn
      dane<- uzupelnianie_danych(metoda,dane)
      
      #usuwanie wierszy empty i buffer
      
      empty_buffer<-c(which(RG$genes$Name=='empty'),which(RG$genes$Name=='buffer'))
      dane<-dane[-empty_buffer,]
         }

# dane - zmienna, ktora mozna wykorzystac do dalszych analiz tj. selekcja

#########################################SELEKCJA#######################################################

# ff - parametr mowiacy o tym, ile chcemy miec genow roznicujacych np. 50,100,200??
# L - klasy factor  - wektor opisujacy przynaleznosci do klas

# wczytanie potrzebnych zmniennych z przykladu  dla 50 genow roznicujacych

L<-as.factor(c(rep(0,4),rep(1,4)))
ff<-50

#########################
load('wczytane_funkcje.RData') # musza zostac wczesniej wgrane funkcje: testy, alg_kor,istotnosc_korelacji,p_wartosci, alg

dane_stale<-dane
dane_do_korelacji<-t(dane_stale)
cor_matrix<-cor(dane_do_korelacji,)
cor_genes<-abs(cor_matrix)
selekcja<-function(dane,L,ff){
    for(f in 1:(length(ff))){
	    dane<-dane_stale
          klasy<-unique(L)
          dane<-as.matrix(dane)
          X<- t(dane)
          bledy_100_200<-list()
          
          # 1 metoda t-test  albo F
          
          stat<-testy(dane,L)
          indeksy_t<-order(stat,decreasing=T)[1:ff[f]]
          
          # 2 metoda Kruskall-Wallis dla wiecej niz 2 klas
          
          if(length(unique(L))>2){
              p_k<-numeric()
              proba<-list()
              for(i in 1:(dim(dane)[2])){
                  for(j in 1:(length(klasy))){
                        proba[[j]]<-dane[i,which(L==klasy[j])]
                                              }
              a <-kruskal.test(proba)
              p_k[i]<-a$statistic
                                         }
          indeksy_Kruskall<-order(p_k,decreasing=TRUE)[1:ff[f]]
                                    }
                                    
          # 3 metoda - limma eBayes moderated t statistic
          
          design<-matrix(,dim(dane)[2],length(klasy))
          for(j in 1: (length(klasy))){
            klasa<- rep(0,dim(dane)[2])
            klasa[which(L==klasy[j])]<-1
            design[,j]<-klasa
                                      }
          fit <- eBayes(lmFit(dane,design))
          geny_naj<-topTable(fit,adjust="BY",number=ff[f])
          indeksy_limma<-as.numeric(rownames(geny_naj))
          
          # 4 metoda - sam statystyka  Tusher i in.(2001)
          
          x<-dane
          y<-as.numeric(L)
          d=list(x=x,y=y,logged2=TRUE)                             #x geny w wierszach, proby w kolumnach, y - wektor etykiet L
          samr.obj <- samr(d, resp.type="Two class unpaired",nperms=100)
          pv=samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
          indeksy_sam<-order(pv)[1:ff[f]]
          
          # 5 metoda - funkcja liczaca t-score dla dwoch grup wg. Tibshirani i Wasserman (2006)
          
          if(length(unique(L))==2){
          score = cst.stat(X, L)
          idx = order(abs(score), decreasing=TRUE)
          indeksy_TW<-idx[1:ff[f]]
          
          # 6 metoda Zuber i Strimmer (2009)  dla dwoch klas
          
          score = shrinkcat.stat(X, L)
          idx = order(abs(score), decreasing=TRUE)
          indeksy_st<-idx[1:ff[f]]
                                  }
                                  
          # 7 metoda (Ahdesmaki i Strimmer (2009))  - dla wielu klas

          ranking.DDA = sda.ranking(X, L, diagonal=TRUE)
          ranking<-ranking.DDA[1:ff[f],]
          indeksy_DDA<-ranking[,1]
          
          # metoda  1    - tFK1
          
          calosc<-list()
          stat_alg<-list()
          w<-seq(0,1,by=0.1)
          for(k in 1:11){
              stat_alg[[k]]<-alg(dim(dane)[1]-1,dim(dane)[1],w[k],stat)
              calosc[[k]]<-order(stat_alg[[k]],decreasing=TRUE)[1:ff[f]]
                        }
                        
          # moja metoda 2   - tFK2
          
          cor_matrix<-abs(cor_matrix)
          macierz<-istotnosc_korelacji(cor_matrix,0.8,dim(dane)[1])         #cor_matrix ma wartosci bezwzgledne korelacji
          macierz<-abs(macierz)                 # macierz z istotnoscia korelacji
          p_value<-p_wartosci(dane,stat,macierz,50)          #50 permutacji
          indexALG1<-list()
          indexALG<-list()
          p_mean1<-list()
          for(d in 1:11){
              p_mean1[[d]]<-p.adjust(p_value$ p_wartosci[[d]],method=c("BY"))
              indexALG1[[d]]<-order(p_mean1[[d]])[1:ff[f]]
              indexALG[[d]]<- p_value$indeksy[[d]][1:ff[f]]
               }
               
#lista indeksow 1

if(length(unique(L))<=2){
  lista<-list(F=indeksy_t,L=indeksy_limma,SAM=indeksy_sam,TW=indeksy_TW,ST=indeksy_st,DDA=indeksy_DDA,tFK11=calosc[[1]],tFK12=calosc[[2]],tFK13=calosc[[3]],tFK14=calosc[[4]],tFK15=calosc[[5]],tFK16=calosc[[6]],tFK17=calosc[[7]],tFK18=calosc[[8]],tFK19=calosc[[9]],tFK110=calosc[[10]],tFK111=calosc[[11]]
  ,tFK21=indexALG[[1]],tFK22=indexALG[[2]],tFK23=indexALG[[3]],tFK24=indexALG[[4]],tFK25=indexALG[[5]],tFK26=indexALG[[6]],tFK27=indexALG[[7]],tFK28=indexALG[[8]],tFK29=indexALG[[9]],tFK210=indexALG[[10]],tFK211=indexALG[[11]])
                        } else {
                        
#lista indeksow  2

  lista<-list(F=indeksy_t,K=indeksy_Kruskall,L=indeksy_limma,SAM=indeksy_sam,DDA=indeksy_DDA,tFK11=calosc[[1]],tFK12=calosc[[2]],tFK13=calosc[[3]],tFK14=calosc[[4]]
  ,tFK15=calosc[[5]],tFK16=calosc[[6]],tFK17=calosc[[7]],tFK18=calosc[[8]],tFK19=calosc[[9]],tFK110=calosc[[10]],tFK111=calosc[[11]],tFK21=indexALG1[[1]],tFK22=indexALG1[[2]],tFK23=indexALG[[3]],tFK24=indexALG[[4]],tFK25=indexALG[[5]],tFK26=indexALG[[6]],tFK27=indexALG[[7]],tFK28=indexALG[[8]],tFK29=indexALG[[9]],tFK210=indexALG[[10]],tFK211=indexALG[[11]])
                                }
                                
 # uczenie maszyn
 
 mis_classified1<-numeric()
 mis_classified2<-numeric()
 mis_classified3<-numeric()
 mis_classified4<-numeric()
 mis_classified5<-numeric()
 for(k in 1:(length(lista))){
    dane<-dane_stale
    dane_k<-dane[lista[[k]],]
    dane_k<-t(dane_k)
    L<-as.factor(L)
    
    # Naive Bayes
    
    dane<-data.frame(dane_k,L)
    True<-as.numeric(L)-1
    poprawne<-0
    for (i in 1:(dim(dane)[1])){
        training_set = c(1:(dim(dane)[1]))[-i]
        klasyfikator = MLearn(L~ ., data = dane, naiveBayesI, training_set )
        if(testPredictions(klasyfikator)==True[i])
            poprawne<-poprawne+1;
                                }
    mis_classified1[k]<-(dim(dane)[1])-poprawne
    
    # knn
    
    poprawne<-0
    for (i in 1:(dim(dane)[1])){
        training_set= c(1:(dim(dane)[1]))[-i]
        klasyfikator = MLearn(L~ ., data=dane, knnI(k=3,l=0), training_set)
        if(testPredictions(klasyfikator)==True[i])
            poprawne<-poprawne+1;
                                }
    mis_classified2[k]<-(dim(dane)[1])-poprawne
    
    # SVM
    
    poprawne<-0
    for (i in 1:(dim(dane)[1])){
        training_set= c(1:(dim(dane)[1]))[-i]
        klasyfikator = MLearn(L~ ., data=dane, svmI, training_set)
        if(testPredictions(klasyfikator)==True[i])
            poprawne<-poprawne+1;
                                 }
    mis_classified3[k]<-(dim(dane)[1])-poprawne
    
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
    mis_classified4[k]<-mean(mis_classified4a)
    
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
    mis_classified5[k]<-mean(mis_classified5a)
                                      }
bledy<-matrix(c(mis_classified1,mis_classified2,mis_classified3,mis_classified4,mis_classified5),5,length(lista),byrow=T)
if(length(unique(L))<=2){
    colnames(bledy)<-c('t','limma','sam','TW','st','DDA','0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6' ,'0.7' ,'0.8', '0.9', '1.0','0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6' ,'0.7' ,'0.8', '0.9', '1.0')
                        } else {
    colnames(bledy)<-c('F','Kruskall','limma','sam','DDA','0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6' ,'0.7' ,'0.8', '0.9', '1.0','0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6' ,'0.7' ,'0.8', '0.9', '1.0')
                                }
rownames(bledy)<-c('NB','knn','svm','nnet','rf')
bledy_100_200[[f]]<-bledy
                        }
wyniki_selekcji<-return(list(bledy_100_200,lista))
attach(wyniki_selekcji)
                          }


wyniki_selekcji<-selekcja(dane,L,ff)                          
#wywolanie funkcji selekcja
bledy<-wyniki_selekcji[[1]]
lista<-wyniki_selekcji[[2]]

save(bledy,lista,file='wyniki_selekcji.RData')
write.xls(bledy,file='wyniki_selekcji.xls')


#################################################REGRESJA STOSOWA###################################
# przykladowy wybor najlepszych metod

najlepsze<-c(1:6,14,19)
index<-lista[najlepsze]			# indeksy genow roznicujacych dla roznych najlepszych metod

library(Hmisc)

pole_auc<-numeric(length(index))
krzywe_roc<-list()
names(pole_auc)<-names(index)
dane<-dane_stale
ilosc.prob <- dim(dane)[2]
ilosc.klas<-length(unique(L))
aposteriori<-list()
for(r in 1: length(index)){
    cat(paste(r,'... \n'))
    dane<-dane_stale
    dane_k<-dane[index[[r]],]                 #indeksy genow roznicujacych to index
    dane_k<-t(dane_k)           
    L<-as.factor(L)
    dane<-data.frame(dane_k,L)             # wiersze to proby, a kolumny to geny + ostatnia kolumna to labels
    	
	#Svm
      cat(paste('svm...'))
    x.svm.prob_LOO<-list()
    x.svm.prob_LOO_nowe<-matrix(,ilosc.prob,ilosc.klas)
    for(i in 1:ilosc.prob)
    {
        x.svm_LOO <- svm(L~., data = dane[-i,], cost=1, gamma=0.015625, probability = TRUE)
        x.svm.prob_LOO[[i]] <- predict(x.svm_LOO, type="prob", newdata=dane[i,], probability = TRUE)
        for (k in 1:ilosc.klas)
        {
            x.svm.prob_LOO_nowe[i,k]<-attr(x.svm.prob_LOO[[i]],"probabilities")[,which(colnames(attr(x.svm.prob_LOO[[i]],"probabilities"))==k-1)]
        }
    }
        
     #knn
    cat(paste('knn...'))
    u<-list()
    knn_prob.pomoc<-matrix(,ilosc.prob,ilosc.klas)
    knn_pred <-matrix(,ilosc.prob,ilosc.klas)
    for (k in 1:ilosc.klas)
    {
        u[[k]] <- as.numeric(as.numeric(L)==k)
        knn_pred[,k] <- c(knn.cv(dane,u[[k]],k=3, prob=TRUE))-1
        knn_prob.pomoc[,k] <- attr((knn.cv(dane,u[[k]],k=3, prob=TRUE)),"prob")
    }

    p_knn<-matrix(0,ilosc.prob,ilosc.klas)
    for (i in 1:ilosc.prob)
    {
        for (k in 1:ilosc.klas)
        {
            if (knn_pred[i,k]==1){
                p_knn[i,k]<-knn_prob.pomoc[i,k]
            } else
            {
                p_knn[i,k]<-1-knn_prob.pomoc[i,k]
            }
        }
    }
          
      #nnet
      
      cat(paste('nnet...'))
    if (ilosc.klas==2)
    {
        p_nnet<-nnet(L~., data=dane, size=3, decay=.01,trace=F)$fitted.values
        p_nnet<-cbind(p_nnet,1-p_nnet)
        rownames(p_nnet)<-NULL
    } else
    {
        p_nnet<-nnet(L~., data=dane, size=3, decay=.01,trace=F)$fitted.values
        rownames(p_nnet)<-NULL
    }

    p_svm<-x.svm.prob_LOO_nowe

      
#tworzenie macierzy P

    P<-matrix(0,ilosc.prob,3*ilosc.klas)
    for (i in 1:ilosc.prob)
    {
        P[i,]<-c(p_svm[i,],p_knn[i,],p_nnet[i,])
    }

 Moore_Penrose<-ginv(t(P)%*%P)
u1<- abs(1-(as.numeric(L)-1))                         #dla grupy oznaczonej jako 0
u2<-as.numeric(L)-1                                   #dla grupy oznaczonej jako 1
Estymator_Beta<-list()
u<-list(u1,u2)
for(k in 1:2){   #k-ilosc klas
Estymator_Beta[[k]]<-Moore_Penrose%*%t(P)%*%u[[k]]
              }
# Estymator u_k
estymator_u<-list()
estymator_u[[1]]<-0
estymator_u[[2]]<-0
for(k in 1:2){
  for(i in 1:dim(dane)[1]){
      estymator_u[[k]]<-c(estymator_u[[k]],sum(P[i,]*Estymator_Beta[[k]]))               #iloczyn skalarny wektorow
                }
              }
estymator_u[[1]]<-estymator_u[[1]][-1]    #dla grupy oznaczonej jako 0
estymator_u[[2]]<-estymator_u[[2]][-1]    #dla grupy oznaczonej jako 1
#which(estymator_u1>estymator_u2)
klasa0<-which(estymator_u[[1]]>estymator_u[[2]])
klasa1<-which(estymator_u[[1]]<estymator_u[[2]])
# liczenie prawdopodobienstw a posteriori dla klasyfikatora regresji stosowej
#estymator_u1[which(estymator_u1<0)]<-0
#estymator_u2[which(estymator_u2<0)]<-0
estymator_u[[1]][which(estymator_u[[1]]<0)]<-0
estymator_u[[2]][which(estymator_u[[2]]<0)]<-0
estymator_u[[1]][which(estymator_u[[1]]>1)]<-1
estymator_u[[2]][which(estymator_u[[2]]>1)]<-1
p_aposteriori<-list()
for(k in 1:2){
  p_aposteriori[[k]]<-numeric(40)
  for(i in 1:dim(dane)[1]){
    p_aposteriori[[k]][i]<-estymator_u[[k]][i]/(estymator_u[[1]][i]+estymator_u[[2]][i])
                }
                }

##################KRZYWA ROC DLA REGRESJI STOSOWEJ############################################                

if(r==1){
a<-roc(dane[,'L'],p_aposteriori[[2]] )
krzywe_roc[[r]]<-a
plsmo(1-a$specificities,a$sensitivities,method="supsmu",type='l',col=1,xlab='1 - specyficznoœæ', ylab='czu³osc', main='Porównianie krzywych ROC dla ró¿nych metod selekcji')
pole_auc[r]<- roc(dane[,'L'],p_aposteriori[[2]])$ auc
legend(0.75, 0.9, names(index), 1:length(index))
        } else {
a<-roc(dane[,'L'],p_aposteriori[[2]])
krzywe_roc[[r]]<-a
plsmo(1-a$specificities,a$sensitivities,method="supsmu",type='l',col=r,add=T)
pole_auc[r]<-roc(dane[,'L'],p_aposteriori[[2]])$ auc
}
}

save(krzywe_roc,pole_auc,file='wczytanie_krzywe_roc_pole_auc.RData')

#######test dla hipotezy H0: AUC1=AUC2


np<- 8       # np - liczba porownan - nalezy wstawic liczbe elementów listy o nazwie index
p.values.auc.test<-matrix(,np,np)

for (l in 1:(np-1))
    {
    for (m in (l+1):np)
    {
        test.proba<-roc.test(krzywe_roc[[l]],krzywe_roc[[m]],method="bootstrap")
        p.values.auc.test[l,m]<-test.proba$p.value
        }
}

#koniec programu


