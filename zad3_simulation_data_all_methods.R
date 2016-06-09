
library('animation')
library('DESeq')
library('edgeR')
library('samr')
setwd('/Users/joannazyprych/Desktop/praca porownania')

#wczytanie danych
fileNames<-c('B1250.csv','B625.csv','B4000.csv','B2000.csv','P0.csv','P625.csv','S0.csv','S625.csv')
Names<-c('B1250','B625','B4000','B2000','P0','P625','S0','S625')
# readFiles<- function(fileNames) { 
#   for(i in seq_along(fileNames)) { 
#     data1<-read.csv2(fileNames[i] ,sep=",")[,-1]
#   } 
# }

for(i in seq_along(fileNames)) { 
  data<-read.csv2(fileNames[i] ,sep=",")[,-1]
  
#   #DESeq
#   cds <- newCountDataSet(data, c(rep(1,10),rep(2,10)))
#   cds <- estimateSizeFactors( cds )
#   cds <- estimateDispersions( cds) #Błąd w parametricDispersionFit(means, disps) : 
#   #Parametric dispersion fit failed. Try a local fit and/or a pooled estimation.
#   res.deseq <- nbinomTest( cds, '1','2')
#   o<-order(res.deseq$padj)
#   write.csv(o,'kolejnosc_genow_dla_calego_eksperymentu.csv')
#   o<-read.csv('kolejnosc_genow_dla_calego_eksperymentu.csv')[,2]
#   granica<-length(which(res.deseq$padj<.1))
# 

plotDE <- function( res,con )
{
  if (dim(res)[1]<granica) {
    plot(res$baseMean,res$log2FoldChange,log="x",main=con,xlab='baseMean',ylab='log2FoldChange', pch=20, cex=.3,col = 'blue')
    lines(res$baseMean,res$log2FoldChange,main=con,xlab='baseMean',ylab='log2FoldChange', type='p',pch=20, cex=.3,col = ifelse( res$padj < .1, "green", "red" ) )
    legend("bottomright", legend=c('fp','tp','fn','tn'), fill =c('purple','green','red','black') ,col=c('purple','green','red','black'),bty='n')
  }else{
    fng<-which(res$padj[1:(min(a*b,granica))] > .1)
    fpg<-which(res$padj[((min(a*b,granica))+1):dim(res)[1]] < .1)
    if(length(fng)==0){
      plot(res$baseMean[c(1:(min(a*b,granica)))],res$log2FoldChange[c(1:(min(a*b,granica)))],main=con,xlab='baseMean',
           ylab='log2FoldChange', type='p',pch=20, cex=.3,col = ifelse( res$padj < .1, "green", "black" ) )
    } else {
      plot(res$baseMean[c(1:(min(a*b,granica)))[-fng]],res$log2FoldChange[c(1:(min(a*b,granica)))[-fng]],main=con,xlab='baseMean',
           ylab='log2FoldChange', type='p',pch=20, cex=.3,col = ifelse( res$padj < .1, "green", "black" ) ) 
    }
    lines(res$baseMean[fng],res$log2FoldChange[fng],log="x",main=con,xlab='baseMean',ylab='log2FoldChange', 
         pch=20, cex=.3,col = 'red')
    legend("bottomright", legend=c('fp','tp','fn','tn'), fill =c('purple','green','red','black') ,col=c('purple','green','red','black'),bty='n')
    lines(res$baseMean[c(((min(a*b,granica))+1):dim(res)[1])[fpg]],res$log2FoldChange[c(((min(a*b,granica))+1):dim(res)[1])[fpg]],main=con,xlab='baseMean',
          ylab='log2FoldChange', type='p',pch=20, cex=.3,col = 'purple' )
    
    lines(res$baseMean[c(((min(a*b,granica))+1):dim(res)[1])[-fpg]],res$log2FoldChange[c(((min(a*b,granica))+1):dim(res)[1])[-fpg]],main=con,xlab='baseMean',
          ylab='log2FoldChange', type='p',pch=20, cex=.3,col = ifelse( res$padj < .1, "green", "black" ) )}
  #text(10000,par("usr")[3]+2,paste(length(which(res$log2FoldChange<0)),'\ndownregulated\ngenes'),cex=1)
  text(1000,par("usr")[4]-.5,paste(length(which(res$padj<0.1)),' found of ',(min(a*b,granica)),' DEG'),cex=1.5)
}
  
# koniec<-length(o)
# za ilosc podstawiamy ilosci genow, dla ktorych chcemy spr jak dziala deseq, dla przykladu zrobimy dla 50, 150,...,20050 genow 
ilosc<-seq(50,900,100)
# disp<-matrix(,2,length(ilosc))
# wspolne<-numeric(length(ilosc))
# w=1
# saveHTML(
# {
#   for (con in ilosc)
#   {
#     b<-con
#     for(a in seq(0.1,1,0.1))
#     {
#       cc.small<-cc[c(o[1:(a*con)],o[(koniec-(1-a)*con):koniec]),]
#       if(a==1){
#         cc.small<-data[o[1:(a*con)],] 
#       }
#       cds <- newCountDataSet(cc.small,  c(rep(1,10),rep(2,10)))
#       cds <- estimateSizeFactors( cds )
#       cds <- estimateDispersions( cds)
#       res.deseq <- nbinomTest( cds, '1','2')
#       plotDE(res.deseq, paste('MAplot for ',con,' genes',sep=''))
#       w<-w+1
#       cat(w,"\n")
#     }
#   }
# } , img.name = paste("deseq for data ",fileNames[i]), autobrowse = FALSE,  imgdir = paste("deseq_dir ",fileNames[i]), 
# htmlfile = paste("deseq for data ",fileNames[i],".html"), outdir = getwd(), title = "Deseqi_all", 
# description = c("Na razie deseqi tak szaleja..."))


######################################################################
y <- DGEList(data, group=as.character(c(rep(1,10),rep(2,10))))
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
et <- et$table
padj<-p.adjust(et[,3],method='BH')
et<-data.frame(et,padj)
o<-order(et$padj)

granica<-length(which(et$padj<.1))
koniec<-length(o)
baseMean<-rowMeans(data/y$samples[,3])
w=1

saveHTML(
{
  for (con in ilosc)
  {
    for(a in seq(0.1,1,0.1))
    {
      b<-koniec-((1-a)*con-1)
      cc.small<-data[c(o[1:(a*con)],o[b:koniec]),]
      if(a==1){
        cc.small<-data[o[1:(a*con)],] 
      }
      y <- DGEList(cc.small, group=as.character(c(rep(1,10),rep(2,10))))
      y <- calcNormFactors(y)
      y <- estimateCommonDisp(y)
      y <- estimateTagwiseDisp(y)
      et <- exactTest(y)
      et <- et$table
      padj<-p.adjust(et[,3],method='BH')
      et<-data.frame(et,padj,baseMean[o[1:con]])
      colnames(et)[c(2,5)]<-c('log2FoldChange','baseMean')
      plotDE(et, paste('MAplot for ',con,' genes',sep=''))
      w<-w+1
      cat(w,"\n")
    }
  }
} , img.name = paste("edgeR for data ",Names[i]), autobrowse = FALSE,  imgdir = paste("edgeR_dir ",Names[i]), 
htmlfile = paste('edgeR for data',Names[i],'.html'), outdir = getwd(), title = paste("edgeR ",Names[i]), 
description = c("EdgeR tak szaleje..."))
}


##################################################################################
ilosc<-seq(50,20650,100)
y<-c(rep(1,10),rep(2,10))
cc1<-data[rowMeans(data)>5,]
sam.result<-SAMseq(cc1,y,resp.type="Two class unpaired",nperm=100,fdr.output = 1)
sam_comb <- as.data.frame(rbind(sam.result$siggenes.table$genes.up,
                                sam.result$siggenes.table$genes.lo))
sam_comb[,3]<-as.numeric(levels(sam_comb[,3]))[sam_comb[,3]]
sam_comb[,4]<-as.numeric(levels(sam_comb[,4]))[sam_comb[,4]]
sam_comb[,5]<-as.numeric(levels(sam_comb[,5]))[sam_comb[,5]]
baseMean<-rowMeans(cc1[sam_comb[,2],]/sam.result[[1]]$depth)
sam_comb<-data.frame(sam_comb,log(sam_comb[,4],2),baseMean)

colnames(sam_comb)[c(6)]<-c('log2FoldChange')
sam_comb<-sam_comb[order(sam_comb[,5],-abs(sam_comb[,6])),]

o<-sam_comb[,2]
write.csv(o,'kolejnosc_genow_dla_calego_eksperymentu_SAMseq.csv')
granica<-length(which(sam_comb[,5]<10))
w=1
saveHTML(
{
  for (con in ilosc)
  {
    b<-con
    for(a in seq(0.1,1,0.1))
    {
      cc.small<-cc[c(o[1:(a*con)],o[(koniec-(1-a)*con):koniec]),]
      if(a==1){
        cc.small<-data[o[1:(a*con)],] 
      }
      sam.result<-SAMseq(cc.small,y,resp.type="Two class unpaired",nperm=100,fdr.output = 1)
      sam_comb <- as.data.frame(rbind(sam.result$siggenes.table$genes.up,
                                      sam.result$siggenes.table$genes.lo))
      sam_comb[,3]<-as.numeric(levels(sam_comb[,3]))[sam_comb[,3]]
      sam_comb[,4]<-as.numeric(levels(sam_comb[,4]))[sam_comb[,4]]
      sam_comb[,5]<-as.numeric(levels(sam_comb[,5]))[sam_comb[,5]]/100
      baseMean<-rowMeans(cc1[sam_comb[,2],]/sam.result[[1]]$depth)
      sam_comb<-data.frame(sam_comb,log(sam_comb[,4],2),baseMean)
      
      colnames(sam_comb)[c(6)]<-c('log2FoldChange')
      sam_comb<-sam_comb[order(sam_comb[,5],-abs(sam_comb[,6])),]
      plotDE(sam_comb, paste('MAplot for ',con,' genes',sep=''))
      w<-w+1
      cat(w,"\n")
    }
  }
} , img.name = "SAMseq_proba", autobrowse = FALSE,  imgdir = "SAMseq_dir", 
htmlfile = "SAMseq_proba.html", outdir = getwd(), title = "SAMseq", 
description = c("SAMseq tak szaleje..."))

owd = setwd('/Users/joannazyprych/Desktop/praca porownania')
ani.options(convert = shQuote('/opt/ImageMagick/bin/convert'), interval=.1,outdir='/Users/joannazyprych/Desktop/praca porownania')
im.convert(paste('SAMseq_proba',1:length(ilosc),'.png',sep=''), output = "animation_SAMseq.gif")

%%%%%%%%%%%%%%%%
library(NBPSeq)  #nie dziala pakiet qvalue
NBPSeq.dgelist = DGEList(counts = count.matrix, group = factor(class))
NBPSeq.dgelist = calcNormFactors(NBPSeq.dgelist, method = "TMM")
NBPSeq.norm.factors = as.vector(NBPSeq.dgelist$samples$norm.factors)
NBPSeq.test = nbp.test(counts = count.matrix, grp.ids = class,
                         + grp1 = 1, grp2 = 2, norm.factors = NBPSeq.norm.factors,
                         + method.disp = "NBP")
NBPSeq.pvalues = NBPSeq.test$p.values
NBPSeq.adjpvalues = NBPSeq.test$q.values
  ##########
library(baySeq)
class<-c(rep('A1',10),rep('A2',10))
baySeq.cd = new("countData", data = as.matrix(data), replicates = class,groups = list(NDE = rep(1, length(class)), DE = class))
baySeq.cd@libsizes = getLibsizes(baySeq.cd, estimationType = "edgeR")
baySeq.cd = getPriors.NB(baySeq.cd, samplesize = 50,equalDispersions = TRUE, estimation = "QL", cl = NULL)
baySeq.cd = getLikelihoods.NB(baySeq.cd, prs = c(0.5,0.5), pET = "BIC", cl = NULL)
baySeq.posteriors.DE = exp(baySeq.cd@posteriors)[, 2]
baySeq.table = topCounts(baySeq.cd, group = "DE", FDR = 1)
baySeq.FDR = baySeq.table$FDR[match(rownames(data),rownames(baySeq.table))]
##############
source("http://bioconductor.org/biocLite.R")
devel = "http://bioconductor.org/packages/2.13/bioc"
biocLite("EBSeq", siteRepos = devel)
library(EBSeq)
sizes = MedianNorm(data)
EBSeq.test = EBTest(Data = data, Conditions = factor(class),sizeFactors = sizes, maxround = 10)
EBSeq.ppmat = GetPPMat(EBSeq.test)
EBSeq.probabilities.DE = EBSeq.ppmat[, "PPDE"]
EBSeq.lFDR = 1 - EBSeq.ppmat[, "PPDE"]
EBSeq.FDR = rep(NA, length(EBSeq.lFDR))
for (i in 1:length(EBSeq.lFDR)) {
  EBSeq.FDR[i] = mean(EBSeq.lFDR[which(EBSeq.lFDR <=EBSeq.lFDR[i])])
                                  }
############

source("TSPM.R")
TSPM.dgelist = DGEList(counts = count.matrix, group = factor(class))
TSPM.dgelist = calcNormFactors(TSPM.dgelist, method = "TMM")
norm.lib.sizes = as.vector(TSPM.dgelist$samples$norm.factors) *as.vector(TSPM.dgelist$samples$lib.size)x0 = rep(1, length(class)), lib.size = norm.lib.sizes)
TSPM.pvalues = TSPM.test$pvalues
TSPM.adjpvalues = TSPM.test$padj
################
source('http://bioconductor.org/biocLite.R')
biocLite('NOISeq')
library('NOISeq')
source("noiseq.r")
class<-c(rep(1,10),rep(2,10))
nf = calcNormFactors(data)
libsizes = apply(data, 2, sum)
common.libsize = prod(libsizes^(1/length(libsizes)))
normfactors = nf * libsizes/common.libsize
norm.matrix = sweep(data, 2, normfactors, "/")
NOISeq.test = noiseq(norm.matrix[, class == 1], norm.matrix[, class == 2], repl = "bio", k = 0.5, norm = "n")
NOISeq.probabilities = NOISeq.test$probab