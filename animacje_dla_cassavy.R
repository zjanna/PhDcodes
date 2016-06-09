library('animation')
library('DESeq2')
library('edgeR')
library('samr')

setwd("~/Desktop/praca porownania/dane rzeczywiste")

#wczytanie danych
source("plotDE.R")

data<-read.csv2("cc.csv")
rownames(data)<-data[,1]
data<-data[,-1]
data<-data[rowMeans(data)>1,]

### edgeR
y <- DGEList(data, group=as.character(c(rep(1,6),rep(2,6))))
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
et <- et$table
padj<-p.adjust(et[,3],method='BH')
baseMean<-rowMeans(t(t(data)/y$samples$norm.factors))
et<-data.frame(et,padj,baseMean)
colnames(et)[c(1,5)]<-c('log2FoldChange','baseMean')
o<-order(et$padj)
et<-et[o,]
save(et,file="wyniki_edgeR_cassava.RData")

### DESeq
colData<-data.frame(condition=as.factor(c(rep(1,6),rep(2,6))))
rownames(colData)<-colnames(data)
dds <- DESeqDataSetFromMatrix(countData = data,colData = colData,design = ~ condition) 
dds <- DESeq(dds)
res <- results(dds)
o<-order(res$padj)
res<-res[o,]
save(res,file="wyniki_DESeq_cassava.RData")

### SAMseq
y<-c(rep(1,6),rep(2,6))
sam.result<-SAMseq(data,y,resp.type="Two class unpaired",nperm=100,geneid =rownames(data),
                   fdr.output = 1,random.seed =10)
sam_comb <- as.data.frame(rbind(sam.result$siggenes.table$genes.up,
                                sam.result$siggenes.table$genes.lo))
sam_comb[,3]<-as.numeric(levels(sam_comb[,3]))[sam_comb[,3]]
sam_comb[,4]<-as.numeric(levels(sam_comb[,4]))[sam_comb[,4]]
sam_comb[,5]<-as.numeric(levels(sam_comb[,5]))[sam_comb[,5]]/100
baseMean<-rowMeans(t(t(data[sam_comb[,2],])/sam.result[[1]]$depth))
sam_comb<-data.frame(sam_comb,log2(sam_comb[,4]),baseMean)

colnames(sam_comb)[c(6)]<-c('log2FoldChange')
sam_comb<-sam_comb[order(sam_comb[,5],-abs(sam_comb[,6])),]
save(sam_comb,file="wyniki_SAMseq_cassava.RData")

######################################################################
###### animacje dla edgeR ##########
n<-seq(50,2000,50)
deg<-seq(.1,.9,.1)

kolejnosc<-match(rownames(et),rownames(data))
tpg<-length(which(et$padj < .1))
tpg<-kolejnosc[1:tpg]
w<-0
saveHTML(
{
  for (j in deg)
  {
    for(i in n) 
    { 
      n.deg<-min(i*j,length(tpg))
      if(n.deg!=length(tpg)){
      a<-c(head(kolejnosc,n.deg),sample(kolejnosc[-tpg],i-n.deg))
      dane<-data[a,]    
      y <- DGEList(dane, group=as.character(c(rep(1,6),rep(2,6))))
      y <- calcNormFactors(y)
      y <- estimateCommonDisp(y)
      y <- estimateTagwiseDisp(y)
      et <- exactTest(y)
      et <- et$table
      padj<-p.adjust(et[,3],method='BH')
      baseMean<-rowMeans(t(t(dane)/y$samples$norm.factors))
      et<-data.frame(et,padj,baseMean)
      colnames(et)[c(1,5)]<-c('log2FoldChange','baseMean')
      dobre<-tpg[1:n.deg]
      plotDE(et, dobre)
      w<-w+1
      cat(w,"\n")
    }
    }
  }
} , img.name = "EdgeR_cassava_small_data", autobrowse = FALSE,  imgdir = "edgeR_cassava_dir", 
htmlfile = "EdgeR_cassava_small_data.html", outdir = getwd(), title = "EdgeR", 
description = c("EdgeR tak szaleje..."))

###### animacje dla edgeR ##########
n<-seq(50,2000,50)
deg<-seq(.1,.9,.1)

kolejnosc<-match(rownames(res),rownames(data))
tpg<-length(which(res$padj < .1))
tpg<-kolejnosc[1:tpg]

w<-0
saveHTML(
{
  for (j in deg)
  {
    for(i in n) 
    { 
      n.deg<-i*j
      a<-c(head(kolejnosc,n.deg),tail(kolejnosc,i-n.deg))
      dane<-data[a,]    
      colData<-data.frame(condition=as.factor(c(rep(1,6),rep(2,6))))
      rownames(colData)<-colnames(dane)
      dds <- DESeqDataSetFromMatrix(countData = dane,colData = colData,design = ~ condition) 
      dds <- DESeq(dds)
      res <- results(dds)
      dobre<-ifelse(length(tpg)==min(length(tpg),n.deg),tpg,tpg[1:n.deg])
      plotDE(et, dobre,data)
      w<-w+1
      cat(w,"\n")
    }
  }
} , img.name = "DESeq_cassava_small_data", autobrowse = FALSE,  imgdir = "DESeq_cassava_dir", 
htmlfile = "DESeq_cassava_small_data.html", outdir = getwd(), title = "DESeq", 
description = c("DESeq tak szaleje..."))

##################################################
n<-seq(50,2000,50)
deg<-seq(.1,.9,.1)

kolejnosc<-match(rownames(sam_comb),rownames(data))
tpg<-length(which(res$padj < .1))
tpg<-kolejnosc[1:tpg]


w=0
saveHTML(
{
  for (j in deg)
  {
    for(i in n) 
    { 
      n.deg<-i*j
      a<-c(head(kolejnosc,n.deg),tail(kolejnosc,i-n.deg))
      dane<-data[a,]    
      sam.result<-SAMseq(dane,c(rep(1,6),rep(2,6)),resp.type="Two class unpaired",nperm=100,
                         genenames =rownames(dane),fdr.output = .1,random.seed =10)
      sam_comb <- as.data.frame(rbind(sam.result$siggenes.table$genes.up,
                                      sam.result$siggenes.table$genes.lo))
      sam_comb[,3]<-as.numeric(levels(sam_comb[,3]))[sam_comb[,3]]
      sam_comb[,4]<-as.numeric(levels(sam_comb[,4]))[sam_comb[,4]]
      sam_comb[,5]<-as.numeric(levels(sam_comb[,5]))[sam_comb[,5]]/100
      baseMean<-rowMeans(t(t(dane[sam_comb[,2],])/sam.result[[1]]$depth))
      sam_comb<-data.frame(sam_comb[,-c(2,3,4)],log2(sam_comb[,4]),baseMean)
      rownames(sam_comb)<-NULL
      colnames(sam_comb)[c(2,3)]<-c('fdr','log2FoldChange')
      
      b<-setdiff(rownames(dane),as.character(sam_comb[,1]))
      log2FoldChange<-rowMeans(t(t(dane[b,7:12])/sam.result[[1]]$depth[7:12]))/rowMeans(t(t(dane[b,1:6])/sam.result[[1]]$depth[1:6]))
      baseMean<-rowMeans(t(t(dane[b,])/sam.result[[1]]$depth))
      sam_comb<-rbind(sam_comb,data.frame(Gene.ID=b,fdr=1,log2FoldChange=log2FoldChange,baseMean=baseMean))

      sam_comb<-sam_comb[order(sam_comb[,2],-abs(sam_comb[,3])),]
      rownames(sam_comb)<-sam_comb[,1]
      sam_comb<-sam_comb[rownames(dane),]
      dobre<-ifelse(length(tpg)==min(length(tpg),n.deg),tpg,tpg[1:n.deg])
      plotDE(sam_comb, dobre,data)
      
      w<-w+1
      cat(w,"\n")
    }
  }
} , img.name = "SAMseq_cassava_small_data", autobrowse = FALSE,  imgdir = "SAMseq_cassava_dir", 
htmlfile = "SAMseq_cassava_small_data.html", outdir = getwd(), title = "SAMseq", 
description = c("DESeq tak szaleje..."))

#owd = setwd('/Users/joannazyprych/Desktop/praca porownania')
#ani.options(convert = shQuote('/opt/ImageMagick/bin/convert'), interval=.1,outdir='/Users/joannazyprych/Desktop/praca porownania')
#im.convert(paste('SAMseq_proba',1:length(ilosc),'.png',sep=''), output = "animation_SAMseq.gif")

# library(NBPSeq)  #nie dziala pakiet qvalue
# NBPSeq.dgelist = DGEList(counts = count.matrix, group = factor(class))
# NBPSeq.dgelist = calcNormFactors(NBPSeq.dgelist, method = "TMM")
# NBPSeq.norm.factors = as.vector(NBPSeq.dgelist$samples$norm.factors)
# NBPSeq.test = nbp.test(counts = count.matrix, grp.ids = class,
#                          + grp1 = 1, grp2 = 2, norm.factors = NBPSeq.norm.factors,
#                          + method.disp = "NBP")
# NBPSeq.pvalues = NBPSeq.test$p.values
# NBPSeq.adjpvalues = NBPSeq.test$q.values
#   ##########
# library(baySeq)
# class<-c(rep('A1',10),rep('A2',10))
# baySeq.cd = new("countData", data = as.matrix(data), replicates = class,groups = list(NDE = rep(1, length(class)), DE = class))
# baySeq.cd@libsizes = getLibsizes(baySeq.cd, estimationType = "edgeR")
# baySeq.cd = getPriors.NB(baySeq.cd, samplesize = 50,equalDispersions = TRUE, estimation = "QL", cl = NULL)
# baySeq.cd = getLikelihoods.NB(baySeq.cd, prs = c(0.5,0.5), pET = "BIC", cl = NULL)
# baySeq.posteriors.DE = exp(baySeq.cd@posteriors)[, 2]
# baySeq.table = topCounts(baySeq.cd, group = "DE", FDR = 1)
# baySeq.FDR = baySeq.table$FDR[match(rownames(data),rownames(baySeq.table))]
# ##############
# source("http://bioconductor.org/biocLite.R")
# devel = "http://bioconductor.org/packages/2.13/bioc"
# biocLite("EBSeq", siteRepos = devel)
# library(EBSeq)
# sizes = MedianNorm(data)
# EBSeq.test = EBTest(Data = data, Conditions = factor(class),sizeFactors = sizes, maxround = 10)
# EBSeq.ppmat = GetPPMat(EBSeq.test)
# EBSeq.probabilities.DE = EBSeq.ppmat[, "PPDE"]
# EBSeq.lFDR = 1 - EBSeq.ppmat[, "PPDE"]
# EBSeq.FDR = rep(NA, length(EBSeq.lFDR))
# for (i in 1:length(EBSeq.lFDR)) {
#   EBSeq.FDR[i] = mean(EBSeq.lFDR[which(EBSeq.lFDR <=EBSeq.lFDR[i])])
#                                   }
# ############
# 
# source("TSPM.R")
# TSPM.dgelist = DGEList(counts = count.matrix, group = factor(class))
# TSPM.dgelist = calcNormFactors(TSPM.dgelist, method = "TMM")
# norm.lib.sizes = as.vector(TSPM.dgelist$samples$norm.factors) *as.vector(TSPM.dgelist$samples$lib.size)x0 = rep(1, length(class)), lib.size = norm.lib.sizes)
# TSPM.pvalues = TSPM.test$pvalues
# TSPM.adjpvalues = TSPM.test$padj
# ################
# source('http://bioconductor.org/biocLite.R')
# biocLite('NOISeq')
# library('NOISeq')
# source("noiseq.r")
# class<-c(rep(1,10),rep(2,10))
# nf = calcNormFactors(data)
# libsizes = apply(data, 2, sum)
# common.libsize = prod(libsizes^(1/length(libsizes)))
# normfactors = nf * libsizes/common.libsize
# norm.matrix = sweep(data, 2, normfactors, "/")
# NOISeq.test = noiseq(norm.matrix[, class == 1], norm.matrix[, class == 2], repl = "bio", k = 0.5, norm = "n")
# NOISeq.probabilities = NOISeq.test$probab