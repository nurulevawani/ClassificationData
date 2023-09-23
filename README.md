# ClassificationData
CClassification data using R Studio

library(affy)
library(GEOquery)
library(Biobase)
library(AnnotationDbi)
#library(simpleaffy)
library(affyPLM)
library(hgu133plus2.db)
library(hgu133acdf)
library(hgu133a.db)
library(hgu133plus2cdf)
library(genefilter)
#library(clariomsmousecdf)
BiocManager::install("pd.clariom.s.mouse.ht")
BiocManager::install("oligo")
library(oligo)


## Call the data in R
gse<- list.celfiles("E:/BIOINFORMATIKA/GSE179717_RAW", full.names=T)
class(gse)
gse
#gse <- read.celfiles(list.celfiles())

affy.data = ReadAffy(filenames=gse)

affy.data

#get pheno data
gset <- getGEO(GEO="GSE179717",GSEMatrix =TRUE)
gset
data.gse <- exprs(gset[[1]])

#BOXPLOT before pre-processing
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE179717", '/', annotation(gset[[1]]), " selected samples", sep ='')
boxplot(exprs(gset[[1]]), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

pheno <- pData(phenoData(gset[[1]]))
varLabels(phenoData(gset[[1]]))
dim(pheno)

#Take alook the Pheno
table(pheno$`age:ch1`)
table(pheno$`gender:ch1`)
table(pheno$`strain:ch1`)
table(pheno$`tissue:ch1`)
table(pheno$`molecule_ch1`)
#table(pheno$`data_row_count`)
#table(pheno$`taxid_ch1`)
table(pheno$`label_ch1`)
#table(pheno$`taxid_ch1`)
pheno$`characteristics_ch1.1`
#table(pheno$`age:ch1`)
#table(pheno$`Sex:ch1`)
#table(pheno$`m stage:ch1`)
#pheno$`m stage:ch1`
#table(pheno$`subgroup:ch1`)
#table(pheno$`ethnic:ch1`)
table(pheno$`characteristics_ch1`)
table(pheno$`characteristics_ch1.1`)
table(pheno$`characteristics_ch1.2`)
table(pheno$`characteristics_ch1.3`)
table(pheno$`characteristics_ch1.4`)
table(pheno$`characteristics_ch1.5`)

#remove subgroup SHH OUTLIER and U
phenob<-pheno[-c(41,51,67),]
phenob
dim(phenob)

#table(phenob$`geo_accession`)

#pie chart characteristics_ch1.4
des <-table(pheno$`characteristics_ch1.4`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c('green','blue',"pink"))

#pie chart Age
des <-table(pheno$`age:ch1`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c('blue',"pink"))

#pie chart gender
des <-table(pheno$`gender:ch1`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c('blue',"pink"))

#pie chart strain
des <-table(pheno$`strain:ch1`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c(1:4))


#pie chart molecule
des <-table(pheno$`molecule_ch1`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c(1:7))

#pie chart label
des <-table(pheno$`label_ch1`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c(1:7))

#pie chart characteristics_ch1
des <-table(pheno$`characteristics_ch1`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c('green','blue',"pink"))

#pie chart characteristics_ch1.1
des <-table(pheno$`characteristics_ch1.1`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c('green','blue',"pink"))

#pie chart characteristics_ch1.2
des <-table(pheno$`characteristics_ch1.2`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c('green','blue',"pink"))

#pie chart characteristics_ch1.3
des <-table(pheno$`characteristics_ch1.3`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c('green','blue',"pink"))

#pie chart characteristics_ch1.5
des <-table(pheno$`characteristics_ch1.5`)
des
persen <- round(des/sum(des)*100)
des <- as.data.frame(des)
lbls <- paste(des$Var1,'-',persen, '%', sep='')
#quartz()
pie(des$Freq, label= lbls, col=c('green','blue',"pink"))

# New data based on phenob
gseb<- gse[-c(41,51,67)]

affy.datab = ReadAffy(filenames=gseb)

affy.datab
#pre processing data
#library(affyPLM)
eset.dChip=threestep(affy.data,background.method = "RMA.2", normalize.method="quantile",summary.method="median.polish")
class(eset.dChip)
dim(eset.dChip)

Ekspres <- exprs(eset.dChip)
class(Ekspres)
dim(Ekspres)

#BOXPLOT after pre-processing
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(pheno$geo_accession))/2),4,2,1))
title <- paste ("GSE179717", '/',annotation(gset[[1]]), " selected samples", sep ='')
boxplot(Ekspres, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)


#filtering
filterdataa <- nsFilter(eset.dChip, require.entrez =T,var.func = IQR, 
                        remove.dupEntrez = T,var.cutoff = 0.5, feature.exclude = "^AFFX")
log <-filterdataa$filter.log
eset <- filterdataa$eset
featureNames(eset) <- make.names(featureNames(eset))

#Filtering results (eset)
dim(eset) 
head(eset)
class(eset)
eset

# create the matrix/data frame
databaru <- exprs(eset)
dim(databaru)
class(databaru)
head(databaru)

#View(databaru)

#There are 4 Subgroup classes: G3, G4, SHH and WNT
phenob$`subgroup:ch1`
datacl <- c(1,1,4,1,3,1,1,1,2,2,2,2,1,1,1,4,2,1,1,
            1,1,1,2,2,1,4,1,1,3,3,3,3,1,2,2,1,2,1,
            1,1,3,3,3,1,1,2,4,1,1,2,1,1,1,1,1,1,2,
            2,4,2,2,4,4,4,1,1,1,1,3,1,3,1,1)
length(datacl)
class(datacl)
datafac <- factor(datacl,levels=1:4, labels= c("G4","G3","SHH","WNT"))
datafac

#Filtering with multtest
#BiocManager::install("multtest")
library(multtest) 
datattest <- mt.teststat(databaru,datafac,test="f") 
class(datattest)
length(datattest)
qqnorm(datattest)
qqline(datattest)

#Adjusted p-value  
rawp = 2 * (1 - pnorm(abs(datattest))) 
prosedur = c("Bonferroni", "Holm", "Hochberg", "BH", "BY") 
adjusted = mt.rawp2adjp(rawp, prosedur) 
data <- adjusted$adjp[,] 
data1 <- data[order(adjusted$index), ] 
head(data1)
dim(data1)

#Take the rawp data column
ffs <- data1[,1] 
class(ffs)
length(ffs)
ffs[1 : 10]
View(ffs)

#Adjusted rawp 
datarawp <- data.frame(databaru, ffs) 
row.names(datarawp) <- row.names(databaru)
class(datarawp)
head(datarawp)
dim(datarawp) 

library(dplyr) 
#datadatarawpfilter <- filter(datarawp, ffs < 0.0005) 
#class(datarawpfilter) dimensi = 96 39
datarawpfilterfinal <- subset(datarawp, ffs < 0.0005)
rownames(datarawpfilterfinal)
class(datarawpfilterfinal)
dim(datarawpfilterfinal) 
head(datarawpfilterfinal)

## Define a new data after filtering
datadef <- datarawpfilterfinal[,1:73]
head(datadef)
dim(datadef)
summary(datadef)
colnames(datadef)

#Preparing data for classification
library(e1071)
library(pROC)
data   = as.data.frame (t((datadef)))
dim(data)
head(data)
#dataY = as.factor(datacl) 
dataY = datafac
dataY

#if it is possible add the demographical variables from phenob data: `age:ch1`, `Sex:ch1`, `m stage:ch1`, `subgroup:ch1`, `ethnic:ch1`
datause = as.data.frame(cbind(data,dataY)) 
datause
dim(datause)
head(datause)
dim(data)

#Classification with SVM kernel Linier
set.seed(978)
rasio = 8/10
train = sample(length(dataY), size = floor(rasio*length(dataY)))
datatrain<- datause[train,]
dim(datatrain)
datatrain
datatest <- datause[-train,]
dim(datatest)

#tuning (best parameter)
tunaslin <- tune(svm, dataY~. ,data = datatrain, kernel="linear",types = "C-clasification",ranges= list( cost = c(0.1,0.01, 0.001,1 , 10 , 100)))
summary(tunaslin)

#Build the model with the best tuning parameter
model <- svm(dataY~. , datatrain, kernel = "linear", cost=0.1,scale=F, types = "C-clasification",decision.value=T)
predictions <- predict(model, datatest)
table(predictions,datatest$dataY)
mean(predictions == datatest$dataY)
predictions <- predict(model, datatrain) 
table(predictions,datatrain$dataY) 
mean(predictions == datatrain$dataY)

#plot(model,datatrain, datatrain$gender~datatrain$ethnic, slice = list(ethnic = 4, gender = 4))
#weighted
w <- t(model$coefs) %*% model$SV        
w <- apply(w, 2, function(v){sqrt(sum(v^2))})   
w <- sort(w, decreasing = T)
x<- as.data.frame(w)
View(x)
all(colnames(w)==colnames(datatrain))

#See the genes
id <- substring(as.character(head(rownames(x), n=10)),2) 
View(id)
