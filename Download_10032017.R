
library(GenomicDataCommons)
library(magrittr)
library(TCGAbiolinks)
library(RTCGAToolbox)
library(SummarizedExperiment)
library(parallel)
library(TCGAWorkflowData)

GenomicDataCommons::status()
##$status
##[1] "OK"

#########################################################
##Downloading TCGA LIHC data aligned to hg19(Legacy data)
#########################################################
##Creat data directory
# data_dir="/home/xiaodian/Dropbox/Shirley/UM_Works/TCGA_DATA/LIHC/"
# oldTCGA_dir="/home/xiaodian/Dropbox/Shirley/UM_Works/TCGA_DATA/LIHC/GDC_Legacy_GRCh37(hg19)"

data_dir="/media/H_driver/Shirley/TCGA/LIHC/"
oldTCGA_dir="/media/H_driver/Shirley/TCGA/LIHC/GDC_Legacy/"

##Number of each file type
res = files() %>% facet('type') %>% aggregations()
res$type
write.table(res$type[,1],file=paste0(oldTCGA_dir,"datatype.txt"))

##Project ID
projResults = projects() %>% results_all()
id=projResults$project_id

##LIVER CANCER(Liver Hepatocellular Carcinoma)
ID=id[1]

##(1)GISTIC results from GDAC firehose
##====================================
library(RTCGAToolbox)

# Download GISTIC results
lastAnalyseDate = getFirehoseAnalyzeDates(1)
getFirehoseDatasets()
gistic = getFirehoseData(dataset="LIHC",gistic2Date = lastAnalyseDate)

# get GISTIC results
gistic.allbygene = getData(gistic,type = "GISTIC", CN = "All")
#gistic.allbygene[1:10,]
gistic.thresholedbygene = getData(gistic,type = "GISTIC", CN = "Thresholded")
#gistic.thresholedbygene = gistic.thresholedbygene[1:10,]
save(gistic.allbygene,gistic.thresholedbygene,file = paste0(oldTCGA_dir,"GISTIC/LIHCGistic.rda"), compress = "xz")

##(2)Copy number variations (CNVs)
##================================
library(TCGAbiolinks)

data_categories=(TCGAbiolinks:::getProjectSummary(ID)$data_categories)$data_category
# [1] "DNA Methylation"             "Clinical"                   
# [3] "Biospecimen"                 "Transcriptome Profiling"    
# [5] "Simple Nucleotide Variation" "Copy Number Variation"      
# [7] "Raw Sequencing Data

query.nocnv= GDCquery(project = ID,
                      data.category = "Copy number variation",
                      platform="Affymetrix SNP Array 6.0",
                      file.type = "nocnv_hg19.seg",
                      legacy = TRUE)
GDCdownload(query.nocnv)
LIHC.nocnvhg19 = GDCprepare(query.nocnv)
save(LIHC.nocnvhg19,file = paste0(oldTCGA_dir,"CNVs/LIHC.nocnvhg19.rda"))


query.cnv= GDCquery(project = ID,
                    data.category = "Copy number variation",
                    platform="Affymetrix SNP Array 6.0",
                    file.type = "hg19.seg",
                    legacy = TRUE)
GDCdownload(query.cnv)
LIHC.cnvhg19 = GDCprepare(query.cnv)
save(LIHC.cnvhg19,file = paste0(oldTCGA_dir,"CNVs/LIHC.cnvhg19.rda"))


##(3)DNA methylation 450 data
##---------------------------
query.met450 = GDCquery(project = ID,
                        data.category = "DNA methylation",
                        platform = "Illumina Human Methylation 450",
                        legacy = TRUE)
GDCdownload(query.met450)
LIHC.met450 = GDCprepare(query.met450)
save(LIHC.met450,file = paste0(oldTCGA_dir,"met450/LIHC.met450.rda"))
LIHC.met450.matrix = assay(LIHC.met450)
save(LIHC.met450.matrix,file = paste0(oldTCGA_dir,"met450/LIHC.met450.matrix.rda"))

##(4)miRNA Expression data
##------------------------
##mirna
query.mirna = GDCquery(project = ID, 
                       data.category = "Gene expression",
                       data.type = "miRNA gene quantification",
                       file.type = "mirna",
                       legacy = TRUE)
GDCdownload(query.mirna)
LIHC.mirna = GDCprepare(query.mirna)
save(LIHC.mirna,file = paste0(oldTCGA_dir,"mirna/LIHC.mirna.rda"))


##hg19.mirbase20
query.hg19.mirbase20 = GDCquery(project = ID, 
                                data.category = "Gene expression",
                                data.type = "miRNA gene quantification",
                                file.type = "hg19.mirbase20",
                                legacy = TRUE)
GDCdownload(query.hg19.mirbase20)
LIHC.hg19.mirbase20 = GDCprepare(query.hg19.mirbase20)
save(LIHC.hg19.mirbase20,file = paste0(oldTCGA_dir,"mirna/LIHC.hg19.mirbase20.rda"))

##(5)Simple nucleotide variation(not work for TCGAbiolinks)
##---------------------------------------------------------
##Check maf availables
query.maf.hg19 = GDCquery(project = ID,
                          data.category = "Simple nucleotide variation",
                          data.type = "Simple somatic mutation",
                          file.type = "bcgsc.ca_LIHC.IlluminaHiSeq_DNASeq.1.somatic.maf",
                          legacy = TRUE)

GDCdownload(query.maf.hg19)
LIHC.mutect.maf = GDCprepare(query.maf.hg19)
save(LIHC.mutect.maf,file = paste0(oldTCGA_dir,"maf/GDC/LIHC.mutect.maf.rda"))


##(6)Clinical from ucsc xena
##--------------------------
element = "LIHC"
link ="https://tcga.xenahubs.net/download/"
pathname = paste0("TCGA.",element,".sampleMap")

base_url = paste0(link,pathname,"/")

names = paste0(element,'_clinicalMatrix.gz')

# Download data
library(curl)
library(RCurl)

des_dir =paste0(oldTCGA_dir,"clin/")
url = paste0(base_url, names)
options('download.file.method'='auto')
download.file(url, destfile = paste0(des_dir,names), method="auto")

##read in and process the RFS clinical data
LIHC.clin.ucsc = read.delim(paste0(des_dir,names),header=TRUE)

time = LIHC.clin.ucsc[,"X_RFS"]
event = LIHC.clin.ucsc[,"X_RFS_IND"]                                                
sampleID = LIHC.clin.ucsc[,"sampleID"]
sampleID = gsub("-", ".",sampleID)
DFS_survial = data.frame(sampleID, time,event)
colnames(DFS_survial) = c("sampleID","DFS.time","DFS.event")

## remove NA
DFS_survial =DFS_survial[!is.na(DFS_survial$DFS.time),]

## Tumor and normal
barcode = DFS_survial[,"sampleID"]

n_idx = which(substr(barcode,14,14) == "1")
t_idx = which(substr(barcode,14,14) == "0")

normal =substr(barcode[n_idx],1,12)
tumor =substr(barcode[t_idx],1,12) 

Normal_or_Tumor = matrix("T",nrow=length(DFS_survial[,"sampleID"]),ncol=1)
colnames(Normal_or_Tumor) = "Normal_or_Tumor"
rownames(Normal_or_Tumor) = DFS_survial[,"sampleID"]

idx = match(intersect(normal,substr(rownames(Normal_or_Tumor),1,12)),substr(rownames(Normal_or_Tumor),1,12))  
Normal_or_Tumor[idx,] = "N"
DFS_survial = cbind(DFS_survial,Normal_or_Tumor)

write.csv (DFS_survial, paste0(des_dir,element,"_DFS_survial.csv"),row.names = F)
write.csv (LIHC.clin.ucsc, paste0(des_dir,element,"_clinical.ucsc.xena.csv"),row.names = F)


##(7)Gene expression 
##------------------
##RSEM
query.rsem = GDCquery(project = ID,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type  = "results", 
                      legacy = TRUE)

GDCdownload(query.rsem)
LIHC.rsem = GDCprepare(query.rsem)
save(LIHC.rsem,file = paste0(oldTCGA_dir,"exp/LIHC.rsem.Illumina_HiSeq.rda"))

##Abstract matrix and gene and sampel information
library(SummarizedExperiment)

## get expression matrix
LIHC.rsem.matrix = assay(LIHC.rsem)

## get genes information
genes.info = rowRanges(LIHC.rsem)

## get sample information
sample.info = colData(LIHC.rsem)
save(LIHC.rsem.matrix,genes.info,sample.info,file = paste0(oldTCGA_dir,"exp/LIHC.rsem.matrix.gene.sample.info.Illumina_HiSeq.rda"))


##RSEM NORMALIZED
query.rsem.normalized = GDCquery(project = ID,
                                 data.category = "Gene expression",
                                 data.type = "Gene expression quantification",
                                 platform = "Illumina HiSeq",
                                 file.type  = "normalized_results",
                                 legacy = TRUE)
GDCdownload(query.rsem.normalized)
LIHC.rsem.normalized = GDCprepare(query.rsem.normalized)
save(LIHC.rsem.normalized,file = paste0(oldTCGA_dir,"exp/LIHC.rsem.normalized.Illumina_HiSeq.rda"))

## get expression matrix
LIHC.rsem.normalized.matrix = assay(LIHC.rsem.normalized)

## get genes information
genes.info = rowRanges(LIHC.rsem.normalized)

## get sample information
sample.info = colData(LIHC.rsem.normalized)
save(LIHC.rsem.normalized.matrix,genes.info,sample.info,file = paste0(oldTCGA_dir,"exp/LIHC.rsem.normalized.matrix.gene.sample.info.Illumina_HiSeq.rda"))

##(8)Protein expression 
##---------------------
query.pro = GDCquery(project = ID,
                     data.category = "Protein expression",
                     platform="MDA_RPPA_Core",
                     legacy = TRUE)
GDCdownload(query.pro)
LIHC.Protein = GDCprepare(query.pro)
save(LIHC.Protein,file = paste0(oldTCGA_dir,"protein/LIHC.ProteinExpression.rda"))

##(9)Clinical and Biospecian
query.clin = GDCquery(project = ID, 
                      data.category = "Clinical",
                      legacy = TRUE)
GDCdownload(query.clin,method="client")
clin = GDCprepare_clinic(query.clin, clinical.info = "patient")
clin.drug = GDCprepare_clinic(query.clin, clinical.info = "drug")
clin.radiation = GDCprepare_clinic(query.clin, clinical.info = "radiation")
clin.admin = GDCprepare_clinic(query.clin, clinical.info = "admin")
save(clin,clin.drug,clin.radiation,clin.admin,file = paste0(newTCGA_dir,"clin/LIHC.clin.GDC.rda"))

query.bios = GDCquery(project = ID, 
                      data.category = "Biospecimen",
                      legacy = TRUE)
GDCdownload(query.bios,method="client")
bios.admin = GDCprepare_clinic(query.bios, clinical.info = "admin")
bios.sample = GDCprepare_clinic(query.bios, clinical.info = "sample")
bios.slide = GDCprepare_clinic(query.bios, clinical.info = "slide")
bios.portion = GDCprepare_clinic(query.bios, clinical.info = "portion")
save(bios.admin,bios.sample,bios.slide,bios.portion,file = paste0(newTCGA_dir,"clin/LIHC.bios.GDC.rda"))


# ##(9)Methylation microarray data(not work by TCGAbiolinks)
# ##---------------------------------------------------------
# # Searching idat file for DNA methylation
# query.met450.microarray = GDCquery(project = ID,
#                                    data.category = "Raw microarray data",
#                                    legacy = TRUE)
# GDCdownload(query.met450.microarray)
# LIHC.met450.microarray = GDCprepare(query.met450.microarray)
# save(LIHC.met450.microarray,file = paste0(oldTCGA_dir,"met450.microarray/LIHC.met450.microarray.rda"))


#############################################################
##Downloading TCGA LIHC data aligned to hg38(Harmonized data)
#############################################################
##Creat data directory
data_dir="/media/H_driver/Shirley/TCGA/LIHC/"
newTCGA_dir="/media/H_driver/Shirley/TCGA/LIHC/GDC_Harmonized/"
setwd("/media/H_driver/Shirley/TCGA/LIHC/Junk/")

##(1)Clinical and Biospecimen from GDC
##------------------------------------
query = GDCquery(project = ID, 
                 data.category = "Clinical")
GDCdownload(query)
clin = GDCprepare_clinic(query, clinical.info = "patient")
clin.drug = GDCprepare_clinic(query, clinical.info = "drug")
clin.radiation = GDCprepare_clinic(query, clinical.info = "radiation")
clin.admin = GDCprepare_clinic(query, clinical.info = "admin")
save(clin,clin.drug,clin.radiation,clin.admin,file = paste0(newTCGA_dir,"clin/LIHC.clin.GDC.rda"))

query = GDCquery(project = ID, 
                 data.category = "Biospecimen")
GDCdownload(query)
bios.admin = GDCprepare_clinic(query, clinical.info = "admin")
bios.sample = GDCprepare_clinic(query, clinical.info = "sample")
bios.slide = GDCprepare_clinic(query, clinical.info = "slide")
bios.portion = GDCprepare_clinic(query, clinical.info = "portion")
save(bios.admin,bios.sample,bios.slide,bios.portion,file = paste0(newTCGA_dir,"clin/LIHC.bios.GDC.rda"))

##(2)Copy number variation 
##------------------------
query = GDCquery(project = ID,
                 data.category = "Copy Number Variation",
                 data.type = "Copy Number Segment")
GDCdownload(query)
LIHC.cns = GDCprepare(query)
save(LIHC.cns,file = paste0(newTCGA_dir,"CNVs/LIHC.cns.rda"))

query = GDCquery(ID,
                 data.category = "Copy Number Variation",
                 data.type = "Masked Copy Number Segment") 
GDCdownload(query)
LIHC.mcns = GDCprepare(query)
save(LIHC.mcns,file = paste0(newTCGA_dir,"CNVs/LIHC.mcns.rda"))

##(3)Simple nucleotide variation 
##------------------------------
LIHC.muse.maf = GDCquery_Maf(element, pipelines = "muse")
LIHC.varscan2.maf = GDCquery_Maf(element, pipelines = "varscan2")
LIHC.somaticsniper.maf = GDCquery_Maf(element, pipelines = "somaticsniper")
LIHC.mutect2.maf = GDCquery_Maf(element, pipelines = "mutect2")
save(LIHC.muse.maf,LIHC.varscan2.maf,LIHC.somaticsniper.maf,LIHC.mutect2.maf,file = paste0(newTCGA_dir,"maf/LIHC.all.maf.rda"))

##(4)Transcriptome Profiling 
##--------------------------
workflow.type = c("HTSeq - Counts", "HTSeq - FPKM","HTSeq - FPKM-UQ")
for(i in workflow.type){
  print(i)
  query = GDCquery(project =ID,
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = i)
  GDCdownload(query)
  LIHC.HTSeq = GDCprepare(query)
  save(LIHC.HTSeq,file = paste0(newTCGA_dir,"Transcriptome_Profiling/","LIHC.",i,".rda"))
}

##(5)miRNA Expression Quantification 
##----------------------------------
query.mirna = GDCquery(project = ID, 
                       data.category = "Transcriptome Profiling", 
                       data.type = "miRNA Expression Quantification")
GDCdownload(query.mirna )
LIHC.mirna = GDCprepare(query.mirna)
save(LIHC.mirna,file = paste0(newTCGA_dir,"Transcriptome_Profiling/LIHC.mirna.rda"))
     
##(6)Isoform Expression Quantification 
##------------------------------------
query.isoform = GDCquery(project = ID, 
                       data.category = "Transcriptome Profiling", 
                       data.type = "Isoform Expression Quantification")
GDCdownload(query.isoform)
LIHC.isoform = GDCprepare(query.isoform)
save(LIHC.isoform,file = paste0(newTCGA_dir,"Transcriptome_Profiling/LIHC.isoform.rda"))


