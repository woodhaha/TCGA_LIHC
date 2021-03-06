
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
                               
##(5)Simple nucleotide variation
##------------------------------
# Function to get two samples to test the function
legacyPipeline = function(project, data.category, platform){
  query = GDCquery(project = project,
                    data.category = data.category,
                    platform = platform,
                    legacy = TRUE)
  #cases = query$results[[1]]$cases[1:2]
  query = GDCquery(project = project,
                    data.category = data.category,
                    platform = platform,
                    legacy = TRUE)#,
                    #barcode = cases)
  GDCdownload(query)
  data = GDCprepare(query)
  return(data)
}

LIHC.maf = legacyPipeline(ID,"Simple nucleotide variation")


##Check maf availables
query.maf.hg19 = GDCquery(project = ID,
                          data.category = "Simple nucleotide variation",
                          data.type = "Simple somatic mutation",
                          access = "open",
                          legacy = TRUE)

GDCdownload(query.maf.hg19)
LIHC.mutect.maf = GDCprepare(query.maf.hg19)
save(LIHC.mutect.maf,file = paste0(oldTCGA_dir,"maf/LIHC.mutect.maf.rda"))


##(6)Clinical and Biospecimen from GDC
##------------------------------------
query = GDCquery(project = ID, 
                  data.category = "Clinical",
                  legacy = TRUE)
GDCdownload(query)
clinical = GDCprepare_clinic(query, clinical.info = "patient")
clinical.drug = GDCprepare_clinic(query, clinical.info = "drug")
clinical.radiation = GDCprepare_clinic(query, clinical.info = "radiation")
clinical.admin = GDCprepare_clinic(query, clinical.info = "admin")

query = GDCquery(project = ID, 
                  data.category = "Biospecimen",
                  legacy = TRUE)
GDCdownload(query)
clinical.admin = GDCprepare_clinic(query, clinical.info = "admin")
clinical.sample = GDCprepare_clinic(query, clinical.info = "sample")
clinical.slide = GDCprepare_clinic(query, clinical.info = "slide")
clinical.portion = GDCprepare_clinic(query, clinical.info = "portion")

LIHC.clin.GDC = clinical
save(LIHC.clin.GDC,file = paste0(oldTCGA_dir,"clin/LIHC.clin.GDC.rda"))


##(7)Clinical from ucsc xena
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


##(8)Gene expression 
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

##RAW COUNT
releaseDate = "2016-01-28"
setwd(oldTCGA_dir)
downloadTCGA( cancerTypes = ID,
              dataSet = "mRNAseq_Preprocess.Level_3",
              destDir = "./exp/rawcount",
              date = releaseDate )

folder = grep( paste0("_", element,"\\.", "|","_",element,"-FFPE)", ".*mRNA_Preprocess"),
                list.files("./exp/rawcount"),value = TRUE)
file = grep( paste0(element, ".uncv2.mRNAseq_raw_counts"), list.files( file.path( "./exp/rawcount",folder ) ),
              value = TRUE)
path = file.path( "./exp/rawcount", folder, file )
assign( value = path, x = paste0(element, ".rnaseq.path"), envir = .GlobalEnv)

LIHC.rawcount = readTCGA(get(paste0(element,".rnaseq.path"),envir = .GlobalEnv),"rnaseq"),
                         x = paste0(element, ".rnaseq"),
                         envir = .GlobalEnv )
save(LIHC.rawcount,file = paste0(oldTCGA_dir,"exp/LIHC.rawcount.Illumina_HiSeq.rda"))

##(9)Protein expression 
##---------------------
query.pro = GDCquery(project = ID,
                     data.category = "Protein expression",
                     platform="MDA_RPPA_Core",
                     legacy = TRUE)
GDCdownload(query.pro)
LIHC.Protein = GDCprepare(query.pro)
save(LIHC.Protein,file = paste0(oldTCGA_dir,"protein/LIHC.ProteinExpression.rda"))


##(10)Methylation microarray data
##-------------------------------
# Searching idat file for DNA methylation
query.met450.microarray = GDCquery(project = ID,
                                   data.category = "Raw microarray data",
                                   legacy = TRUE)
GDCdownload(query.met450.microarray)
LIHC.met450.microarray = GDCprepare(query.met450.microarray)
save(LIHC.met450.microarray,file = paste0(oldTCGA_dir,"met450.microarray/LIHC.met450.microarray.rda"))


#############################################################
##Downloading TCGA LIHC data aligned to hg38(Harmonized data)
#############################################################
##Creat data directory
data_dir="/media/H_driver/Shirley/TCGA/LIHC/"
newTCGA_dir="/media/H_driver/Shirley/TCGA/LIHC/GDC_Harmonized/"


##(1)Copy number variation 
##------------------------
query = GDCquery(project = ID,
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query)
LIHC.cns = GDCprepare(query)
save(LIHC.cns,file = paste0(oldTCGA_dir,"CNVs/LIHC.cns.rda"))

query = GDCquery(ID,
                  "Copy Number Variation",
                  data.type = "Masked Copy Number Segment") 
GDCdownload(query)
LIHC.mcns = GDCprepare(query)
save(LIHC.mcns,file = paste0(oldTCGA_dir,"CNVs/LIHC.mcns.rda"))

##(2)Simple nucleotide variation 
##------------------------------
LIHC.muse.maf = GDCquery_Maf(element, pipelines = "muse")
LIHC.varscan2.maf = GDCquery_Maf(element, pipelines = "varscan2")
LIHC.somaticsniper.maf = GDCquery_Maf(element, pipelines = "somaticsniper")
LIHC.mutect2.maf = GDCquery_Maf(element, pipelines = "mutect2")
save(LIHC.muse.maf,LIHC.varscan2.maf,LIHC.somaticsniper.maf,LIHC.mutect2.maf,file = paste0(oldTCGA_dir,"CNVs/LIHC.maf.rda"))

##(2)Transcriptome Profiling 
##--------------------------
workflow.type = c("HTSeq - Counts", "HTSeq - FPKM","HTSeq - FPKM-UQ")
for(i in workflow.type){
  print(i)
  query = GDCquery(project =ID,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = i)
  GDCdownload(query)
  LIHC.exp.HTSeq = GDCprepare(query)
}

save(LIHC.exp.HTSeq,file = paste0(oldTCGA_dir,"exp/LIHC.exp.HTSeq.counts.FPKM.FPKM-UQ.rda"))

data.type = c("miRNA Expression Quantification","Isoform Expression Quantification")
#data.type = c("miRNA Expression Quantification")
for(i in data.type){
  print(i)
  query = GDCquery(project = ID,
                    data.category = "Transcriptome Profiling",
                    data.type = i,
                    workflow.type = "BCGSC miRNA Profiling")
  GDCdownload(query)
  LIHC.exp.miRNA.Isoform = GDCprepare(query)
  print(head(LIHC.miRNA))
}
save(LIHC.exp.miRNA.Isoform,file = paste0(oldTCGA_dir,"exp/LIHC.exp.miRNA.Isoform.rda"))


##(3)miRNA Expression Quantification 
##----------------------------------
query.mirna = GDCquery(project = ID, 
                        data.category = "Transcriptome Profiling", 
                        data.type = "miRNA Expression Quantification")


##(3)Clinical and Biospecimen from GDC
##------------------------------------
clin.GDC = GDCquery_clinic(ID, type = "clinical", save.csv = TRUE)
biospecimen.GDC = GDCquery_clinic(ID, type = "biospecimen", save.csv = TRUE)
save(clin.GDC,biospecimen.GDC,file = paste0(oldTCGA_dir,"clin/LIHC.clin.biospecimen.rda"))


