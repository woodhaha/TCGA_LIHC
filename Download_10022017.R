##Copy number variation aligned to hg19
##-------------------------------------
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

##Creat data directory
# data_dir="/home/xiaodian/Dropbox/Shirley/UM_Works/TCGA_DATA/LIHC/"
# oldTCGA_dir="/home/xiaodian/Dropbox/Shirley/UM_Works/TCGA_DATA/LIHC/GDC_Legacy_GRCh37(hg19)"

data_dir="/media/H_driver/Shirley/TCGA/LIHC/"
oldTCGA_dir="/media/H_driver/Shirley/TCGA/LIHC/GDC_Legacy_GRCh37(hg19)/"

##Number of each file type
res = files() %>% facet('type') %>% aggregations()
res$type
write.table(res$type[,1],file=paste0(data_dir,"datatype.txt"))

##Project ID
projResults = projects() %>% results_all()
count(projects())
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
LIHC.nocnv = GDCprepare(query.nocnv, save = TRUE, save.filename = paste0(oldTCGA_dir,"CNVs/LIHCnocnvhg19.rda"))


query.cnv= GDCquery(project = ID,
                      data.category = "Copy number variation",
                      platform="Affymetrix SNP Array 6.0",
                      file.type = "hg19.seg",
                      legacy = TRUE)
GDCdownload(query.cnv)
LIHC.cnv = GDCprepare(query.cnv, save = TRUE, save.filename = paste0(oldTCGA_dir,"CNVs/LIHCcnvhg19.rda"))


##(3)DNA methylation 450 data
##---------------------------
query.met450 = GDCquery(project = ID,
                        data.category = "DNA methylation",
                        platform = "Illumina Human Methylation 450",
                        legacy = TRUE)
GDCdownload(query.met450)
LIHC.met450 = GDCprepare(query.met450, save = TRUE, save.filename = paste0(oldTCGA_dir,"met450/LIHCmet450.rda"))

##(4)miRNA Expression data
##------------------------
##mirna
query.mirna = GDCquery(project = ID, 
                       data.category = "Gene expression",
                       data.type = "miRNA gene quantification",
                       file.type = "mirna",
                       legacy = TRUE)
GDCdownload(query.mirna)
LIHC.mirna = GDCprepare(query.mirna, save = TRUE, save.filename = paste0(oldTCGA_dir,"mirna/LIHCmirna.rda"))

##hg19.mirbase20
query.hg19.mirbase20 = GDCquery(project = ID, 
                       data.category = "Gene expression",
                       data.type = "miRNA gene quantification",
                       file.type = "hg19.mirbase20",
                       legacy = TRUE)
GDCdownload(query.hg19.mirbase20)
LIHC.hg19.mirbase20 = GDCprepare(query.hg19.mirbase20, save = TRUE, save.filename = paste0(oldTCGA_dir,"mirna/LIHChg19.mirbase20.rda"))


##(5)Simple nucleotide variation
##------------------------------
##Check maf availables
query.maf.hg19 = GDCquery(project = ID,
                           data.category = "Simple nucleotide variation",
                           data.type = "Simple somatic mutation",
                           access = "open",
                           legacy = TRUE)
# Check maf availables
knitr::kable(getResults(query.maf.hg19)[,c("created_datetime","file_name")])

query.maf.hg19 = GDCquery(project = ID,
                           data.category = "Simple nucleotide variation",
                           data.type = "Simple somatic mutation",
                           access = "open",
                           file.type = "gsc_LIHC_pairs.aggregated.capture.tcga.uuid.automated.somatic.maf",
                           legacy = TRUE)
GDCdownload(query.maf.hg19)
LIHC.mutect.maf = GDCprepare(query.maf.hg19,save = TRUE, save.filename = paste0(oldTCGA_dir,"maf/LIHCmafhg19.rda"))
# LIHC.mutect.maf = GDCquery_Maf("LIHC", pipelines = "mutect2")
# save(LIHC.mutect.maf,file=paste0(oldTCGA_dir,"maf/LIHCmafhg19.rda"))

##(6)Clinical 
##-----------
query.clin = GDCquery(project = ID, 
                  data.category = "Clinical", 
                  data.type = "Clinical data",
                  legacy = TRUE,
                  file.type = "txt")
GDCdownload(query.clin )
LIHC.clin = GDCprepare(query.clin, save = TRUE, save.filename = paste0(oldTCGA_dir,"clin/LIHCclin.rda"))
names(LIHC.clin)

# ##(7)Biospecimen 
# ##--------------
# query.bio = GDCquery(project = ID, 
#                      data.category = "Biospecimen", 
#                      legacy = TRUE)
# GDCdownload(query.bio )
# LIHC.bio = GDCprepare(query.bio, save = TRUE, save.filename = paste0(oldTCGA_dir,"clin/LIHCbio.rda"))
# names(LIHC.bio)

##(7)Protein expression 
##---------------------
query.pro = GDCquery(project = ID,
                  data.category = "Protein expression",
                  platform="MDA_RPPA_Core",
                  legacy = TRUE)
GDCdownload(query.pro)
LIHC.pro = GDCprepare(query.pro, save = TRUE, 
                   save.filename = paste0(oldTCGA_dir,"Protein/LIHCProteinExpression.rda"))

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
LIHC.rsem = GDCprepare(query.rsem, save = TRUE, summarizedExperiment = TRUE, save.filename = paste0(oldTCGA_dir,"exp/rsem/LIHCrsemIllumina_HiSeq.rda"))

# ##RSEM NORMALIZED
# query.rsem.normalized = GDCquery(project = ID,
#                                  data.category = "Gene expression",
#                                  data.type = "Gene expression quantification",
#                                  platform = "Illumina HiSeq", 
#                                  file.type  = "normalized_results", 
#                                  legacy = TRUE)
# GDCdownload(query.rsem.normalized)
# LIHC.rsem.normalized = GDCprepare(query.rsem.normalized, save = TRUE, summarizedExperiment = TRUE, save.filename = paste0(oldTCGA_dir,"exp/rsem.normalized/LIHC.rsem.normalizedIllumina_HiSeq.rda"))

# ##RAW COUNT
# releaseDate <- "2016-01-28"
# setwd(oldTCGA_dir)
# downloadTCGA( cancerTypes = ID, 
#               dataSet = "mRNAseq_Preprocess.Level_3",
#               destDir = "./exp/rawcount", 
#               date = releaseDate )
# 
# folder <- grep( paste0("_", element,"\\.", "|","_",element,"-FFPE)", ".*mRNA_Preprocess"),
#                 list.files("./exp/rawcount"),value = TRUE)
# file <- grep( paste0(element, ".uncv2.mRNAseq_raw_counts"), list.files( file.path( "./exp/rawcount",folder ) ), 
#               value = TRUE)
# path <- file.path( "./exp/rawcount", folder, file )
# assign( value = path, x = paste0(element, ".rnaseq.path"), envir = .GlobalEnv)
# 
# LIHC.rawcount = readTCGA(get(paste0(element,".rnaseq.path"),
#                              envir = .GlobalEnv),
#                          "rnaseq"),
#         x = paste0(element, ".rnaseq"),
#         envir = .GlobalEnv )

# ##(9)Methylation array
# ##---------------------
# # Searching idat file for DNA methylation
# query.met450.microarray = GDCquery(project = ID,
#                                    data.category = "Raw microarray data",
#                                    legacy = TRUE)
# GDCdownload(query.met450.microarray)
# LIHC.met450.microarray = GDCprepare(query.met450.microarray, save = TRUE, 
#                       save.filename = paste0(oldTCGA_dir,"met450.microarray/LIHCmet450.microarray.rda"))

##(9)Clinical from ucsc xena
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
# download.file(url, des_dir)

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

# dir.create(paste0(base_dir, element))
# save_dir = paste0(base_dir,element,"/")

write.csv (DFS_survial, paste0(des_dir,element,"_DFS_survial.csv"),row.names = F)
write.csv (LIHC.clin.ucsc, paste0(des_dir,element,"_clinical.ucsc.xena.csv"),row.names = F)
