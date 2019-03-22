###### GEO
GEO_mRNA_assay<-read.csv('GEO/mRNA/CAV_IIandIII_All126wNormPla_PartekRMA_rawLog2Int_2016Dec09_forBRI.csv')
GEO_mRNA_assay %>% dim
# GEO_mRNA_assay %<>% dplyr::rename(CAIII_c06_911285.H02 = CAIII_c06_911285H02)
GEO_sample_names<-GEO_mRNA_assay %>% names %>% str_sub(-10, -1)
GEO_mRNA_assay_selected<-GEO_mRNA_assay[GEO_mRNA_assay$Probeset.ID %in% mRNA.data_filt$X, c('Probeset.ID', names(GEO_mRNA_assay)[GEO_sample_names %in% Box_sample_names])]
GEO_mRNA_assay_selected[1:5, 1:5]
dim(GEO_mRNA_assay_selected)

Box_sample_names<-mRNA.data_filt %>% names %>% substr(2, 20)
setdiff(Box_sample_names, GEO_sample_names)






############## mRNA AFTER recreated using Elizabeth's code
mRNA_design_GEO<-read.csv(file=paste0(Folder_location,"designAtBaseline.csv"))
mRNA_design_GEO %>% dim # VVV

mRNA_assay_GEO<-read.csv(file=paste0(Folder_location,"filteredAffymetrixData_RemovedDateSent.csv"))
mRNA_assay_GEO %>% dim
mRNA_assay_GEO[1:3, 1:4]

######### Box: Used at Dror's analysis

mRNA.design.file.box<-("CAV Data Sharing/Affy arrays from serum exposure expt - Hessner/Cohort 2/DataForShinyApp/designAtBaseline.csv")
mRNA.design<-read.csv(paste(Active.Box.folder, mRNA.design.file.box, sep=''))
mRNA.design %>% dim # VVV

mRNA.data_filt.file.box<-"CAV Data Sharing/Affy arrays from serum exposure expt - Hessner/Cohort 2/DataForShinyApp/filteredAffymetrixData_RemovedDateSent.csv"
mRNA.data_filt<-read.csv(paste(Active.Box.folder, mRNA.data_filt.file.box, sep=''))
mRNA.data_filt %>% dim
mRNA.data_filt[1:3, 1:4]



setdiff(mRNA.data_filt$X, mRNA_assay_GEO$X)







##### Compare probe id overlap

mRNA.data_filt.file.box<-"CAV Data Sharing/Affy arrays from serum exposure expt - Hessner/Cohort 2/DataForShinyApp/filteredAffymetrixData_RemovedDateSent.csv"
mRNA.data_filt<-read.csv(paste(Active.Box.folder, mRNA.data_filt.file.box, sep=''))
dim(mRNA.data_filt)
probes_Box<-mRNA.data_filt$X %>% as.character # 21,168

### GEO:
affyEset # [54,567, 53]
affyEsetFilter<-nsFilter(affyEset, remove.dupEntrez=FALSE, var.func=IQR, var.cutoff=0.5)$eset
dim(affyEsetFilter)	# now have 21,168 genes to look at


probes_GEO<-featureNames(affyEsetFilter)

setdiff(probes_GEO, probes_Box)
setdiff(probes_Box, probes_GEO)

list(probes_GEO = probes_GEO, probes_Box = probes_Box) %>% gplots::venn()


