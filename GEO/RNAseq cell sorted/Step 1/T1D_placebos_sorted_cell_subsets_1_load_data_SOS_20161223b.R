## scripts to load and QC data for T1D placebos sorted cell subsets project ##

##### set up environment: load/save data #####

# rm(list=ls())
# opar <- par()
setwd("/Users/samskinner/Box Sync/Projects/CAV-Cohort-2-RNAseq-CellSubsets/DufortAnalysis/")

## load data for use in these scripts
# load(file="data/T1D_placebos_sorted_cell_subsets_data_2016-06-30.Rdata")  # load R objects from saved environment
# save(file="data/T1D_placebos_sorted_cell_subsets_data_2016-06-30.Rdata", list=ls())  # save R objects in environment


##### set up environment: load packages #####

## load general packages
library(xlsx)
library(dplyr)
library(stringr)
library(ggplot2)
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1)))
library(ggthemes)

## load analysis-specific packages
library(dendextend)
library(limma)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(PCGSE)
library(readr)

# load my packages needed for these analyses
library(RNAseQC)
library(countSubsetNorm)
library(miscHelpers)
library(geneSetTools)


##### load data from other sources (usually only do this once, when first running the analyses) #####

## load data from whole blood placebos data processing
load(file="T1D_placebos_data_for_cell_subset_analyses_2016-07-21.Rdata")


##### Read in count values #####

# read counts files
counts_C8BBUANXX.tmp <- read.csv("../data_AC8BBUANXX/P133-1_C8BBUANXX_160610_combined_counts.csv")
counts_C8BUUANXX.tmp <- read.csv("../data_BC8BUUANXX/P133-1_C8BUUANXX_160610_combined_counts.csv")
counts_C8HAWANXX.tmp <- read.csv("../data_AC8HAWANXX/P133-1_C8HAWANXX_160610_combined_counts.csv")
# str(counts_C8BBUANXX.tmp)

# Merge counts files
counts.merged <- counts_C8BBUANXX.tmp %>%
  merge(counts_C8BUUANXX.tmp, by ="geneName") %>%
  merge(counts_C8HAWANXX.tmp, by = "geneName")
rownames(counts.merged) <- counts.merged$geneName
counts.merged <- counts.merged[, -which(colnames(counts.merged) == "geneName")]

# Trim Lib ID's to lib#### (lib plus numbers)
colnames(counts.merged) <- colnames(counts.merged) %>%
  str_extract(pattern="lib[0-9]+") %>%
  make.unique(sep="_")

# Keep protein coding genes with HGNC symbols, and drop non-protein-coding genes
counts.merged.tmp <- counts.merged
grch38.tmp <- grch38[grch38[["biotype"]] %in% "protein_coding",]
counts.merged.tmp <- counts.merged.tmp[which(rownames(counts.merged.tmp) %in% grch38.tmp$ensgene),]
counts.merged.tmp$ensgene <- rownames(counts.merged.tmp)

## sum counts.merged for duplicated HGNC symbols (why are there duplicates, and why were they being averaged rather than summed?)
# this also drops rows with HGNC.symbols==NA, which should include any non-protein-coding genes, as we populated the names
counts.merged.aggregated <-
  aggregate(
    counts.merged.tmp[,1:(ncol(counts.merged.tmp)-1)],
    by=list(counts.merged.tmp$ensgene), sum)
rownames(counts.merged.aggregated) <- counts.merged.aggregated$Group.1
counts.merged.aggregated <-
  counts.merged.aggregated[,-which(colnames(counts.merged.aggregated) == "Group.1")]
dim(counts.merged.aggregated)
# 19751 genes, 420 libraries

rm_tmp(ask=FALSE)  # clean up unneeded variables


##### plot saturation curves of read counts #####

saturation.merged.aggregated <-
  estimate_saturation(counts.merged.aggregated, method="division", ndepths=30, min_counts=1)
plot_saturation_curve(saturation.merged.aggregated, plot_lines=TRUE, plot_terminal_points=TRUE)
# looks good! reaching saturation at around 13,000 genes at 3,000,000 reads in all libraries
# a few libraries with very few reads, which will likely get dropped during QC


##### Read in library prep data #####

sample_annotation.full <-
  read.xlsx(file="../data_BC8BUUANXX/P133-1 Final Annotation w-correction.xlsx",
            sheetIndex=1) %>%
  dplyr::select(-starts_with("NA")) # drop NA columns

# drop rows with all NAs
sample_annotation.full <- sample_annotation.full %>%
  filter(rowSums(!is.na(sample_annotation.full)) > 0)

# fix column names
colnames(sample_annotation.full) <- standardize_names(colnames(sample_annotation.full))

sample_annotation.full <- sample_annotation.full %>%
  plyr::rename(replace=c("group"="study",
                         "week"="weeks",
                         "visit"="visit_number"))

# check that duplicated columns have identical values, then drop them
for (i in grep("_1", colnames(sample_annotation.full), value=TRUE)) {
  print(all(sample_annotation.full[[i]]==
              sample_annotation.full[[str_replace(i, "_1", "")]]))
}
sample_annotation.full <-
  sample_annotation.full[, grep("_1", colnames(sample_annotation.full), invert=TRUE)]

# add a "participant_id" column; need to pull in information from other object
sample_annotation.full$participant_id <- as.character(NA)
sample_annotation.full$participant_id[sample_annotation.full$study=="START"] <-
  with(sample_annotation.full[sample_annotation.full$study=="START",],
       paste(study, str_replace_all(subject_id, "-", ""), sep="_"))
sample_annotation.full$participant_id[sample_annotation.full$study=="T1DAL"] <-
  patient_data.by_study[["T1DAL"]]$participant_id[
    match(sample_annotation.full$subject_id[sample_annotation.full$study=="T1DAL"],
          patient_data.by_study[["T1DAL"]]$id)]
sample_annotation.full[,c("study","subject_id","participant_id")]

# pull in additional data from patient_data
sample_annotation.full <-
  merge(sample_annotation.full, patient_data.merged[,2:ncol(patient_data.merged)],
        by="participant_id", all.x=TRUE)

# add visit_id column for pulling data from master
sample_annotation.full$visit_id <-
  with(sample_annotation.full,
       paste(participant_id, visit_number, sep="_"))

# merge in cpeptide data
sample_annotation.full <-
  merge(sample_annotation.full,
        master[, setdiff(colnames(master), colnames(sample_annotation.full))],
        by.x="visit_id", by.y="rnaseq_visit_id", all.x=TRUE)
sample_annotation.full <- sample_annotation.full %>%
  plyr::rename(replace=c("libid"="library_id_whole_blood"))

# might also want to bring in CBC data?


##### Read in library metrics #####

metrics_C8BBUANXX.tmp <- read.csv("../data_AC8BBUANXX/P133-1_C8BBUANXX_160610_combined_metrics.csv")
metrics_C8BUUANXX.tmp <- read.csv("../data_BC8BUUANXX/P133-1_C8BUUANXX_160610_combined_metrics.csv")
metrics_C8HAWANXX.tmp <- read.csv("../data_AC8HAWANXX/P133-1_C8HAWANXX_160610_combined_metrics.csv")
# str(metrics_C8BBUANXX.tmp)

metrics.merged <-
  rbind(metrics_C8BBUANXX.tmp, metrics_C8BUUANXX.tmp, metrics_C8HAWANXX.tmp)
# glimpse(metrics.merged)
colnames(metrics.merged) <- str_to_lower(colnames(metrics.merged))
for (i in 2:ncol(metrics.merged)) {
  if (!is.numeric(metrics.merged[,i])) {
    if (sum(str_detect(metrics.merged[,i], "%")) > 0) {
      metrics.merged[,i] <- as.numeric(str_replace(metrics.merged[,i], "%", "")) / 100
    } else metrics.merged[,i] <- as.numeric(metrics.merged[,i])
  }
}
metrics.merged <- dplyr::rename(metrics.merged, lib_id = libid)

# Trim Lib ID's to lib#### (lib plus 4 digit)
metrics.merged$lib_id <- metrics.merged$lib_id %>%
  str_extract("lib[0-9]+") %>%
  make.unique(sep="_")
# which(duplicated(metrics.merged$lib_id))

rm_tmp(ask=FALSE)


##### Remove problematic libraries #####

#Plot total counts
metrics.merged <- arrange(metrics.merged, fastq_total_reads)

plot_read_counts(metrics.merged, file_prefix="plots/T1D_sorted_cell_subsets",
                 n_lowcount=30, id_col="lib_id")
# about 15 low-count libraries

## plot total reads, mapped_reads_w_dups, median_cv_coverage all against each other
plot_metrics(
  metrics.merged, metrics.libID_col="lib_id",
  design=sample_annotation.full, design.libID_col="library_id",
  threshold.perc_aligned=0.7,
  point_size=2,
  by_var="cell_type",
  file_prefix="plots/T1D_sorted_cell_subsets_all_libs", plotdims=c(11,9),
  by_var_levels=c("Bcell", "CD4", "CD8", "monocyte", "PBMC"),
  my_cols=colorblind_pal()(5))
# a bunch of libraries that have relatively low percent alignment, but otherwise look good
# it turns out that the ones between 0.45 and 0.7 are more problematic, and only a few in number
# lowered percent alignment threshold to 0.7 to keep those libraries
# most libraries look good; some obviously bad, appear to be a mix of cell types

# ## plot total reads, mapped_reads_w_dups, median_cv_coverage all against each other by different variables
# for (i in colnames(sample_annotation.full)[2:19])
#   plot_metrics(metrics.merged, metrics.libID_col="lib_id",
#                design=sample_annotation.full, design.libID_col="library_id",
#                threshold.percent_aligned=0.45,
#                by_var=i,
#                file_prefix="plots/T1D_sorted_cell_subsets_all_libs", plotdims=c(11,9))
# # nothing apparent with any of these variables; for some reason I cannot figure out, some libraries have lots of junk reads

# in very low alignment (~20%) libraries, over-represented sequences look like reagents (they have lots of perfect hits in GenBank, from lots of different organisms)
# in medium alignment (~60%) libraries, over-represented sequences look like reagents (they have lots of perfect hits in GenBank, from lots of different organisms)

# problematic libraries: lib6373, lib6446, lib6405, lib6383, lib6645, lib6310, lib6428
lib.bad.tmp <-
  metrics.merged$lib_id[
    with(metrics.merged,
         fastq_total_reads < 5e6 | mapped_reads_w_dups < 0.7 | median_cv_coverage > 1.0)]

# inspect metrics and sample data for problematic libraries
metrics.merged[metrics.merged$lib_id %in% lib.bad.tmp,]
sample_annotation.full[sample_annotation.full$library_id %in% lib.bad.tmp,1:10]
# nothing obvious that they have in common; no ABATE because it wasn't included!!!

# check sample percent duplication
metrics.merged[order(metrics.merged$percent_duplication, decreasing=F),
               c("lib_id", "percent_duplication")]
# some with really high duplication, but look good otherwise; weird


## Make quality control cuts based on metrics
metrics.merged.qc <-
  metrics.merged[
    with(metrics.merged,
         fastq_total_reads > 5e6 & mapped_reads_w_dups > 0.7 & median_cv_coverage < 1.0),]
nrow(metrics.merged.qc)  # 380 libraries (QC removed 40, or 9.5% of libraries)

# Remove from counts data libraries that fail QC cuts
counts.merged.aggregated.qc <-
  counts.merged.aggregated[
    , colnames(counts.merged.aggregated) %in% metrics.merged.qc$lib_id]
counts.merged.aggregated.qc <-
  counts.merged.aggregated.qc[,order(colnames(counts.merged.aggregated.qc))]

# Remove from sample annotation data libraries that fail QC cuts
sample_annotation.full.qc <-
  sample_annotation.full[sample_annotation.full$library_id %in%
                           colnames(counts.merged.aggregated.qc),]
sample_annotation.full.qc <-
  sample_annotation.full.qc[order(sample_annotation.full.qc$library_id),]


##### some basic checks for problematic data #####

## calculate total number of reads mapping to X and to Y chromosome, for each library
logXY.tmp <- logXYratio(counts.merged.aggregated.qc, gene_ID="symbol")


hist(logXY.tmp, breaks=40, main="log-ratio of X reads to Y reads")
sort(logXY.tmp) # obvious break between 5 and 9

pdf("plots/logXYratio.pdf", width=8.5, height=6)
ggplot(data=data.frame(logXYratio=logXY.tmp), mapping=aes(x=logXYratio) )+
  geom_histogram(colour="black", fill="gray60")
dev.off()


# examine list of libraries by assigned sex and logXYratio
cbind(sample_annotation.full.qc$sex,
      logXY.tmp[match(sample_annotation.full.qc$library_id, names(logXY.tmp))])
# looks like libraries 11600-11606 are supposed to be male, but actually female
# looks like libraries 11577-11582 are supposed to be female, but are actually male
# perhaps a row of 8 swapped for another row of 8 somewhere in the process?

# assign sex based on total number of reads mapping to X and Y chromosomes, for each library
sample_annotation.full.qc$sex_by_rna <-
  ifelse(logXY.tmp[sample_annotation.full.qc$library_id] > 7, "F", "M")
sample_annotation.full.qc[,c("library_id","sex","sex_by_rna")]

# check that inferred sex is same for all libraries for each patient
table(sample_annotation.full.qc[,c("participant_id", "sex_by_rna")])
which(rowSums(table(sample_annotation.full.qc[,c("participant_id", "sex_by_rna")])==0)!=1)
# 4 patients have libraries called male and female

# plot them
data.tmp <- data.frame(
  sex=sample_annotation.full.qc$sex[
    match(names(logXY.tmp), sample_annotation.full.qc$library_id)],
  logXYratio=logXY.tmp)

pdf("plots/logXYratio_by_sex_expected_inferred.pdf", width=9, height=6)
ggplot(data=data.tmp, mapping=aes(x=logXYratio, fill=sex)) +
  geom_histogram(position="dodge") +
  scale_fill_manual(values=c("red", "blue"))
dev.off()

# check that inferred sex matches sex from database
all(sample_annotation.full.qc$sex == sample_annotation.full.qc$sex_by_rna)
# nope; identify the ones that look wrong
sample_annotation.full.qc[
  which(sample_annotation.full.qc$sex != sample_annotation.full.qc$sex_by_rna),]


## check logXYratio for all libraries, including those removed due to low quality
logXY.tmp <- logXYratio(counts.merged.aggregated, gene_ID="symbol")
hist(logXY.tmp, breaks=15, main="log-ratio of X reads to Y reads")
sort(logXY.tmp) # obvious break between 6.1 and 8.8

pdf("plots/logXYratio.all_libs.pdf", width=8.5, height=6)
ggplot(data=data.frame(logXYratio=logXY.tmp), mapping=aes(x=logXYratio) )+
  geom_histogram(colour="black", fill="gray60")
dev.off()

# assign sex based on total number of reads mapping to X and Y chromosomes, for each library
sample_annotation.full$sex_by_rna <-
  ifelse(logXY.tmp[sample_annotation.full$library_id] > 7, "F", "M")
sample_annotation.full[,c("library_id","sex","sex_by_rna")]

# histogram by sex
data.tmp <- data.frame(
  sex=sample_annotation.full$sex[
    match(names(logXY.tmp), sample_annotation.full$library_id)],
  logXYratio=logXY.tmp)

pdf("plots/logXYratio_by_sex_expected_inferred.all_libs.pdf", width=9, height=6)
ggplot(data=data.tmp, mapping=aes(x=logXYratio, fill=sex)) +
  geom_histogram(position="dodge") +
  scale_fill_manual(values=c("red", "blue"))
dev.off()

glimpse(sample_annotation.full)
cbind(sample_annotation.full$library_id, sample_annotation.full$sex,
      logXY.tmp[sample_annotation.full$library_id])
# some additional mismatches (lib11731, lib11716, lib11635, lib11597, lib11410), but they're all really low-quality libraries; set them to NA
sample_annotation.full$sex_by_rna[
  sample_annotation.full$library_id %in%
    c("lib11731", "lib11716", "lib11635", "lib11597", "lib11410")] <- NA


# check UTY expression for those libraries with mis-matched sex
counts.merged.aggregated.qc["UTY",
                            sample_annotation.full.qc$library_id[
  which(sample_annotation.full.qc$sex != sample_annotation.full.qc$sex_by_rna)],]

table(sample_annotation.full.qc[,c("study","sex_by_rna")])
# some imbalance in sex by study, which could mess up the sex-linked genes if I correct by study and not Sex

## extract library names to include in tests
libs.test_sex.tmp <- colnames(counts.merged.aggregated)

counts.test_sex.tmp <-
  data.frame(lib_id=libs.test_sex.tmp,
             sex=factor(sample_annotation.full$sex[
               match(libs.test_sex.tmp, sample_annotation.full$library_id)]),
             UTY=log(unlist(counts.merged[match("ENSG00000183878", rownames(counts.merged)), libs.test_sex.tmp])+0.5),
             XIST=log(unlist(counts.merged[match("ENSG00000229807", rownames(counts.merged)), libs.test_sex.tmp])+0.5))

counts.test_sex.tmp <- counts.test_sex.tmp %>%
  group_by(sex) %>%
  mutate(UTY_outlier = ifelse(is_outlier(UTY), lib_id, NA)) %>%
  mutate(XIST_outlier = ifelse(is_outlier(XIST), lib_id, NA))


## scatterplot of points on UTY and XIST to identify sex (with points colored by sex)
pdf(file="plots/sex_by_UTY_XIST.pdf", h=9,w=9)
ggplot(data=counts.test_sex.tmp, aes(x=XIST, y=UTY, colour=sex)) +
  geom_point(position=position_jitter(w=0.02,h=0.02), size=3) +
  geom_text(aes(label=UTY_outlier), colour="black", vjust=-1.3) +
  geom_text(aes(label=XIST_outlier), colour="black", vjust=-1.3) +
  scale_color_manual(values=c("red", "blue")) +
  labs(x="log10 XIST expression", y="log10 UTY expression")
dev.off()
# interesting; lib11577 has high XIST and UTY counts, suggesting it may be a mixed sample

## boxplots of X and Y chromosome genes to check sex
ggplot(data=counts.test_sex.tmp, aes(x = sex, y=UTY)) +
  # geom_boxplot(position="jitter") +
  geom_boxplot() +
  labs(x="Sex", y="log10 UTY expression")

ggplot(data=counts.test_sex.tmp, aes(x = sex, y=XIST)) +
  # geom_boxplot(position="jitter") +
  geom_boxplot() +
  labs(x="Sex", y="log10 XIST expression")


## drop problematic libraries, and further filter quality-controlled master and counts objects
# skipped this step 2016-07-01, because I want to be able to identify the libraries where the sex is wrong!

# sample_annotation.full.qc <-
#   sample_annotation.full.qc[(sample_annotation.full.qc$sex_by_rna %in% c("F","M")) &
#                        ((sample_annotation.full.qc$sex_by_rna == sample_annotation.full.qc$sex) |
#                           (sample_annotation.full.qc$sex == "")),]
# sample_annotation.full.qc <- sample_annotation.full.qc[sample_annotation.full.qc$libID %nin% bad.libs,]

# remove problematic libraries from counts file
# counts.merged.aggregated.qc <-
#   counts.merged.aggregated.qc[,match(sample_annotation.full.qc$library_id, colnames(counts.merged.aggregated.qc))]

rm_tmp(ask=FALSE)  # remove unneeded variables


# remove problematic libraries from counts file


##### output QC'd count data as csv files (should be normalized later in the process) #####

write.csv(counts.merged.aggregated.qc, file="output/counts_postQC.csv")


##### run PCA on (normalized) full data set to check for bias and batch effects #####

# Remove libraries not in sample annotation object,
# filter out genes that have a count of at least one in < 15% of libraries,
# and normalize counts across libraries
counts.merged.aggregated.qc.normalized <-
  calcNormCounts(counts=counts.merged.aggregated.qc,
                 min_count=1,
                 design=sample_annotation.full,
                 libID_col="library_id")
# glimpse(counts.merged.aggregated.qc.normalized)

## run PCA on filtered, normalized counts

# run PCA
# only need pcaAll for the scores, the sds, and the plot
pcaAll.full <-
  calc_PCAs(counts.merged.aggregated.qc.normalized, log2_transform=TRUE)

# scree plot, the number of informative PCs = elbow
plot(pcaAll.full, type="l")
# 3 informative PCs

# attach sample info to PC scores  
scores_sample_annotation_pcaAll.full <-
  merge(sample_annotation.full.qc,
        as.data.frame(pcaAll.full$x), by.x = "library_id", by.y="row.names")


### Color plots of PCs by different variables

## plots of PCs1 1-3, with no color
plot_PCAs(scores_sample_annotation_pcaAll.full,
          pvars.labs=pcaAll.full$pvars.labs, PCs = 1:3,
          file_prefix="plots/T1D_placebos_sorted_cells.all_libs", plotdims=c(7.3,6))

## plots of PC1s 1-3, colored by cell type
plot_PCAs(scores_sample_annotation_pcaAll.full,
          pvars.labs=pcaAll.full$pvars.labs, PCs = 1:3,
          color_by_var="cell_type", color_by_var_levels=c("CD4","CD8","Bcell","monocyte","PBMC"),
          my_cols=colorblind_pal()(5),
          file_prefix="plots/T1D_placebos_sorted_cells.all_libs", plotdims=c(9,6))

# plot PC1 vs. PC2, colored by cell type, with the problem libraries labeled in their color
pdf("plots/T1D_placebos_sorted_cells.all_libs.PC1_vs_PC2.color_by_cell_type.text_by_library_id.pdf",
    width=9, height=6)
ggplot(data=scores_sample_annotation_pcaAll.full,
       mapping=aes(x=PC1, y=PC2, color=cell_type)) +
  geom_point(alpha=0.8, size=3) +
  geom_text(data=scores_sample_annotation_pcaAll.full[
    scores_sample_annotation_pcaAll.full$library_id %in%
      paste0("lib", c(11577:11583,11600:11606)),],
    mapping=aes(label=library_id),
    nudge_y=5, size=5,
    show.legend=FALSE) +
  geom_text(data=scores_sample_annotation_pcaAll.full[
    scores_sample_annotation_pcaAll.full$PC2 > 50 &
      scores_sample_annotation_pcaAll.full$PC2 < 75,],
    mapping=aes(label=library_id),
    nudge_y=5, size=5,
    show.legend=FALSE) +
  scale_color_manual(values=setNames(colorblind_pal()(5), c("CD4","CD8","Bcell","monocyte","PBMC")))
dev.off()

# cluster libraries on PC1 and PC2, so that I can flag cell_type outliers
hclust.PC12.full <- hclust(dist(scores_sample_annotation_pcaAll.full[,c("PC1","PC2")]))

pdf("plots/cluster_dendrogram_PC12.label_by_cell_type.pdf", width=10, height=6)
plot(hclust.PC12.full, labels=scores_sample_annotation_pcaAll.full$cell_type, cex=0.5)
dev.off()

# plot(hclust.PC12.full, labels=scores_sample_annotation_pcaAll.full$library_id, cex=0.5)

clusters.PC12.full <- cutree(hclust.PC12.full, h=125)
plot(hclust.PC12.full, labels=clusters.PC12.full, cex=0.5)

# check the cell types that clusters correspond to
table(data.frame(scores_sample_annotation_pcaAll.full$cell_type, clusters.PC12.full))

cluster.names <- c("CD4/CD8","Bcell","monocytes","PBMC")

scores_sample_annotation_pcaAll.full$cell_type_by_rna <-
  cluster.names[clusters.PC12.full]


##### construct an object with library information for deconvoluting problem libraries #####

library_qc_summary <-
  data.frame(library_id=sample_annotation.full$library_id,
             trna_sampleid=sample_annotation.full$trna_sampleid,
             sex_by_annotation=sample_annotation.full$sex,
             sex_by_rna=sample_annotation.full$sex_by_rna,
             cell_type_by_annotation=sample_annotation.full$cell_type)
library_qc_summary$cell_type_by_rna <-
  scores_sample_annotation_pcaAll.full$cell_type_by_rna[match(
    library_qc_summary$library_id, scores_sample_annotation_pcaAll.full$library_id)]
library_qc_summary$sex_match <-
  library_qc_summary$sex_by_annotation == library_qc_summary$sex_by_rna
library_qc_summary$cell_type_match <-
  (library_qc_summary$cell_type_by_annotation == library_qc_summary$cell_type_by_rna) |
  (library_qc_summary$cell_type_by_annotation %in% c("CD4","CD8") &
     library_qc_summary$cell_type_by_rna == "CD4/CD8")

library_qc_summary <- library_qc_summary[order(library_qc_summary$trna_sampleid),]

write.csv(library_qc_summary,
          file="output/T1D_placebos_cell_subsets_match_annotation_RNA.csv")


##### create a corrected version of the sample annotation object, and rerun the stuff above #####

sample_annotation.full.corrected <-
  sample_annotation.full[,-which(colnames(sample_annotation.full)=="sex_by_rna")]
sample_annotation.full.corrected <-
  sample_annotation.full.corrected[order(sample_annotation.full.corrected$trna_sampleid),]

# flip the library_id and trna_sampleid for the 24 problem samples (tRNAs 26465-26488)
rows.tmp <-
  which(sample_annotation.full.corrected$trna_sampleid %in%
          paste0("tRNA", 26465:26488))
sample_annotation.full.corrected[rows.tmp, c("library_id", "trna_sampleid")] <-
  sample_annotation.full.corrected[rev(rows.tmp), c("library_id", "trna_sampleid")]

# redo the sex_by_rna, normalization, filtering, and PCA with the corrected version

# Remove from sample annotation data libraries that fail QC cuts
sample_annotation.full.corrected.qc <-
  sample_annotation.full.corrected[sample_annotation.full.corrected$library_id %in%
                           colnames(counts.merged.aggregated.qc),]
sample_annotation.full.corrected.qc <-
  sample_annotation.full.corrected.qc[order(sample_annotation.full.corrected.qc$library_id),]

# Remove from sample annotation data libraries that fail QC cuts
sample_annotation.full.uncorrected.qc <-
  sample_annotation.full[sample_annotation.full$library_id %in%
                           colnames(counts.merged.aggregated.qc),]
sample_annotation.full.uncorrected.qc <-
  sample_annotation.full.uncorrected.qc[
    order(sample_annotation.full.uncorrected.qc$library_id),]


##### repeat basic checks for problematic data #####

## calculate total number of reads mapping to X and to Y chromosome, for each library
logXY.tmp <- logXYratio(counts.merged.aggregated.qc, gene_ID="symbol")
hist(logXY.tmp, breaks=40, main="log-ratio of X reads to Y reads")
sort(logXY.tmp) # obvious break between 5 and 9

# examine list of libraries by assigned sex and logXYratio
cbind(sample_annotation.full.corrected.qc$sex,
      logXY.tmp[match(sample_annotation.full.corrected.qc$library_id, names(logXY.tmp))])
# libraries look right now!

# assign sex based on total number of reads mapping to X and Y chromosomes, for each library
sample_annotation.full.corrected.qc$sex_by_rna <-
  ifelse(logXY.tmp[sample_annotation.full.corrected.qc$library_id] > 7, "F", "M")
sample_annotation.full.corrected.qc[,c("library_id","sex","sex_by_rna")]

# check that inferred sex is same for all libraries for each patient
table(sample_annotation.full.corrected.qc[,c("participant_id", "sex_by_rna")])
which(rowSums(table(sample_annotation.full.corrected.qc[,c("participant_id", "sex_by_rna")])==0)!=1)
# now they all look right!!!

# plot them
data.tmp <- data.frame(
  sex=sample_annotation.full.corrected.qc$sex[
    match(names(logXY.tmp), sample_annotation.full.corrected.qc$library_id)],
  logXYratio=logXY.tmp)

pdf("plots/logXYratio_by_sex_expected_inferred.corrected.pdf", width=9, height=6)
ggplot(data=data.tmp, mapping=aes(x=logXYratio, fill=sex)) +
  geom_histogram(position="dodge") +
  scale_fill_manual(values=c("red", "blue"))
dev.off()

# check that inferred sex matches sex from database
all(sample_annotation.full.corrected.qc$sex == sample_annotation.full.corrected.qc$sex_by_rna)
# looks good now!

## check logXYratio for all libraries, including those removed due to low quality
logXY.tmp <- logXYratio(counts.merged.aggregated, gene_ID="symbol")
hist(logXY.tmp, breaks=15, main="log-ratio of X reads to Y reads")
sort(logXY.tmp) # obvious break between 6.1 and 8.8

# assign sex based on total number of reads mapping to X and Y chromosomes, for each library
sample_annotation.full.corrected$sex_by_rna <-
  ifelse(logXY.tmp[sample_annotation.full.corrected$library_id] > 7, "F", "M")
sample_annotation.full.corrected[,c("library_id","sex","sex_by_rna")]

# histogram by sex
data.tmp <- data.frame(
  sex=sample_annotation.full.corrected$sex[
    match(names(logXY.tmp), sample_annotation.full.corrected$library_id)],
  logXYratio=logXY.tmp)

pdf("plots/logXYratio_by_sex_expected_inferred.all_libs.corrected.pdf", width=9, height=6)
ggplot(data=data.tmp, mapping=aes(x=logXYratio, fill=sex)) +
  geom_histogram(position="dodge") +
  scale_fill_manual(values=c("red", "blue"))
dev.off()


glimpse(sample_annotation.full.corrected)
cbind(sample_annotation.full.corrected$library_id, sample_annotation.full.corrected$sex,
      logXY.tmp[sample_annotation.full.corrected$library_id])
# some additional mismatches (lib11731, lib11716, lib11635, lib11597, lib11410), but they're all really low-quality libraries; set them to NA
sample_annotation.full.corrected$sex_by_rna[
  sample_annotation.full.corrected$library_id %in%
    c("lib11731", "lib11716", "lib11635", "lib11597", "lib11410")] <- NA

table(sample_annotation.full.corrected.qc[,c("study","sex_by_rna")])
# some imbalance in sex by study, which could mess up the sex-linked genes if I correct by Study and not Sex

## plot XIST and UTY for libraries that pass QC
counts.test_sex.tmp <-
  data.frame(
    lib_id=sample_annotation.full.corrected.qc$library_id,
    sex=factor(sample_annotation.full.corrected.qc$sex, levels=c("F", "M")),
    UTY=log2(unlist(cpm(counts.merged)[
      match("ENSG00000183878", rownames(counts.merged)),
      match(sample_annotation.full.corrected.qc$library_id,
            colnames(counts.merged))])+1),
    XIST=log2(unlist(cpm(counts.merged)[
      match("ENSG00000229807", rownames(counts.merged)),
      match(sample_annotation.full.corrected.qc$library_id,
            colnames(counts.merged))])+1))

counts.test_sex.tmp <- counts.test_sex.tmp %>%
  group_by(sex) %>%
  mutate(UTY_outlier = ifelse(is_outlier(UTY), lib_id, NA)) %>%
  mutate(XIST_outlier = ifelse(is_outlier(XIST), lib_id, NA)) %>%
  as.data.frame()


## scatterplot of points on UTY and XIST to identify sex (with points colored by sex)
pdf(file="plots/sex_by_UTY_XIST.corrected.pdf", h=8,w=9)
ggplot(data=counts.test_sex.tmp, aes(x=XIST, y=UTY, colour=sex)) +
  geom_point(position=position_jitter(w=0.02,h=0.02), size=3) +
  geom_text(aes(label=UTY_outlier), colour="black", vjust=-1.3) +
  geom_text(aes(label=XIST_outlier), colour="black", vjust=-1.3) +
  scale_color_manual(values=c("red", "blue")) +
  labs(x="XIST expression\nlog2 (CPM + 1)", y="UTY expression\nlog2 (CPM + 1)")
dev.off()
# interesting; lib11577 has high XIST and UTY counts, suggesting it may be a mixed sample

## scatterplot of points on UTY and XIST to identify sex (with points colored by sex, and no point labels)
pdf(file="plots/sex_by_UTY_XIST.corrected.no_text.pdf", h=8,w=9)
ggplot(data=counts.test_sex.tmp, aes(x=XIST, y=UTY, colour=sex)) +
  geom_point(position=position_jitter(w=0.02,h=0.02), size=3) +
  scale_color_manual(values=c("red", "blue")) +
  labs(x="XIST expression\nlog2 (CPM + 1)", y="UTY expression\nlog2 (CPM + 1)")
dev.off()
# interesting; lib11577 has high XIST and UTY counts, suggesting it may be a mixed sample


## boxplots of X and Y chromosome genes to check sex
ggplot(data=counts.test_sex.tmp, aes(x = sex, y=UTY)) +
  # geom_boxplot(position="jitter") +
  geom_boxplot() +
  labs(x="Sex", y="log10 UTY expression")

ggplot(data=counts.test_sex.tmp, aes(x = sex, y=XIST)) +
  # geom_boxplot(position="jitter") +
  geom_boxplot() +
  labs(x="Sex", y="log10 XIST expression")


## drop problematic libraries, and further filter quality-controlled master and counts objects
# skipped this step 2016-07-01, because I want to be able to identify the libraries where the sex is wrong!
# 
# sample_annotation.full.corrected.qc <-
#   sample_annotation.full.corrected.qc[(sample_annotation.full.corrected.qc$sex_by_rna %in% c("F","M")) &
#                        ((sample_annotation.full.corrected.qc$sex_by_rna == sample_annotation.full.corrected.qc$sex) |
#                           (sample_annotation.full.corrected.qc$sex == "")),]
# sample_annotation.full.corrected.qc <- sample_annotation.full.corrected.qc[sample_annotation.full.corrected.qc$libID %nin% bad.libs,]
# 
# remove problematic libraries from counts file
# counts.merged.aggregated.qc <-
#   counts.merged.aggregated.qc[,match(sample_annotation.full.corrected.qc$library_id, colnames(counts.merged.aggregated.qc))]

rm_tmp(ask=FALSE)  # remove unneeded variables


##### repeat plotting PCA on (normalized) full data set to check for bias and batch effects #####

# PCA doesn't change, because counts don't change, so I can use the existing versions with the updated sample annotation

# attach sample info to PC scores  
scores_sample_annotation_pcaAll.corrected.full <-
  merge(sample_annotation.full.corrected.qc,
        as.data.frame(pcaAll.full$x), by.x = "library_id", by.y="row.names")


### Color plots of PCs by different variables

## plots of PC1s 1-3, colored by cell type

plot_PCAs(scores_sample_annotation_pcaAll.corrected.full,
          pvars.labs=pcaAll.full$pvars.labs, PCs = 1:3,
          color_by_var="cell_type", color_by_var_levels=c("CD4","CD8","Bcell","monocyte","PBMC"),
          my_cols=colorblind_pal()(5),
          file_prefix="plots/T1D_placebos_sorted_cells.all_libs.corrected", plotdims=c(9,6))

# plot PC1 vs. PC2, colored by cell type, with the problem libraries labeled
pdf("plots/T1D_placebos_sorted_cells.all_libs.corrected.PC1_vs_PC2.color_by_cell_type.text_by_library_id.pdf",
    width=9, height=6)
ggplot(data=scores_sample_annotation_pcaAll.corrected.full,
       mapping=aes(x=PC1, y=PC2, color=cell_type)) +
  geom_point(alpha=0.8, size=3) +
  geom_text(data=scores_sample_annotation_pcaAll.corrected.full[
    scores_sample_annotation_pcaAll.corrected.full$library_id %in%
      paste0("lib", c(11577:11583,11600:11606)),],
    mapping=aes(label=library_id),
    nudge_y=5, size=5,
    show.legend=FALSE) +
  geom_text(data=scores_sample_annotation_pcaAll.corrected.full[
    scores_sample_annotation_pcaAll.corrected.full$PC2 > 50 &
      scores_sample_annotation_pcaAll.corrected.full$PC2 < 75,],
    mapping=aes(label=library_id),
    nudge_y=5, size=5,
    show.legend=FALSE) +
  scale_color_manual(values=setNames(colorblind_pal()(5), c("CD4","CD8","Bcell","monocyte","PBMC")))
dev.off()

## plots of PCs 1-3, colored by contaminated samples
scores_sample_annotation_pcaAll.corrected.full$contaminated <-
  ifelse(scores_sample_annotation_pcaAll.corrected.full$library_id == "lib11577",
         "Y", "N") %>%
  factor(levels=c("N", "Y"))
plot_PCAs(scores_sample_annotation_pcaAll.corrected.full,
          pvars.labs=pcaAll.full$pvars.labs, PCs = 1:3,
          color_by_var="contaminated",
          pch_by_var="contaminated",
          my_pch=c(4,16),
          my_cols=c("black", "red"),
          file_prefix="plots/T1D_placebos_sorted_cells.all_libs.corrected",
          plotdims=c(9,6))


##### plot cell-type specific markers, with libraries colored by cell type #####
## among other things, this could help to determine the cell type of lib11577

master.tmp <- sample_annotation.full.corrected.qc
counts.tmp <- log2(cpm(counts.merged.aggregated.qc)+1)

genes.tmp <-
  c("CD4", "CD8A", "CD8B", "CD19", "CD79A", "BLK", "BLNK", "CSPG2", "CXCL3", "CTSL", "CASP1")
for (i in genes.tmp) {
  if (i %nin% rownames(counts.tmp)) {
    print(paste(i, "not found in counts"))
    next
  }
  master.tmp[[paste0(i, "_cpm")]] <-
    counts.tmp[
      i, match(master.tmp$library_id, colnames(counts.tmp))]
}


# ggplot(data=master.tmp,
#        mapping=aes(x=CD4_cpm, y=CD8B_cpm, color=cell_type)) +
#   geom_point(size=2, alpha=0.8) +
#   scale_color_colorblind() +
#   geom_text(data=master.tmp[master.tmp$library_id=="lib11577",],
#             color="black",
#             mapping=aes(label=library_id))

pdf("plots/cell_type_plot_CD4_vs_CD8B.T_cells_only.pdf", w=8, h=6)
ggplot(data=master.tmp[master.tmp$cell_type %in% c("CD4", "CD8"),],
       mapping=aes(x=CD4_cpm, y=CD8B_cpm, color=cell_type)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual("Cell Type", values=setNames(colorblind_pal()(5)[c(1,2)], c("CD4", "CD8"))) +
  geom_text(data=master.tmp[master.tmp$library_id=="lib11577",],
            color="black",
            mapping=aes(label=library_id)) +
  labs(x="CD4 (log2 CPM)", y=("CD8B (log2 CPM)"))
dev.off()

# ggplot(data=master.tmp,
#        mapping=aes(x=CD4_cpm, y=CXCL3_cpm, color=cell_type)) +
#   geom_point(size=3, alpha=0.7) +
#   scale_color_colorblind("Cell Type") +
#   geom_text(data=master.tmp[master.tmp$library_id=="lib11577",],
#             color="black",
#             mapping=aes(label=library_id))

pdf("plots/cell_type_plot_CTSL_vs_CD8B.pdf", w=8, h=6)
ggplot(data=master.tmp,
       mapping=aes(x=CTSL_cpm, y=CD8B_cpm, color=cell_type)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(
    "Cell Type",
    values=setNames(colorblind_pal()(5), c("CD4", "CD8", "Bcell", "monocyte", "PBMC"))) +
  geom_text(data=master.tmp[master.tmp$library_id=="lib11577",],
            color="black",
            mapping=aes(label=library_id)) +
  labs(x="CTSL (log2 CPM)", y=("CD8B (log2 CPM)"))
dev.off()

# ggplot(data=master.tmp,
#        mapping=aes(x=CD4_cpm, y=CXCL3_cpm, color=cell_type)) +
#   geom_point(size=3, alpha=0.7) +
#   scale_color_colorblind("Cell Type") +
#   geom_text(data=master.tmp[master.tmp$library_id=="lib11577",],
#             color="black",
#             mapping=aes(label=library_id))

pdf("plots/cell_type_plot_CD19_vs_CD8B.pdf", w=8, h=6)
ggplot(data=master.tmp,
       mapping=aes(x=CD19_cpm, y=CD8B_cpm, color=cell_type)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(
    "Cell Type",
    values=setNames(colorblind_pal()(5), c("CD4", "CD8", "Bcell", "monocyte", "PBMC"))) +
  geom_text(data=master.tmp[master.tmp$library_id=="lib11577",],
            color="black",
            mapping=aes(label=library_id)) +
  labs(x="CD19 (log2 CPM)", y="CD8B (log2 CPM)")
dev.off()

# ggplot(data=master.tmp,
#        mapping=aes(x=CD8B_cpm, y=BLNK_cpm, color=cell_type)) +
#   geom_point(size=3, alpha=0.7) +
#   scale_color_colorblind("Cell Type") +
#   geom_text(data=master.tmp[master.tmp$library_id=="lib11577",],
#             color="black",
#             mapping=aes(label=library_id)) +
#   labs(x="CD8B (log2 CPM)", y="BLNK (log2 CPM)")

### check the one weird-looking CD4 library
master.tmp %>%
  filter(cell_type=="CD4") %>%
  arrange(desc(CD19_cpm)) %>%
  dplyr::select(one_of(c("library_id", "cell_type", "CD19_cpm")))
# lib11465 has unusually high CD19_cpm

master.tmp %>%
  filter(cell_type=="CD4") %>%
  arrange(desc(CTSL_cpm)) %>%
  dplyr::select(one_of(c("library_id", "cell_type", "CTSL_cpm")))
# lib11655 has unusually high CTSL


##### Drop outliers, re-normalize, and re-run PCA #####
libs.bad.tmp <- c("lib11766", "lib11577")
# lib11766 is not an outlier in PCA using log2(cpm(counts)+1)
# but it is really different from all the other libraries

# remove outlier from sample annotation file
sample_annotation.full.corrected.qc.no_outliers <-
  sample_annotation.full.corrected.qc[
    sample_annotation.full.corrected.qc$library_id %nin% libs.bad.tmp,]
# glimpse(sample_annotation.full.corrected.qc.no_outliers)
# this leaves 378 libraries of the original 420

# and the uncorrected version
sample_annotation.full.uncorrected.qc.no_outliers <-
  sample_annotation.full.uncorrected.qc[
    sample_annotation.full.uncorrected.qc$library_id %nin% libs.bad.tmp,]
# glimpse(sample_annotation.full.uncorrected.qc.no_outliers)

# trim metrics object to include only desired libraries
metrics.merged.qc.no_outliers <-
  metrics.merged.qc[
    match(sample_annotation.full.corrected.qc.no_outliers$library_id,
          metrics.merged.qc$lib_id),]

# trim counts object to include only desired libraries
keepCounts.tmp <-
  colnames(counts.merged.aggregated.qc) %in%
  sample_annotation.full.corrected.qc.no_outliers$library_id
counts.merged.aggregated.qc.no_outliers <- counts.merged.aggregated.qc[,keepCounts.tmp]
# glimpse(counts.merged.aggregated.qc.no_outliers)

# order design by order of counts libraries, and drop samples from design if not found in counts object
sample_annotation.full.corrected.qc.no_outliers <-
  sample_annotation.full.corrected.qc.no_outliers[
    match(colnames(counts.merged.aggregated.qc.no_outliers),
          sample_annotation.full.corrected.qc.no_outliers$library_id),]
# table(sample_annotation.full.corrected.qc.no_outliers$cell_type)  # check numbers of cell types

# Remove libraries not in sample annotation object,
# filter out genes that have a cpm of at least one in < 15% of libraries,
# and normalize counts across libraries
counts.merged.aggregated.qc.no_outliers.normalized <-
  calcNormCounts(counts=counts.merged.aggregated.qc.no_outliers,
                 min_cpm=1,
                 design=sample_annotation.full.corrected.qc.no_outliers,
                 libID_col="library_id")
# glimpse(counts.merged.aggregated.qc.no_outliers.normalized)
# retain 12,755 of 19,751 genes

rm_tmp(ask=FALSE)


## run PCA on filtered, normalized counts

# run PCA
# only need pcaAll for the scores, the sds, and the plot
pcaAll.no_outliers <-
  calc_PCAs(counts.merged.aggregated.qc.no_outliers.normalized, log2_transform=TRUE)

# scree plot, the number of informative PCs = elbow
plot(pcaAll.no_outliers, type="l")
# 3 informative PCs

# attach sample info to PC scores  
scores_sample_annotation_pcaAll.corrected.qc.no_outliers <-
  merge(sample_annotation.full.corrected.qc.no_outliers,
        as.data.frame(pcaAll.no_outliers$x), by.x = "library_id", by.y="row.names")


### Color plots of PCs by different variables

## plots of PCs1 1-3, with no color
plot_PCAs(scores_sample_annotation_pcaAll.corrected.qc.no_outliers,
          pvars.labs=pcaAll.no_outliers$pvars.labs, PCs = 1:3,
          file_prefix="plots/T1D_placebos_sorted_cells.corrected.qc.no_outliers",
          plotdims=c(7.3,6))

## plots of PC1s 1-3, colored by cell_type
plot_PCAs(scores_sample_annotation_pcaAll.corrected.qc.no_outliers,
          pvars.labs=pcaAll.no_outliers$pvars.labs, PCs = 1:3,
          color_by_var="cell_type",
          color_by_var_levels=c("CD4","CD8","Bcell","monocyte","PBMC"),
          my_cols=colorblind_pal()(5),
          file_prefix="plots/T1D_placebos_sorted_cells.corrected.qc.no_outliers",
          plotdims=c(9,6))

## plots of PC1s 1-3, colored by participant_id
plot_PCAs(scores_sample_annotation_pcaAll.corrected.qc.no_outliers,
          pvars.labs=pcaAll.no_outliers$pvars.labs, PCs = 1:3,
          color_by_var="participant_id",
          add_legend=FALSE,
          my_cols=colorblind_pal()(8),
          file_prefix="plots/T1D_placebos_sorted_cells.corrected.qc.no_outliers",
          plotdims=c(9,6))

## plots of PC1s 1-3, colored by abc_timepoint
plot_PCAs(scores_sample_annotation_pcaAll.corrected.qc.no_outliers,
          pvars.labs=pcaAll.no_outliers$pvars.labs, PCs = 1:3,
          color_by_var="abc_timepoint",
          color_by_var_levels=c("A","B","C"),
          my_cols=colorblind_pal()(3),
          file_prefix="plots/T1D_placebos_sorted_cells.corrected.qc.no_outliers",
          plotdims=c(9,6))

## plots of PC1s 1-3, colored by flowcell
plot_PCAs(scores_sample_annotation_pcaAll.corrected.qc.no_outliers,
          pvars.labs=pcaAll.no_outliers$pvars.labs, PCs = 1:3,
          pch_by_var="flowcell",
          file_prefix="plots/T1D_placebos_sorted_cells.corrected.qc.no_outliers",
          plotdims=c(9,6))

## plots of PCs 1-3, colored by flowcell/lane
scores_sample_annotation_pcaAll.corrected.qc.no_outliers$flowcell_lane <-
  with(scores_sample_annotation_pcaAll.corrected.qc.no_outliers,
       paste(flowcell, lane, sep="_"))
plot_PCAs(scores_sample_annotation_pcaAll.corrected.qc.no_outliers,
          pvars.labs=pcaAll.no_outliers$pvars.labs, PCs = 1:3,
          color_by_var="flowcell_lane",
          add_legend=FALSE,
          my_cols=colorblind_pal()(8),
          file_prefix="plots/T1D_placebos_sorted_cells.corrected.qc.no_outliers",
          plotdims=c(9,6))

## plots of PCs 1-3, colored by swapped / unswapped samples (to see if anything else is weird with them)
scores_sample_annotation_pcaAll.corrected.qc.no_outliers$swapped <-
  ifelse(scores_sample_annotation_pcaAll.corrected.qc.no_outliers$trna_sampleid %in%
           paste0("tRNA", 26465:26488), "Y", "N") %>%
  factor(levels=c("N", "Y"))
plot_PCAs(scores_sample_annotation_pcaAll.corrected.qc.no_outliers,
          pvars.labs=pcaAll.no_outliers$pvars.labs, PCs = 1:3,
          color_by_var="swapped",
          my_cols=colorblind_pal()(8),
          file_prefix="plots/T1D_placebos_sorted_cells.corrected.qc.no_outliers",
          plotdims=c(9,6))


##### Calculate correlations between PCs and metrics/clinical variables #####

# merge metrics, sample annotation, and PCA scores
scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers <-
  merge(metrics.merged.qc.no_outliers,
        scores_sample_annotation_pcaAll.corrected.qc.no_outliers[
          , 1:max(which(str_detect(colnames(scores_sample_annotation_pcaAll.corrected.qc.no_outliers), "PC[0-9]")))],
        by.x="lib_id", by.y="library_id")

# generate empty matrix to store correlation data; drop "sample" columns, which mostly relate to replication of sequencing
cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers <-
  matrix(
    NA,
    nrow=sum(str_detect(colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers), "PC[0-9]") &
               !str_detect(colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers), "sample")),
    ncol=sum(!str_detect(colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers), "PC[0-9]") &
               !str_detect(colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers), "sample"))-1,
    dimnames=list(
      colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers)[
        str_detect(colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers), "PC[0-9]") &
          !str_detect(colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers), "sample")],
      colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers)[
        !str_detect(colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers), "PC[0-9]") &
          !str_detect(colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers), "sample")][-1]))

min_libs.tmp <- 10 # set minimum number of libraries with data for calculations to be done

## Loading ICC library SOS 12/23/16
# library(ICC)
# calculate correlations of all variables with PCs, using pearson for numeric variables, and ICC for non-numeric
for (i in colnames(cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers)[-1]) {
  if (sum(!is.na(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers[[i]])) <= min_libs.tmp) next
  if (is.numeric(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers[[i]])) {
    cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers[,i] <-
      cor(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers[[i]],
          scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers[
            , match(rownames(cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers),
                    colnames(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers))],
          method="pearson", use="pairwise")
  } else if (any(duplicated(na.omit(scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers[[i]])))) { # skip columns with all unique values
    for (j in rownames(cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers)) {
      cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers[j,i] <-
        ICC::ICCbare(
          data=scores_metrics_sample_annotation_pcaAll.full.corrected.qc.no_outliers,
          x=i, y=j)
    }
  }
}

# drop all-NA columns from correlation matrix
cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers <-
  cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers[
    ,(apply(cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers,
            MARGIN=2, function(x) sum(!is.na(x))) > 2)]
sort(cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers[1,])

pdf("plots/heatmap_correlations_clinvars_metrics_pcaAll.full.corrected.qc.no_outliers.pdf",
    w=15, h=6)
heatmap.2(x=cor_all_clinvars_metrics_vs_pcaAll.full.corrected.qc.no_outliers[1:10,],
          Rowv=FALSE, Colv=FALSE, dendrogram="none",
          col=colorRampPalette(c("blue", "white", "red"))(100),
          trace="none", margins=c(8,15),
          keysize=1.1, density.info="none")
dev.off()
# looks like cell_type is strongly correlated with PCs 1,2,3,5
# looks like a bunch of quality metrics are correlated with PC1 and PC3
# clinical variables not at all correlated with those low-number PCs

rm_tmp(ask=FALSE)


##### Calculate correlations or distances among samples, and compare similarity by participant_id (corrected and uncorrected) #####

# vwts.merged.aggregated.qc.no_outliers <-
#   voomWithQualityWeights(counts.merged.aggregated.qc.no_outliers.normalized, plot=TRUE, span=0.1)

cor.spearman.merged.aggregated.qc.no_outliers <-
  cor(counts.merged.aggregated.qc.no_outliers.normalized, method="spearman")

cor.pearson.merged.aggregated.qc.no_outliers <-
  cor(log2(counts.merged.aggregated.qc.no_outliers.normalized+1), method="pearson")

## plot heatmap of all correlations

# generate a vector for the order of objects in the heatmap
sample_annotation.full.corrected.qc.no_outliers$heatmap_order <-
  with(sample_annotation.full.corrected.qc.no_outliers,
       order(cell_type, participant_id, weeks))

sample_annotation.full.uncorrected.qc.no_outliers$heatmap_order <-
  with(sample_annotation.full.uncorrected.qc.no_outliers,
       order(cell_type, participant_id, weeks))

cor.spearman.merged.aggregated.qc.no_outliers.corrected_order <-
  cor.spearman.merged.aggregated.qc.no_outliers[
    sample_annotation.full.corrected.qc.no_outliers$heatmap_order,][
      ,sample_annotation.full.corrected.qc.no_outliers$heatmap_order]

cor.spearman.merged.aggregated.qc.no_outliers.uncorrected_order <-
  cor.spearman.merged.aggregated.qc.no_outliers[
    sample_annotation.full.uncorrected.qc.no_outliers$heatmap_order,][
      ,sample_annotation.full.uncorrected.qc.no_outliers$heatmap_order]

cor.pearson.merged.aggregated.qc.no_outliers.corrected_order <-
  cor.pearson.merged.aggregated.qc.no_outliers[
    sample_annotation.full.corrected.qc.no_outliers$heatmap_order,][
      ,sample_annotation.full.corrected.qc.no_outliers$heatmap_order]

cor.pearson.merged.aggregated.qc.no_outliers.uncorrected_order <-
  cor.pearson.merged.aggregated.qc.no_outliers[
    sample_annotation.full.uncorrected.qc.no_outliers$heatmap_order,][
      ,sample_annotation.full.uncorrected.qc.no_outliers$heatmap_order]


# spearman correlations with corrected annotation
pdf(file="plots/heatmap.cor_counts.spearman.merged.aggregated.qc.no_outliers.corrected_order.pdf",
    w=12,h=9)
heatmap.2(cor.spearman.merged.aggregated.qc.no_outliers.corrected_order,
          scale="none", Rowv=FALSE, Colv=FALSE, symm=TRUE, trace="none",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
          ColSideColors=
            setNames(colorblind_pal()(5),c("CD4","CD8","Bcell","monocyte","PBMC"))[
              sample_annotation.full.corrected.qc.no_outliers$cell_type[
                sample_annotation.full.corrected.qc.no_outliers$heatmap_order]],
          RowSideColors=
            setNames(colorblind_pal()(5),c("CD4","CD8","Bcell","monocyte","PBMC"))[
              sample_annotation.full.corrected.qc.no_outliers$cell_type[
                sample_annotation.full.corrected.qc.no_outliers$heatmap_order]])
dev.off()

# spearman correlations with uncorrected annotation
pdf(file="plots/heatmap.cor_counts.spearman.merged.aggregated.qc.no_outliers.uncorrected_order.pdf",
    w=12,h=9)
heatmap.2(cor.spearman.merged.aggregated.qc.no_outliers.uncorrected_order,
          scale="none", Rowv=FALSE, Colv=FALSE, symm=TRUE, trace="none",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
          ColSideColors=
            setNames(colorblind_pal()(5),c("CD4","CD8","Bcell","monocyte","PBMC"))[
              sample_annotation.full.uncorrected.qc.no_outliers$cell_type[
                sample_annotation.full.uncorrected.qc.no_outliers$heatmap_order]],
          RowSideColors=
            setNames(colorblind_pal()(5),c("CD4","CD8","Bcell","monocyte","PBMC"))[
              sample_annotation.full.uncorrected.qc.no_outliers$cell_type[
                sample_annotation.full.uncorrected.qc.no_outliers$heatmap_order]])
dev.off()

# pearson correlations with corrected annotation
pdf(file="plots/heatmap.cor_counts.pearson.merged.aggregated.qc.no_outliers.corrected_order.pdf",
    w=12,h=9)
heatmap.2(cor.pearson.merged.aggregated.qc.no_outliers.corrected_order,
          scale="none", Rowv=FALSE, Colv=FALSE, symm=TRUE, trace="none",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
          ColSideColors=
            setNames(colorblind_pal()(5),c("CD4","CD8","Bcell","monocyte","PBMC"))[
              sample_annotation.full.corrected.qc.no_outliers$cell_type[
                sample_annotation.full.corrected.qc.no_outliers$heatmap_order]],
          RowSideColors=
            setNames(colorblind_pal()(5),c("CD4","CD8","Bcell","monocyte","PBMC"))[
              sample_annotation.full.corrected.qc.no_outliers$cell_type[
                sample_annotation.full.corrected.qc.no_outliers$heatmap_order]])
dev.off()

# pearson correlations with uncorrected annotation
pdf(file="plots/heatmap.cor_counts.pearson.merged.aggregated.qc.no_outliers.uncorrected_order.pdf",
    w=12,h=9)
heatmap.2(cor.pearson.merged.aggregated.qc.no_outliers.uncorrected_order,
          scale="none", Rowv=FALSE, Colv=FALSE, symm=TRUE, trace="none",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
          ColSideColors=
            setNames(colorblind_pal()(5),c("CD4","CD8","Bcell","monocyte","PBMC"))[
              sample_annotation.full.uncorrected.qc.no_outliers$cell_type[
                sample_annotation.full.uncorrected.qc.no_outliers$heatmap_order]],
          RowSideColors=
            setNames(colorblind_pal()(5),c("CD4","CD8","Bcell","monocyte","PBMC"))[
              sample_annotation.full.uncorrected.qc.no_outliers$cell_type[
                sample_annotation.full.uncorrected.qc.no_outliers$heatmap_order]])
dev.off()


##### check on calculation of progressor status (done in a separate set of scripts and imported) #####

# count of progressors and non-progressors using threshold method
colSums(table(unique(sample_annotation.full.corrected.qc.no_outliers[,c("participant_id","progressor_min_percent_pre104weeks")])))
# 17 progressors, 11 non-progressors

# count of progressors and non-progressors using split at ~100 weeks
colSums(table(unique(sample_annotation.full.corrected.qc.no_outliers[,c("participant_id","progressor_split_at100weeks")])))
# 10 progressors, 15 non-progressors

table(sample_annotation.full.corrected.qc.no_outliers[is.na(sample_annotation.full.corrected.qc.no_outliers$progressor_split_at100weeks),"weeks"]) # sanity check that all NAs are not at ~100 weeks

rm_tmp(ask=FALSE)


##### save Rdata objects for downstream use #####

## these objects have had bad libraries and outliers removed
## they include all samples
## the data are not normalized, and include all treatment groups

# generate simply-named combined objects, with all libraries passing qc and outlier cuts
counts.final <- counts.merged.aggregated.qc.no_outliers
metrics.final <- metrics.merged.qc.no_outliers
sample_annotation.final <- sample_annotation.full.corrected.qc.no_outliers

save(file="output/T1D_placebos_sorted_cells_data_for_analysis_2016-07-11.Rdata",
     list=c("counts.final", "metrics.final", "sample_annotation.final", "pcaAll.no_outliers"))
