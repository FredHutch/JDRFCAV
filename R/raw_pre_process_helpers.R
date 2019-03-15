#' raw data pre-processing functions
#' @param all_kind all kind...
#'
#' @name raw_pre_process_helpers
NULL

#' @return all kind...
#' @details TBA







#' @export
#' @rdname raw_pre_process_helpers
# store All list into MAE -------------------------------------------------

store_All_list_into_MAE<-function(All.list.raw.only, pData.in){
  # All.list.raw.only=Assay.raw.list[2:17]
  # pData.in=Assay.raw.list[1][[1]]  #similar to Cohort2_Meta_Design, but with all Yis (will still require later on to select the Yi of interest, and calculate the y.cont and y.categ)

  ## 1. raw exprs
  # All.list.raw.only %>% map(dim)
  # names(All.list.raw.only)
  # All.list.raw.only %>% map(function(x) rownames(x))
  All.list.raw.only.t<-All.list.raw.only %>% map(function(x){
    # x=All.list.raw.only[[1]]
    x %>%
      data.frame %>%
      rownames_to_column %>%
      gather(var, value, -rowname) %>% ## t coerced the numeric values into character. use gather()+spread() instead.
      spread(rowname, value) %>%
      remove_rownames %>%
      column_to_rownames('var')
  }) %>%
    map(as.matrix) %>%
    map(SummarizedExperiment)

  # assays(All.list.raw.only.t[[1]])[[1]]
  # All.list.raw.only.t %>% map(dim)

  ## 2. pData
  # pData.in<-distinct(Reduce( list(Merged_Any_AssayeSet2$CohortMembership.Uniqe.Length,...) ), function(x,y) full_join(x,y, by='??'))

  ## 3. sampleMap
  # listmap1<-data.frame(primary=rownames(All.list.raw.only[[1]]), assay=names(All.list.raw.only[[1]])
  # listmap.i<-map2(.x=All.list.raw.only, .y=All.list.raw.design.only, function(x,y) data.frame(primary=rownames(x), assay=rownames(y) ) )
  listmap.i<-All.list.raw.only %>% map(function(x) data.frame(primary=rownames(x), assay=rownames(x) ) )
  # names(listmap.i) is automatically inherited from the All.list.raw.only (list)
  mySampleMap<-listToMap(listmap.i)

  ## 4. MAE.all
  # PrepMultiAssay(All.list.raw.only.t, pData.in, mySampleMap)
  # All.list.raw.only %>% map(function(x) rownames(x) %in% rownames(pData.in)) # Pro.Insulin ABATE_003018
  # rownames(pData.in)[which( !rownames(pData.in) %in% (All.list.raw.only %>% map(function(x) rownames(x) ) %>% unlist %>% unique) )]
  MultiAssay.i<-MultiAssayExperiment(experiments=All.list.raw.only.t, pData.in, mySampleMap)
  MultiAssay.i
}



#' @export
#' @rdname raw_pre_process_helpers
# setup + MAE ---------------------------------------------------------------

### metadata.raw.table need to be added as parameter, and be checked internally!!!
Setup.plue.MAE.to.vertical.function<-function(MAE, param.Y.name, param.assays.vector, param.subjects.study, param.cohort, param.covariates, param.rv144.mRNA.collapse.multiple.probes, parame.gene.or.module, param.Selected.Meta.sets, Assay.Analyte.sep){
  # MAE=RNAseq.Wholeblood.MAE; RNAseq.Wholeblood.MAE=1, param.cohort
  # MAE=MAE.CAV; param.assays.vector=1
  select<-dplyr::select

  ### 2.1 test MAE level:
  if (! param.Y.name %in% names(colData(MAE))) cat('param.Y.name is not in colData(MAE)')

  ### 2.2 pre-screening for issues with specific features. use only in rare cases. ML.function should know how to address such cases, and report a non-valid results
  # MAE<-pre.screening.for.issues.Function(MAE)


  ## 2.3 update assay type (store in metadata slot)
  experiments(MAE)@listData<-pmap(list(experiments(MAE)@listData, metadata.raw.table$Assay %>% as.character, metadata.raw.table$Assay.type %>% as.character), function(x, y, z) Insert.Metadata.Assay.type.FUNCTION(x, y, z))


  ### 2.4 Filter out non-selected assays
  MAE.selected.assays<-MAE[,,param.assays.vector]



  ### 2.5 filter out subjects:

  ### 2.5.1 NAs in Y / outcome
  # Filter out NA in Y (and their associated slots):
  # if Y has missing values, they are removed ASAP at the MAE, which will automatically remove it from the raw/assays data
  y.count.missing.values<-sum(is.na(colData(MAE.selected.assays)[param.Y.name]))
  remove.subject.location<-which( is.na(colData(MAE.selected.assays)[param.Y.name])==TRUE )
  if( y.count.missing.values!=0 ) cat(paste('outcome has', y.count.missing.values, 'missisng values. these subjects will be removed'))
  if(length(remove.subject.location)>0) MAE.selected.assays<-MAE.selected.assays[,-remove.subject.location, ]
  ### check:
  if( is.factor(colData(MAE.selected.assays)[,param.Y.name]) ) {
    cat(paste('Y should be numeric for both binary and continuous outcomes','\n'))
    colData(MAE.selected.assays)[,param.Y.name]<-as.numeric(as.character(colData(MAE.selected.assays)[,param.Y.name]))
  }
  if( is.numeric(colData(MAE.selected.assays)[,param.Y.name]) ) cat(paste('Y was sucessfully converted into numeric format','\n'))
  cat('Y=', str(colData(MAE.selected.assays)[,param.Y.name]))



  ### 2.5.2 filter by Study (Cohort 2)
  if( ! is.null(param.subjects.study)){
    if( param.cohort=='cohort2' ){
      subject.ids.keep<-colData(MAE.selected.assays) %>% data.frame %>% filter(study %in% param.subjects.study) %>% pull(participant_id)
      MAE.selected.assays<-MAE.selected.assays[ , subject.ids.keep, , drop=FALSE]
    }
  }


  # assay(MAE.selected.assays[,,5][[1]]) %>% dim [1:10,c('START_013003', 'START_013005')]

  # ff<-assay(experiments(MAE.selected.assays)[[1]]) %>% rownames %>% table
  # ff[ff>=2]

  ### 2.6 Covariates
  # 2.6.1 check if covariates in colData:

  if( length(param.covariates)>=1 ){
    if(!all(param.covariates %in% colnames(colData(MAE.selected.assays)) ) ) cat('not all covariates names are valid in MAE')
    Covariates.raw.MAE<-colData(MAE.selected.assays)[,param.covariates, drop=FALSE] %>% data.frame %>% as.tibble ## this is at the MAE level. nor relevant to Join_SumEset

    ## 2.6.2 filter out subjects with NAs
    Covariates.MAE.count.missing.values<-sum(is.na(Covariates.raw.MAE))
    which.covariate.has.missing<-Covariates.raw.MAE %>% map(function(x) which( is.na(x)==TRUE ) )
    remove.subject.location.all.covariates<-which.covariate.has.missing[which(which.covariate.has.missing %>% map(function(x) length(x)!=0) %>% unlist)] %>% unlist

    if( !is.null(remove.subject.location.all.covariates) ) {
      keep.subject.location.all.covariates.ID<-rownames(colData(MAE.selected.assays))[-remove.subject.location.all.covariates]
      if( length(remove.subject.location.all.covariates)>=1 ) cat('Subjects removed due to missing values in covariates:', rownames(colData(MAE.selected.assays))[remove.subject.location.all.covariates])
      if(length(remove.subject.location.all.covariates)>=1) MAE.selected.assays<-MAE.selected.assays[,keep.subject.location.all.covariates.ID, ]
    }
    # colData(MAE.selected.assays[,,]) %>% nrow
  }
  print(MAE.selected.assays)



  ### 2.7 Modules / Pathways:

  ## rv144 mRNA only: collapse multiple probes of same gene
  if(param.rv144.mRNA.collapse.multiple.probes){
    MAE.mRNA.collapsed<-Collapse.Multiple.Probes.to.Symbol.Function(MAE.selected.assays[,,'mRNA'][[1]], Prefix='mRNA')
    experiments(MAE.selected.assays)@listData$mRNA<-MAE.mRNA.collapsed
  }


  ## filter out non-mapped assays
  mod_all.sets.tib.selected_merged_meta_sets<-NA # default parameter

  if(parame.gene.or.module=='module'){
    ### 3.1. load GMTs, and keep only relevant ones, and merge all (relevant) sets
    load(file='mod_all.sets.tib.Rdata') # mod_all.sets.tib
    # mod_all.sets.tib$GMT.list %>% map(function(x) x[1] %>% map(function(x) x[1:4]))
    mod_all.sets.tib.selected<-mod_all.sets.tib %>% filter( Gmt.name %in% param.Selected.Meta.sets) # Filter OUT meta-sets

    # Before unlist, assure there is no overlap pathway names across all sets, if there is, make the name uniqe, e.g. for each set: unite('set_name', pathawy_name)
    #if( all(unlist(BMT$GMT.list ,recursive = FALSE) %>% names %>% table !=1)) '!!! some overlap pathway names across all sets'
    mod_all.sets.tib.selected_merged_meta_sets<-unlist(mod_all.sets.tib.selected$GMT.list, recursive = FALSE) # assume list names at each of 4 sublists is unique, and no overlap across 4 lists
    # mod_all.sets.tib.selected_merged_meta_sets[1:4]

    some.Analyte.names.at.assay<-experiments(MAE.selected.assays)@listData %>% map(function(x) assay(x) %>% rownames %>% head) ### enough to test only with few analytes
    modules.all.SYMBOL<-mod_all.sets.tib.selected_merged_meta_sets %>% unlist %>% unique
    crit.SYMBOL.in.assay<-some.Analyte.names.at.assay %>% map(function(x) any(x %in% modules.all.SYMBOL)) %>% unlist
    if (any(crit.SYMBOL.in.assay)) MAE.selected.assays<-MAE.selected.assays[,,crit.SYMBOL.in.assay]
  }



  # assay(MAE.selected.assays[,,1][[1]]) %>% is.na %>% any
  # assay(Join_SumEset.C2.complete.Y[,1:10]) %>% is.na %>% any
  # DF.long[1:10,1:30]




  # 3. MAE -> summarizedExperiment -> c(y, DF long) ----------------------------------------------------------

  ## 3.1 SummarizedExperiment
  Join_SumEset.C2.complete.Y<-MAE.to.Join_SumEset.FUNCTION(MAE.selected.assays, Assay.Analyte.sep) ### missing values in Y were ALRAEDY removed in previous MAE level
  assay.vec.orig<-rowData(Join_SumEset.C2.complete.Y) %>% data.frame %>% pull(Assay) # to be used later at mlr screening customed made feature selection with clustering

  Y.raw<-colData(Join_SumEset.C2.complete.Y) %>% data.frame %>% select(param.Y.name)
  DF.long<-assay(Join_SumEset.C2.complete.Y) %>% t
  if(!all(colData(Join_SumEset.C2.complete.Y)$participant_id==rownames(DF.long) )) cat('subjects order is not the same')

  Covariates.raw<-colData(Join_SumEset.C2.complete.Y)[,param.covariates, drop=FALSE] %>% data.frame  # should inherit the entire colData from MAE
  rownames(Covariates.raw)<-rownames(colData(Join_SumEset.C2.complete.Y))

  ### 3.2 scaling
  ## not part of screening, and should NOT be fused into CV
  DF.long.scaled<-Screen.DF.Standarize.by.assay.FUNCTION(DF.long, Assay.Analyte.sep)
  # DF.long.scaled %>% dim
  # str(DF.long)
  My.family<-ifelse( length(table(Y.raw))==2 & all(names(table(Y.raw))[1:2]==c(0,1)) , 'binomial', 'gaussian')
  cat(My.family)

  # boxplot(t(DF.long[, ])); DF.long %>% t %>% data.frame %>% map(sd)
  # boxplot(t(DF.long.scaled[, ]));  DF.long.scaled %>% t %>% data.frame %>% map(sd)
  # plot(DF.long %>% t %>% data.frame %>% map(sd), DF.long.scaled %>% t %>% data.frame %>% map(sd), xlab='orig', ylab='scaled')
  # plot(DF.long %>% t %>% data.frame %>% map(mean), DF.long.scaled %>% t %>% data.frame %>% map(mean), xlab='orig', ylab='scaled')


  ## 3.2.1 scaling continous covariates only
  if( ncol(Covariates.raw)==0 ) Covariates.scaled<-NULL
  if( ncol(Covariates.raw)>=1 ) {
    Covariates.scaled<-data.frame(Covariates.raw[Covariates.raw %>% map(is.factor) %>% unlist],
                                  Covariates.raw[Covariates.raw %>% map(is.numeric) %>% unlist] %>% scale)
    # apply(Covariates.scaled, 2, mean)
    # apply(Covariates.scaled, 2, sd)
    Covariates.scaled<-Covariates.scaled %>% map(function(x) as.numeric(x) ) %>% bind_cols %>% data.frame # can be also done before scaling. should not affect results!
    rownames(Covariates.scaled)<-rownames(Covariates.raw)
  }

  out<-list(
    MAE.selected.assays                        = MAE.selected.assays, ## !! Not scaled, but after subjects filtering
    Y.raw                                      = Y.raw,
    DF.long.scaled                             = DF.long.scaled,
    Covariates.scaled                          = Covariates.scaled,
    My.family                                  = My.family,
    mod_all.sets.tib.selected_merged_meta_sets = mod_all.sets.tib.selected_merged_meta_sets,
    Join_SumEset.C2.complete.Y                 = Join_SumEset.C2.complete.Y
  )
  return(out)

} # end function










#' @export
#' @rdname raw_pre_process_helpers

# 0.Set MultiAssay.cohort.2 with assay type !!! (MAE)------------------------------------------------------------------------
Insert.Metadata.Assay.type.FUNCTION<-function(SumEset.i, Assay, Assay.type){
  # SumEset.i<-experiments(MAE)[[1]]
  # Assay=metadata.raw.table$Assay[[1]]
  # Assay.type=metadata.raw.table$Assay.type[[1]]
  metadata(SumEset.i)$Name=Assay
  metadata(SumEset.i)$Assay.type=Assay.type
  SumEset.i
}




#' @export
#' @rdname raw_pre_process_helpers

# MAE to Join_SumEset ----------------------------------------------------------------------
MAE.to.Join_SumEset.FUNCTION<-function(MAE, Assay.Analyte.sep){
  # MAE=MAE.selected.assays
  ## MAE=MAE.univ.filtered

  ## 3.1 assays:colData (same as MAE)
  MAE.joint.ColData<-colData(MAE) ### Before univariate filtering
  # dim(MAE.joint.ColData)

  ## 3.2 assays:
  MAE.joint.assay<-assays(MAE)@listData %>% map(data.frame) %>% map(rownames_to_column) %>% bind_rows(.id='Assay') ### This is replacing the reduce(join()) accross all assay, and create the NAs!!!
  # head(MAE.joint.assay)
  MAE.joint.assay.ord.subjects<-MAE.joint.assay[,rownames(MAE.joint.ColData)] %>% as.matrix # reorder subjects same order as in ColData. should not affect row order
  rownames(MAE.joint.assay.ord.subjects)<-paste(MAE.joint.assay$Assay, MAE.joint.assay$rowname, sep=Assay.Analyte.sep) ## rows/analyte order has NOT been changed

  # dim(MAE.joint.union.assay.ord.subjects)
  # MAE.joint.union.assay %>% data.frame %>% map_chr(function(x) all(!is.na(x))) %>% table
  # MAE.joint.union.assay2<-longFormat(MAE.univ.filtered) %>% as.tibble %>% spread(primary, value) # doesn't work!!, but not crucial. done above with bind_rows


  ## 3.3 rowData: join (bind_rows)
  MAE.joint.rowData<-experiments(MAE)@listData %>%
    map(function(x) rowData(x)) %>% map(data.frame) %>%
    bind_rows(.id = "Assay") %>% as.tibble %>%
    mutate(Analyte=rownames(MAE.joint.assay.ord.subjects)) # rownames are the same, not oredered

  # dim(MAE.joint.intersect.rowData)
  # MAE.joint.rowData$Assay %>% table
  # if(!all(MAE.joint.rowData$Feature.name==MAE.joint.assay$rowname)) cat('rows order of assays and rowData is not the same')
  #if(!all(MAE.joint.rowData$Analyte==rownames(MAE.joint.union.assay.ord.subjects)) ) cat('rows order of assays and rowData is not the same')

  #rownames(MAE.joint.assay.ord.subjects)<-MAE.joint.rowData %>% unite(col='Assay.Analyte', Assay, Analyte, sep=Assay.Analyte.sep) %>% pull(Assay.Analyte)



  ## 3.4 SummarizedExperiment
  if( !identical( rownames(MAE.joint.ColData), colnames(MAE.joint.assay.ord.subjects)) ) cat('SummarizedExperiment will not work')
  Join_SumEset<-SummarizedExperiment(assays =MAE.joint.assay.ord.subjects,
                                     rowData=MAE.joint.rowData,
                                     colData=MAE.joint.ColData)
  # cat(paste('MAE to SummarizedExperiment completed', "\n"))
  # cat(paste(' ',ncol(Join_SumEset), 'subjects;  ', nrow(Join_SumEset), ' analytes (total across all assays)'))
  print(Join_SumEset)
  return(Join_SumEset)
}






#' @export
#' @rdname raw_pre_process_helpers
# Screen.scale  ----------------------------------------------------------------------

Screen.DF.Standarize.by.assay.FUNCTION<-function(DF, Assay.Analyte.sep){ # no additional parameters
  # DF=DF.long; DF[1:3, 1:5]

  # DF=DF.long2
  DF<-DF %>% t

  DF.Assay.id<-rownames(DF) %>% data.frame %>% separate(1, into=c('Assay','Analyte'), sep=Assay.Analyte.sep)
  DF.tibble.by.assay<-data.frame(Assay=DF.Assay.id$Assay, DF) %>% as.tibble %>% nest(-Assay, .key = DF.by.assay)# %>% unnest(DF.by.assay)

  DF.tibble.by.assay<-DF.tibble.by.assay %>%
    mutate(DF.scaled=DF.by.assay %>% map(function(x) x.t.scaled<-scale(t(x), center = TRUE, scale = TRUE) %>% t %>% as.tibble) ) # scale features, require t, and then undo(t)
  #boxplot(t(DF.tibble.by.assay$DF.by.assay[[1]]))
  #boxplot(t(DF.tibble.by.assay$DF.scaled[[1]]))
  DF.long.scaled<-DF.tibble.by.assay %>% select(-DF.by.assay) %>% unnest(DF.scaled) %>% select(-Assay) %>% data.frame
  rownames(DF.long.scaled)<-rownames(DF)
  # DF.long.scaled[1:3, 1:5]

  cat( paste('Scale completed. ', "\n", nrow(DF.long.scaled), 'Analytes;  ',ncol(DF.long.scaled), 'subjects'))
  return(DF.long.scaled %>% t)
}
