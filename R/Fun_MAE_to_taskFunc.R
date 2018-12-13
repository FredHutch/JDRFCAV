#' Convert MAE to mlr task
#'
#' @param MAE_obj MAE class
#' @param param.Y.name Vector of dependent variable name
#' @param param.covariates Vector of coaraiate variable(s) name
#' @param param_positive_y_level if ClassifTask, value (character or numeric) to be considered as the positive factor outcome
#' @return mlr's \code{ClassifTask} or \code{RegrTask}
#' @details In case of individual MAE assay (omic) with multiple sub-assays, only first sub-assay will be used.
#'
#'Either ClassifTask or RegrTask will be returned, based on the type of the param.Y.name variable
#' @examples
#' data(miniACC, package = 'MultiAssayExperiment') # ExpressionSet
#' miniACC
#' Fun_MAE_to_taskFunc(miniACC, param.Y.name = 'vital_status', param.covariates = c('gender','days_to_death'), param_positive_y_level = '1')
#'
#' @export



Fun_MAE_to_taskFunc<-function(MAE_obj, param.Y.name, param.covariates, param_positive_y_level){
  # i=4; MAE = MAE_DF$MAE_selected[[i]]; param.Y.name = MAE_DF$param.Y.name[[i]]; param.covariates = MAE_DF$param.covariates[[i]]
  # MAE_obj<-miniACC; param.Y.name = 'vital_status'; param.covariates = c('gender','days_to_death'); param_positive_y_level = '1'
  # sampleMap(MAE) %>% data.frame %>% pull(primary) %>% table

  DF_ColData<-MAE_obj %>% colData %>% data.frame# %>% dplyr::slice(1:10)# %>% select(-!!param.Y.name)
  DF_functionals<-DF_ColData[,c(param.covariates, param.Y.name), drop = FALSE]
  # DF_functionals$target<-DF_ColData$cpep_model_decayrate

  ## all non-numeric coariates must be factors!
  DF_functionals %<>% mutate_if(is.character, as.factor)
  # DF_functionals %>% str

  DF_exprsS<-assays(MAE_obj) %>% as.list %>% map(t) ## assume single assay within each SE
  # str(DF_exprsS)
  # DF_exprsS %>% map(dim)
  # DF_exprsS %>% map(rownames) %>% gplots::venn()

  ## !! makeFunctionalData() require only complete cases.
  ## !! instead of removing non-overlap subjects, fill in dummy NAs: (rowbind (fill) dropouts with NAs)
  complete_subject_id<-DF_ColData %>% rownames

  DF_exprsS_completeNA<-DF_exprsS %>% map(function(x){
    # x = DF_exprsS[[5]]
    DF_exprsS_i<-x %>% data.frame %>% rownames_to_column('colname')

    ## replace colname (assay-unique id) with primary (shared) id. ## should tolerate biological replicates within assays
    DF_exprsS_i_Joint<-DF_exprsS_i %>% left_join(sampleMap(MAE_obj) %>% data.frame %>% dplyr::select(-assay), by='colname') %>%
      dplyr::select(-colname) %>%
      column_to_rownames('primary')

    DF_dropouts_i<-data.frame(rowname = setdiff(complete_subject_id, DF_exprsS_i_Joint %>% rownames))
    DF_dropouts_i_NA<-data.frame(DF_dropouts_i, data.frame(matrix(ncol = ncol(DF_exprsS_i_Joint), nrow = nrow(DF_dropouts_i), NA)) %>% set_colnames(names(DF_exprsS_i_Joint)))
    #fF<-bind_rows(DF_exprsS_i_Joint %>% rownames_to_column, DF_dropouts_i_NA) %>% column_to_rownames('rowname') %>% as.matrix
    rbind(DF_exprsS_i_Joint %>% rownames_to_column, DF_dropouts_i_NA) %>% column_to_rownames('rowname') %>% as.matrix
  })
  # DF_exprsS_completeNA %>% str
  # DF_exprsS_completeNA %>% map(dim)

  ## Add assays to DF_functionals (currently it has only the covariates)
  for (i in 1: length(MAE_obj) ){ # for each assay
    # i=1
    DF_functionals[,names(DF_exprsS_completeNA)[[i]]]<-DF_exprsS_completeNA[[i]]
  }
  # str(DF_functionals)

  ## Classif:
  if( ( DF_functionals[,param.Y.name] %>% factor %>% levels  %>% length )==2 ){
    DF_functionals[,param.Y.name]<-DF_functionals[,param.Y.name] %>% as.factor
    task_MAE_functionals<-makeClassifTask(data = DF_functionals, target = param.Y.name, positive = param_positive_y_level) # with covariates!
    ## Regr
  } else task_MAE_functionals<-makeRegrTask(data = DF_functionals, target = param.Y.name) # with covariates!

  return(task_MAE_functionals)
}
