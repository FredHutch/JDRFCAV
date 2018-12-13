#' Dataset, JDRF modules sets tibble.
#'
#' A tibble containing 10 meta-lists of module / gene sets.
#' use function Fun_Load_module_sets() to unlist into a single long list that include all gene sets
#'
#' @format A tibble with 10 rows and 2 variables:
#' \describe{
#'   \item{Gmt.name}{meta gene set name}
#'   \item{GMT.list}{nested lists}
#'   ...
#' }
#' @source \url{http://software.broadinstitute.org/gsea/msigdb}
#' JDRF internal curation, and Broad's GSEA / MSigDB
"JDRF_modules_sets_tibble"



