#' validate.gene.methy.data
#'
#' Check whether `gene.methy.data` contains all genes required by `models` and that there is an acceptable level of missingness for each required gene.
#' Note that genes with acceptable levels of missing values are later imputed using KNN imputation when calling \link{predict.features}.
#' If you'd rather use a different imputation method, then make sure to impute missing values before calling \link{predict.features}.
#'
#' @param gene.methy.data A data frame with gene-level methylation data, created by \link{gene.methylation}.  Patients are rows and columns are genes.
#' @param models A list of models used to predict features from gene-level methylation data.  The models should come from `data('all.models')`.
#' @param prop.missing.cutoff The maximum proportion of missing values allowed for each required gene.
#' @export
#' @return
#' * `val.passed` a logical indicating whether the data passed validation
#' * `required.genes` a list of genes required by each model
#' * `missing.genes` a list of genes that are required but completely missing in the data
#' * `required.genes.with.high.missing` a list of genes that are required and have a proportion of missing values greater than `prop.missing.cutoff`

#' @examples
#'data('example.models');
#'
#'### example gene-level methylation data
#'data('example.data.gene.methy');
#'# note this dataset is derived from the following commands:
#'# data('example.data');
#'# example.data.gene.methy <- gene.methylation(example.data);
#'check <- validate.gene.methy.data(example.data.gene.methy, example.models);
#'stopifnot(check$val.passed);
#'
#' #check$required.genes; # genes required to fit each model
#' #check$missing.genes; # genes that are required but completely missing in your data
#' #check$required.genes.with.high.missing; # genes that are required and have a proportion of missing values greater than `prop.missing.cutoff`
validate.gene.methy.data <- function(gene.methy.data, models, prop.missing.cutoff = 0.2) {
    required.genes <- lapply(models, function(x) colnames(x$ptype));
    required.genes.prop.missing <- lapply(
        X = required.genes,
        FUN = function(x) {
            sapply(
                X = x,
                FUN = function(y) {
                    if (y %in% colnames(gene.methy.data)) {
                        prop.missing <- mean(is.na(gene.methy.data[,y]));
                    } else {
                        prop.missing <- 1;
                        }
                    }
                );
            }
        );
    required.genes.with.high.missing <- lapply(required.genes.prop.missing, function(x) x[x > prop.missing.cutoff]);
    missing.genes <- lapply(required.genes, function(x) setdiff(x, colnames(gene.methy.data)));
    val.passed <- length(unlist(required.genes.with.high.missing )) == 0 & length(unlist(missing.genes)) == 0;
    return(list(
        val.passed = val.passed,
        required.genes = required.genes,
        missing.genes = missing.genes,
        required.genes.with.high.missing = required.genes.with.high.missing
        ));
    }
