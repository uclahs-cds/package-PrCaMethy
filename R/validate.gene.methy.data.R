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
#' * `features.you.can.predict` logical vector indicating which features you can predict (i.e. you have the required genes with missing data rates < prop.missing.cutoff)
#' * `required.genes` a list of genes required by each model
#' * `missing.genes` a list of genes that are required but completely missing in the data
#' * `required.genes.with.high.missing` a list of genes that are required and have a proportion of missing values greater than `prop.missing.cutoff`

#' @examples
#'data('all.models');
#'
#'### example gene-level methylation data
#'data('example.data.gene.methy');
#'# note this dataset is derived from the following commands:
#'# data('example.data');
#'# example.data.gene.methy <- gene.methylation(example.data);
#'check <- validate.gene.methy.data(example.data.gene.methy, all.models);
#'stopifnot(check$val.passed);
#'
#'# genes required to fit each model:
#' #check$required.genes;
#'
#'# genes that are required but completely missing in your data:
#' #check$missing.genes;
#'
#'# genes that are required and have a proportion of missing values greater than `prop.missing.cutoff`
#' #check$required.genes.with.high.missing;
validate.gene.methy.data <- function(gene.methy.data, models, prop.missing.cutoff = 0.3) {
    # check if all columns of gene.methy.data are numeric class
    stopifnot('all columns in gene.methy.data should be numeric class' = all(sapply(gene.methy.data, is.numeric)));

    required.genes <- lapply(
        X = models,
        FUN = function(x) {
            if (inherits(x, 'randomForest')) {
                return(x$xNames);
            } else if (inherits(x, 'glmnet')) {
                return(x$xNames);
            } else {
                stop('each model in the models list object must class randomForest or glmnet');
                }
            }
        );
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
    features.you.can.predict <- sapply(required.genes.with.high.missing, function(x) length(x) == 0);
    return(list(
        val.passed = val.passed,
        features.you.can.predict = features.you.can.predict,
        required.genes = required.genes,
        missing.genes = missing.genes,
        required.genes.with.high.missing = required.genes.with.high.missing
        ));
    }
