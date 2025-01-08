#' validate.gene.methy.data
#'
#' Check whether `gene.methy.data` contains all genes required by `models`.
#'
#' @param gene.methy.data A data frame with gene-level methylation data, created by \link{gene.methylation}.  Patients are rows and columns are genes.
#' @param models A list of models used to predict features from gene-level methylation data.  The models should come from `data('all.models')`.
#' @export
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
#' check$required.genes; # genes required to fit each model
#' check$missing.genes; # genes that are required but missing in your data
validate.gene.methy.data <- function(gene.methy.data, models) {
    required.genes <- lapply(models, function(x) colnames(x$ptype));
    missing.genes <- lapply(required.genes, function(x) setdiff(x, colnames(gene.methy.data)));
    val.passed <- all(sapply(missing.genes, function(x) length(x) == 0));
    return(list(
        val.passed = val.passed,
        required.genes = required.genes,
        missing.genes = missing.genes
        ));
    }
