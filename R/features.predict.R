#' features.predict
#'
#' Main function used to predict various clinical and molecular features from gene-level methylation data in prostate cancer patients.
#' @inheritParams validate.gene.methy.data
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
#'features <- features.predict(example.data.gene.methy, example.models);
#'features;
features.predict <- function(gene.methy.data, models) {
    check <- validate.gene.methy.data(gene.methy.data, models);
    if (!check$val.passed) {
        print('Error: gene.methy.data is missing genes that are required for predicting features.  See the returned results for a list of missing genes per feature.')
        return(check);
        }
    features <- lapply(models, function(x) caret:::predict.train(x, newdata = gene.methy.data));
    features <- do.call(data.frame, features);
    stopifnot(all(rownames(features) == rownames(gene.methy.data)));
    return(features);
    }
