#' predict.features
#'
#' Main function used to predict various clinical and molecular features from gene-level methylation data in prostate cancer patients.
#' @inheritParams validate.gene.methy.data
#' @export "predict.features"
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
#'features <- predict.features(example.data.gene.methy, all.models);
#'features;
predict.features <- function(gene.methy.data, models, prop.missing.cutoff = 0.3) {
    check <- validate.gene.methy.data(gene.methy.data, models, prop.missing.cutoff);
    if (!check$val.passed) {
        print('Error: gene.methy.data has genes with high missingness that are required for predicting features.  See the returned results for more details.')
        return(check);
        }
    # impute missing values
    gene.methy.data.imp <- impute::impute.knn(t(gene.methy.data))$data;
    gene.methy.data.imp <- data.frame(t(gene.methy.data.imp), check.names = FALSE);
    features <- lapply(
        X = seq_along(models),
        FUN = function(x) {
            m <- models[[x]];
            print(names(models)[x]);
            temp <- gene.methy.data.imp[,m$xNames];
            if (inherits(m, 'randomForest')) {
                pred <- data.frame(predict(m, newdata = temp), check.names = FALSE);
            } else if (inherits(m, 'glmnet')) {
                pred <- predict(
                    object = m,
                    newx = as.matrix(temp),
                    s = m$lambdaOpt,
                    type = 'response'
                    );
            } else {
                stop('each model in the models list object must class randomForest or glmnet');
                }
            colnames(pred)[1] <- names(models)[x];
            return(pred);
            }
        );
    check <- sapply(features, function(x) all(rownames(x) == rownames(features[[1]])));
    stopifnot(all(check));
    features <- do.call(data.frame, features);
    stopifnot(all(rownames(features) == rownames(gene.methy.data.imp)));
    return(features);
    }
