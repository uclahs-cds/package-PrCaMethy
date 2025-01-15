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
#'str(features);
predict.features <- function(gene.methy.data, models, prop.missing.cutoff = 0.3) {
    check <- validate.gene.methy.data(gene.methy.data, models, prop.missing.cutoff);
    if (!check$val.passed) {
        print('Error: gene.methy.data has genes with high missingness that are required for predicting features.  See the returned results for more details.')
        return(check);
        }
    # impute missing values
    gene.methy.data.imp <- impute::impute.knn(t(gene.methy.data))$data;
    gene.methy.data.imp <- data.frame(t(gene.methy.data.imp), check.names = FALSE);
    #pred.rf <- utils::getS3method('predict', 'randomForest');
    #pred.glmnet <- utils::getS3method('predict', 'glmnet');
    #pred.lognet <- utils::getS3method('predict', 'lognet');

    # requireNamespace in order to get predict() S3 methods to work correctly
    requireNamespace('glmnet', quietly = TRUE);
    requireNamespace('randomForest', quietly = TRUE);

    features <- lapply(
        X = seq_along(models),
        FUN = function(x) {
            m <- models[[x]];
            print(names(models)[x]);
            temp <- gene.methy.data.imp[,m$xNames];
            if (inherits(m, 'randomForest')) {
                pred <- data.frame(predict(m, newdata = temp, type = 'response'), check.names = FALSE);
                #pred <- data.frame(randomForest:::predict.randomForest(m, newdata = temp, type = 'response'), check.names = FALSE);
                #pred <- data.frame(pred.rf(m, newdata = temp, type = 'response'), check.names = FALSE);
            } else if (inherits(m, 'multnet')) {
                #pred <- pred.glmnet(
                pred <- predict(
                #pred <- glmnet:::predict.multnet(
                    object = m,
                    newx = as.matrix(temp),
                    s = m$best.lambda,
                    type = 'class'
                    );
            } else if (inherits(m, 'lognet')) {
                pred <- predict(
                #pred <- glmnet::predict.glmnet(
                    object = m,
                    newx = as.matrix(temp),
                    s = m$best.lambda,
                    type = 'class'
                    );
            } else if (inherits(m, 'elnet')) {
                #pred <- pred.glmnet(
                pred <- predict(
                #pred <- glmnet::predict.glmnet(
                    object = m,
                    newx = as.matrix(temp),
                    s = m$best.lambda,
                    type = 'response'
                    );
            } else {
                stop('each model in the models list object must class randomForest, elnet, lognet, or multnet (the latter 3 are glmnet models for gaussian, binomial or multinomial families, respectively)');
                }
            if (!is.null(m$classnames)) {
                pred <- data.frame(pred, check.names = FALSE);
                pred[,1] <- factor(pred[,1], levels = m$classnames);
                rownames(pred) <- rownames(temp);
                }
            colnames(pred)[1] <- names(models)[x];
            return(pred);
            }
        );
    check <- sapply(features, function(x) all(rownames(x) == rownames(gene.methy.data.imp)));
    stopifnot(all(check));
    features <- do.call(data.frame, features);
    stopifnot(all(rownames(features) == rownames(gene.methy.data.imp)));
    return(features);
    }
