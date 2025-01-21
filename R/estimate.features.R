#' Predict clinical and molecular features
#'
#' Main function used to predict various clinical and molecular features from gene-level methylation data in prostate cancer patients.
#' @inheritParams validate.gene.methy.data
#' @validate.data TRUE/FALSE, whether to validate input data.  This should generally always be TRUE, but developers may set to FALSE to speedup testing/development.
#' @export
#' @return
#' A list with the following objects:
#' * `features` a data frame with predicted features as columns and patients as rows
#' * `validation` results from validation checks telling you which features you were able to predict.  See \link{validate.gene.methy.data} for more details.
#' @examples
#'data('all.models');
#'
#'### example gene-level methylation data
#'data('example.data.gene.methy');
#'# note this dataset is derived from the following commands:
#'# data('example.data');
#'# example.data.gene.methy <- gene.methylation(example.data);
#'
#'features <- estimate.features(example.data.gene.methy, all.models);
#'str(features);
estimate.features <- function(gene.methy.data, models, prop.missing.cutoff = 0.3, validate.data = TRUE) {
    if (validate.data) {
        check <- validate.gene.methy.data(gene.methy.data, models, prop.missing.cutoff);
        if (all(check$features.you.can.predict == FALSE)) {
            print('Error: gene.methy.data does not have the required genes (with an an acceptable level of missing data, i.e. proportion missing < prop.missing.cutoff) for any of the features. See returned output for more info.')
            return(check);
            }
        if (any(check$features.you.can.predict == FALSE)) {
            stopifnot(all(names(check$features.you.can.predict) == names(models)));
            models <- models[names(check$features.you.can.predict)[which(check$features.you.can.predict)]];
            print('Warning: gene.methy.data does not have the required genes (with an an acceptable level of missing data, i.e. proportion missing < prop.missing.cutoff) for some of the features. see $validation output for more info.')
            }
    } else {
        check <- NULL;
    }
    # impute missing values
    if (sum(is.na(gene.methy.data)) > 0) {
        print('Starting imputation...');
        base::invisible(utils::capture.output(gene.methy.data.imp <- impute::impute.knn(t(gene.methy.data))$data));
        gene.methy.data.imp <- data.frame(t(gene.methy.data.imp), check.names = FALSE);
        print('Finished imputation.')
    } else {
        gene.methy.data.imp <- gene.methy.data;
        }

    # requireNamespace in order to get predict() S3 methods to work correctly
    requireNamespace('glmnet', quietly = TRUE);
    requireNamespace('randomForest', quietly = TRUE);

    features <- lapply(
        X = seq_along(models),
        FUN = function(x) {
            m <- models[[x]];
            #print(names(models)[x]);
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
    check2 <- sapply(features, function(x) all(rownames(x) == rownames(gene.methy.data.imp)));
    stopifnot(all(check2));
    features <- do.call(data.frame, features);
    stopifnot(all(rownames(features) == rownames(gene.methy.data.imp)));
    return(list(features = features, validation = check));
    }
