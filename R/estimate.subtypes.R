#' @title Predict methylation subtype
#' @description Assign patients to four prostate cancer DNA methylation subtypes
#' @inheritParams validate.subtype.model.cpgs
#' @param impute.using.all.cpgs TRUE/FALSE indicating whether to impute missing values using all CpGs in `methy.data` or only the CpGs required by \link{subtype.model}.  When TRUE, imputation will be slower and use more memory, but should be more accurate.
#' @export
#' @examples
#'data('subtype.model');
#'
#'### example CpG data
#'data('example.data');
#'
#'subtypes <- estimate.subtypes(example.data);
#'head(subtypes);
estimate.subtypes <- function(methy.data, prop.missing.cutoff = 0.3, impute.using.all.cpgs = TRUE) {
    check <- validate.subtype.model.cpgs(methy.data, prop.missing.cutoff);
    if (!check$val.passed) {
        print('Error: methy.data has CpGs with high missingness that are required for predicting subtypes.  See the returned results for more details.')
        return(check);
        }
    # impute missing values
    if (sum(is.na(methy.data)) == 0) {
        methy.data.imp <- methy.data;
    } else {
        print('Starting imputation...');
        if (!impute.using.all.cpgs) {
            methy.data <- methy.data[,check$required.cpgs, drop = FALSE];
            }
        base::invisible(utils::capture.output(methy.data.imp <- impute::impute.knn(t(methy.data))$data));
        methy.data.imp <- data.frame(t(methy.data.imp), check.names = FALSE);
        print('Finished imputation.');
        }

    # requireNamespace in order to get predict() S3 methods to work correctly
    #requireNamespace('pamr', quietly = TRUE);
    data(subtype.model.pamr, envir = environment());
    methy.data.imp.sub <- methy.data.imp[,check$required.cpgs];
    methy.data.imp.sub <- t(methy.data.imp.sub);
    subtypes <- pamr::pamr.predict(
        fit = subtype.model,
        newx = methy.data.imp.sub, # CpGs in rows, samples in columns
        type = 'class',
        threshold = 0
        );
    stopifnot(length(subtypes) == ncol(methy.data.imp.sub));
    subtypes <- data.frame(
        subtype = subtypes,
        check.names = FALSE
        );
    rownames(subtypes) <- colnames(methy.data.imp.sub);

    return(subtypes);
    }
