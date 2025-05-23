#' @title Predict methylation subtype
#' @description Assign patients to four prostate cancer DNA methylation subtypes
#' @inheritParams validate.subtype.model.cpgs
#' @param subtype.model Which subtype model to use ('PAMR' or 'RF' for random forest).  Although slower, we recommend 'RF' for its increased accuracy and intrinsic imputation for missing values.  Further, if some of the required CpGs are completely missing, then you must use 'RF'.
#' @param pamr.impute.using.all.cpgs If using `subtype.model = 'PAMR'`, should imputation be done using all CpGs in `methy.data` (TRUE) or only the CpGs required by \link{subtype.model.pamr} (FALSE).  When TRUE, imputation will be slower and use more memory, but should be more accurate.
#' @param seed integer seed used for imputation.
#' @export
#' @return
#' * `subtypes`: data.frame with the estimated subtypes and sample IDs (rownames of `methy.data`)
#' * `validation`: output from \link{validate.subtype.model.cpgs} to check if `methy.data` contains the required CpGs and whether any CpGs have high missingness.
#' @examples
#'### example CpG data
#'data('example.data');
#'
#'subtypes <- estimate.subtypes(example.data);
#'
#'# estimated subtypes
#'head(subtypes$subtypes);
#'
#'# validation results:
#'# length(subtypes$validation$required.cpgs)
#'# length(subtypes$validation$required.cpgs.with.high.missing)
#'# length(subtypes$validation$missing.cpgs)
estimate.subtypes <- function(
    methy.data,
    subtype.model = 'RF',
    prop.missing.cutoff = 0.3,
    pamr.impute.using.all.cpgs = TRUE,
    seed = 123
    ) {
    set.seed(seed);
    stopifnot('subtype.model must be RF or PAMR' = subtype.model %in% c('PAMR', 'RF'));
    stopifnot('prop.missing.cutoff must be between 0 and 1' = prop.missing.cutoff >= 0 & prop.missing.cutoff <= 1);
    stopifnot('PAMR cannot handle CpGs that are 100% missing; consider using RF if some of the required CpGs are completely missing' = !(prop.missing.cutoff == 1 & subtype.model == 'PAMR'));

    check <- validate.subtype.model.cpgs(methy.data, prop.missing.cutoff);
    if (!check$val.passed & prop.missing.cutoff < 1) {
        message('Error: methy.data has CpGs with high missingness that are required for predicting subtypes.  See the returned $validation results for more details.  If you insist on predicting the subtypes despite the high missingness (which will decrease the accuracy of subtype assignment), consider using subtype.model = \'RF\' with prop.missing.cutoff = 1.');
        return(list(
            subtypes = NULL,
            validation = check
            ));
        }

    ### PAMR
    if (subtype.model == 'PAMR') {
        # Impute missing values
        if (sum(is.na(methy.data)) == 0) {
            methy.data.imp <- methy.data;
        } else {
            print('Starting imputation...');
            if (!pamr.impute.using.all.cpgs) {
                methy.data <- methy.data[,check$required.cpgs, drop = FALSE];
                }
            base::invisible(utils::capture.output(methy.data.imp <- impute::impute.knn(t(methy.data))$data));
            methy.data.imp <- data.frame(t(methy.data.imp), check.names = FALSE);
            print('Finished imputation.');
            }
        data(subtype.model.pamr, envir = environment());
        methy.data.imp.sub <- methy.data.imp[,check$required.cpgs];
        methy.data.imp.sub <- t(methy.data.imp.sub);

        subtypes <- pamr::pamr.predict(
            fit = subtype.model.pamr,
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
        }
    ### RF
    if (subtype.model == 'RF') {
        methy.data <- data.frame(methy.data, check.names = FALSE);
        data(subtype.model.rf, envir = environment());
        subtype.cpgs <- colnames(subtype.model.rf$xvar);

        # for cpgs in subtype.cpgs that are not present in methy.data column names, add them as a new column of NAs.  RF will impute them.
        missing.cpgs <- setdiff(subtype.cpgs, colnames(methy.data));

        if (length(missing.cpgs) > 0) {
            message(sprintf(
                'Warning: %d of %d required CpGs are missing from the data. See the $validation outcome for more details.  Although random forest imputes missing values, having many CpGs that are missing may decrease accuracy of subtype assignment.',
                length(missing.cpgs),
                length(subtype.cpgs)
                ));
            for (cpg in missing.cpgs) {
                methy.data[,cpg] <- NA;
                };
            }
        stopifnot(all(subtype.cpgs %in% colnames(methy.data)));
        methy.data <- methy.data[,subtype.cpgs, drop = FALSE];
        stopifnot(all(colnames(methy.data) == subtype.cpgs));

        subtypes <- predict(
            object = subtype.model.rf,
            newdata = methy.data,
            na.action = 'na.impute'
            );
        stopifnot(length(subtypes$class) == nrow(methy.data));
        subtypes <- data.frame(
            sample.id = rownames(methy.data),
            subtype = subtypes$class
            );
        rownames(subtypes) <- subtypes$sample.id;
        }

    return(list(
        subtypes = subtypes,
        validation = check
        ));
    }
