#' Validate input data for estimate.subtypes()
#'
#' Check whether `methy.data` contains all CpGs required by \link{subtype.model} for assigning patients to four prostate cancer DNA methylation subtypes.
#'
#' @param methy.data A data.frame with patients as rows (rownames give patient ids) and column names give CpG ids.
#' @param prop.missing.cutoff The maximum proportion of missing values allowed for each required CpG. KNN imputation is used to impute missing values.
#' @export
#' @return
#' * `val.passed` a logical indicating whether the data passed validation
#' * `check$required.cpgs` a vector of CpG ids that are required for predicting the subtypes
#' * `missing.cpgs` a vector of CpG ids that are required but completely missing in the data
#' * `required.cpgs.with.high.missing` a vector of CpG ids that are required and have a proportion of missing values greater than `prop.missing.cutoff`

#' @examples
#'data('example.data');
#'check <- validate.subtype.model.cpgs(example.data);
#'stopifnot(check$val.passed);
#'
#'# CpGs required to fit each model:
#' #check$required.cpgs;
#'
#'# CpGs that are required but completely missing in your data:
#' #check$missing.cpgs;
#'
#'# CpGs that are required and have a proportion of missing values greater than `prop.missing.cutoff`
#' #check$required.cpgs.with.high.missing;
validate.subtype.model.cpgs <- function(methy.data, prop.missing.cutoff = 0.3) {
    # check whether all values of methy.data data.frame are between 0 and 1:
    methy.data.nomiss <- na.omit(methy.data);
    stopifnot('All values of methy.data should be between 0 and 1' = all(methy.data.nomiss >= 0 & methy.data.nomiss <= 1));

    data(subtype.model, envir = environment());
    required.cpgs <- rownames(subtype.model$centroids);
    missing.cpgs <- setdiff(required.cpgs, colnames(methy.data));
    nonmissing.cpgs <- setdiff(required.cpgs, missing.cpgs);
    if (length(nonmissing.cpgs) > 0) {
        methy.sub <- methy.data[,nonmissing.cpgs, drop = FALSE];
        }

    required.cpgs.prop.missing <- sapply(
        X = required.cpgs,
        FUN = function(y) {
            if (y %in% colnames(methy.sub)) {
                prop.missing <- mean(is.na(methy.sub[,y]));
            } else {
                prop.missing <- 1;
                }
            }
        );
    required.cpgs.with.high.missing <- lapply(required.cpgs.prop.missing, function(x) x[x > prop.missing.cutoff]);
    val.passed <- length(unlist(required.cpgs.with.high.missing)) == 0 & length(unlist(missing.cpgs)) == 0;
    val.passed;
    unlist(required.cpgs.with.high.missing)
    return(list(
        val.passed = val.passed,
        required.cpgs = required.cpgs,
        missing.cpgs = missing.cpgs,
        required.cpgs.with.high.missing = required.cpgs.with.high.missing
        ));
    }
