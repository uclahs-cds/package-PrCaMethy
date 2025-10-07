#' Validate input data for estimate.subtypes()
#'
#' Check whether `methy.data` contains all CpGs required by \link{subtype.model.pamr} or \link{subtype.model.rf} for assigning patients to four prostate cancer DNA methylation subtypes.
#'
#' @param methy.data A data.frame with patients as rows (rownames give patient ids) and column names give CpG ids.
#' @param prop.missing.cutoff The maximum proportion of missing values allowed for each required CpG.
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

    data(subtype.model.pamr, envir = environment());
    required.cpgs <- rownames(subtype.model.pamr$centroids);
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

    # regardless of what the user specifies for prop.missing.cutoff, we should
    # print a warning if some CpGs have high missing.
    cpgs.high.miss.warn <- sum(required.cpgs.prop.missing > 0.5);
    if (cpgs.high.miss.warn > 0) {
        message('Warning: ', cpgs.high.miss.warn, ' out of ', length(required.cpgs) ,' required CpGs have > 50% missing values. Having many CpGs with high missing data may decrease accuracy of subtype assignment.');
        }

    ### check missingness in samples
    miss <- apply(methy.sub, 1, function(x) mean(is.na(x)));
    n.sample.high.miss <- sum(miss > 0.5);
    if (n.sample.high.miss > 0) {
        message('Warning: ', n.sample.high.miss, ' out of ', nrow(methy.sub) ,' samples have > 50% missing values. Samples with high rates of missing values may have inaccurate subtype assignment.');
        }

    required.cpgs.with.high.missing <- lapply(required.cpgs.prop.missing, function(x) x[x > prop.missing.cutoff]);
    val.passed <- length(unlist(required.cpgs.with.high.missing)) == 0 & length(unlist(missing.cpgs)) == 0;

    return(list(
        val.passed = val.passed,
        required.cpgs = required.cpgs,
        missing.cpgs = missing.cpgs,
        required.cpgs.with.high.missing = required.cpgs.with.high.missing
        ));
    }
