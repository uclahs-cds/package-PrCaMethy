#' gene.methylation
#'
#' Calculate gene-level methylation for a dataset containing CpGs from the Illumina 450k methylation array.
#' Gene-level methylation is calculated as the median beta-value among CpG islands in the gene promoter region.
#'
#' @param methy methylation dataset where rownames give patient ids and columns use CpG ids
#' @param print.progress TRUE/FALSE to show progress bar
#' @export
#' @importFrom stats median predict
#' @importFrom utils data
#' @examples
#'data(example.data);
#'
#'### full example (may take a few mins to run)
#'#example.data.gene.methy <- gene.methylation(example.data);
#'
#'### fast example
#'set.seed(123);
#'example.data <- example.data[,sample(1:ncol(example.data), 10)];
#'example.data.gene.methy <- gene.methylation(example.data);
gene.methylation <- function(methy, print.progress = TRUE) {
    data(gene.promoter.cpgi, envir = environment());

    methy <- methy[,colnames(methy) %in% unlist(gene.promoter.cpgi)];

    # only keep elements in the list gene.promoter.cpgi that contain CpGs in the methy dataset
    gene.promoter.cpgi <- gene.promoter.cpgi[sapply(gene.promoter.cpgi, function(x) any(x %in% colnames(methy)))];
    if (length(gene.promoter.cpgi) == 0) {
        stop('Error: no CpG islands in the gene promoter regions are present in the methy dataset');
        }

    gene.promoter.methy <- data.frame(matrix(NA, nrow = nrow(methy), ncol = length(gene.promoter.cpgi)));
    rownames(gene.promoter.methy) <- rownames(methy);
    colnames(gene.promoter.methy) <- names(gene.promoter.cpgi);

    all(unlist(gene.promoter.cpgi) %in% colnames(methy));

    progressr::handlers(list(
        progressr::handler_progress(
            format   = ':spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta',
            width    = 60,
            complete = '+'
            )
        ));
    progressr::with_progress({
        p <- progressr::progressor(along = 1:length(gene.promoter.cpgi));
        gene.promoter.median.methy <- lapply(
            X = seq_along(gene.promoter.cpgi),
            FUN = function(x) {
                p(); # update progress bar
                temp <- methy[,colnames(methy) %in% gene.promoter.cpgi[[x]], drop = FALSE];
                med <- apply(temp, 1, median, na.rm = TRUE);
                med <- data.frame(med);
                colnames(med)[1] <- names(gene.promoter.cpgi)[x];
                return(med);
                }
            );
        },
        enable = print.progress
        );

    check <- sapply(
        X = gene.promoter.median.methy,
        FUN = function(x) {
            all(rownames(x) == rownames(gene.promoter.median.methy[[1]]))
            }
        );

    gene.promoter.median.methy.final <- do.call(data.frame, gene.promoter.median.methy);

    check <- which(apply(gene.promoter.median.methy.final, 2, function(x) all(is.na(x))));
    if (length(check) > 0) {
        gene.promoter.median.methy.final <- gene.promoter.median.methy.final[,-check];
        }
    return(gene.promoter.median.methy.final);
    }
