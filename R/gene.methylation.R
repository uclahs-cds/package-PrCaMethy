#' gene.methylation
#'
#' Calculate gene-level methylation for a dataset containing CpGs from the Illumina 450k methylation array.
#' Gene-level methylation is calculated as the median beta-value among CpG islands in the gene promoter region.
#'
#' @param methy methylation dataset where rownames give patient ids and columns use CpG ids
#' @param print.progress TRUE/FALSE to show progress bar
#' @examples
#'data(example.data);
#'example.data.gene.methy <- gene.methylation(example.data);
gene.methylation <- function(methy, print.progress = TRUE) {
    data(gene.promoter.cpgi);

    methy <- methy[,colnames(methy) %in% unlist(gene.promoter.cpgi)];

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
    stopifnot(all(check));

    gene.promoter.median.methy.final <- do.call(data.frame, gene.promoter.median.methy);
    return(gene.promoter.median.methy.final);
    }
