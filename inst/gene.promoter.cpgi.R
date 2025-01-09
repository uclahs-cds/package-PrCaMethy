data(cpg.annotation);
data(example.data);
methy <- example.data; # methylation data with patient ids as rownames, CpGs as columns
print.progress <- TRUE;

# subset to CpG islands in gene promoter region
cpg.annotation <- subset(cpg.annotation, subset = genomic.location == 'Promoter' & !is.na(UCSC_RefGene_Name) & Relation_to_Island == 'Island');
genes <- unique(unlist(strsplit(cpg.annotation$UCSC_RefGene_Name, ';')));
names(genes) <- genes;


# the gene column can contain multiple genes separated by ;
progressr::handlers(list(
    progressr::handler_progress(
        format   = ':spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta',
        width    = 60,
        complete = '+'
        )
    ));
# progressr::with_progress({
#     p <- progressr::progressor(along = 1:length(genes));
    gene.promoter.cpgi <- lapply(
    X = seq_along(genes),
    FUN = function(x) {
        #p(); # update progress bar
        #print(x);
        gene <- genes[x];
        index <- grep(paste0('^', gene, '|;', gene, '$|;', gene, ';'), cpg.annotation$UCSC_RefGene_Name);
        cpgs <- cpg.annotation$cpg[index];

        if (x %% 100 == 0) {
            print(x);
            }
        return(cpgs[cpgs %in% colnames(methy)]);
        }
    );
    # },
    # enable = print.progress
    # );
stopifnot(length(gene.promoter.cpgi) == length(genes));
names(gene.promoter.cpgi) <- genes;

usethis::use_data(gene.promoter.cpgi, overwrite = TRUE);
