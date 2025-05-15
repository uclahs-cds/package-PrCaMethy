devtools::load_all();
source('config.R') # see project PRAD-000101-MethySubtypes/PrCaMethy/config.R

cpg.annotation <- readRDS(arg$path.cpg.annotation);
#data(example.data);
#methy <- example.data; # methylation data with patient ids as rownames, CpGs as columns
methy <- readRDS(arg$path.example.methy.data);
print.progress <- TRUE;
methy.cpgs <- colnames(methy)[grep('^cg', colnames(methy))];
methy <- NULL;
gc();

# subset to CpG islands in gene promoter region
cpg.annotation <- subset(cpg.annotation, subset = promoter == 'yes' & !is.na(UCSC_RefGene_Name) & cpg.type == 'Island');
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
        return(cpgs[cpgs %in% methy.cpgs]);
        }
    );
    # },
    # enable = print.progress
    # );
stopifnot(length(gene.promoter.cpgi) == length(genes));
names(gene.promoter.cpgi) <- genes;

check.genes <- c('MIR572', 'WSB2', 'GSTM1', 'TIAM2');
stopifnot(all(check.genes %in% names(gene.promoter.cpgi)));
check.genes.cpgs <- sapply(gene.promoter.cpgi[check.genes], function(x) length(x) > 0);
stopifnot(all(check.genes.cpgs));

stopifnot(!'.' %in% names(gene.promoter.cpgi));

# replace all '-' with '.' in gene names in order to match with gene names used in models.
check <- grepl('-', names(gene.promoter.cpgi));
table(check);
names(gene.promoter.cpgi)[check];
names(gene.promoter.cpgi) <- gsub('-', '.', names(gene.promoter.cpgi));

# remove genes with 0 cpgs
check <- sapply(gene.promoter.cpgi, function(x) length(x) == 0);
if (any(check)) {
    gene.promoter.cpgi <- gene.promoter.cpgi[!check];
    }

usethis::use_data(gene.promoter.cpgi, overwrite = TRUE, compress = 'xz');

# also save a data.frame long version
gene.promoter.cpgi.df <- lapply(
    X = names(gene.promoter.cpgi),
    FUN = function(x) {
        d <- data.frame(
            gene = x,
            cpg = gene.promoter.cpgi[[x]],
            stringsAsFactors = FALSE,
            check.names = FALSE
            );
        return(d);
        }
    );
gene.promoter.cpgi.df <- do.call(rbind, gene.promoter.cpgi.df);
data.table::fwrite(
    x = gene.promoter.cpgi.df,
    file = file.path(arg$path.save.annot, paste0(Sys.Date(), '_gene_promoter_cpgi.csv')),
    row.names = FALSE,
    col.names = TRUE
    );
gene.promoter.cpgi.list <- gene.promoter.cpgi;
save(
    gene.promoter.cpgi.list,
    gene.promoter.cpgi.df,
    file = file.path(arg$path.save.annot, paste0(Sys.Date(), '_gene_promoter_cpgi.RData')),
    compress = 'xz'
    );
