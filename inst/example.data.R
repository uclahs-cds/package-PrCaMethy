devtools::load_all();
path.data <- '/hot/project/disease/ProstateTumor/PRAD-000101-MethySubtypes/data/2024-07-16_pooled_tumour_normal_all_cpgs_all_cohorts.rds';

example.data <- readRDS(path.data);
rownames(example.data) <- example.data$patient.id;
example.data <- subset(example.data, subset = cohort == 'TCGA', select = - c(patient.id, cohort));

example.data.gene.methy <- gene.methylation(example.data);

### required genes
data(all.models);
required.genes <- lapply(
        X = all.models,
        FUN = function(x) {
            if (inherits(x, 'randomForest')) {
                return(x$xNames);
            } else if (inherits(x, 'glmnet')) {
                return(x$xNames);
            } else {
                stop('each model in the models list object must class randomForest or glmnet');
                }
            }
        );
genes <- unique(unlist(required.genes));
stopifnot(all(genes %in% colnames(example.data.gene.methy)));

# check whether '-' is in genes
check <- grepl('-', genes);
table(check);
names(genes)[check];

miss <- apply(example.data.gene.methy[,genes], 1, function(x) mean(is.na(x)));
sum(miss == 0);
nomiss <- names(miss)[miss == 0];
example.data.no.miss <- example.data[nomiss,];

set.seed(1234);
example.data <- example.data.no.miss[sample(1:nrow(example.data.no.miss), 10),];
dim(example.data);

data('subtype.model');
required.cpgs <- rownames(subtype.model$centroids);
check <- apply(example.data[,required.cpgs, drop = FALSE], 2, function(x) mean(is.na(x)));
stopifnot(max(check) <= 0.3);

format(object.size(example.data), 'Mb');

check <- gene.methylation(example.data);
stopifnot(all(genes %in% colnames(check)));

# only keep CpGs necessary for outcomes
data('gene.promoter.cpgi');
stopifnot(all(genes %in% names(gene.promoter.cpgi)));
genes[!genes %in% names(gene.promoter.cpgi)];

reqired.genes.cpgs <- unique(unlist(gene.promoter.cpgi[genes]));
cpgs.keep <- colnames(example.data)[colnames(example.data) %in% c(required.cpgs, reqired.genes.cpgs)];
example.data <- example.data[,cpgs.keep, drop = FALSE];

check2 <- gene.methylation(example.data);
stopifnot(all(genes %in% colnames(check2)));

usethis::use_data(example.data, overwrite = TRUE, compress = 'xz');
