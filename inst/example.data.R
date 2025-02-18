devtools::load_all();
source('config.R') # see project PRAD-000101-MethySubtypes/PrCaMethy/config.R
library(impute);

example.data <- readRDS(arg$path.example.methy.data);
rownames(example.data) <- example.data$patient.id;
example.data <- subset(example.data, subset = cohort == 'TCGA', select = - c(patient.id, cohort));

# remove features with 100% missing
example.data <- example.data[, apply(example.data, 2, function(x) mean(is.na(x)) < 1), drop = FALSE];

# put in long format
example.data <- t(example.data);

# impute missing values
example.data <- data.frame(impute.knn(as.matrix(example.data))$data, check.names = FALSE);
stopifnot(sum(is.na(example.data)) == 0);

example.data <- t(example.data);

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

# check whether '-' is in genes
check <- grepl('-', genes);
table(check);
names(genes)[check];


check <- gene.methylation(example.data);
stopifnot(all(genes %in% colnames(check)));
stopifnot(sum(is.na(check)) == 0);

# only keep CpGs necessary for outcomes
data('gene.promoter.cpgi');
stopifnot(all(genes %in% names(gene.promoter.cpgi)));
genes[!genes %in% names(gene.promoter.cpgi)];

reqired.genes.cpgs <- unique(unlist(gene.promoter.cpgi[genes]));
cpgs.keep <- colnames(example.data)[colnames(example.data) %in% c(required.cpgs, reqired.genes.cpgs)];
example.data <- example.data[,cpgs.keep, drop = FALSE];


format(object.size(example.data), 'Mb');

# randomly keep only 2 patients:
set.seed(123);
example.data <- example.data[sample(rownames(example.data), 2), ];
format(object.size(example.data), 'Mb');


usethis::use_data(example.data, overwrite = TRUE, compress = 'xz');
