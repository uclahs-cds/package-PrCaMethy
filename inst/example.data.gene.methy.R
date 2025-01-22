### save example data with gene methylation
devtools::load_all();
data(example.data);

example.data.gene.methy <- gene.methylation(example.data);

usethis::use_data(example.data.gene.methy, overwrite = TRUE, compress = 'xz');

stopifnot(max(apply(example.data.gene.methy, 2, function(x) mean(is.na(x)))) <= 0.3);

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
