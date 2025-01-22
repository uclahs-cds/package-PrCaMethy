test_that(
    desc = 'check model format',
    code = {
        data(all.models);
        model.class <- lapply(all.models, class);
        allowed.classes <- c('randomForest', 'elnet', 'lognet', 'multnet', 'glmnet');

        expect_true(all(sapply(all.models, function(x) inherits(x, allowed.classes))));

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
        required.genes <- unique(unlist(required.genes));
        data(example.data.gene.methy);
        expect_true(all(required.genes %in% colnames(example.data.gene.methy)));

        }
    );
