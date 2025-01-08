#res.date <- '2024-12-13';
res.date <- '2025-01-07'; # includes knn imputation
test.mode <- TRUE;

path.ml.res <- paste0('/hot/project/disease/ProstateTumor/PRAD-000101-MethySubtypes/output/prediction/', res.date, '_F72-predict-clinical-and-drivers_gene-methy_association-filter_discrete-methyFALSE.RData');
tolerance <- 0.03 # use smallest model within __ of best model

load(path.ml.res);

outcomes <- unique(ml.res.params$outcome);

if (test.mode) {
    outcomes <- c('log2.psa.continuous', 't.stage');
    example.models.index <- 1:2
} else {
    example.models.idnex <- c(1,4);
    }

models <- lapply(
    X = seq_along(outcomes),
    FUN = function(x) {
        print(x);
        outcome <- outcomes[x];
        if (outcome %in% res.regr.test$outcome) {
            res <- res.regr.test[res.regr.test$outcome == outcome,];
        } else {
            res <- res.classif.test[res.classif.test$outcome == outcome,];
            }
        res <- res[res$metric %in% c('pearson', 'mcc'),];
        max.est <- res$estimate[which.max(res$estimate)];
        #res <- res[which.max(res$estimate),];
        res.tol <- res[res$estimate >= max.est - tolerance,];

        # use smallest model within tolerance of best model
        res <- res.tol[which(res.tol$top.features == min(res.tol$top.features, na.rm = TRUE)),];
        res <- res[which.max(res$estimate),];

        mod.id <- which(ml.res.params$outcome == outcome & ml.res.params$top.features == res$top.features);

        file <- file.path(dirname(path.ml.res), paste0(res.date, '_F72-predict-clinical-and-drivers_discrete-methyFALSE_models-', mod.id, '-', outcome, '-', res$top.features, '.RData'));

        if (!file.exists(file)) {
            res.date2 <- as.Date(res.date) - 1;
            file <- file.path(dirname(path.ml.res), paste0(res.date2, '_F72-predict-clinical-and-drivers_discrete-methyFALSE_models-', mod.id, '-', outcome, '-', res$top.features, '.RData'));
            }
        stopifnot(file.exists(file));

        load(file);


        if (res$model == 'glmnet') {
            model <- fit.glmnet;
        } else {
            model <- fit.rf
            }
        #print(format(object.size(model), 'Mb'));
        #lapply(model, function(x) format(object.size(x), 'Mb'))
        #model$terms <- NULL;
        #model$trainingData <- NULL;
        print(format(object.size(model), 'Mb'));

        return(model);
        }
    );

# fix outcome names
outcomes[outcomes == 'log2.psa.continuous'] <- 'log2p1.psa.continuous';
outcomes[outcomes == 'log2.snvs.per.mbps'] <- 'log2p1.snvs.per.mbps';
outcomes[outcomes == 't.stage'] <- 't.category';
outcomes[outcomes == 'MYC_cna'] <- 'MYC.cna.gain';
outcomes[grepl('_cna', outcomes)] <- paste0(outcomes[grepl('_cna', outcomes)], '.loss');
outcomes[outcomes == 'NKX3.1_cna.loss'] <- 'NKX3-1_cna.loss';
outcomes <- gsub('_', '.', outcomes);

names(models) <- outcomes;

print(format(object.size(models), 'Mb'));

# save all models
all.models <- models;
usethis::use_data(all.models, overwrite = TRUE);

# smaller example models for examples/testing
lapply(all.models, function(x) format(object.size(x), 'Mb'));

example.models <- all.models[example.models.index];
format(object.size(example.models), 'Mb');
usethis::use_data(example.models, overwrite = TRUE);
