# To do:
# fix line 26 since there could be 2 min top.features (i.e. both glmnet and RF), and pick model with better est
res.date <- '2024-12-13';
path.ml.res <- paste0('/hot/project/disease/ProstateTumor/PRAD-000101-MethySubtypes/output/prediction/', res.date, '_F72-predict-clinical-and-drivers_gene-methy_association-filter_discrete-methyFALSE.RData');
tolerance <- 0 # use smallest model within __ of best model

load(path.ml.res);

outcomes <- unique(ml.res.params$outcome);

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

        mem <- data.frame(
            outcome = outcome,
            top.features = res$top.features,
            rf.mb = as.numeric(gsub(' Mb', '', format(object.size(fit.rf), 'Mb'))),
            glmnet.mb = as.numeric(gsub(' Mb', '', format(object.size(fit.glmnet), 'Mb')))
            );

        if (res$model == 'glmnet') {
            model <- fit.glmnet;
        } else {
            model <- fit.rf
            }
        #print(format(object.size(model), 'Mb'));
        lapply(model, function(x) format(object.size(x), 'Mb'))
        model$terms <- NULL;
        model$trainingData <- NULL;
        print(format(object.size(model), 'Mb'));


        return(list(
            model = model,
            mem = mem
            ));
        }
    );
names(models) <- outcomes;
# save(
#     models,
#     file = 'models.RData'
#     );

print(format(object.size(models), 'Mb'));

usethis::use_data(models, overwrite = TRUE);
