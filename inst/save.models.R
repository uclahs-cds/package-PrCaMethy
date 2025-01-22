res.date <- '2025-01-13';

discrete.methy <- FALSE;
compress <- 'xz';
test.mode <- FALSE;

path.ml.res <- paste0('/hot/project/disease/ProstateTumor/PRAD-000101-MethySubtypes/output/prediction/', res.date, '_F72-predict-clinical-and-drivers_gene-methy_association-filter_discrete-methy', discrete.methy, '.RData');
tolerance <- 0.02 # use smallest model within __ of best model

load(path.ml.res);

outcomes <- unique(ml.res.params$outcome);

if (test.mode) {
    outcomes <- c('log2.psa.continuous', 't.stage');
    example.models.index <- 1:2
} else {
    example.models.index <- c(1,4);
    }

####
reduce.glmnet.memory <- function(glmnet.fit, lambda.opt) {
    # Identify the index of the optimal lambda
    lambda.idx <- which(glmnet.fit$lambda == lambda.opt);
    if (length(lambda.idx) == 0) {
        stop('The specified lambda value is not in the glmnet fit object.')
        }

    # Check if the model is multinomial
    is.multinomial <- is.list(glmnet.fit$beta);

    if (is.multinomial) {
    # Multinomial case: beta is a list of matrices (one for each class)
    glmnet.fit$beta <- lapply(
        X = glmnet.fit$beta,
        FUN = function(beta.matrix) {
            # Zero out coefficients for non-optimal lambdas
            beta.matrix[, -lambda.idx] <- 0
            beta.matrix
            }
        );
    glmnet.fit$a0 <- glmnet.fit$a0;  # Keep intercepts intact for all lambdas
    } else {
        # Standard case: beta is a single matrix
        beta.matrix <- glmnet.fit$beta;
        beta.matrix[, -lambda.idx] <- 0;  # Zero out non-optimal lambda columns
        glmnet.fit$beta <- beta.matrix;
        glmnet.fit$a0 <- glmnet.fit$a0;  # Keep intercepts intact for all lambdas
        }

    return(glmnet.fit)
    }
####

# devtools::load_all();
# data(example.data.gene.methy);
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
            model <- fit.glmnet$finalModel2;
            lam <- fit.glmnet$bestTune$lambda;
            xnames <- fit.glmnet$finalModel$xNames;

            ### reduce glmnet model size
            model.red <- reduce.glmnet.memory(model, lam);
            # stopifnot(all(xnames %in% colnames(example.data.gene.methy)));

            # newx <- as.matrix(example.data.gene.methy[, xnames]);
            # stopifnot(sum(is.na(newx)) == 0)

            # pred.red <- predict(
            #     object = model.red,
            #     newx = newx,
            #     s = lam,
            #     type = 'response'
            #     );
            # pred.full <- predict(
            #     object = model,
            #     newx = newx,
            #     s = lam,
            #     type = 'response'
            #     );
            # stopifnot(identical(pred.red, pred.full))
            # format(object.size(model), 'Mb');
            # format(object.size(model.red), 'Mb');
            model <- model.red;
            model$best.lambda <- lam;
            model$xNames <- xnames;
        } else {
            model <- fit.rf$finalModel;
            }
        #print(format(object.size(model), 'Mb'));
        #lapply(model, function(x) format(object.size(x), 'Mb'))
        #model$terms <- NULL; # unncessary high memory object
        #model$trainingData <- NULL;
        print(format(object.size(model), 'Mb'));

        return(model);
        }
    );

model.size <- sapply(models, function(x) format(object.size(x), 'Mb'));
model.size <- as.numeric(gsub(' Mb', '', model.size));
max(model.size);
median(model.size);


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

# Remove log2.psa.continuous due to potential cohort bias
# Remove age.categorical since unnecessary, can just use continuous age.
models <- models[!names(models) %in% c('log2p1.psa.continuous', 'age.categorical')];
length(models);
names(models);

# save all models
all.models <- models;
usethis::use_data(all.models, overwrite = TRUE, compress = compress);
