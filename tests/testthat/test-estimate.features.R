test_that(
    desc = 'check predicted features',
    code = {
        data(all.models);
        data(example.data.gene.methy);
        example.data.imp <- example.data.gene.methy;

        # # save time in test by imputing beforehand
        # example.data.imp <- apply(
        #     X = example.data.gene.methy,
        #     MARGIN = 2,
        #     FUN = function(column) {
        #         column[is.na(column)] <- median(column, na.rm = TRUE)
        #         return(column);
        #         }
        #     );

        features <- estimate.features(
            gene.methy.data = example.data.imp,
            models = all.models,
            validate.data = FALSE
            )$features;
        expect_true(ncol(features) == 14);

        # validate CNA variables
        cna.vars <- colnames(features)[grepl('\\.cna\\.', colnames(features))];
        expect_true(length(cna.vars) == 8);
        check.cna <- sapply(
            X = cna.vars,
            FUN = function(x) {
                expect_true(inherits(features[,x], 'factor'));
                expect_true(all(levels(features[,x]) == c('Yes', 'No')))
                }
            );
        # age.continuous
        expect_true(is.numeric(features$age.continuous));
        expect_true(min(features$age.continuous) > 40);
        expect_true(max(features$age.continuous) < 100);

        # ISUP.grade
        expect_true(inherits(features$ISUP.grade, 'factor'));
        expect_true(all(levels(features$ISUP.grade) == paste0('ISUP', 1:5)));

        # t.category
        expect_true(inherits(features$t.category, 'factor'));
        expect_true(all(levels(features$t.category) == paste0('T', 1:4)));

        # psa.categorical
        expect_true(inherits(features$psa.categorical, 'factor'));
        expect_true(all(levels(features$psa.categorical) == c('psa_lt_10', 'psa_10_19.9', 'psa_gte_20')));

        # pga
        expect_true(inherits(features$pga, 'numeric'));
        expect_true(min(features$pga) >= 0);
        expect_true(max(features$pga) <= 100);

        # log2p1.snvs.per.mbps
        expect_true(inherits(features$log2p1.snvs.per.mbps, 'numeric'));
        expect_true(min(features$log2p1.snvs.per.mbps) >= 0);
        }
    );
