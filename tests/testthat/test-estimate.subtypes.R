test_that(
    desc = 'check PAMR predicted subtypes',
    code = {
        data(subtype.model.pamr);
        required.cpgs <- rownames(subtype.model.pamr$centroids);


        # example data
        data(example.data);
        example.data <- example.data[,required.cpgs];

        # estimate.subtypes should be able to handle missing data
        set.seed(123);
        example.data[sample(1:prod(dim(example.data)), 10)] <- NA;

        subtypes <- estimate.subtypes(
            methy.data = example.data,
            subtype.model = 'PAMR',
            prop.missing.cutoff = 0.6
            );
        subtypes <- subtypes$subtypes;

        expect_true(inherits(subtypes$subtype, 'factor'));
        expect_true(all(levels(subtypes$subtype) == paste0('MS-', 1:4)));
        expect_true(all(subtypes$sample.id == rownames(example.data)));
        }
    );

test_that(
    desc = 'When CpG missingness exceeds prop.missing.cutoff, error msg should be printed and validation results returned.',
    code = {
        data(subtype.model.pamr);
        data(example.data);
        required.cpgs <- rownames(subtype.model.pamr$centroids);
        example.data <- example.data[,required.cpgs];

        # estimate.subtypes should be able to handle missing data
        set.seed(123);
        example.data[sample(1:prod(dim(example.data)), 10)] <- NA;

        # Error msg should be print
        expect_message(
            estimate.subtypes(
                methy.data = example.data,
                subtype.model = 'PAMR',
                prop.missing.cutoff = 0
                ),
            'high missingness'
            );

        # Subtypes should be NULL but validation results should be returned
        check <- estimate.subtypes(
            methy.data = example.data,
            subtype.model = 'PAMR',
            prop.missing.cutoff = 0
            );
        expect_null(check$subtypes);
        expect_true(inherits(check$validation, 'list'));
        }
    );

test_that(
    desc = 'Random forest can handle CpGs with missing values',
    code = {
        data(subtype.model.rf);
        data(example.data);
        required.cpgs <- colnames(subtype.model.rf$xvar);
        example.data <- example.data[,required.cpgs];

        set.seed(123);
        example.data[sample(1:prod(dim(example.data)), 10)] <- NA;

        # slow
        subtypes <- estimate.subtypes(
            methy.data = example.data,
            subtype.model = 'RF',
            prop.missing.cutoff = 0.6
            );
        subtypes <- subtypes$subtypes;
        expect_true(inherits(subtypes$subtype, 'factor'));
        expect_true(all(levels(subtypes$subtype) == paste0('MS-', 1:4)));
        expect_true(all(subtypes$sample.id == rownames(example.data)));
        }
    );

test_that(
    desc = 'Random forest with prop.missing.cutoff = 1 can handle required CpGs that are completely missing',
    code = {
        data(subtype.model.rf);
        data(example.data);
        required.cpgs <- colnames(subtype.model.rf$xvar);

        set.seed(123);
        cpgs.rm <- sample(required.cpgs, 10);
        example.data <- example.data[,!(colnames(example.data) %in% cpgs.rm)];

        # slow
        subtypes <- estimate.subtypes(
            methy.data = example.data,
            subtype.model = 'RF',
            prop.missing.cutoff = 1
            );
        subtypes <- subtypes$subtypes;
        expect_true(inherits(subtypes$subtype, 'factor'));
        expect_true(all(levels(subtypes$subtype) == paste0('MS-', 1:4)));
        expect_true(all(subtypes$sample.id == rownames(example.data)));
        }
    );

test_that(
    desc = 'PAMR should error if prop.missing.cutoff = 1 since it cannot handle CpGs that are completely missing',
    code = {
        data(subtype.model.pamr);
        data(example.data);
        required.cpgs <- rownames(subtype.model.pamr$centroids);
        example.data <- example.data[,required.cpgs];

        # Error msg should be printed
        expect_error(
            estimate.subtypes(
                methy.data = example.data,
                subtype.model = 'PAMR',
                prop.missing.cutoff = 1
                ),
            'PAMR cannot handle CpGs that are 100% missing'
            );
        }
    );
test_that(
    desc = 'When using PAMR, previously impute.knn would throw an error if any CpGs had >50% missing or any samples had >80% missing.
    Now we want the program to run without error  in these cases, instead errors/warnings should be caught by validate.subtype.model.cpgs()',
    code = {
        data(subtype.model.pamr);
        data(example.data);
        required.cpgs <- rownames(subtype.model.pamr$centroids);
        example.data <- example.data[,required.cpgs];

        # add a third patient
        example.data <- rbind(example.data, patient = example.data[1,]);

        # CpG 1 has 66% missing
        example.data[1:2,1] <- NA;
        check <- estimate.subtypes(
            methy.data = example.data,
            subtype.model = 'PAMR',
            prop.missing.cutoff = 0.7
            );
        expect_true(check$validation$val.passed);

        # make sure that some samples have >80% missing
        example.data[1, 1:5000] <- NA;
        check <- estimate.subtypes(
            methy.data = example.data,
            subtype.model = 'PAMR',
            prop.missing.cutoff = 0.7
            );
        expect_true(check$validation$val.passed);
        }
    );
