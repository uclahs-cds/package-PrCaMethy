test_that(
    desc = 'check predicted subtypes',
    code = {
        data(subtype.model.pamr);
        required.cpgs <- rownames(subtype.model.pamr$centroids);


        # example predictions
        data(example.data);
        example.data <- example.data[,required.cpgs];

        # estimate.subtypes should be able to handle missing data
        set.seed(123);
        example.data[sample(1:prod(dim(example.data)), 10)] <- NA;

        subtypes <- estimate.subtypes(example.data, prop.missing.cutoff = 0.6);

        expect_true(inherits(subtypes$subtype, 'factor'));
        expect_true(all(levels(subtypes$subtype) == paste0('MS-', 1:4)));
        expect_true(all(rownames(subtypes) == rownames(example.data)));
        }
    );
