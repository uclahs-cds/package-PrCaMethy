test_that(
    desc = 'check predictions',
    code = {
        data(subtype.model);
        required.cpgs <- rownames(subtype.model$centroids);


        # example predictions
        data(example.data);
        example.data <- example.data[,required.cpgs];

        # estimate.subtypes should be able to handle missing data
        expect_true(sum(is.na(example.data)) > 0);


        subtypes <- estimate.subtypes(example.data);

        expect_true(inherits(subtypes$subtype, 'factor'));
        expect_true(all(levels(subtypes$subtype) == paste0('MS-', 1:4)));
        expect_true(all(rownames(subtypes) == rownames(example.data)));
        }
    );
