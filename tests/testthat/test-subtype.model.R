test_that(
    desc = 'check model format',
    code = {
        data(subtype.model);
        expect_true(inherits(subtype.model, 'pamrtrained'));
        required.cpgs <- rownames(subtype.model$centroids);
        expect_true(all(startsWith(required.cpgs, 'cg')));
        expect_true(length(required.cpgs) == 5486);
        
        # example predictions
        data(example.data);
        example.data <- example.data[,required.cpgs];
        # median imputation
        example.data.imp <- apply(
            X = example.data,
            MARGIN = 2,
            FUN = function(column) {
                column[is.na(column)] <- median(column, na.rm = TRUE)
                return(column);
                }
            );
        example.data.imp <- t(example.data.imp);
        expect_true(sum(is.na(example.data.imp)) == 0);
        subtypes <- pamr::pamr.predict(
            fit = subtype.model,
            newx = example.data.imp,
            type = 'class',
            threshold = 0
            );
        expect_true(inherits(subtypes, 'factor'));
        expect_true(all(levels(subtypes) == paste0('MS-', 1:4)));
        }
    );
