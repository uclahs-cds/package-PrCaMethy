test_that(
    desc = 'check data format',
    code = {
        data(example.data);
        check.num <- sapply(example.data, is.numeric);
        check.range <- sapply(
            X = example.data,
            FUN = function(x) {
                x <- na.omit(x);
                all(x >= 0 & x <= 1);
                }
            );
        expect_true(all(check.num));
        expect_true(all(check.range));
        }
    );

test_that(
    desc = 'required CpGs',
    code = {
        data(example.data);
        data(subtype.model.pamr);
        required.cpgs <- rownames(subtype.model.pamr$centroids);

        expect_true(all(startsWith(required.cpgs, 'cg')));

        expect_true(all(required.cpgs %in% colnames(example.data)));

        example.data.sub <- example.data[,required.cpgs];
        check.miss <- apply(example.data.sub, 2, function(x) mean(is.na(x)));
        expect_true(all(check.miss <= 0.3));
        }
    );
