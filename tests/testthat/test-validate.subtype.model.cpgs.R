test_that(
    'catch missing CpGs', {
        data(example.data);
        data(subtype.model.pamr);
        required.cpgs <- rownames(subtype.model.pamr$centroids);
        example.data <- example.data[, required.cpgs, drop = FALSE];
        check <- validate.subtype.model.cpgs(example.data);
        expect_true(check$val.passed);

        set.seed(123);
        example.data.miss <- example.data[, sample(1:ncol(example.data), 1000)];
        check2 <- validate.subtype.model.cpgs(example.data.miss);
        expect_false(check2$val.passed);

        miss.cpgs <- required.cpgs[!required.cpgs %in% colnames(example.data.miss)];
        expect_equal(check2$missing.cpgs[order(check2$missing.cpgs)], miss.cpgs[order(miss.cpgs)]);


        # high missing is caught
        cpgs <- colnames(example.data.miss)[c(10, 100)];
        for (i in cpgs) {
            example.data.miss[sample(1:nrow(example.data.miss), nrow(example.data.miss) / 2),i] <- NA;
            }
        check3 <- validate.subtype.model.cpgs(example.data.miss);
        expect_true(all(cpgs %in% names(check3$required.cpgs.with.high.missing)));
        }
    );
