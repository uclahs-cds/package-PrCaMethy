test_that(
    desc = 'check data format',
    code = {
        data(example.data.gene.methy);
        check.num <- sapply(example.data, is.numeric);
        check.range <- sapply(
            X = example.data,
            FUN = function(x) {
                x <- na.omit(x);
                all(x >=0 & x <= 1);
                }
            );
        expect_true(all(check.num));
        expect_true(all(check.range));
        }
    );

test_that(
    desc = 'required genes',
    code = {
        data(all.models);
        data(example.data.gene.methy);
        required.genes <- unique(unlist(sapply(all.models, function(x) x$xNames)));
       
        expect_true(all(required.genes %in% colnames(example.data.gene.methy)));
        
        example.data.sub <- example.data.gene.methy[,required.genes];
        check.miss <- apply(example.data.sub, 2, function(x) mean(is.na(x)));
        expect_true(all(check.miss <= 0.3));
        }
    );