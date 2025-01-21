test_that(
    'catch missing genes', {
        data(example.data.gene.methy);
        data(all.models);

        model <- all.models['psa.categorical'];
        required.genes <- model[[1]]$xNames
        example.data.gene.methy <- example.data.gene.methy[,required.genes];

        check <- validate.gene.methy.data(
            gene.methy.data = example.data.gene.methy,
            models = model
            );
        expect_true(length(required.genes) == length(check$required.genes[[1]]) & all(required.genes %in% check$required.genes[[1]]));

        gene <- colnames(example.data.gene.methy)[10];

        ex.data <- example.data.gene.methy[, colnames(example.data.gene.methy) != gene];
        check2 <- validate.gene.methy.data(
            gene.methy.data = ex.data,
            models = model
            );
        expect_true(gene %in% check2$missing.genes);

        # high missing is caught
        set.seed(123);
        genes <- colnames(example.data.gene.methy)[sample(1:ncol(example.data.gene.methy), 10)];
        for(i in genes) {
            example.data.gene.methy[sample(1:nrow(example.data.gene.methy), nrow(example.data.gene.methy) / 2),i] <- NA;
            }
        check3 <- validate.gene.methy.data(
            gene.methy.data = example.data.gene.methy,
            models = model
            );
        expect_true(all(genes %in% names(check3$required.genes.with.high.missing[[1]])));
        }
    );