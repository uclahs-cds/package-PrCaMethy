test_that(
    'required genes', {
        data(all.models);
        data(example.data.gene.methy);
        data(gene.promoter.cpgi);
        required.genes <- unique(unlist(sapply(all.models, function(x) x$xNames)));

        expect_true(all(required.genes %in% names(gene.promoter.cpgi)));
        gene.promoter.cpgi.sub <- gene.promoter.cpgi[required.genes];
        check <- sapply(
            X = gene.promoter.cpgi.sub,
            FUN = function(x) {
                length(x) > 0 & all(startsWith(x, 'cg'))
                }
            );
        expect_true(all(check));
        }
    );
