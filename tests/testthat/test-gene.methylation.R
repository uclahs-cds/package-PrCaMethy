test_that(
    desc = 'check data format',
    code = {
        data(example.data);
        set.seed(123);
        example.data <- example.data[,sample(1:ncol(example.data), 10)];
        example.data.gene.methy <- gene.methylation(example.data);
        expect_true(all(example.data.gene.methy >= 0))
        expect_true(all(example.data.gene.methy <= 1))

        data(gene.promoter.cpgi);
        gene.promoter.cpgi <- gene.promoter.cpgi[colnames(example.data.gene.methy)];
        for (gene in colnames(example.data.gene.methy)) {
            print(gene);
            gene.cpgs <- gene.promoter.cpgi[[gene]];
            gene.cpgs <- gene.cpgs[gene.cpgs %in% colnames(example.data)];
            gene.methy <- apply(example.data[,gene.cpgs, drop = FALSE], 1, function(x) median(x, na.rm = TRUE));
            expect_true(all(gene.methy == example.data.gene.methy[,gene]));
            }
        }
    );
