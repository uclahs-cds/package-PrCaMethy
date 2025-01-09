# https://www.bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylation450kanno.ilmn12.hg19.html
# v 0.6.1
library(IlluminaHumanMethylation450kanno.ilmn12.hg19);

data(Other);
cpg.annotation <- as.data.frame(Other);
cpg.annotation$genomic.location <- factor(ifelse(
    test = cpg.annotation$UCSC_RefGene_Name == '',
    yes = 'Intergenic',
    no = ifelse(
        test = grepl('TSS200|TSS1500|5\'UTR', cpg.annotation$UCSC_RefGene_Group) | cpg.annotation$Regulatory_Feature_Group == 'Promoter_Associated',
        yes = 'Promoter',
        no = 'Body'
        )
    ));
table(cpg.annotation$genomic.location);

cpg.annotation$cpg <- rownames(cpg.annotation);

data(Locations);
locations <- as.data.frame(Locations);
locations$chr <- gsub('chr', '',locations$chr);
stopifnot(all(as.character(c(1:22, 'X', 'Y')) %in% locations$chr));
locations$chr <- factor(
    x = locations$chr,
    levels = as.character(c(1:22, 'X', 'Y'))
    );
locations$cpg <- rownames(locations);
stopifnot(all(locations$cpg %in% cpg.annotation$cpg));
cpg.annotation <- merge(
    x = cpg.annotation,
    y = locations,
    by = 'cpg',
    all.x = TRUE
    );

# island, open sea, shelf, shore
data(Islands.UCSC);
islands <- as.data.frame(Islands.UCSC);
islands$cpg.type <- ifelse(
    test = grepl('Shelf', islands$Relation_to_Island),
    yes = 'Shelf',
    no = ifelse(
        test = grepl('Shore', islands$Relation_to_Island),
        yes = 'Shore',
        no = ifelse(
            test = grepl('OpenSea', islands$Relation_to_Island),
            yes = 'Open sea',
            no = islands$Relation_to_Island
            )
        )
    );
islands$cpg.type <- factor(islands$cpg.type);
islands$cpg <- rownames(islands);
stopifnot(all(cpg.annotation$cpg %in% rownames(islands)));
cpg.annotation <- merge(
    x = cpg.annotation,
    y = islands,
    by = 'cpg',
    all.x = TRUE
    );

usethis::use_data(cpg.annotation, overwrite = TRUE, compress = 'xz');
