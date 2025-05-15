devtools::load_all();
library(IlluminaHumanMethylation450kanno.ilmn12.hg19);
source('config.R') # see project PRAD-000101-MethySubtypes/PrCaMethy/config.R

data(Other);
cpg.annotation <- as.data.frame(Other);
cpg.annotation$promoter <- ifelse(
        test = grepl('TSS200|TSS1500|5\'UTR', cpg.annotation$UCSC_RefGene_Group) | cpg.annotation$Regulatory_Feature_Group == 'Promoter_Associated',
        yes = 'yes',
        no = 'no'
        );
table(cpg.annotation$promoter);

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

cpg.annotation.450k <- cpg.annotation;
#usethis::use_data(cpg.annotation.450k, overwrite = TRUE, compress = 'xz');
saveRDS(
    cpg.annotation.450k,
    file = file.path(arg$path.save.annot, paste0(Sys.Date(), '_cpg_annotation_450k.rds')),
    );
