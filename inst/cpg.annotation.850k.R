devtools::load_all();
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38);
source('config.R') # see project PRAD-000101-MethySubtypes/PrCaMethy/config.R

data(Other);
cpg.annotation <- as.data.frame(Other);
cpg.annotation$promoter <- ifelse(
        test = grepl('TSS200|TSS1500|5\'UTR', cpg.annotation$UCSC_RefGene_Group) | cpg.annotation$Regulatory_Feature_Group == 'Promoter_Associated',
        yes = 'yes',
        no = 'no'
        );

cpg.annotation$cpg.id.full <- rownames(cpg.annotation);
cpg.annotation$cpg <- sapply(strsplit(cpg.annotation$cpg.id.full, '_'), function(x) x[1]);

data(Locations);
locations <- as.data.frame(Locations);
locations$chr <- gsub('chr', '',locations$chr);
stopifnot(all(as.character(c(1:22, 'X', 'Y')) %in% locations$chr));
locations$chr <- factor(
    x = locations$chr,
    levels = as.character(c(1:22, 'X', 'Y'))
    );
locations$cpg.id.full <- rownames(locations);
locations$cpg <- sapply(strsplit(locations$cpg.id.full, '_'), function(x) x[1]);
stopifnot(all(locations$cpg %in% cpg.annotation$cpg));
cpg.annotation <- merge(
    x = cpg.annotation,
    y = locations,
    by = 'cpg.id.full'
    );
stopifnot(nrow(cpg.annotation) == nrow(locations));

# island, open sea, shelf, shore
data(Islands.UCSC);
islands <- as.data.frame(Islands.UCSC);
islands$cpg.type <- islands$Relation_to_Island
islands$cpg.type <- factor(islands$cpg.type);
islands$cpg.id.full <- rownames(islands);
islands$cpg <- sapply(strsplit(islands$cpg.id.full, '_'), function(x) x[1]);
cpg.annotation <- merge(
    x = cpg.annotation,
    y = islands,
    by = 'cpg.id.full',
    all.x = TRUE
    );
stopifnot(nrow(cpg.annotation) == nrow(islands));

cpg.annotation.850k <- cpg.annotation;
#usethis::use_data(cpg.annotation.450k, overwrite = TRUE, compress = 'xz');
saveRDS(
    cpg.annotation.850k,
    file = file.path(arg$path.save.annot, paste0(Sys.Date(), '_cpg_annotation_850k.rds')),
    );
