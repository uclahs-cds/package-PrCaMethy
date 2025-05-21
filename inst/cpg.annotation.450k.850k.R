library(rtracklayer);
library(GenomicRanges);
library(R.utils);

devtools::load_all();
source('config.R') # see project PRAD-000101-MethySubtypes/PrCaMethy/config.R

cpg.450k <- readRDS(arg$path.cpg.annotation.450k);
cpg.850k <- readRDS(arg$path.cpg.annotation.850k);
cpg.850k <- cpg.850k[!cpg.850k$cpg %in% cpg.450k$cpg, ];
cpg.850k <- cpg.850k[!duplicated(cpg.850k$cpg), ];

cpg.450k$array <- '450k';
cpg.850k$array <- '850k';

colnames(cpg.450k)[which(colnames(cpg.450k) == 'chr')] <- 'chr.hg19';
colnames(cpg.450k)[which(colnames(cpg.450k) == 'pos')] <- 'pos.hg19';
colnames(cpg.850k)[which(colnames(cpg.850k) == 'chr')] <- 'chr.hg38';
colnames(cpg.850k)[which(colnames(cpg.850k) == 'pos')] <- 'pos.hg38';


cpg.450k.sub <- cpg.450k[, c('cpg', 'chr.hg19', 'pos.hg19', 'strand', 'cpg.type', 'promoter', 'array','UCSC_RefGene_Name')];
cpg.850k.sub <- cpg.850k[, c('cpg', 'chr.hg38', 'pos.hg38', 'strand', 'cpg.type', 'promoter', 'array','UCSC_RefGene_Name')];

levels(cpg.850k.sub$cpg.type)[levels(cpg.850k.sub$cpg.type) == 'OpenSea'] <- 'Open sea';
stopifnot(all(levels(cpg.450k.sub$cpg.type) == levels(cpg.850k.sub$cpg.type)));
for (i in 1:ncol(cpg.450k.sub)) {
    if (is.factor(cpg.450k.sub[, i])) {
        #print(i);
        stopifnot(all(levels(cpg.450k.sub[, i]) == levels(cpg.850k.sub[, i])));
        }
    }


########## convert 450k coordinates to hg38
# Create GRanges object from hg19 coordinates
gr.hg19 <- GRanges(
    seqnames = paste0('chr', as.character(cpg.450k.sub$chr)),
    ranges = IRanges(start = cpg.450k.sub$pos, end = cpg.450k.sub$pos),
    strand = cpg.450k.sub$strand,
    cpg = cpg.450k.sub$cpg
    );

# Download the liftover chain file
# path.file <- file.path(arg$path.save.annot, paste0(chain.file, '.gz'));
# download.file(
#     'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',
#     destfile = path.file
#     );
# gunzip(path.file, remove = FALSE);


# Load the chain file and perform liftOver
chain <- import.chain(arg$path.hg19.to.hg38.chain);
gr.hg38 <- liftOver(gr.hg19, chain)
gr.hg38 <- unlist(gr.hg38)

# Convert back to data frame and merge with original
lifted.df <- data.frame(
    cpg = gr.hg38$cpg,
    chr.hg38 = gsub('chr', '', as.character(seqnames(gr.hg38))),
    pos.hg38 = start(gr.hg38)
    );

# Merge hg38 coordinates back into original data
cpg.450k.hg38 <- merge(cpg.450k.sub, lifted.df, by = 'cpg', all.x = TRUE)
mean(is.na(cpg.450k.hg38$pos.hg38));



# pool 450k and 850k annotations
cpg.annotation <- data.frame(data.table::rbindlist(list(cpg.450k.hg38, cpg.850k.sub), fill = TRUE), check.names = FALSE)
cpg.annotation <- cpg.annotation[, c('cpg', 'chr.hg38', 'pos.hg38', 'chr.hg19', 'pos.hg19', 'cpg.type', 'promoter', 'array','UCSC_RefGene_Name')];
mean(cpg.annotation$chr.hg38 == cpg.annotation$chr.hg19, na.rm = TRUE);

saveRDS(
    cpg.annotation,
    file = file.path(arg$path.save.annot, paste0(Sys.Date(), '_cpg_annotation_450k-850k.rds')),
    );

data.table::fwrite(
    cpg.annotation,
    file = file.path(arg$path.save.annot, paste0(Sys.Date(), '_cpg_annotation_450k-850k.csv.gz')),
    row.names = FALSE,
    col.names = TRUE
    );
