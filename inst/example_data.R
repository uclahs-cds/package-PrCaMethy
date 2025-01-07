path.data <- '/hot/project/disease/ProstateTumor/PRAD-000101-MethySubtypes/data/2024-07-16_pooled_tumour_normal_all_cpgs_all_cohorts.rds';

ex.data <- readRDS(path.data);
rownames(ex.data) <- ex.data$patient.id;
ex.data <- subset(ex.data, subset = cohort == 'TCGA', select = - c(patient.id, cohort));
check.miss <- apply(ex.data, 2, function(x) mean(is.na(x)));
ex.data <- ex.data[, which(check.miss < 0.3)];

set.seed(123);
ex.data <- ex.data[sample(1:nrow(ex.data), 5),];

format(object.size(ex.data), 'Mb');

usethis::use_data(ex.data, overwrite = TRUE);
