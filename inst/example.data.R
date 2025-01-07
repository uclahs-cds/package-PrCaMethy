path.data <- '/hot/project/disease/ProstateTumor/PRAD-000101-MethySubtypes/data/2024-07-16_pooled_tumour_normal_all_cpgs_all_cohorts.rds';

example.data <- readRDS(path.data);
rownames(example.data) <- example.data$patient.id;
example.data <- subset(example.data, subset = cohort == 'TCGA', select = - c(patient.id, cohort));
check.miss <- apply(example.data, 2, function(x) mean(is.na(x)));
example.data <- example.data[, which(check.miss < 0.3)];

set.seed(123);
example.data <- example.data[sample(1:nrow(example.data), 5),];

format(object.size(example.data), 'Mb');

usethis::use_data(example.data, overwrite = TRUE);
