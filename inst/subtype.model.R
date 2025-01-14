path.subtype.model <- '/hot/project/disease/ProstateTumor/PRAD-000101-MethySubtypes/data/2024-11-13_example-data_for-predicting-subtypes.Rdata';

load(path.subtype.model);

subtype.model <- model.for.predicting.subtypes;
subtype.model.required.cpgs <- rownames(subtype.model$centroids);

usethis::use_data(subtype.model, overwrite = TRUE, compress = 'xz');
usethis::use_package('pamr');

