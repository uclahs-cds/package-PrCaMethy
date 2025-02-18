source('config.R') # see project PRAD-000101-MethySubtypes/PrCaMethy/config.R
load(arg$path.subtype.model);

subtype.model <- model.for.predicting.subtypes;
subtype.model.required.cpgs <- rownames(subtype.model$centroids);

usethis::use_data(subtype.model, overwrite = TRUE, compress = 'xz');
usethis::use_package('pamr');
