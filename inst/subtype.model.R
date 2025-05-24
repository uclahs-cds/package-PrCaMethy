# there are 2 models for assigning the 4 methylation subtypes to new samples: PAMR and random forest

source('config.R') # see project PRAD-000101-MethySubtypes/PrCaMethy/config.R


### PAMR
load(arg$path.subtype.model.pamr);

subtype.model.pamr <- model.for.predicting.subtypes;
required.cpgs <- rownames(subtype.model.pamr$centroids);

usethis::use_data(subtype.model.pamr, overwrite = TRUE, compress = 'xz');
usethis::use_package('pamr');

### Random Forest
load(arg$path.subtype.model.rf);
subtype.model.rf <- rf;
usethis::use_data(subtype.model.rf, overwrite = TRUE, compress = 'xz');
usethis::use_package('randomForestSRC');
