#' Random Forest model for predicting methylation subtype
#'
#' [randomForestSRC](https://cran.r-project.org/web/packages/randomForestSRC/index.html) model used for assigning new patients to four prostate cancer DNA methylation subtypes. Note [`subtype.model.pamr`] requires all 5,486 subtype-defining CpGs to be measured, whereas `subtype.model.rf` can handle missing CpGs through imputation (although ideally you should have as many of the CpGs as possible).
'subtype.model.rf'
