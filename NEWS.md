## PrCaMethy 1.0.1 (2025-06-06)

### Bug Fixes

- Previously `estimate.subtypes()` with `subtype.model = 'PAMR'` would error if any samples had >50% missing CpGs or any CpGs had >80% missing values (these are the default cutoffs used by `impute::impute.knn` for imputing missing values). This has been fixed such that now error handling for high missingness in CpGs is strictly handled by argument `prop.missing.cutoff` and `validate.subtype.model.cpgs()`. Further, if any CpGs or samples have >50% missing values, a warning will be printed to alert the user that the subtype assignment may be inaccurate.

## PrCaMethy 1.0.0 (2025-05-22)

### New Features

- Added a random forest model for assigning the 4 methylation subtypes (`subtype.model.rf`).  Unlike `subtype.model.pamr` which requires all 5,486 subtype-defining CpGs to be measured, `subtype.model.rf` can handle missing CpGs through imputation (although ideally you should have as many of the CpGs as possible).

### Changed

- Renamed `subtype.model` to `subtype.model.pamr` since now the package has 2 models for assigning the methylation subtypes (`subtype.model.pamr` and `subtype.model.rf`).  This will not cause any breaking changes to the user.

### Breaking changes

- `estimate.subtypes()` now returns a list with 2 elements: `subtypes` and `validation` where the latter checks the validity of the input methylation data.  Previously it only returned the subtypes data.frame.  This is a minor breaking change.

## PrCaMethy 0.2.0 (2025-05-15)

### New Features

* Support CpGs from Illumina 850K array when using `gene.methylation` to calculate gene-level methylation.  Previously only 450K array was supported.

### Enhancements

* Documentation updated to explain both Illumina 450K and 850K arrays are supported (previously only 450K was supported)

## PrCaMethy 0.1.0 (2025-02-12)

* First release of the package.
