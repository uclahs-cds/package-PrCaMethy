## PrCaMethy 0.3.0 (2025-05-22)

### New Features

- Added a random forest model for assigning the 4 methylation subtypes (`subtype.model.rf`).  Unlike `subtype.model.pamr` which requires all 5,486 subtype-defining CpGs to be measured, `subtype.model.rf` can handle missing CpGs through imputation (although ideally you should have as many of the CpGs as possible).

### Changed

- Renamed `subtype.model` to `subtype.model.pamr` since now the package has 2 models for assigning the methylation subtypes (`subtype.model.pamr` and `subtype.model.rf`).  This will not cause any breaking changes to the user.

## PrCaMethy 0.2.0 (2025-05-15)

### New Features

* Support CpGs from Illumina 850K array when using `gene.methylation` to calculate gene-level methylation.  Previously only 450K array was supported.

### Enhancements

* Documentation updated to explain both Illumina 450K and 850K arrays are supported (previously only 450K was supported)


## PrCaMethy 0.1.0 (2025-02-12)

* First release of the package.
