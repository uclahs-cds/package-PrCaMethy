## PrCaMethy 0.2.0 (2025-05-20)

### New Features

- Added a random forest model for assigning the 4 methylation subtypes (`subtype.model.rf`).  Unlike `subtype.model.pamr` which requires all 5,486 subtype-defining CpGs to be measured, `subtype.model.rf` can handle missing CpGs through imputation (although ideally you should have as many of the CpGs as possible).

### Changed

- Renamed `subtype.model` to `subtype.model.pamr` since now the package has 2 models for assigning the methylation subtypes (`subtype.model.pamr` and `subtype.model.rf`).  This will not cause any breaking changes to the user.

## PrCaMethy 0.1.0 (2025-02-12)

* First release of the package.
