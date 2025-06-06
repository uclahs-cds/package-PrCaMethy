% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate.gene.methy.data.R
\name{validate.gene.methy.data}
\alias{validate.gene.methy.data}
\title{Validate input data for gene.methylation()}
\usage{
validate.gene.methy.data(gene.methy.data, models, prop.missing.cutoff = 0.3)
}
\arguments{
\item{gene.methy.data}{A data frame with gene-level methylation data, created by \link{gene.methylation}.  Patients are rows and columns are genes.}

\item{models}{A list of models used to predict features from gene-level methylation data.  The models should come from \code{data('all.models')}.}

\item{prop.missing.cutoff}{The maximum proportion of missing values allowed for each required gene.}
}
\value{
\itemize{
\item \code{val.passed} a logical indicating whether the data passed validation
\item \code{features.you.can.predict} logical vector indicating which features you can predict (i.e. you have the required genes with missing data rates < prop.missing.cutoff)
\item \code{required.genes} a list of genes required by each model
\item \code{missing.genes} a list of genes that are required but completely missing in the data
\item \code{required.genes.with.high.missing} a list of genes that are required and have a proportion of missing values greater than \code{prop.missing.cutoff}
}
}
\description{
Check whether \code{gene.methy.data} contains all genes required by \code{models} and that there is an acceptable level of missingness for each required gene.
Note that genes with acceptable levels of missing values are later imputed using KNN imputation when calling \link{estimate.features}.
If you'd rather use a different imputation method, then make sure to impute missing values before calling \link{estimate.features}.
}
\examples{
data('all.models');

### example gene-level methylation data
data('example.data.gene.methy');
# note this dataset is derived from the following commands:
# data('example.data');
# example.data.gene.methy <- gene.methylation(example.data);
check <- validate.gene.methy.data(example.data.gene.methy, all.models);
stopifnot(check$val.passed);

# genes required to fit each model:
#check$required.genes;

# genes that are required but completely missing in your data:
#check$missing.genes;

# genes that are required and have a proportion of missing values greater than `prop.missing.cutoff`
#check$required.genes.with.high.missing;
}
