#' Models for predicting clinical and molecular features
#'
#' [Random forest](https://cran.r-project.org/web/packages/randomForest/index.html) and [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) models for predicting various clinical and molecular features of patients diagnosed with prostate cancer.
#' The predictors are gene-level methylation estimated by \link{gene.methylation}.
#' 
#' Models are available for predicting the following features:
#' 
#' * `age.continuous`: patient age in years
#' * `ISUP.grade`: International Society of Urological Pathology (ISUP) grade risk group (1-5). See [here](https://www.prostate.org.au/testing-and-diagnosis/grading-genetics/your-gleason-score-isup-grade/#:~:text=The%20ISUP%20Grade&text=It%20is%20done%20using%20prostate,the%20corresponding%20ISUP%20Grade%20Group.) for further details.
#' * `t.category`: TNM tumour category (1-4), measures the size and extent of the primary tumour
#' * `psa.categorical`: Prostate-specific antigen (PSA) category: <= 10, [10, 20), and >= 20 ng/mL.
#' * `pga`: percentage of the genome altered was defined as PGA = (base-pair length of all genome regions with gain or loss) / 3.2 billion bases \code{*} 100
#' * `<gene>.cna.<loss/gain>`: features with `.cna.` in their name give the gene name and then identify whether there is a copy number loss or gain event.  See the Examples section for the full list of cna features.
#' * `log2p1.snvs.per.mbps`: single nucleotide variants (SNVs) per mega-base pairs (Mbps) with a log2(x + 1) transformation.
#' @examples
#' data(all.models);
#' 
#' # Models for predicting the following features:
#' names(all.models);
#' 
#' # Model class per feature, e.g. randomForest or glmnet:
#' lapply(all.models, class);
#' 
#' # Required genes for predicting each feature:
#' # lapply(all.models, function(x) x$xNames)
'all.models'
