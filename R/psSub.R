#' Simulated \code{"phyloseq"} object.
#'
#' A \code{"phyloseq"} class object with 1644 taxa from 268 SIRS, healthy, and negative control samples.
#'
#' @format An object of class \code{"phyloseq"}.
#' \describe{
#'   \item{SubjectID}{IDs of individuals}
#'   \item{SampleType}{A factor with three levels: SIRS, Control, and Healthy}
#'   \item{Sample_Type}{A factor with two levels: Plasma, Control}
#'   \item{SampleCode}{Sample ID for sample}
#'   }
#' @source A subset of samples from \href{https://academic.oup.com/ofid/article/4/suppl_1/S71/4294264}{where}
#' @examples
#' \dontrun{sample_data(psSub)
#' otu_table(psSub)}
"psSub"
