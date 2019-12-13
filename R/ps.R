#' Simulated \code{'phyloseq'} object.
#'
#' A \code{'phyloseq'} class object with 195 taxa in 8 samples that are from seven rounds of six-fold dilutions (1:1 up to 1:279,936) from standard ZymoBIOMICS microbial community and 10 negative extraction controls.
#'
#' @format An object of class \code{'phyloseq'}.
#' \describe{
#'   \item{SubjectID}{IDs for each sample}
#'   \item{sampleID}{IDs for each sample}
#'   \item{SampleType}{A factor variable with two levels: Standard and Negative}
#'   \item{Names}{Extended names of the sampleIDs}
#'   }
#' @source A samples from \href{https://academic.oup.com/ofid/article/4/suppl_1/S71/4294264}{where}
#' @examples
#' \dontrun{sample_data(ps)
#' otu_table(ps)}
"ps"
