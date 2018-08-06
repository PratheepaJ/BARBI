#' Simulated \code{"phyloseq"} object.
#'
#' A \code{"phyloseq"} class object with 1638 taxa from 48 patient plasma, healthy plasma, and negative control samples.
#'
#' @format An object of class \code{"phyloseq"}.
#' \describe{
#'   \item{SubjectCode}{IDs of individuals}
#'   \item{Sample_Type}{A factor with three levels: Patient_Plasma, Control, and Healthy_Plasma}
#'   \item{Extraction_Number}{Factor to define blocks/batches}
#'   }
#' @source A subset of samples from \href{https://academic.oup.com/ofid/article/4/suppl_1/S71/4294264}{where}
#' @examples
#' \dontrun{sample_data(ps)
#' otu_table(ps)}
"ps"
