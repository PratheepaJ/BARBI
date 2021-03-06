% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps.R
\docType{data}
\name{ps}
\alias{ps}
\title{ZymoBIOMICS microbial community \code{'phyloseq'} object.}
\format{An object of class \code{'phyloseq'}.
\describe{
  \item{SubjectID}{IDs for each sample}
  \item{sampleID}{IDs for each sample}
  \item{SampleType}{A factor variable with two levels: Standard and Negative}
  \item{Names}{Extended names of the sampleIDs}
  }}
\source{
Valdiation data for the BARBI method (manuscript is under review)
}
\usage{
ps
}
\description{
A \code{'phyloseq'} class object with 195 taxa in 8 samples that are from seven rounds of six-fold dilutions (1:1 up to 1:279,936) from standard ZymoBIOMICS microbial community and 10 negative extraction controls.
}
\examples{
\dontrun{sample_data(ps)
otu_table(ps)}
}
\keyword{datasets}
