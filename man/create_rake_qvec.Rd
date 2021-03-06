% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multical.R
\name{create_rake_qvec}
\alias{create_rake_qvec}
\title{Create linear term vector in QP}
\usage{
create_rake_qvec(D, sample_counts, target_counts, lambda)
}
\arguments{
\item{D}{Design matrix, output from \code{\link{create_design_matrix}}}

\item{sample_counts}{Vector of sample counts for each cell}

\item{target_counts}{Vector of target counts for each cell}
}
\description{
Create linear term vector in QP
}
\keyword{internal}
