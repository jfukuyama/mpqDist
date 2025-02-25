% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpqDist.R
\name{get_mpq_distances}
\alias{get_mpq_distances}
\title{Compute MPQr distances for a range of values of r}
\usage{
get_mpq_distances(X, tr, rvec = r_transform(0:100/100))
}
\arguments{
\item{X}{A \eqn{n \times p} data matrix.}

\item{tr}{A tree (of class \code{phylo}) with p leaves.}
}
\value{
A list containing distances between the rows of \code{X},
and the value of \eqn{r} that was used (\code{r}).
}
\description{
Compute MPQr distances for a range of values of r
}
\examples{
data(AntibioticSmall)
out.agpca = adaptivegpca(AntibioticSmall$X, AntibioticSmall$Q, k = 2)
}
