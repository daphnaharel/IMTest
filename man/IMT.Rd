% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMT.R
\name{IMT}
\alias{IMT}
\title{Runs information matrix test for an information matrix test GPCM model.}
\usage{
IMT(mod, constraint, R, ind_vec)
}
\arguments{
\item{mod}{An IMT GPCM model.}

\item{constraint}{Constraint is either "rasch" or "gpcm".}

\item{R}{number of iterations for simulation of the variance-covariance matrix.}

\item{ind_vec}{Vector of 0's and 1's for item-level parameters to be tested in the information matrix test.}
}
\value{
A list containing the information matrix test statistic and the associated degrees of freedom.
}
\description{
Runs information matrix test for an information matrix test GPCM model.
}
\examples{
data(dataset)
collapse = split(rep(c(1:4), 10), rep(1:10, each = 4))
my_data = collapse_data(dataset, collapse, "rasch")
model = gpcm_IMT(my_data$data, constraint = "rasch")
\donttest{
test_fit = IMT(model, "rasch", R = 5000, my_data$ind)
#This line of code takes longer than 10 seconds to run
pvalue = pchisq(test_fit$Tstat, test_fit$df, lower.tail = FALSE)
}
# See vignette("IMT-vignette") for more examples
}
