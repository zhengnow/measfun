% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tbaccu.R
\name{tbaccu}
\alias{tbaccu}
\title{Table the measures of accuracy}
\usage{
tbaccu(var, reference, thres = NULL, ci = TRUE, cluster = FALSE, id = NULL)
}
\arguments{
\item{var}{a vector, test variable}

\item{reference}{a binary reference variable}

\item{thres}{a numeric value, the threshold value when the test variable is continuous}

\item{ci}{TRUE/FALSE, whether to include the confidence interval in the output}

\item{cluster}{TRUE/FALSE, an indicator of whether the data is clustered. The data was assume to adopt simple random sampling with equal weights}

\item{id}{subject ID variable, when data is clustered}
}
\value{
a data frame for the estimates, with the 95\%CI (optional), and the counts
}
\description{
The function estimate the measures of accuracy (without adjustment)
for a test variable and a binary reference variable, independent or clustered data.
}
\details{
tbaccu outputs the table including estimates of sensitivity, specificity, PPV and NPV.
The test variable, when applicable, are dichotomized at the threshold value.
When the data is clustered, the survey package is used to estimate the confidence interval,
assumed simple random sampling.
}
\examples{
y <- rbinom(n=50, size=1, prob=0.3) # reference variable
x <- rnorm(n=50, mean=20, sd=5)
tbaccu(x,y, thres = 20)

}
