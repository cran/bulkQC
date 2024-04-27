\name{getOutliers}
\alias{getOutliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to obtain values outside of whiskers on boxplot
}
\description{
Within vector of continuous data, identifies and outputs those values sufficiently smaller than the first quartile (Q1) or larger than the third quartile (Q3).  Used in consort with \code{\link{ind_uni}} function to identify individual univariate outliers.
}
\usage{
getOutliers(x, mult = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of continuous data}
  \item{mult}{A multiplier on the interquartile range (IQR) to flag outliers that are mult*IQR less than Q1 or mult*IQR greater than Q3; default is 1.5}
}

\value{
Returns the subset of observations in x that are mult*IQR less than Q1 or mult*IQR greater than Q3
}
\references{
Tukey J. Exploratory Data Analysis. 1st edition. Reading, Mass: Pearson; 1977. 712 p.	
}

\seealso{\code{\link{ind_uni}}}
\examples{
ex = c(rnorm(95), -10, -8, 7, 9, 11)
getOutliers(ex)
getOutliers(ex, mult=3) 
}
