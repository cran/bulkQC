\name{ind_uni}
\alias{ind_uni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identifies individual univariate outliers}
\description{
Discovers potential individual univariate outliers by identifying and returning those observations outside of the whiskers on a boxplot
}
\usage{
ind_uni(d0, exclude = c("pid", "site"), n_uniq = 10, m = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d0}{A data frame with columns as variables and rows as observations}
  \item{exclude}{A vector of names of variables to exclude in outlier identification}
  \item{n_uniq}{Number of unique observations of a variable needed for outlier identification to be performed}
  \item{m}{A multiplier on the interquartile range (IQR) to flag outliers that are mult*IQR less than Q1 or mult*IQR greater than Q3; default is 1.5}
}

\value{
  \item{nID}{The number of observations evaluated}
  \item{nVar}{The number of variables evaluated}
  \item{data}{A data frame containing those observations deemed to be potential outliers that appends the outliers with the excluded variables to aid in interpretation}	
}

\references{
Tukey J. Exploratory Data Analysis. 1st edition. Reading, Mass: Pearson; 1977. 712 p. 
}

\seealso{\code{\link{getOutliers}}}

\examples{
data(iris)
iris2 = iris
iris2$pid = 1:dim(iris2)[1]
ind_uni(iris2, exclude=c("pid", "Species"), m=1.5)
ind_uni(iris2, exclude=c("pid", "Species"), m=3)
}

