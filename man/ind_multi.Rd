\name{ind_multi}
\alias{ind_multi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identifies individual multivariate outliers}
\description{
Discovers potential individual multivariate outliers by identifying and returning those observations with outlier score greater than a threshold. The outlier score is calculated using single isolation forests.
}
\usage{
ind_multi(d0, exclude = c("pid", "site"), thresh = 0.7, n_uniq = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d0}{A data frame with columns as variables and rows as observations}
  \item{exclude}{A vector of names of variables to exclude in outlier identification}
  \item{thresh}{Threshold (0-1) that an outlier score must exceed to be flagged for further investigation}
  \item{n_uniq}{Number of unique observations of a variable needed for outlier identification to be performed}
}
\details{
The function evaluates multivariate observations from each row consisting of those variables not excluded by the 'exclude' argument above. For each multivariate observation, an outlier score is calculated using single isolation forests. Those multivariate observations that are isolated earliest in a decision tree have a lower tree depth, in turn have higher outlier scores, and are thought more likely to be outliers. 

}
\value{
  \item{nID}{The number of observations evaluated}
  \item{nVar}{The number of variables evaluated}
  \item{data}{A data frame containing those observations deemed to be potential outliers that appends the outliers with the excluded variables to aid in interpretation, and includes an outlier score for each row}
}

\references{
Cortes D. Explainable outlier detection through decision tree conditioning. arXiv:200100636 [cs, stat] [Internet]. 2020 Jan 2 [cited 2021 Nov 12]; Available from: http://arxiv.org/abs/2001.00636 
}

\examples{
data(iris)
iris2 = iris
iris2$pid = 1:dim(iris2)[1]
ind_multi(iris2, exclude=c("pid", "Species"), thresh=0.7, n_uniq=10)
ind_multi(iris2, exclude=c("pid", "Species"), thresh=0.6, n_uniq=10)
}

