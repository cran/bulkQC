\name{site_outliers}
\alias{site_outliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identifies site level outliers}
\description{
Discovers potential site level outliers by using unadjusted and adjusted regression models and standardized difference calculations. 
}
\usage{
site_outliers(d0, exclude = c("pid"), siteID = "site", covs = c("age"), threshG = 0.001, 
thresh2 = 0.05, threshS = 0.5, n_uniq = 10, n_dec = 4, n_decS = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d0}{A data frame with columns as variables and rows as observations}
  \item{exclude}{A vector of names of variables to exclude in outlier identification}
  \item{siteID}{The name of the variable in the data frame that identifies sites}
  \item{covs}{A vector of covariates to adjust for in the adjusted regression models}
  \item{threshG}{P-value threshold for global test equal means across sites}
  \item{thresh2}{P-value threshold for comparison of reference site vs. all other sites}
  \item{threshS}{Standardized difference threshold above which a site difference is deemed meaningfully large}
  \item{n_uniq}{Number of unique observations of a variable needed for outlier identification to be performed}
  \item{n_dec}{Number of decimals to display for p-values in output}
  \item{n_decS}{Number of decimals to display for standardized differences in output}
}
\details{
The function compares the distribution of a given variable across sites by first conducting a global test of equal means (without and with adjustment for covariates of interest). Among those variables where the null hypothesis of equal means across sites is rejected, the function then compares each site vs. all other sites using unadjusted and adjusted comparisons. The unadjusted comparisons include a two-sample t-test with equal variance and a standardized difference calculation. The adjusted comparisons include a linear regression model with an indicator variable for reference site and user-specified covariates, and an adjusted standardized difference calculated as the model coefficient for site divided by the model estimated root mean squared error.
}

\value{
  \item{overall}{A matrix with rows as variables where global test of equal means is rejected and columns as the corresponding p-values from the unadjusted and adjusted statistical tests}
  \item{sitewise_P}{For the variables identified by the global tests (columns), the unadjusted p-values (from two-sample t-test) comparing each site to all other sites (rows). Values above threshold printed as missing.} 
  \item{sitewise_P_adj}{For the variables identified by the global tests (columns), the adjusted p-values (from linear regression model) comparing each site to all other sites (rows). Values above threshold printed as missing.}  
  \item{sitewise_StDf}{For the variables identified by the global tests (columns), the unadjusted standardized differences comparing each site to all other sites (rows). Values below threshold printed as missing.}    
  \item{sitewise_StDf_adj}{For the variables identified by the global tests (columns), the adjusted standardized differences comparing each site to all other sites (rows). Values below threshold printed as missing.}     
}

\references{
Yang D, Dalton JE. A unified approach to measuring the effect size between two groups using SAS. 2012;6 
}

\examples{
data(iris)
iris2 = iris
iris2$temp = rnorm(dim(iris2)[1]) #for covariate adjustment
site_outliers(iris2, site="Species", covs=c("temp"))
}

