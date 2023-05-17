\name{MRprollim-package}
\alias{MRprollim-package}
\alias{MRprollim}
\docType{package}
\title{
  MR-PROLLIM: Mendelian Randomization with Pleiotropy-Robust Log-Linear Model
}
\description{
  This package implements the MR-PROLLIM algorithms.
}
\details{
Mendelian randomization (MR) is a useful statistical technique to infer
causality between traits. In dealing with a binary outcome, there are two challenging
barriers on the way towards a valid MR analysis. The first one is the inconsistency of
the traditional ratio estimator. Briefly, this estimator is consistent for neither
the conditional logarithmic odds ratio (OR) defined by a logistic regression nor the
population-averaged one given by an ideal randomized controlled trial (RCT). The second
barrier, which may cause substantial distortions to MR results, is the potential
horizontal pleiotropy of SNPs. In particular, the correlated pleiotropy is more harmful
and more difficult to handle. MR-PROLLIM aims to surmount both the barriers. It adopts
Pleiotropy-RObust Log-LInear Model (PROLLIM) to parameterize the associations between
traits and uses risk ratio (RR) to describe the causality. MR-PROLLIM relies on maximum
likelihood estimation (MLE) and requires individual-level data. See our paper
\url{https://github.com/limintao-pku/MR-PROLLIM} for more details.
}
\author{
Mintao Li, Peking University, \email{limintao@pku.edu.cn}.\cr
Jinzhu Jia, Peking University, \email{jzjia@pku.edu.cn}.\cr
Maintainer: Mintao Li.
}
\references{
\url{https://github.com/limintao-pku/MR-PROLLIM}
}
\keyword{ package }
\seealso{
\code{\link[MRprollim]{est_proc_cont}} and \code{\link[MRprollim]{est_proc_bi}}.
}
\examples{
#see examples at https://github.com/limintao-pku/MR-PROLLIM/wiki.
}
