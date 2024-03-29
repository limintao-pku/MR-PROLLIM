\name{summary.mrp}
\alias{summary.mrp}
\title{
  Summaries
}
\description{
  \code{summary.mrp} is a function that produces result summaries for MR-PROLLIM.
}

\usage{
summary.mrp(est_out,sd,se_u1_method=c("hessian","sandwich","b1"),
            s1_cut_k=0.01,p12_hat_cut=0.5)
}
\arguments{
\item{est_out}{
Output of \code{est_proc_bi} or \code{est_proc_cont}. It should be an 'MR-PROLLIM
output' object.
}
\item{sd}{
The sample standard deviation (sd) of the continuous exposure. This value is required
to produce result summaries corresponding to 1-unit increase of the exposure.
\code{sd = 1} means 1-sd increase of the exposure. For output from \code{est_proc_bi},
this argument will be fixed at 1.
}
\item{se_u1_method}{
The method used to compute the standard error for u1. "hessian" indicates the inverse Hessian method.
"sandwich" indicates the sandwich method. "b1" means the same method as that used for b1. "hessian" is
chosen as the default.
}
\item{s1_cut_k}{
If the estimate of s1 for random-effects MR-PROLLIM is < \code{median(sqrt(est_out$data$m_sigma[,1]))*s1_cut_k}, a
message will be generated. Users are recommended to fit the intercept model or Egger model if the current model is
neither of them.
}
\item{p12_hat_cut}{
If the estimate of p1*p2 for random-effects MR-PROLLIM is > \code{p12_hat_cut}, a
message will be generated.
}
}

\details{
Produce result summaries for MR-PROLLIM.
}

\value{
A list or data frame will be returned depending on which MR-PROLLIM method is used.
}

\references{
\url{https://github.com/limintao-pku/MR-PROLLIM}.
}

\seealso{
\code{\link[MRprollim]{est_proc_bi}} and \code{\link[MRprollim]{est_proc_cont}}.
}
\examples{
#see examples at https://github.com/limintao-pku/MR-PROLLIM/wiki.
}