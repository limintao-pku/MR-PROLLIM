\name{data_p32data_me_mo_q_re}
\alias{data_p32data_me_mo_q_re}
\title{
  Transform data_p3 to data_me_mo_q_re
}
\description{
  \code{data_p32data_me_mo_q_re} is a function that transforms data_p3 to data_me_mo_q_re
  for \code{est_proc_bi}. The output of random-effects MR-PROLLIM with  
  \code{control_p3$inspect_data_p3 = T} can be transformed to data_me_mo_q_re to 
  accelerate the classical algorithms.
}

\usage{
data_p32data_me_mo_q_re(data_p3,length_all,p_cut=5e-8,adj_m="none",trace=F)
}
\arguments{
\item{data_p3}{
The output of random-effects MR-PROLLIM with \code{control_p3$inspect_data_p3 = T}.
}
\item{length_all}{
The number of individuals.
}
\item{p_cut}{
Cutoff Wald \emph{P} value. It is recommended to apply a stricter \emph{P} threshold
for fixed-effects MR-PROLLIM or classical algorithms than random-effects MR-PROLLIM.
}
\item{adj_m}{
Adjustment method that will be passed to \code{\link[stats]{p.adjust}}.
}
\item{trace}{
A logical argument indicating whether to track the progress.
}
}

\details{
Transform data_p3 to data_me_mo_q_re for \code{est_proc_bi}.
}

\value{
A list containing necessary elements that define the data_me_mo_q_re will be returned.
}

\references{
\url{https://github.com/limintao-pku/MR-PROLLIM}.
}

\note{
Note that \code{data_p3} should be output with \code{control_p3$inspect_data_p3 = T}. The
final results from \code{est_proc_bi} with \code{control_p3$inspect_data_p3 = F} do not
contain enough data for transformation. This function should not be used for
\code{est_proc_cont}.
}

\seealso{
\code{\link[MRprollim]{est_proc_bi}}.
}
\examples{
#see examples at https://github.com/limintao-pku/MR-PROLLIM/wiki.
}