\name{data_p32data_p12}
\alias{data_p32data_p12}
\title{
  Transform data_p3 to data_p12
}
\description{
  \code{data_p32data_p12} is a function that transforms data_p3 to data_p12
  for \code{est_proc_bi}. The output of random-effects MR-PROLLIM with  
  \code{control_p3$inspect_data_p3 = T} can be transformed to data_p12 to 
  accelerate the fixed-effects MR-PROLLIM.
}

\usage{
data_p32data_p12(data_p3,p_cut=5e-8,adj_m="none")
}
\arguments{
\item{data_p3}{
The output of random-effects MR-PROLLIM with \code{control_p3$inspect_data_p3 = T}.
}
\item{p_cut}{
Cutoff Wald \emph{P} value. It is recommended to apply a stricter \emph{P} threshold
for fixed-effects MR-PROLLIM or classical algorithms than random-effects MR-PROLLIM.
}
\item{adj_m}{
Adjustment method that will be passed to \code{\link[stats]{p.adjust}}.
}
}

\details{
Transform data_p3 to data_p12 for \code{est_proc_bi}.
}

\value{
A list containing necessary elements that define the data_p12 will be returned.
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