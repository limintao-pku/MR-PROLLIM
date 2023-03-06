\name{data_remp2data_p3_opt}
\alias{data_remp2data_p3_opt}
\title{
  Transform Random-Effects MR-PROLLIM Output to data_p3_opt
}
\description{
  \code{data_remp2data_p3_opt} is a function that transforms random-effects
  MR-PROLLIM (Procedure 3) output to data_p3_opt for \code{est_proc_bi} and
  \code{est_proc_cont}. This function may be useful if users want to re-run
  the final optimization with changes of some random-effects MR-PROLLIM arguments.
}

\usage{
data_remp2data_p3_opt(data_remp)
}
\arguments{
\item{data_remp}{
The final output of random-effects MR-PROLLIM.
}
}

\details{
Transform random-effects MR-PROLLIM (Procedure 3) output to data_p3_opt for 
\code{est_proc_bi} and \code{est_proc_cont}.
}

\value{
A list containing necessary elements that define the data_p3_opt will be returned.
}

\references{
\url{https://github.com/limintao-pku/MR-PROLLIM}.
}

\note{
The function will automatically judge where \code{data_remp} is derived from and
whether or not \code{control_p3$nome} is set to True. 
}

\seealso{
\code{\link[MRprollim]{est_proc_bi}} and \code{\link[MRprollim]{est_proc_cont}}.
}
\examples{
#see examples at https://github.com/limintao-pku/MR-PROLLIM/wiki.
}