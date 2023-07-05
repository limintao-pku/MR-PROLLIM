\name{get_sig_snp}
\alias{get_sig_snp}
\title{
  Estimate the SNP-Exposure Effects and Select Significant SNPs
}
\description{
  \code{get_sig_snp} is a function that estimates the SNP-exposure effects and calculates
  Wald \emph{P} values with individual-level data. Parallel computation is supported.
}

\usage{
get_sig_snp(x,g,c,c_inherit=T,start=NULL,type=c("c","b"),bi_type=c("log","logit"),
                      p_cut=c(1e-3,1e-5,5e-8),return_dt=T,
                      cd=T,max_unique=10,cd_g_code=T,
                      trinary_only=T,n_min_limit=100,
                      control_limit_c=NULL,scale=T,
                      mc.cores=1,PSOCK=F,parallel_trace=F,dt=T,
                      nlm=T,
                      nlm_control=list(gradtol=1e-8,stepmax=5,steptol=1e-8,iterlim=300),
                      nlminb_control=list(
                        nlminb_control=list(scale=1,eval.max=300,iter.max=300),
                        nlm_control=list(gradtol=1e-10,stepmax=2,steptol=1e-10)
                      )
            )
}

\arguments{
\item{x}{The exposure vector. The trait can be binary or continuous.}
\item{g}{
The SNP matrix. Typical genotypes can be coded with numeric or integer (recommended)
0, 1, and 2. NAs are accepted. Columns of the matrix represent different SNPs.
Column names are required and should be unique as SNP identifiers.
}
\item{c}{
A list containing matrixes of the control variables. If all SNPs share the same
control variables, please set \code{c_inherit = T} and add one control variable
matrix into the list. It is also possible to assign different control variables
to different SNPs by providing a list with multiple control variable matrixes.
In this case, set \code{c_inherit = F} and ensure \code{length(c)==ncol(g)}. NULL
indicates no control variable.
}
\item{c_inherit}{
A logical argument indicating whether or not all SNPs share the same control variables.
}
\item{start}{
A list with \code{length(start)==ncol(g)} containing the starting point vectors for
all SNPs. This argument works only when \code{type = 'b'}.
}
\item{type}{
What is the type of the exposure.\cr'c' indicates continuous.\cr'b' indicates binary.
}
\item{bi_type}{
Use which model to parameterize the SNP-exposure association.\cr'log' indicates
log-linear model.\cr'logit' indicates logistic model.\cr Default is the log-linear model.
}
\item{p_cut}{
A vector or number indicating the cutoff Wald \emph{P} value(s).
}
\item{return_dt}{
A logical argument indicating whether to return detailed data that are used to
calculate the Wald \emph{P} values.
}
\item{cd}{
A logical argument (check data) indicating whether to check the input data.
}
\item{max_unique}{
MR-PROLLIM currently only accepts categorical variants as potential IVs. This
argument indicates the limit of the number of unique values for each input variant.
}
\item{cd_g_code}{
A logical argument indicating whether to check the coding of trinary SNPs. We
recommend to use 2 to represent two copies of the minor allele.
}
\item{trinary_only}{
A logical argument indicating whether to use trinary SNPs only. Both random-effects
and fixed-effects MR-PROLLIMs currently rely on trinary SNPs. 
}
\item{n_min_limit}{
The lower bound of the number of the least common genotypes.
}
\item{control_limit_c}{
A list containing arguments controlling the \code{limit_c} procedure.
See the Control Lists section of \code{\link[MRprollim]{est_proc_bi}} for details.
}
\item{scale}{
A logical argument indicating whether to standardize the continuous exposure and control variables. It
is recommended to set \code{scale = T} to avoid failures in optimizing procedures.
}
\item{mc.cores}{
Number of cores used for parallel computing. Parallel computation is implemented
with the \code{\link[parallel]{parallel}} package. By default, this function will
use forking on Unix-alike platforms and PSOCK clusters on Windows.
}
\item{PSOCK}{
A logical argument indicating whether to use PSOCK clusters for parallelization.
This argument only works on Unix-alike platforms.
}
\item{parallel_trace}{
A logical argument indicating whether to track the progress of the time-consuming steps.
}
\item{dt}{A logical argument indicating whether to print the detailed information.}
\item{nlm}{
A logical argument indicating which optimizer to be used for binary exposures. If
\code{nlm = T}, the \code{\link[stats]{nlm}} optimizer will be used. Otherwise, the
\code{\link[stats]{nlminb}} optimizer (combined with \code{\link[stats]{nlm}}) will be used. According to our experience,
\code{\link[stats]{nlm}} is faster than \code{\link[stats]{nlminb}}, but the former
may sometimes get stuck in regions producing abnormal probabilities (probabilities > 1)
due to oversized \code{stepmax}.
}
\item{nlm_control}{
See \code{\link[stats]{nlm}} for details.
}
\item{nlminb_control}{
\code{\link[stats]{nlminb}} will be used to derive initial estimates, and \code{\link[stats]{nlm}} will be used to improve the precision. See Note below.
}
}

\details{
This function can be used to examine the candidate SNPs. Those satisfying \code{p_cut}
can then be put into \code{est_proc_bi} or \code{est_proc_cont}. This step may save time
and memory if one exposure will be tested with multiple outcomes. This function can also
be used to see whether there are enough significant SNPs for a certain method.
}

\value{
A list containing effect estimates and Wald \emph{P} values will be returned.
}

\note{
The log-linear model has an innate limitation. It cannot guarantee the probabilities to
lie between 0 and 1 if control variables are included (especially those strongly
correlated with the exposure). This may result in NAs for a large number of SNPs. See \url{https://github.com/limintao-pku/MR-PROLLIM/wiki} for ways to handle this case.\cr\cr
This function with \code{nlm=F} uses two optimizers for log-linear regressions. According to our simulations, the \code{\link[stats]{nlminb}} optimizer is less likely to get stuck in regions that produce >1 probabilities
than the \code{\link[stats]{nlm}} optimizer. But it does not have a gradient control as good as \code{\link[stats]{nlm}}, thus leading to less precise estimates. By default, this function first uses
\code{\link[stats]{nlminb}} to derive initial estimates and refine them with \code{\link[stats]{nlm}} (the initial results are put into \code{\link[stats]{nlm}} as the starting point).
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