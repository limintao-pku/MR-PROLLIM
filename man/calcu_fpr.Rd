\name{calcu_fpr}
\alias{calcu_fpr}
\title{
  Calculate the False-Positive Rate of Two Sequential Tests
}
\description{
  \code{calcu_fpr} is a function that implements a parametric bootstrap algorithm to
  estimate the false-positive rate of two sequential tests with both summary data
  and individual-level data of zero-effect SNPs. The two tests may be correlated with
  each other. Currently, this function assumes Z-test for the stage 1 test and Wald
  test for the stage 2 test.
}

\usage{
calcu_fpr(p1,p2,summary_data,x,g,c,c_inherit=T,start=NULL,type=c("c","b"),
                    cd=T,max_unique=10,cd_g_code=T,n_min_limit=100,
                    control_limit_c=NULL,scale=T,
                    mc.cores=1,PSOCK=F,parallel_trace=F,dt=T,
                    nlm=T,
                    nlm_control=list(
                      gradtol=1e-8,steptol=1e-8,stepmax=5,terlim=100
                    ),
                    nlminb_control=list(
                      nlminb_control=list(scale=1,eval.max=300,iter.max=300),
                      nlm_control=list(gradtol=1e-10,stepmax=2,steptol=1e-10)
                    ),
                    rm_sig_p_cut=0.05,rm_sig_adj_m="fdr",cover=0.9999,
                    n_boot=10000,get_sig_snp_data=NULL)
}

\arguments{
\item{p1}{The cutoff \emph{P} value for stage 1 test. Vector is accepted.}
\item{p2}{The cutoff Wald \emph{P} value for stage 2 test. Vector is accepted.}
\item{summary_data}{
A data frame containing GWAS summary statistics for a group of zero-effect SNPs.
\code{summary_data} should include 3 named columns ('snp', 'eff', and 'se'), of which
'snp' should match exactly with \code{colnames(g)}, and 'eff' and 'se' will be used
to calculate \emph{Z} statistics and \emph{P} values. Note that users should ensure SNPs
in \code{summary_data} and \code{g} have been aligned (i.e., coded in the same way).
}
\item{x}{The exposure vector. The trait can be binary or continuous.}
\item{g}{
The SNP matrix. Typical genotypes can be coded with numeric or integer (recommended)
0, 1, and 2. NAs are accepted. Columns of the matrix represent different SNPs.
Column names are required and should be unique as SNP identifiers. Note that these
SNPs are expected to have no effect on the exposure \code{x}.
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
\item{rm_sig_p_cut}{
Cutoff \emph{P} value defining the statistical significance. If SNPs are significant,
they will be removed. 
}
\item{rm_sig_adj_m}{
Adjustment method that will be passed to \code{\link[stats]{p.adjust}}.
}
\item{cover}{
The coverage of the sampling region for estimating the conditional probability.
\code{1 - cover} should be small (e.g., < 0.001 * the conditional probability) but not zero.
}
\item{n_boot}{
Number of bootstrap repeats.
}
\item{get_sig_snp_data}{
If provided, the \code{get_sig_snp} procedure will be skipped.
}
}

\details{
Unlike traditional MR analyses, it is recommended to adopt a two-step SNP selection to prepare the input SNPs. This function may help users decide the cutoff
Wald \emph{P} value for \code{est_proc_bi} and \code{est_proc_cont}. See our paper
or \url{https://github.com/limintao-pku/MR-PROLLIM/wiki} for more details about the two-step SNP selection.
}

\value{
A list containing estimates of the false-positive rate and conditional probability will
be returned.
}

\note{
The input SNPs are expected to be mutually independent and have no effect on the exposure.
A sample of a few thousand such SNPs may be appropriate.
}

\references{
\url{https://github.com/limintao-pku/MR-PROLLIM}.
}

\seealso{
\code{\link[MRprollim]{get_sig_snp}}.
}
\examples{
#see examples at https://github.com/limintao-pku/MR-PROLLIM/wiki.
}