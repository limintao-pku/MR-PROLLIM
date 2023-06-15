\name{est_proc_bi}
\alias{est_proc_bi}
\title{
  Run MR-PROLLIM Estimating Procedures with a Binary Exposure
}
\description{
  \code{est_proc_bi} is a function that implements the designated MR-PROLLIM
  algorithms (see argument \code{est_type}). This function only accepts
  individual-level data measured in one sample. Parallel computation is supported.
}

\usage{
est_proc_bi(x,y,g,c=NULL,c_inherit=T,dum_loc_list="auto",
                      dt=T,mc.cores=1,PSOCK=F,parallel_trace=F,
                      est_type=c("p3","p2","p1","me_mo_q_re"),
                      start=NULL,
                      snp_exp_check=T,p_snp="ask",control_snp_exp_check=NULL,
                      cd=T,max_unique=10,cd_g_code=T,
                      n_min_limit=100,
                      control_limit_c=NULL,
                      control_p12=NULL,data_p12=NULL,
                      control_p3=NULL,data_p3_k=NULL,data_p3=NULL,data_p3_opt=NULL,
                      control_est_k_prior=NULL,control_est_k_post=NULL,
                      control_global_search=NULL,
                      control_me_mo_q_re=NULL,data_me_mo_q_re=NULL,...)
}
\arguments{
\item{x}{The exposure vector. The trait should be coded with numeric or integer 0 and 1.}
\item{y}{The outcome vector. The trait should be coded with numeric or integer 0 and 1.}
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
\item{dum_loc_list}{
An argument indicating the locations of dummy control variables. It can be "auto",
in which case MR-PROLLIM will automatically mark control variables with only two unique
values as dummy variables. Users may also provide a list containing numeric vectors
with \code{length(dum_loc_list)==length(c)}. This argument affects the estimates of 
iterative expectations (see our paper). Note that categorical control variables
should be manually transformed to dummy variables by users. Otherwise, they may
be treated as continuous variables.
}
\item{dt}{A logical argument indicating whether to print the detailed information.}
\item{mc.cores}{
Number of cores used for parallel computing. Parallel computation is implemented
with the \code{\link[parallel]{parallel}} package. By default, MR-PROLLIM will
use forking on Unix-alike platforms and PSOCK clusters on Windows.
}
\item{PSOCK}{
A logical argument indicating whether to use PSOCK clusters for parallelization.
This argument only works on Unix-alike platforms.
}
\item{parallel_trace}{
A logical argument indicating whether to track the progress of time-consuming
estimating steps.
}
\item{est_type}{
A character argument indicating which estimating procedure to run.\cr
'p3': random-effects MR-PROLLIM.\cr
'p2': robust fixed-effects MR-PROLLIM.\cr
'p1': ordinary fixed-effects MR-PROLLIM.\cr
'me_mo_q_re': extensions of classical MR methods.
}
\item{start}{
A list with \code{length(start)==ncol(g)} containing the starting point vectors for
all SNPs. This argument may be passed to different optimizing procedures. So, it is
safer to leave this argument as NULL and let MR-PROLLIM use the default starting points.
}
\item{snp_exp_check}{
A logical argument indicating whether to check the SNP-exposure effects according
to \code{p_snp}. Sometimes setting \code{snp_exp_check=F} may avoid duplicates if
users have provided SNPs that are known to satisfy \code{p_snp}.
}
\item{p_snp}{
The cutoff Wald \emph{P} value for SNP-exposure effects.
}
\item{control_snp_exp_check}{
A list containing additional arguments controlling the \code{snp_exp_check} procedure.
See more details below.
}
\item{cd}{
A logical argument (check data) indicating whether to check the input data.
}
\item{max_unique}{
MR-PROLLIM currently only accepts categorical variants as potential IVs. This
argument indicates the limit of the number of unique values for each input variant.
}
\item{cd_g_code}{
A logical argument indicating whether to check the coding of SNPs. We recommend
to use 2 to represent two copies of the minor allele.
}
\item{n_min_limit}{
The lower bound of the number of the least common genotype for each SNP.
}
\item{control_limit_c}{
A list containing arguments controlling the \code{limit_c} procedure.
See more details below.
}
\item{control_p12}{
A list containing arguments controlling the fixed-effects MR-PROLLIM procedures.
See more details below.
}
\item{data_p12}{
Data for fixed-effects MR-PROLLIM procedures. If provided, the procedure can be
substantially accelerated.
}
\item{control_p3}{
A list containing arguments controlling the random-effects MR-PROLLIM procedures.
See more details below.
}
\item{data_p3_k}{
Stage 1 data for random-effects MR-PROLLIM procedures.
}
\item{data_p3}{
Data for random-effects MR-PROLLIM procedures. If provided, the procedure can be
substantially accelerated.
}
\item{data_p3_opt}{
Data for random-effects MR-PROLLIM final optimization. If provided, the procedure can be
substantially accelerated. See \code{\link[MRprollim]{data_remp2data_p3_opt}}.
}
\item{control_est_k_prior}{
A list containing arguments controlling the \code{est_k_prior} procedure.
See more details below.
}
\item{control_est_k_post}{
A list containing arguments controlling the \code{est_k_post} procedure.
See more details below.
}
\item{control_global_search}{
A list containing arguments controlling the \code{global_search} procedure.
See more details below.
}
\item{control_me_mo_q_re}{
A list containing arguments controlling the classical algorithms.
See more details below.
}
\item{data_me_mo_q_re}{
Data for classical algorithms. If provided, the procedure can be
substantially accelerated.
}
\item{...}{
Other arguments that will be passed to \code{\link[rgenoud]{genoud}}. Note that
these arguments should be different from those in
\code{control_global_search$genoud_control}.
}
}

\details{Mendelian randomization (MR) is a useful statistical technique to infer
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
likelihood estimation (MLE) and requires individual-level data. See our paper and
\url{https://github.com/limintao-pku/MR-PROLLIM/wiki} for more details.}

\value{
Different methods will yield different outputs. In general, a list containing what
MR-PROLLIM has obtained will be returned. Users may be interested in the effect
estimates (in \code{eff}, \code{eff_m}, or \code{beta}) and the standard errors
(in \code{se}, \code{se_m}, or \code{vcov}). Users may use the
\code{\link[MRprollim]{summary.mrp}} function to view the results.
}

\note{
For MR-PROLLIM methods assuming double log-linear models (i.e., all methods in
\code{est_proc_bi}), it is possible that the point estimate or upper bound of the
confidence interval given by \code{\link[MRprollim]{summary.mrp}} may be positive
infinity (see our paper for reasons). To avoid this, users are recommended to code
the exposure inversely. This step will shift the causal effect
\code{beta1} to \code{-beta1}, and users can shift the point and interval estimates
backwards manually.\cr\cr
The log-linear model has an innate limitation. It cannot guarantee the probabilities to
lie between 0 and 1 if control variables are included (especially those strongly
correlated with the exposure or outcome). This may trigger large SNP removals through the
built-in quality control mechanisms of MR-PROLLIM. See \url{https://github.com/limintao-pku/MR-PROLLIM/wiki} for ways to handle this case.\cr\cr
MR-PROLLIM uses two optimizers for log-linear regressions. According to our simulations, the \code{\link[stats]{nlminb}} optimizer is less likely to get stuck in regions that produce >1 probabilities
than the \code{\link[stats]{nlm}} optimizer. But it does not have a gradient control as good as \code{\link[stats]{nlm}}, thus leading to less precise estimates. By default, MR-PROLLIM first uses
\code{\link[stats]{nlminb}} to derive initial estimates and refine them with \code{\link[stats]{nlm}} (the initial results are put into \code{\link[stats]{nlm}} as the starting point).\cr\cr
The log likelihood function for the full model of random-effects MR-PROLLIM is not globally concave. MR-PROLLIM uses \code{\link[rgenoud]{genoud}} to search for the global optimum. The parallel computation of \code{\link[rgenoud]{genoud}} seems to be inefficient if the codes are evaluated in non-global environment on
Linux platforms. MR-PROLLIM therefore runs \code{\link[rgenoud]{genoud}} codes in the global environment (the user's workspace) to speed up the procedure. To do this, MR-PROLLIM needs to create 5 objects into the global environment (.mrp_p3_cl, .mrp_p3_f1, .mrp_p3_f2, .mrp_p3_par, and .mrp_p3_fit). \emph{If the names have been used, please backup the data you need before running this function}. These objects will be removed when random-effects MR-PROLLIM finishes normally. If \code{mc.cores=1}, \code{\link[rgenoud]{genoud}} will not be run in the global environment.
}

\section{Control Lists}{
This section shows the default settings of all control lists in \code{est_proc_bi}. Users may use them to change the behavior of MR-PROLLIM or leave them as NULLs, in which case the default settings will work. Note that if you want to change a certain argument, you generally don't need to write other arguments into the corresponding control list (except for the nlm, nlminb, and GenSA control lists).\cr\cr

\preformatted{
control_snp_exp_check = list(
      use_nlm=F,start=NULL,
      nlm_control=list(gradtol=1e-8,steptol=1e-8,stepmax=5,iterlim=100),
      nlminb_control=list(nlminb_control=list(scale=1,eval.max=300,iter.max=300),
                          nlm_control=list(gradtol=1e-10,stepmax=2,
                                           steptol=10^c(-10,-12,-8)))
)
}
> \code{use_nlm}, \code{start}, \code{nlm_control}, and \code{nlminb_control}: See \code{\link[MRprollim]{get_sig_snp}} for details.\cr\cr

\preformatted{
control_limit_c = list(limit_c=T,dum_loc="auto",quantile=c(0,1),outlier=T)
}
> \code{limit_c}: Whether to modify the extreme values of the control variables.\cr
> \code{dum_loc}: See \code{dum_loc_list}. Control variables marked by \code{dum_loc} will not be modified.\cr
> \code{quantile}: Values outside the range will be set to the boundary values. c(0,1) means this argument will not work.\cr
> \code{outlier}: Whether to set outliers to boundary values. Outliers are defined as values > Q75+1.5IQR or < Q25-1.5IQR.\cr\cr

\preformatted{
control_p12 = list(k1_td_p=0.1,k1_td_adj_m="fdr",
                   delta_p_cut=0.01,delta_p_cut2=0.1,delta_p2_adj_m="fdr",
                   boot_n=10000,
                   p2_cut=0.05,p2_cut_adj_m="none",
                   p2_cut2=0.05,p2_cut2_adj_m="none",
                   p3_cut=0.05,p3_cut_adj_m="bonferroni",
                   n_random=NULL,n_max=2^16,
                   twosnp_p_cut=1e-5,p_prop_limit=1e-5,
                   nlminb_control=list(
                     nlminb_control=list(scale=1,eval.max=300,iter.max=300),
                     nlm_control=list(gradtol=1e-10,stepmax=2,steptol=10^c(-10,-12,-8))
                   ),
                   stage1_simplification=T,ind_hp=T,outlier_detect=T,
                   outlier_p=0.05,se_cut_k=1.5,
                   hp_p=0.1,adj_m="fdr",
                   c_type=c("linear","median&mode"),mode_phi=1)
}
> \code{k1_td_p}: Cutoff \emph{P} value for k1_td.\cr
> \code{delta_p_cut}: Simulated probability for delta < 0 should be higher than \code{delta_p_cut}.\cr
> \code{delta_p_cut2, delta_p2_adj_m}: Simulated probability for delta < 0 should additionally satisfy these  arguments after adjustment.\cr 
> \code{boot_n}: Number of repeats for calculating the simulated probabilities and the bootstrap standard errors in Procedure 1 and 2.\cr
> \code{p2_cut}: Cutoff \emph{P} value for removing SNPs judged to be with HP before the root selection procedure.\cr
> \code{p2_cut2}: Cutoff \emph{P} value for removing SNPs judged to be with HP after the root selection procedure.\cr
> \code{p3_cut}: Roots with beta1_td being significantly > 1 or exp(h) being significantly < 0 will be removed. Significance is defined according to \code{p3_cut} and \code{p3_cut_adj_m}.\cr
> \code{n_random}, \code{n_max}, and \code{twosnp_p_cut}: Arguments used in the root selection procedure. See our paper for their meanings.\cr
> \code{p_prop_limit}: If an estimated probability is close to 1 or 0, the normal approximation may be poor. If the asymptotic confidence intervals (one-sided) defined by \code{1 - p_prop_limit} include 1 or 0, this SNP will be removed.\cr
> \code{nlminb_control}: \code{\link[stats]{nlminb}} will be used to derive initial estimates, and \code{\link[stats]{nlm}} will be used to improve the precision. See Note below. \code{steptol} here accepts multiple numbers.\cr
> \code{stage1_simplification}: Whether to perform first-stage simplification.\cr
> \code{ind_hp}: Whether to conduct the indirect HP tests.\cr
> \code{outlier_detect}: Whether to remove outliers before \code{ind_hp}.\cr
> \code{outlier_p} and \code{se_cut_k}: Q statistic-based outlier removal. \code{outlier_p} defines significant heterogeneity. SNPs with standard errors > Q3+\code{se_cut_k}*IQR are also considered as outliers.\cr
> \code{hp_p}: Cutoff \emph{P} value for indirect HP tests.\cr
> \code{c_type}: Combination methods.\cr
> \code{mode_phi}: An argument for the weighted mode method.\cr\cr

\preformatted{
control_p3 = list(n_snp_limit=50,p_prop_limit=1e-5,
                  inspect_data_p3_k=F,inspect_data_p3=F,
                  inspect_data_p3_opt=F,nome=F,
                  s_filter=F,s_filter_ask=F,s_k_r_limit=c(0.2,0.2,0.2),
                  s_k_a_limit="auto",auto_k_limit=0.2,
                  beta_start=c(0,0,0,0,-1,2,0),auto_s1=T,ask=F,log_appr=0,
                  p1_sp=NULL,p2_sp=NULL,r_sp=NULL,model_u2=T,Egger="auto",
                  t_b1=F,model_select=T,
                  check_fit_upper=0.999,check_fit_lower=0.001,s1_cut_k=0.01,
                  sandwich=T,hessian=T,
                  boot_se=F,n_boot=3000,n_rep_max=3,
                  nlminb_control=list(
                    nlminb_control=list(rel.tol=1e-12,sing.tol=1e-12,step.min=1,
                                        eval.max=300,iter.max=300),
                    nlm_control=list(gradtol=1e-8,stepmax=2,steptol=1e-8)
                  )
                  )
}
> \code{n_snp_limit}: The minimum SNP number for random-effects MR-PROLLIM.\cr
> \code{p_prop_limit}: See above.\cr
> \code{inspect_data_p3_k}: Return stage 1 data directly.\cr
> \code{inspect_data_p3}: Return stage 1 and 2 data directly.\cr
> \code{inspect_data_p3_opt}: Return the data (stage 1, stage 2, and posterior sampling) before final optimizing.\cr
> \code{nome}: Whether the NOME (NO Measurement Error) assumption holds. If \code{nome=T}, \code{est_proc_bi_p3_f}, instead of \code{est_proc_bi_p3.1_f}, will be used. See our paper for more details.\cr
> \code{s_filter}, \code{s_filter_ask}, \code{s_k_r_limit}, \code{s_k_a_limit}, and \code{auto_k_limit}: These arguments control the procedure of removing SNPs that have large standard errors on k and p_td. By default, this procedure is disabled by \code{s_filter = F}.\cr
> \code{beta_start}: The starting point for the full model of random-effects MR-PROLLIM. If you are interested in what these values mean. See the \code{est_proc_bi_p3.1_f.R} at our website or \code{MRprollim:::est_proc_bi_p3.1_f}.\cr
> \code{auto_s1}: Whether to let MR-PROLLIM automatically sets the starting value for s1_t. This value equals \code{floor(log(var(m_hat[,1])))}. Note that \code{sqrt(exp(s1_t))=s1}.\cr
> \code{ask}: If the starting point is inappropriate, whether to ask the user to input a new one.\cr
> \code{log_appr}: Whether to use Taylor expansions to approximate logarithms in \code{est_proc_bi_p3.1_f}. 0 indicates False, 1 indicates first order, and 2 indicates second order. If \code{est_proc_bi_p3.1_f} outputs NAs many times and causes a failure in the optimization, use Taylor expansions may work. If the NAs do not casual a failure, users may ignore these NA warnings.\cr
> \code{p1_sp}, \code{p2_sp}, \code{r_sp}, and \code{model_u2}: These arguments control the model of random-effects MR-PROLLIM. 'sp' means 'specify'. If given, the corresponding parameters will be fixed. If \code{model_u2=F}, u2 in \code{est_proc_bi_p3.1_f} will be set to 0. See our paper for more details.\cr
> \code{Egger}: If p1_hat is close to 1 or 0, or s1_hat is too small, use the intercept model or the Egger model? If \code{Egger="auto"}, MR-PROLLIM will select the one with a higher log likelihood. \code{Egger} can also be True or False.\cr
> \code{t_b1}: Whether to estimate beta1_td directly.\cr
> \code{model_select}: Whether to conduct model selection according to \code{check_fit_upper}, \code{check_fit_lower}, and \code{s1_cut_k}. See our paper for more details.\cr
> \code{sandwich} and \code{hessian}: If \code{sandwich=T} and \code{hessian=T}, the sandwich variance estimate and inverse Hessian estimate will be compared, and the one with a larger standard error for beta1_td will be adopted. The sandwich method is robust to model misspecifications but tends to underestimate the variance in small samples (i.e., small SNP number).\cr
> \code{boot_se}: Whether to calculate the bootstrap variance estimate.\cr
> \code{n_rep_max}: A maximum of \code{n_boot * n_rep_max} repeats will be conducted.\cr
> \code{nlminb_control}: \code{\link[stats]{nlminb}} is a local optimizer that will be used to implement MLE for random-effects MR-PROLLIM models. \code{\link[stats]{nlm}} will be used to improve the precision of the final estimates (estimates given by nlminb or one global optimizer).\cr\cr

\preformatted{
control_est_k_prior = list(p_snp="previous",start=c(0,0,-1,-1,1),
                           p0_start=c(seq(0.1,0.9,by=0.2),0.99,0.01),auto_s=T,
                           auto_s_offset=1,p0_sp=NULL,p0_cut=1e-8,u1_sp=NULL,
                           nlminb_control=list(rel.tol=1e-10,sing.tol=1e-10,
                                               step.min=1,eval.max=300,
                                               iter.max=300))
}
> \code{p_snp}: Cutoff Wald \emph{P} value for \code{est_k_prior}. 'previous' means being consistent with \code{p_snp} of \code{est_proc_bi}. It can also be 'max', in which case the maximum Wald \emph{P} value will be used. It can also be a number.\cr
> \code{start}: Starting point for \code{est_k_prior}. See our paper for more details.\cr
> \code{p0_start}: Use each of the values as the initial p0 for \code{est_k_prior}. If \code{p0_start} is not NULL, values in \code{p0_start} will be transformed and inserted into \code{control_est_k_prior$start[1]} for the optimization.\cr 
> \code{auto_s} and \code{auto_s_offset}: Whether to let MR-PROLLIM automatically sets starting values for s01_t and s02_t. These values equal floor(log(apply(k_hat[,1:2],2,FUN=var)))+auto_s_offset. Note that sqrt(exp(s01_t))=s01 and sqrt(exp(s02_t))=s02.\cr
> \code{p0_sp} and \code{u1_sp}: Specify some parameters for the model of \code{est_k_prior}. See our paper for more details.\cr
> \code{p0_cut}: If the finally used \code{p_snp} < \code{p0_cut}, \code{p0_sp} will be set to 0.\cr
> \code{nlminb_control}: See \code{\link[stats]{nlminb}}.\cr\cr

\preformatted{
control_est_k_post = list(n_post=3000,p_cover=0.9999,f=NULL)
}
> \code{n_post}: Number of posterior samples for each SNP.\cr
> \code{p_cover}: The coverage of the sampling region. See our paper for more details.\cr
> \code{f}: User-defined function for generating the posterior samples. The first two arguments should accept kp_td and their vcov estimate, respectively.\cr\cr

\preformatted{
control_global_search = list(global_search=T,global_search_EI=T,gs_type="genoud",
                         genoud_control=list(pop.size=min(max(1000,mc.cores*300),10000),
                                             pop.size.EI=100,
                                             max.generations=60,
                                             wait.generations=10,
                                             hard.generation.limit=TRUE,
                                             solution.tolerance=1e-5,
                                             gradient.check=F,print.level=0,
                                             BFGSburnin=10,balance=F,
                                             optim_control=list()),
                         GenSA_control=list(maxit=3000,maxit.EI=300,
                                            max.time=NULL,verbose=F,
                                            trace.mat=F,seed=-100377),
                         lower=c(-10,-10,-10,-10,-30,-10,-5),
                         upper=c(10,10,10,10,10,10,5),auto_s1_k=5)
}
> \code{global_search}: Whether to conduct a global search for random-effects MR-PROLLIM models.\cr
> \code{global_search_EI}: Whether to conduct a global search for the Egger or intercept model.\cr
> \code{gs_type}: Use which method to conduct the global search. 'genoud' (recommended) or 'GenSA'.\cr
> \code{genoud_control}: See \code{\link[rgenoud]{genoud}} for details. \code{pop.size} largely determines the performance of this optimizer, and we generally do not recommend lower \code{pop.size} (\code{pop.size} should be set to keep a high probability that the global optimum is found). \code{pop.size.EI} is used only for the Egger or intercept model in place of \code{pop.size}. Note that MR-PROLLIM will automatically update some parameters if the current model is not the full model.\cr
> \code{GenSA_control}: See the \code{GenSA} package for details. If you want to use this method, be sure the package has been installed.\cr
> \code{lower} and \code{upper}: The parameter space of interest. Note that \code{\link[rgenoud]{genoud}} may go out of the space if a better point is found outside the space.\cr
> \code{auto_s1_k}: Determine the boundary values for s1_t if \code{control_p3$auto_s1=T}. \code{lower(s1_t)=beta_start(s1_t)-auto_s1_k} and \code{upper(s1_t)=beta_start(s1_t)+auto_s1_k}.\cr\cr

\preformatted{
control_me_mo_q_re = list(p_prop_limit=1e-5,mode_phi=1,outlier_p=0.05,
                          se_cut_k=NULL,boot_n=10000,inspect_data_mmqr=F,
                          nlminb_control=list(
                            nlminb_control=list(scale=1,eval.max=300,iter.max=300),
                            nlm_control=list(gradtol=1e-10,stepmax=2,steptol=10^c(-10,-12,-8))
                          )
                          )
}
> \code{p_prop_limit}: See above.\cr
> \code{mode_phi}: An argument for the weighted mode method.\cr
> \code{outlier_p}: \emph{P} value defining significant heterogeneity among individual effect estimators.\cr
> \code{se_cut_k}: See above. NULL means no SNP removal because of large standard errors will be conducted.\cr
> \code{boot_n}: Number of repeats for calculating the bootstrap variance estimates.\cr
> \code{inspect_data_mmqr}: Whether to return the data of individual effect estimators directly.\cr
> \code{nlminb_control}: \code{\link[stats]{nlminb}} will be used to derive initial estimates, and \code{\link[stats]{nlm}} will be used to improve the precision. See Note below. \code{steptol} here accepts multiple numbers.\cr\cr
}

\references{
\url{https://github.com/limintao-pku/MR-PROLLIM}.
}
\author{
Mintao Li, Peking University, \email{limintao@pku.edu.cn}.\cr
Jinzhu Jia, Peking University, \email{jzjia@pku.edu.cn}.\cr
}

\seealso{
\code{\link[MRprollim]{est_proc_cont}}.
}
\examples{
#see examples at https://github.com/limintao-pku/MR-PROLLIM/wiki.
}