\name{plot.mrp.p3}
\alias{plot.mrp.p3}
\title{
  Obtain the Classification Plots for Random-Effects MR-PROLLIM
}
\description{
  \code{plot.mrp.p3} is a function that draws classification plots for random-effects MR-PROLLIM.
}

\usage{
plot.mrp.p3(est_out,ci_cover=0.95,use_trans=T,trans_lower=0,reorder=NULL,
            emph_loc=NULL,emph_col="purple",
            cluster_type=3,cluster_color=c("black","blue","red"),bg_color="white",
            control_plot1=list(),control_plot2=list(),expr1=NULL,expr2=NULL,
            control_lines=list(),control_points=list(),control_abline=list(),
            sub_fig_label=NULL,
            sub_fig_line=3,sub_fig_k=0.08,sub_fig_cex=1,
            control_legend=NULL,
            interactive=F,obj_name=NULL)
}
\arguments{
\item{est_out}{
Final random-effects MR-PROLLIM (Procedure 3) output from \code{est_proc_bi} or \code{est_proc_cont}.
It should be an 'MR-PROLLIM output' object.
}
\item{ci_cover}{
A numeric value denoting the coverage of the two-sided confidence intervals. If \code{ci_cover=NULL},
no confidence interval will be drawn. 
}
\item{use_trans}{
The classification plots use the intensity of colors to indicate the estimated posterior probability (EPP).
For each cluster, a linear transformation is applied to EPPs if \code{use_trans=T}, with the lowest EPP set
to \code{trans_lower} and the highest set to 1. Transformed EPPs will be input to
\code{\link[grDevices]{colorRamp}} functions. If \code{use_trans=F}, no transformation will be conducted.
}
\item{trans_lower}{
See argument \code{use_trans}.
}
\item{reorder}{
A numeric or integer vector with \code{length(reorder)} = the number of SNPs. It will function in this way:
\code{plot(x=x[reorder],y=y[reorder])}. This argument may be useful if users want to emphasize certain points
that are covered by other points.
}
\item{emph_loc,emph_col}{
These two arguments can be used if users want to emphasize certain points. \code{emph_loc} accepts a numeric or integer vector which contains indexes of the SNPs. \code{emph_col} controls the color of the emphasized points.
}
\item{cluster_type}{
A numeric value in c(2,3). If \code{cluster_type=2}, SNPs judged to have uncorrelated or correlated horizontal
pleiotropy (HP) will be merged into one cluster.
}
\item{cluster_color}{
A character vector with 3 elements, corresponding to the colors of SNPs without HP, with uncorrelated HP, and
with correlated HP, respectively.
}
\item{bg_color}{
\code{c(bg_color,cluster_color[i])} will be input to \code{\link[grDevices]{colorRamp}} to generate color
functions. See argument \code{use_trans}.
}
\item{control_plot1,control_plot2}{
Control lists containing arguments that control the classification plots. \code{plot.mrp.p3} first creates a
canvas with only coordinate axes (using \code{\link[base]{plot}}); then draw confidence intervals with
\code{\link[graphics]{lines}} if \code{m_ci_cover} is not NULL; finally draw the points with
\code{\link[graphics]{points}}. Note that \code{plot.mrp.p3} will create two subfigures, as MR-PROLLIM allows
non-additive SNP effects. Additional arguments for \code{\link[base]{plot}} can be passed via
\code{control_plot1} and \code{control_plot2}. They will be used in the first plotting stage.
}
\item{expr1,expr2}{
Expressions that will be evaluated after \code{plot1} and \code{plot2}.
}
\item{control_lines,control_points}{
Control lists that will be used in the second and third plotting stage, respectively. See arguments
\code{control_plot1} and \code{control_plot2}. Both subfigures will share the control lists. 
}
\item{control_abline}{
A control list that will be used to draw line y=x for both subfigures.
}
\item{sub_fig_label}{
Subfigure labels.
}
\item{sub_fig_line,sub_fig_cex}{
Arguments that will be used in \code{\link[graphics]{mtext}} to print subfigure labels.
}
\item{sub_fig_k}{
This argument controls the ordinates of subfigure labels.
}
\item{control_legend}{
A list containing arguments that control the legend. See below.
}
\item{interactive}{
Whether to draw interactive plots. The interactive plots depend on the 'ggplot2' and 'plotly' package. Please ensure that these packages have been installed if you want to set \code{interactive=T}. The interactive plots may help you quickly locate which SNPs are highly likely to have pleiotropic effects or which SNPs are relatively ambiguous.
}
\item{obj_name}{
An object name for saving the plotting data.
}
}

\details{
Two classification subfigures based on maximum estimated posterior probability will be produced.
}

\section{Control Lists}{
\preformatted{
control_legend=list(draw_legend=T,mar=c(2,0,2,2),x_extend=0.05,width=0.4,
                    border_col="gray",border_lwd=1,
                    scale_cex=0.8,scale_offset=0.2,
                    cluster_cex=0.9,c_x_adj=0,c_y_adj=-2,cluster_center=T,
                    point_cex=2,point_lwd=2,p_y_adj=0,
                    title=NULL,title_cex=1,t_x_adj=0,t_y_adj=-1)
}
We will add more information here soon.
}

\value{
A list containing the plotting data will be returned if \code{obj_name} is provided.
}

\references{
\url{https://github.com/limintao-pku/MR-PROLLIM}.
}

\seealso{
\code{\link[MRprollim]{summary.mrp}}.
}
\examples{
#see examples at https://github.com/limintao-pku/MR-PROLLIM/wiki.
}