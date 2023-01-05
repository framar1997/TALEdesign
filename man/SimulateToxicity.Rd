\name{SimulateToxicity}
\alias{SimulateToxicity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SimulateToxicity
}
\description{
This function explores the explorable DLs and it updates the trial information
}
\usage{
SimulateToxicity(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{explorable_DL}{
the list of the DLs explorable at the current stage
}
  \item{toxicity_matrix}{
    matrix containing the toxicity probabilities for each DL
}
  \item{trial_data}{
   data.frame containg the information accrued during the trial
}
  \item{n_cohort}{
    cohort dimension
}
  \item{trial_stage}{
   current trial stage
}
\item{select_lambda}{
   If TRUE you must specify the thresholds values . If FALSE uses the values in Table 2 on the original paper.
}
\item{lambda_e}{
   Escalation threshold.
}
\item{lambda_d}{
  De-escalation threshold.
}
\item{lambda_r}{
  Removotion threshold.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
