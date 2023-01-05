\name{ConductTrial}
\alias{ConductTrial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ConductTrial
}
\description{
This function simulate a single trial and returns a trial_data summary containing the DLs explored during the trial, the observed DLTs, the enrolled patients, etc.
}
\usage{
ConductTrial(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
    total enrolled patients
}
\item{n_cohort}{
    cohort dimension
}
\item{toxicity_matrix}{
    matrix containing the toxicity probabilities for each DL
}
\item{n_max}{
   maximum number of patients that can be assigned at a DL
}
\item{select_lambda}{
    TRUE by default. If FALSE uses the values in Table 2 on the original paper.
}
\item{lambda_e}{
   Escalation threshold. By default = 0.16
}
\item{lambda_d}{
  De-escalation threshold. By default = 0.24
}
\item{lambda_r}{
  Removotion threshold. By default = 0.33
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
tox_table =  rbind(c(0.12,	0.15,	0.22,	0.4,	0.43),
                   c(0.12,	0.15,	0.20,	0.35,	0.42),
                   c(0.1,	0.12,	0.18,	0.33,	0.40),
                   c(0.1,	0.12,	0.15,	0.25,	0.27),
                   c(0,	0.1,	0.1,	0.2,	0.22))
colnames(tox_table) = c(1:ncol(tox_table))
rownames(tox_table) = rev(c(1:nrow(tox_table)))
ConductTrial(N=15, n_cohort=3, toxicity_matrix = tox_table, n_max = 9)
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
