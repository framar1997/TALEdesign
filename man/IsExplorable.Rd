\name{IsExplorable}
\alias{IsExplorable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
IsExplorable
}
\description{
 returns TRUE if a DL is explorable in the next trial stage, FALSE otherwise
}
\usage{
IsExplorable(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DL}{
  a DL
}
  \item{trial_data}{
  data.frame containg the information accrued during the trial
}
  \item{all_DLs}{
  data.frame containing all the DLs of the trial
}
\item{N}{
    total enrolled patients
}
\item{n_cohort}{
    cohort dimension
}
  \item{trial_stage}{
  current stage of the trial
}
 \item{nrow}{
 amount of drug B doses
}
 \item{ncol}{
 amount of drug A doses
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
