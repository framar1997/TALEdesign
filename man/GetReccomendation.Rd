\name{GetReccomendation}
\alias{GetReccomendation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GetReccomendation
}
\description{
 This function returns the reccomendation/suggestion coming from a DL (E: is suggesting escalation; S: is suggesting to stay; D+R: is suggesting de-escalation and revisit; "D" is suggesting remotion)
}
\usage{
GetReccomendation(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{row_explore}{
 An explorable DL
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
