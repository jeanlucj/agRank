\name{rPL}
\alias{rPL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampler for Plackett-Luce Model
}
\description{
Given a vector of scores and an arbitrary positive number (the scale parameter in Gumbell distribution), sample a full ranking from the Plackett-Luce model. The length of the ranking is the same as the length of the score vector.
}
\usage{
rPL(S, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
vector of scores
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
an arbitrary positive number (the scale parameter in Gumbell distribution)
}
}
\details{
The Plackett-Luce model (Luce, 1959; Plackett, 1975) is the generalization of Bradley-Terry model which describes a distribution over full rankings, instead of pairwise comparisons. Let s_i be the score for the i-th item. It has been proved that if each item draws independent values from a Gumbell distribution with location parameter s_i * beta and scale parameter beta, and the drawn values are sorted in the decreasing order, the ranking of those values gives rise to a ranking of items consistent to the Plackett-Luce model (Yellott, 1977).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Return a list with two components:
  \item{ranks}{a vector where the i-th element is the rank assigned to the i-th item.}
  \item{ranking}{a vector where the i-th element is the item ranked in the i-th place}
}
\references{
%% ~put references to the literature/web site here ~
R.. Ducan Luce. Individual Choice Behavior a Theoretical Analysis. John Wiley and sons, 1959.

Plackett, Robin L. "The analysis of permutations." Applied Statistics (1975): 193-202.

Yellott, John I. "The relationship between Luce's choice axiom, Thurstone's theory of comparative judgment, and the double exponential distribution." Journal of Mathematical Psychology 15.2 (1977): 109-144.
}


\examples{
#random score vector of length 10
S = runif(10, 1, 10)
#beta is chosen to be 2 by defause
beta = 2

res = rPL(S, beta)
res$ranks
res$ranking
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
