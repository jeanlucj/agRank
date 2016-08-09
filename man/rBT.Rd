\name{rBT}
\alias{rBT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampler for Bradley-Terry Model
}
\description{
Given a vector of scores, sample a full ranking from the Bradley-Terry model. The length of the ranking is the
same as the length of the score vector.
}
\usage{
rBT(S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
vector of scores
}
}
\details{
Given two items i and j with scores s_i and s_j respectively, according to Bradley Terry model (Bradley and Terry, 1952), the probability of i ranked higher than j is: e^\{s_i\} / (e^\{s_i\} + e^\{s_j\}). Assuming independence among all pairwise comparisons gives rise to the sampler.


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
Bradley, Ralph Allan, and Milton E. Terry. "Rank analysis of incomplete block designs: I. The method of paired comparisons." Biometrika 39.3/4 (1952): 324-345.
}

\examples{
#random scores with length 10
S = runif(10, 1, 10)
res = rBT(S)

res$ranks
res$ranking
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
