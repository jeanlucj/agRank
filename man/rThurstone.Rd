\name{rThurstone}
\alias{rThurstone}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampler for Thurstone Model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a vector of scores and a vector of variances, sample a full ranking from the Thurstone model. The length of the ranking is the same as the length of the score vector.
}
\usage{
rThurstone(S, Svar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
vector of scores
}
  \item{Svar}{
%%     ~~Describe \code{Svar} here~~
vector of variances corresponding to each item
}
}
\details{
Thurstone model describes a ranking process by sampling from normal distributions (Thurstone, 1927). An item i has a score s_i and an associated variance v_i. We sample a value independently for each item from N(s_i, v_i) and sort the resulting values in the decreasing order. This process gives rise to a ranking of all the items.
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
Thurstone, Louis L. "A law of comparative judgment." Psychological review 34.4 (1927): 273.
}


\examples{
#random score vector
S = runif(10, 1, 10)
#random variance vector
Svar = runif(10, 0, 5)

res = rThurstone(S, Svar)
res$ranks
res$ranking
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
