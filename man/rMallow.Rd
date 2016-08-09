\name{rMallow}
\alias{rMallow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampler for Mallows Model
}
\description{
Given a reference ranking and a dispersion parameter, sample a full ranking from Mallows model.
}
\usage{
rMallow(sigma, phi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
vector of the reference ranking, where i-th element is the item ranked in the i-th place by the reference ranking
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
a (0, 1] dispersion parameter indicating the "variance" of the sampled ranking
}
}
\details{
Mallows model describes a distribution over rankings with respect to a distance measure between one ranking and a reference ranking (Mallows, 1957). The sampler is based on the Repeated Insertion Model (Doignon et al., 2004).
}
\value{
Return a list with two components:
  \item{ranks}{a vector where the i-th element is the rank assigned to the i-th item.}
  \item{ranking}{a vector where the i-th element is the item ranked in the i-th place}
}

\references{
%% ~put references to the literature/web site here ~
Doignon, Jean-Paul, Aleksandar Pekec, and Michel Regenwetter. "The repeated insertion model for rankings: Missing link between two subset choice models." Psychometrika 69.1 (2004): 33-54.

Mallows, Colin L. "Non-null ranking models. I." Biometrika 44.1/2 (1957): 114-130.
}


\examples{
#a random reference ranking of length 10
#the first element ranked higher than the second element, etc.
sigma = sample(10)

#a random dispersion parameter
phi = runif(1, 0, 1)

res = rMallow(sigma, phi)
res$ranks
res$ranking
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
