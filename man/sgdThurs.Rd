\name{sgdThurs}
\alias{sgdThurs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate MAP estimators for Thurstone Model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the MAP estimators for the Thurstone model with adherence (or reliability) parameters (Raman and Thorsten, 2014). The prior on scores is a multivariate normal distribution, and the prior on adherence is an improper, flat prior (on the whole real line).
}
\usage{
sgdThurs(data, mu, sigma, rate, maxiter = 1000, tol = 1e-09, start, decay)
}
%- maybe also 'usage' for other objects documented here.
\arguments{   \item{data}{
%%     ~~Describe \code{data} here~~
a n * m matrix, where n is the number of observers and m is the number of items to rank; each row vector is a partial or complete ranking, with i-th element being the rank assined to item i; the entry where that item is not ranked in the partial ranking is replaced by 0
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
the mean vector of the normal prior on scores
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
the covariance matrix of the normal prior on scores
}
  \item{rate}{
%%     ~~Describe \code{rate} here~~
the learning rate/stepsize
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
the maximal number of iterations
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
if the L-2 norm of the change of target function values is less than tol, we assert convergence
}
  \item{start}{
%%     ~~Describe \code{start} here~~
a length m + n vector of the starting value of parameters; the first m elements are the starting values for scores, while the next n elements are the starting values for adherences of observers
}
  \item{decay}{
%%     ~~Describe \code{decay} here~~
the learning rate is divided by decay if the target function value doesn't decrease; this number is set to be 1.1 by default
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
The original Thurstone model is modified to accommodate an extra adherence parameter, which accounts for each observer's degree of disagreement with the real ranking (Ramam and Thorsten, 2014). We assign a normal prior on the score parameter in order to incorporate the information from the additive relationship matrix. The MAP estimators for scores and adherences are calculated using stochastic gradient descent algorithm.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Return a list with 4 components
\item{value}{a vector which stores the value of the target function to optimize over all iterations}
\item{niter}{number of iterations to reach convergence}
\item{score}{the MAP estimator for scores}
\item{adherence}{the MAP estimator for adherences}
}

\references{
%% ~put references to the literature/web site here ~
Raman, Karthik, and Thorsten Joachims. "Methods for ordinal peer grading." Proceedings of the 20th ACM SIGKDD international conference on Knowledge discovery and data mining. ACM, 2014.
}

\examples{
#synthetic data with 3 items to rank and 2 observers
#the first observer ranks item 1 \succ item 2
#the second observer ranks item 1 \succ item 3 \succ item 2
data = rbind(c(1, 2, 0), c(1, 3, 2))

#set up the prior
mu = c(1, 1, 1)
sigma = diag(1, 3)

#set up the starting values
start = rep(1, 5)

res = sgdThurs(data, mu, sigma, 0.1, 1000, 1e-8, start, 1.1)

res$value
res$niter
res$score
res$adherence


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
