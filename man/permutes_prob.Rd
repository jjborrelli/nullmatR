\name{permutes_prob}
\alias{permutes_prob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
permutes_prob(mat, iter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mat, iter = 100) 
{
    cS <- colSums(mat)/nrow(mat)
    rS <- rowSums(mat)/ncol(mat)
    pmat <- matrix(nrow = nrow(mat), ncol = ncol(mat))
    for (i in 1:nrow(mat)) {
        for (j in 1:ncol(mat)) {
            pmat[i, j] <- sum((rS[i] + cS[j])/2)
        }
    }
    mat.list <- list()
    for (q in 1:iter) {
        mat2 <- matrix(nrow = nrow(mat), ncol = ncol(mat))
        for (i in 1:nrow(mat)) {
            for (j in 1:ncol(mat)) {
                mat2[i, j] <- rbinom(1, 1, prob = pmat[i, j])
            }
        }
        mat.list[[q]] <- mat2
    }
    return(permuted.matrices = mat.list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
