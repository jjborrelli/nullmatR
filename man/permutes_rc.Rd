\name{permutes_rc}
\alias{permutes_rc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
permutes_rc(mat, iter = 100)
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
    pattern1 <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
    pattern2 <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
    count <- 0
    mat.list <- list()
    while (count < iter) {
        srow <- sample(1:nrow(mat), 2)
        scol <- sample(1:ncol(mat), 2)
        test <- mat[srow, scol]
        if (sum(test == pattern1) == 4) {
            count <- count + 1
            mat[srow, scol] <- pattern2
            mat.list[[count]] <- mat
            next
        }
        else if (sum(test == pattern2) == 4) {
            count <- count + 1
            mat[srow, scol] <- pattern1
            mat.list[[count]] <- mat
            next
        }
        else {
            next
        }
    }
    matrices <- lapply(mat.list, as.matrix)
    return(permuted.matrices = matrices)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
