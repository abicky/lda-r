# This file is part of LDA-R.

printf <- function(format, ...) {
    cat(sprintf(format, ...))
}

# logsumexp
# cf. http://www.mail-archive.com/r-help@r-project.org/msg126202.html
logSum <- function(x) {
    xmax <- which.max(x)
    log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
}
