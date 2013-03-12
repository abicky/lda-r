#-----------------------------------------------------------------------------
# LDA-R - pure R version of LDA-C (http://www.cs.princeton.edu/~blei/lda-c/)
#
# author    : abicky
# license   : GPL 2 (http://www.gnu.org/licenses/gpl-2.0.html)
# link      : https://github.com/abicky/mathjax_for_pukiwiki
#-----------------------------------------------------------------------------
source("./util.R")
source("./class.R")

LDA <- list()

LDA$estimate <- function(corpus, start = "random", INITIAL.ALPHA = 1, EM.CONVERGED = 1e-4,
                      VAR.CONVERGED = 1e-6, EM.MAX.ITER = 100, ESTIMATE.ALPHA = TRUE,
                      varMaxIter = 20, NUM.TOPICS = 3L) {
    if (class(corpus) == "character") {
        corpus <- Corpus$createFromFile(corpus)
    }

    # allocate variational parameters
    numDocs <- corpus$getNumDocs()
    numTerms <- corpus$getNumTerms()
    maxLength <- corpus$getMaxDocLength()
    gamma <- matrix(0, numDocs, NUM.TOPICS)

    # initialize model
    model = LDAModel$new(numTerms, NUM.TOPICS)
    ss <- SuffStats$new(model)
    if (start == "seeded") {
        ss$initializeWithCorpus(model, corpus)
        LDA$.mle(model, ss, FALSE)
        model$setAlpha(INITIAL.ALPHA)
    } else if (start == "random") {
        ss$initializeParams(model)
        LDA$.mle(model, ss, FALSE)
        model$setAlpha(INITIAL.ALPHA)
    }

    # run expectation maximization
    docs <- corpus$getDocs()
    oldLikelihood <- 0
    converged <- Inf
    i <- 0
    while ((converged < 0 || converged > EM.CONVERGED || i <= 2) && i <= EM.MAX.ITER) {
        i <- i + 1
        printf("**** em iteration %d ****\n", i)
        likelihood <- 0
        ss$reset(model)

        # e-step
        for (d in seq_len(numDocs)) {
            if ((d %% 1000) == 0) {
                printf("document %d\n", d)
            }
            ret <- LDA$.docEStep(docs[[d]], model, ss, varMaxIter, VAR.CONVERGED)
            likelihood <- likelihood + ret$likelihood
            gamma[d, ] <- ret$gamma
        }

        # m-step
        LDA$.mle(model, ss, ESTIMATE.ALPHA)

        # check for convergence
        converged <- (oldLikelihood - likelihood) / (oldLikelihood)
        if (converged < 0) {
            varMaxIter <- varMaxIter * 2
        }
        oldLikelihood <- likelihood

        # output model and likelihood
        printf("%10.10f\t%5.5e\n", likelihood, converged)
    }

    wordAssign <- vector("list", numDocs)
    for (d in seq_len(numDocs)) {
        if ((d %% 1000) == 0) {
            printf("final e step document %d\n", d)
        }
        ret <- LDA$.inference(docs[[d]], model, varMaxIter, VAR.CONVERGED)
        phi <- ret$phi
        wordAssign[[d]] <- LDA$.wordAssignment(docs[[d]], phi)
    }

    return(list(gamma = gamma, corpus = corpus, model = model, wordAssign = wordAssign))
}

LDA$inference <- function(model, corpus, VAR.MAX.ITER = 20, VAR.CONVERGED = 1e-6, NUM.TOPICS = 3L) {
    gamma <- matrix(0, numDocs, NUM.TOPICS)
    docs <- corpus$getDocs()
    for (d in seq_along(docs)) {
        if ((d %% 1000) == 0) {
            printf("document %d\n", d)
        }
        ret <- LDA$.inference(docs[[d]], model, VAR.MAX.ITER, VAR.CONVERGED)
        gamma[d, ] <- ret$gamma
    }

    return(list(gamma = gamma, phi = phi, corpus = corpus, model = model))
}

LDA$.optAlpha <- function(alphaSuffStats, D, K, MAX.ALPHA.ITER = 1000, NEWTON.THRESH = 1e-5) {
    initA <- 100
    logA <- log(initA)  # ??
    a <- initA

    iter <- 0
    df <- Inf
    while (abs(df) > NEWTON.THRESH && iter < MAX.ALPHA.ITER) {
        iter <- iter + 1
        a <- exp(logA)
        if (is.nan(a)) {  # a becomes NaN if a < 0
        #if (a < 0) {
            initA <- initA * 10
            printf("warning : alpha is nan; new init = %5.5f\n", initA)
            a <- initA
            logA <- log(a)  # ??
        }
        f <- LDA$.alhood(a, alphaSuffStats, D, K)   # likelihood related with alpha
        df <- LDA$.dAlhood(a, alphaSuffStats, D, K)
        d2f <- LDA$.d2Alhood(a, D, K)
        logA <- logA - df / (d2f * a + df)  # ??
        #a <- a - df / d2f  # Newton-Raphson method
        printf("alpha maximization : %5.5f   %5.5f\n", f, df)
        #printf("alpha maximization : %5.5f %5.5f   %5.5f %5.5f\n", a, f, df, d2f)
    }
    return(exp(logA))
    #return(a)
}

# cf. A.4.2 (alphaSuffStats = E_q[log(\theta)])
LDA$.alhood <- function(a, alphaSuffStats, D, K) {
    return(D * (lgamma(K * a) - K * lgamma(a)) + (a - 1) * alphaSuffStats)
}
LDA$.dAlhood <- function(a, alphaSuffStats, D, K) {
    return(D * (K * digamma(K * a) - K * digamma(a)) + alphaSuffStats)
}
LDA$.d2Alhood <- function(a, D, K) {
    return(D * (K * K * trigamma(K * a) - K * trigamma(a)))
}

LDA$.inference <- function(doc, model, VAR.MAX.ITER, VAR.CONVERGED) {
    K <- model$getNumTopics()
    V <- doc$getLength()
    N <- doc$getTotal()
    counts <- doc$getCounts()
    words <- doc$getWords()
    logProbW <- model$getLogProbW()

    converged <- Inf
    likelihood <- 0
    oldLikelihood <- 0

    # compute posterior dirichlet

    # initialization
    phi <- matrix(1.0 / K, V, K)
    gamma <- rep(model$getAlpha(), K) +  N / K  # N / K is same as colSums(counts * phi)
    varIter <- 0

    while (converged > VAR.CONVERGED && (varIter < VAR.MAX.ITER || VAR.MAX.ITER == -1)) {
	varIter <- varIter + 1

        phiSum <- 0
        oldPhi <- phi
        phi <- t(digamma(gamma) + logProbW[ , words])  # eq. 6
        phiSum <- apply(phi, 1, logSum)
        phi <- exp(phi - phiSum)  # normalization
        gamma <- gamma + colSums(counts * (phi - oldPhi))  # eq. 7

        likelihood <- LDA$.computeLikelihood(doc, model, phi, gamma)
        if (is.nan(likelihood)) {
            stop("likelihood is NaN")
        }
        converged <- (oldLikelihood - likelihood) / oldLikelihood
        oldLikelihood <- likelihood

        # printf("[LDA INF] %8.5f %1.3e\n", likelihood, converged)
    }
    return(list(likelihood = likelihood, gamma = gamma, phi = phi))
}

LDA$.computeLikelihood <- function(doc, model, phi, gamma) {
    K <- model$getNumTopics()
    words <- doc$getWords()
    counts <- doc$getCounts()
    logProbW <- model$getLogProbW()
    alpha <- model$getAlpha()
    dig <- digamma(gamma)
    gammaSum <- sum(gamma)
    digsum <- digamma(gammaSum)

    # eq. 15
    likelihood <- lgamma(alpha * K) - K * lgamma(alpha) - lgamma(gammaSum)
                + sum((alpha - 1) * (dig - digsum) + lgamma(gamma) - (gamma - 1) * (dig - digsum))
    l <- counts * phi * t((dig - digsum) - t(log(phi)) + logProbW[, words])
    likelihood <- likelihood + sum(l[is.finite(l)])

    return(likelihood)
}

LDA$.docEStep <- function(doc, model, ss, VAR.MAX.ITER, VAR.CONVERGED) {
    # posterior inference

    ret <- LDA$.inference(doc, model, VAR.MAX.ITER, VAR.CONVERGED)
    gamma <- ret$gamma
    phi <- ret$phi

    # update sufficient statistics
    alphaSuffStats <- ss$getAlphaSuffStats() + sum(digamma(gamma)) - model$getNumTopics() * digamma(sum(gamma))  # eq. 8
    words <- doc$getWords()
    counts <- doc$getCounts()
    classWord <- ss$getClassWord()
    classWord[, words] <- classWord[, words] + t(counts * phi)  # eq. 9

    ss$setAlphaSuffStats(alphaSuffStats)
    ss$setClassWord(classWord)
    ss$setClassTotal(rowSums(classWord))

    ss$setNumDocs(ss$getNumDocs() + 1L)

    return(ret)
}

LDA$.mle <- function(model, ss, ESTIMATE.ALPHA) {
    classWord <- ss$getClassWord()
    index <- classWord > 0
    logProvW <- matrix(0, nrow(classWord), ncol(classWord))
    logProvW[index] <- log(classWord) - log(ss$getClassTotal())
    logProvW[!index] <- -100
    model$setLogProbW(logProvW)

    if (ESTIMATE.ALPHA) {
        alpha <- LDA$.optAlpha(ss$getAlphaSuffStats(), ss$getNumDocs(), model$getNumTopics())
        model$setAlpha(alpha)
        printf("new alpha = %5.5f\n", model$getAlpha())
    }
}

LDA$.wordAssignment <- function(doc, phi) {
    words <- doc$getWords()
    return(structure(apply(phi, 1, which.max), names = words))
}
