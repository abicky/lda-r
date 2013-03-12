# This file is part of LDA-R.

Document <- setRefClass(
    "Document",
    fields = list(
      words  = "integer",  # term index of 'counts'
      counts = "integer",  # the frequency of each term
      length = "integer",  # the number of terms used in this document
      total  = "integer"   # the number of words (= N)
      ),
    methods = list(
      initialize = function(words, counts, length, total) {
          words <<- words
          counts <<- counts
          length <<- length
          total <<- total
      })
    )
Document$accessors(names(Document$fields()))


Corpus <- setRefClass(
    "Corpus",
    fields = list(
      docs     = "list",     # list of 'Document'
      numTerms = "integer",  # the number of terms registered in this corpus (= V)
      numDocs  = "integer"   # the number of documents (= M)
      ),
    methods = list(
      initialize = function(docs, numTerms, numDocs) {
          docs <<- docs
          numTerms <<- numTerms
          numDocs <<- numDocs
      },
      getMaxDocLength = function() {
          return(max(sapply(docs, function(doc) doc$getLength())))
      }
      )
    )
Corpus$accessors(names(Corpus$fields()))

Corpus$createFromFile <- function(filename) {
    printf("reading data from %s\n", filename)

    lines <- readLines(filename)
    nd <- length(lines)
    nw <- 1L
    docs <- vector("list", nd)
    for (i in seq_along(lines)) {
        items <- strsplit(lines[i], " ")[[1]][-1]
        total <- 0L
        len <- length(items)
        words <- integer(len)
        counts <- integer(len)
        for (j in seq_along(items)) {
            wordInfo <- strsplit(items[j], ":")[[1]]
            word <- as.integer(wordInfo[1]) + 1L  # R is 1-origin
            count <- as.integer(wordInfo[2])
            words[j] <- word
            counts[j] <- count
            total <- total + count
            if (word > nw) {
                nw <- word
            }
        }
        docs[[i]] <- Document$new(words, counts, len, total)
    }

    printf("number of docs    : %d\n", nd)
    printf("number of terms   : %d\n", nw)
    return(Corpus$new(docs, nw, nd))
}


LDAModel <- setRefClass(
    "LDAModel",
    fields = list(
      alpha     = "numeric",  # word と
      logProbW  = "matrix",   # probability of each term (K x V matrix)
      numTopics = "integer",  # the number of topic (= K)
      numTerms  = "integer"   # the number of terms (= V)
      ),
    methods = list(
      initialize = function(numTerms, numTopics) {
          numTerms <<- numTerms
          numTopics <<- numTopics
          alpha <<- 1
          logProbW <<- matrix(0, numTopics, numTerms)
      }
      )
    )
LDAModel$accessors(names(LDAModel$fields()))


SuffStats <- setRefClass(
    "SuffStats",
    fields = list(
      classWord      = "matrix",   # likelihood of topics for each term (K x V matrix)
      classTotal     = "numeric",  # the total of classWord per topic (K vector)
      alphaSuffStats = "numeric",  # eq. 8
      numDocs        = "integer"   # the number of documents (= M)
      ),
    methods = list(
      initialize = function(model) {
          K <- model$getNumTopics()
          V <- model$getNumTerms()
          classWord <<- matrix(0, K, V)
          classTotal <<- rep(0, K)
      },
      initializeParams = function(model) {
          K <- model$getNumTopics()
          V <- model$getNumTerms()
          classWord <<- matrix(1.0 / V + runif(K * V), K, V)
          classTotal <<- rowSums(classWord)
      },
      initializeWithCorpus = function(model, corpus, numInit) {
          K <- model$getNumTopics()
          D <- corpus$getNumDocs()
          docs <- corpus$getDocs()
          for (k in seq_len(K)) {
              for (i in seq_len(numInit)) {
                  d <- sample(D, 1)
                  printf("initialized with document %d\n", d)
                  doc = docs[[d]]
                  words <- doc$getWords()
                  counts <- doc$getCounts()
                  classWord[k, words] <<- classWord[k, words] + counts
              }
          }
          classWord <<- classWord + 1.0
          #classTotal <<- classTotal + rowSums(classWord)  # ??
          classTotal <<- rowSums(classWord)
      },
      reset = function(model) {
          K <- model$getNumTopics()
          V <- model$getNumTerms()
          classTotal <<- rep(0, K)
          classWord <<- matrix(0, K, V)
          numDocs <<- 0L
          alphaSuffStats <<- 0
      }
      )
    )
SuffStats$accessors(names(SuffStats$fields()))
