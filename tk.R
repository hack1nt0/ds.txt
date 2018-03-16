tk <- function(texts,
               type = "mix",
               dict = jiebaR::DICTPATH,
               hmm = jiebaR::HMMPATH,
               user = jiebaR::USERPATH,
               idf = jiebaR::IDFPATH,
               stop_word = jiebaR::STOPPATH,
               write = T,
               qmax = 20,
               topn = 5,
               encoding = "UTF-8",
               detect = F,
               symbol = T,
               lines = 1e+05,
               output = NULL,
               bylines = T) {
    wk <- jiebaR::worker(type, dict, hmm, user, idf, stop_word, write, qmax, topn, encoding, detect,
                         symbol, lines, output, bylines)
    return(jiebaR::segment(texts, wk))
}
