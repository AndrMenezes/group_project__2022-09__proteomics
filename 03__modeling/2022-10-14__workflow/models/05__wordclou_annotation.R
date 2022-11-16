library(wordcloud)

path_res <- "./03__modeling/2022-10-14__workflow/models/results"


data_annotated <- read.table(file.path(path_res, "Proteins_STRING.tsv"),
                             header = TRUE)
dplyr::glimpse(data_annotated)

x <- data_annotated |> 
  dplyr::filter(method == "limma") |> 
  dplyr::pull(annotation)

stopwords_en <- tm::stopwords("en")
my_stopwords_en <- c("not", "yet", "protein")
stopwords_en <- c(stopwords_en, my_stopwords_en)

x <- tolower(x)
x <- tm::removePunctuation(x)
x <- tm::stripWhitespace(x)
x <- tm::removeWords(x, stopwords_en)
x <- tm::Corpus(tm::VectorSource(x))

tdm_x <- tm::TermDocumentMatrix(x)

mat_x <- as.matrix(tdm_x)
word_freqs <- sort(rowSums(mat_x), decreasing = TRUE)
dm <- data.frame(word = names(word_freqs), freq = word_freqs)
head(dm)
wordcloud::wordcloud(words = dm$word, freq = dm$freq, min.freq = 2,
                     random.order = FALSE, rot.per = 0.05,
                     colors = RColorBrewer::brewer.pal(8, "Dark2"))
