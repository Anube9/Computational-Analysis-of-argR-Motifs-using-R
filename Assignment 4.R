setwd(getwd())
#Question 1
# Read in counts matrix from text file
counts <- as.matrix(read.csv("argR-counts-matrix.txt", header = FALSE))

#Reading it as a dataframe
counts <- as.data.frame(counts, stringsAsFactors = FALSE)


# separate values by tabs and create new data frame with properly named columns
counts <- read.table(text = counts$V1, sep = "\t", header = FALSE, 
                     col.names = c("letter", "col1", "col2", "col3", "col4", "col5", 
                                   "col6", "col7", "col8", "col9", "col10", 
                                   "col11", "col12", "col13", "col14", "col15", 
                                   "col16", "col17", "col18", "col19"), 
                     stringsAsFactors = FALSE, row.names = "letter")

# deleteing th eunused columns from data frame 
counts <- counts[,-1]
print(counts)

#frequency matrix 
colsums_counts <- colSums(counts)
frequency_matrix <- counts/colsums_counts 
frequency_matrix

# Add pseudocounts to counts matrix
pseudo_counts <- (counts+1)/(colsums_counts+4)

#weighted matrix 
W_matrix <- log2(pseudo_counts/0.25)
print(W_matrix)

#QUESTION 2
install.packages("Biostrings")
library(Biostrings)

#reading the unzipped file 
E_coli<- read.csv("E_coli_K12_MG1655.400_50-1", header = FALSE)

#changing it to a dataframe
E_coli <- data.frame(E_coli)

# separate the string and set the first part as the row name
E_coli <- read.table(text = E_coli$V1, sep = "\\",stringsAsFactors = FALSE)


#deleting the third column 
E_coli <- subset( E_coli, select = -c(V3))


# Define motif length and upstream length
motif_length <- ncol(W_matrix)


#loop that takes in the sequences and gene_ids and gets the motifs according to the motif length and calculates the final score by using the weight matrix depending on the base indexes
results <- list()
for (i in c(1:nrow(E_coli))) {
  gene_ids <- E_coli[i,]$V1# append each gene ID to the vector
  sequences <- gsub(' ','',E_coli[i,]$V2)
  final_score <- c()
  for (j in c(1:(nchar(sequences)-(motif_length+1)))) {
    subsequences <- subseq(sequences, start=j, end=j+motif_length-1)
    score = 0
    for (i in c(1:nchar(subsequences))) {
      base_index <- substr(subsequences,i,i)
      score = score + W_matrix[base_index,i]
    }
    final_score <- c(final_score,score)
  }
  max_score <- max(final_score)
  results[[gene_ids]] = max_score
}


#The results from the above step are sorted in the descending order
my_results <- unlist(results)
final_results <- sort(my_results, decreasing = TRUE)
final_results <- data.frame(final_results)

#The top 30 gene_ids with the maximum scores are printed 
top_30 <- head(final_results, n = 30)
print(top_30)



