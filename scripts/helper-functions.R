
make_similarity_matrix_fun <- function(lower_tri_df) {
  n <- max(lower_tri_df$row, lower_tri_df$col)
  similarity_matrix <- matrix(0, nrow = n, ncol = n)
  for (k in 1:nrow(lower_tri_df)) {
    i <- lower_tri_df$row[k]
    j <- lower_tri_df$col[k]
    similarity <- lower_tri_df$similarity[k]
    similarity_matrix[i, j] <- similarity
    similarity_matrix[j, i] <- similarity
  }
  diag(similarity_matrix) <- 1; 
  return(similarity_matrix)}


# Convert to lower triangle with indices (Option 3)
similarity_to_ltr_fun <- function(similarity_matrix){
  
lower_tri_indices <- which(lower.tri(similarity_matrix, diag = TRUE), arr.ind = TRUE)
  lower_tri_df <- data.frame(
  row = lower_tri_indices[, 1],
  col = lower_tri_indices[, 2],
  similarity = similarity_matrix[lower_tri_indices] )

 return(lower_tri_df)}


make_filename_safe <- function(name) {
  name <- iconv(name, from = "latin1", to = "UTF-8")
  name |>
    tolower() |>
    (\(x) gsub("[^a-z0-9]+", "_", x))() |>
    (\(x) gsub("^_|_$", "", x))() |>
    (\(x) substr(x, 1, 100))()
}

to_pylist <- function(vec) {
  if (length(vec) == 0) {
    list()
  } else {
    as.list(vec)
  }
}
