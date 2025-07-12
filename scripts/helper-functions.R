`%notin%` <- Negate(`%in%`)

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

make_distance_matrix_fun <- function(lower_tri_df) {
  n <- max(lower_tri_df$row, lower_tri_df$col)
  distance_matrix <- matrix(0, nrow = n, ncol = n)
  for (k in 1:nrow(lower_tri_df)) {
    i <- lower_tri_df$row[k]
    j <- lower_tri_df$col[k]
    distance <- lower_tri_df$distance[k]
    distance_matrix[i, j] <- distance
    distance_matrix[j, i] <- distance
  }
  diag(distance_matrix) <- 0; 
  return(distance_matrix)}

# Convert to lower triangle with indices
similarity_to_ltr_fun <- function(similarity_matrix){
  
lower_tri_indices <- which(lower.tri(similarity_matrix, diag = TRUE), arr.ind = TRUE)
  lower_tri_df <- data.frame(
  row = lower_tri_indices[, 1],
  col = lower_tri_indices[, 2],
  similarity = similarity_matrix[lower_tri_indices] )

 return(lower_tri_df)}

# Convert to lower triangle with indices
distance_to_ltr_fun <- function(distance_matrix){
  
  lower_tri_indices <- which(lower.tri(distance_matrix, diag = TRUE), arr.ind = TRUE)
  lower_tri_df <- data.frame(
    row = lower_tri_indices[, 1],
    col = lower_tri_indices[, 2],
    distance = distance_matrix[lower_tri_indices] )
  
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

# Function to convert marginal correlations to partial correlations
marginal_to_partial_correlations <- function(marginal_corr, rho = 0.1) {
  
  marginal_corr <- pmax(marginal_corr, t(marginal_corr))
  glasso_fit <- glasso(marginal_corr, rho = rho)
  precision_matrix <- glasso_fit$wi
  
  P <- nrow(precision_matrix)
  partial_corr <- matrix(0, P, P)
  for (i in 1:P) {
    for (j in 1:P) {
      if (i == j) {
        partial_corr[i, j] <- 1.0
      } else {
        partial_corr[i, j] <- -precision_matrix[i, j] / sqrt(precision_matrix[i, i] * precision_matrix[j, j])
      }
    }
  }
  partial_corr <- (partial_corr + t(partial_corr)) / 2
  
  return(list(marginal_corr = marginal_corr, 
              precision_matrix = precision_matrix, 
              partial_corr = partial_corr))
}


extend_fg_hierarchy <- function(fg_hierarchy, new_groups) {
  for (group in new_groups) {
    if (!("name" %in% names(group)) || !("smarts" %in% names(group))) {
      warning("Invalid custom functional group structure.")
      next
    }
    fg_hierarchy <- append(fg_hierarchy, list(group))
  }
  return(fg_hierarchy)
}
