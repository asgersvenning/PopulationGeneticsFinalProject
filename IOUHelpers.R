randomRound <- function(x) {
  decimal_part <- x %% 1
  integer_part <- floor(x)
  integer_part + rbinom(length(x), 1, decimal_part)
}

encode_binary_vector <- function(start, end, length) {
  vec <- bit(length = length)
  
  for (i in 1:length(start)) {
    vec[start[i]:end[i]] <- T
  }
  
  vec
}

bit_vec_iou <- function(v1, v2, r) {
  ######## CAN CAUSE CRASHES - USE WITH CAUTION ###############
  #### Don't call manually - doesn't check correctness of inputs
  if (missing(r)) {
    r <- c(1L, length(v1))
  } else {
    if (!is.integer(r)) {
      r <- as.integer(r)
      warning("r is not an integer! It has been coerced.")
    }
    r <- c(1L, r)
  }
  r2 <- bit(r[2])
  
  
  .Call(bit:::C_R_bit_sum, .Call(bit:::C_R_bit_and, v1, v2, r2), r) / # Intersection
    .Call(bit:::C_R_bit_sum, .Call(bit:::C_R_bit_or, v1, v2, r2), r)  # Union
}

person <- function(start, end, length) {
  list(start = start, end = end, length = length)
}

person_pairwise_dist_df <- function(df, col, ref, ...) {
  # Extract person "objects"
  persons <- df[[col]] 
  # Get person name attribute
  labels <- names(persons)
  # Pre-compute binary bit arrays encodings of persons archaic state
  vecs <- sapply(persons, function(x) do.call(encode_binary_vector, x))
  # Initialize 1 - IoU pairwise distance matrix
  dmat <- matrix(NA, length(persons), length(persons))
  colnames(dmat) <- labels
  rownames(dmat) <- labels
  # Initialize progress bar
  pb <- progress_bar$new(
    total = (length(persons) - 1) * length(persons) / 2,
    clear = T, 
    width = 100, 
    format = "Calculating pairwise IoU [:bar] :percent eta: :eta")
  
  # Main loop
  for (i in 1:length(persons)) {
    person_i <- vecs[[i]]
    for (j in 1:length(persons)) {
      # Don't calculate lower triangle; IoU is symmetric and IoU = 0 for i == j
      if (i <= j) next
      person_j <- vecs[[j]]
      dmat[i,j] <- 1 - bit_vec_iou(person_i, person_j, ...)
      pb$tick()
    }
  }
  
  ## Re-tidy pairwise distance matrix
  # Prepare for appending person meta data
  join1 <- ref
  join2 <- ref
  names(join1) <- paste0(col, "_", 1)
  names(join2) <- paste0(col, "_", 2)
  
  # Re-tidy
  as.dist(dmat) %>% 
    as.matrix %>% 
    as.data.frame %>% 
    rownames_to_column(paste0(col, "_", 1)) %>% 
    pivot_longer(!paste0(col, "_", 1), 
                 names_to = paste0(col, "_", 2), 
                 values_to = "iou") %>% 
    # Append meta data
    full_join(
      df %>%
        select(!all_of(col)) %>% 
        rename_with(function(x) ifelse(x == ref, x, paste0(x, "_", 1))),
      by = join1
    ) %>% 
    full_join(
      df %>% 
        select(!all_of(col)) %>% 
        rename_with(function(x) ifelse(x == ref, x, paste0(x, "_", 2))),
      by = join2
    )
}
