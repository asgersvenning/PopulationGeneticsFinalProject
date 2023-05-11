randomRound <- function(x) {
  decimal_part <- x %% 1
  integer_part <- floor(x)
  integer_part + rbinom(length(x), 1, decimal_part)
}

encode_binary_vector <- function(start, end, length) {
  vec <- bit(length = length)
  if (is.null(start) != is.null(end) || length(start) != length(end)) stop("Number of start and end-points of segments not equal!")
  
  if (!is.null(start) && length(start) > 0) for (i in 1:length(start)) {
    vec[start[i]:end[i]] <- T
  }
  
  vec
}

create_bit_vec <- function(brs, downscale = 1, type = "random") {
  type <- pmatch(type, c("random", "round", "floor", "ceiling", "expand", "contract"))
  if (is.na(type)) stop("type must be either 'random', 'round', 'floor', ceiling', 'expand' or 'contract' or an abbreviation of these.")
  if (!(!is.null(downscale) && is.numeric(downscale) && is.finite(downscale) && downscale >= 1)) stop("Downscale must be greater than or equal to 1 and finite!")
  if (downscale != 0) {
    len <- ceiling(brs$length / downscale)
    if (type == 1) brs <- brs %>% 
        lapply(function(x) randomRound(x / downscale))
    if (type == 2) brs <- brs %>% 
        lapply(function(x) round(x / downscale))
    if (type == 3) brs <- brs %>% 
        lapply(function(x) floor(x / downscale))
    if (type == 4) brs <- brs %>% 
        lapply(function(x) ceiling(x / downscale))
    if (type == 5) {
      brs$start <- floor(brs$start / downscale)
      brs$end   <- ceiling(brs$end / downscale)
    }
    if (type == 6) {
      brs$start   <- ceiling(brs$start / downscale)
      brs$end <- floor(brs$end / downscale)
    }
    brs$length <- len
  }
  invalid_secs <- which(brs$start > brs$end)
  if (length(invalid_secs) > 0) {
    brs$start <- brs$start[-invalid_secs]
    brs$end   <- brs$end[-invalid_secs]   
  }
  
  encode_binary_vector(brs$start, brs$end, brs$length)
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
  merge_sec <- which((start + 1) == lag(end))
  if (length(merge_sec) != 0) {
    start <- start[-merge_sec]
    end <- end[-(merge_sec + 1)]
  }
  list(start = start, end = end, length = length)
}

person_pairwise_dist_df <- function(df, col, ref, implicit = F, ...) {
  # Extract person "objects"
  persons <- df[[col]] 
  # Get person name attribute
  labels <- names(persons)
  # Pre-compute binary bit arrays encodings of persons archaic state
  vecs <- if (!implicit) sapply(persons, function(x) do.call(encode_binary_vector, x)) else persons
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
      dmat[i,j] <- if (!implicit) bit_vec_iou(person_i, person_j, ...) else implicit_iou(person_i, person_j)
      pb$tick()
    }
  }
  
  dmat <- 1 - dmat
  
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


prettify_scientific <- Vectorize(function(x) {
  if (is.na(x)) return(NA)
  if (!is.character(x)) x <- as.character(x)
  if (!str_detect(x, "e\\+|e\\-")) return(x)
  out <- if (str_detect(x, "e\\+")) {
    paste0("$", str_replace(x, "1*e\\+", " \\\\cdot 10^{"),"}$")
  }
  else if (str_detect(x, "e\\-")) {
    paste0("$", str_replace(x, "1*e\\-", " \\\\cdot 10^{-"),"}$")
  }
  else {
    stop("Invalid format!")
  }
  str_remove(out, "(?<=^\\$ )\\\\cdot") %>% 
    latex2exp::TeX()
})
