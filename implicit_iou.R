implicit_iou <- function(p1, p2) {
  intersection <- list(start = c(), end = c(), length = p1$length)
  union <- list(start = c(), end = c(), length = p1$length)
  ## Check if either p1 or p2 contains zero segments
  if (length(p1$start) == 0 | length(p2$start) == 0) {
    intersection <- c()
    if (length(p1$start) == 0 & length(p2$start) == 0) {
      union <- c()
    }
    if (length(p1$start) == 0) {
      union <- p2
    } 
    if (length(p2$start) == 0) {
      union <- p1
    }
  }
  else {
    for (i in 1:length(p1$start)) {
      for (j in 1:length(p2$start)) { 
        if (p2$start[j] > p1$end[i]) { # Case where p2-segment starts after end of p1-segment
          p2$start <- p2$start[j:length(p2$start)]
          p2$end <- p2$end[j:length(p2$end)]
          ## Skip all further segments in p2, since they can no longer intersect with the current
          ## segment in p1 (assuming segments are sorted)
          break
        } 
        else if (p2$end[j] < p1$start[i]) { # Case where p2-segment ends before start of p1-segment
          union$start <- c(union$start, p2$start[j])
          union$end   <- c(union$end, p2$end[j])
          
          ## Remove segment since it is no longer possible to intersect with any further segments in
          ## p1 (assuming segments are sorted)
          p2$end[j]   <- NA
          p2$start[j] <- NA
        }
        else { # Case where the segments intersects
          starts <- c(p1$start[i], p2$start[j])
          ends   <- c(p1$end[i], p2$end[j])
          intersection$start <- c(intersection$start, max(starts))
          intersection$end   <- c(intersection$end, min(ends))
          p2_before <- F
          p2_after  <- F
          if (p2$start[j] < p1$start[i]) { # Case where p2 starts before p1
            union$start <- c(union$start, min(starts))
            union$end   <- c(union$end, max(starts) - 1)
            
          }
          if (p2$end[j] > p1$end[i]) { # Case where p2 ends after p1
            union$start <- c(union$start, min(ends) + 1)
            union$end   <- c(union$end, max(ends))
            
            # Remove the overlapping part of p2, assuming that segments do not overlap within p1 alone
            p2$start[j] <- p1$end[i] + 1
          }
          else {
            p2$end[j]   <- NA
            p2$start[j] <- NA
          }
        }
      }
      union$start <- c(union$start, p1$start[i])
      union$end   <- c(union$end, p1$end[i])
      p2$start <- p2$start[!is.na(p2$start)]
      p2$end   <- p2$end[!is.na(p2$end)]
      if (length(p2$start) == 0) break
    }
    union$start <- c(union$start, p2$start)
    union$end   <- c(union$end, p2$end)
  }
  
  # Calculate length of intersection & union
  union <- person_segment_length(union)
  intersection <- person_segment_length(intersection)
  # Handle case where both union and intersection == 0
  if (union == 0) return(0)
  intersection / union
}

person_segment_length <- function(p) {
  if (length(p) == 0) return(0)
  sum(p$end - p$start + 1)
}

tidy_segments <- function(p) {
  dfs <- p %>% 
    lapply(function(x) tibble(start = x$start, end = x$end))
  
  tibble(
    person = if (is.null(names(p))) factor(1:length(p)) else factor(names(p)),
    breaks = dfs
  ) %>% 
    unnest(breaks)
}

bit_vec_to_person <- function(v) {
  pos <- v %>% 
    as.which() %>% 
    as.vector  
  
  pos_start <- which((pos - lag(pos, default = -10)) != 1)
  pos_end   <- which((lead(pos, default = length(v) + 10) - pos) != 1)
  
  list(
    start = pos[pos_start],
    end = pos[pos_end],
    length = length(v)
  )
}