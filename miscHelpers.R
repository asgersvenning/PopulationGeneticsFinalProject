align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
} 

label_grouped_axis <- function(map, min_size = 5, verbose = F) function(labels) {
  group <- map[labels] %>% 
    unname
  
  group_names <- unique(group)
  group_set   <- c()
  
  groups <- group %>%
    factor(group_names) %>% 
    as.integer %>%
    mod(2) %>%
    as.bit %>%
    decode_binary_vector()
  
  new_start <- c(groups$start, groups$end[1:(length(groups$end) - 1)] + 1)  %>% 
    sort
  new_end <- c(groups$end, groups$start[2:length(groups$start)] - 1) %>% 
    sort
  groups$start <- new_start
  groups$end   <- new_end
  if (groups$end[length(groups$end)] != groups$length) {
    groups$end <- c(groups$end, length(groups$end))
    groups$start <- c(groups$start, groups$start[length(groups$start)] + 1)
  }
  
  if (verbose) {
    oldPar <- par(mfrow = c(2, 1))
    
    group %>%
    factor(group_names) %>% 
    as.integer %>%
    mod(2) %>%
    as.bit %>% 
    plot(type = "s", main = "Before")
  }
  
  groups_length <- groups$end - groups$start + 1
  groups_center_ind <- groups$start + floor(groups_length / 2)
  if (verbose) {
    axis(1, groups_center_ind, group_names)
  }
  groups_remove <- vector("logical", length(groups_length))
  groups_remove <- !groups_remove
  groups_remove[1] <- F
  
  last_center <- groups_center_ind[1]
  this_length <- 0
  for (i in 2:length(groups_length)) {
    this_length <- this_length + (groups_center_ind[i] - last_center)
    if (verbose) print(this_length)
    if (this_length >= min_size) {
      this_length <- 0 
      groups_remove[i] <- F
      group_set <- c(group_set, group_names[i])
      if (verbose) print(paste0("--> ", group_names[i]))
    }
    last_center <- groups_center_ind[i]
  }
  
  
  groups_center_ind <- groups_center_ind[!groups_remove]
  if (verbose) {
    print("##################################")
    print(paste0("group : ", capture.output(str(group))))
    print(paste0("groups : ", capture.output(str(groups))))
    print(paste0("group_center_ind : ", capture.output(str(groups_center_ind))))
    print(paste0("group_length : ", capture.output(str(groups_length))))
    print(paste0("group_remove : ", capture.output(str(groups_remove))))
    print("All groups:")
    print(group_names)
    print("Removed groups:")
    print(unique(group)[groups_remove])
    print("Retained groups: ")
    print(group[groups_center_ind])
    print("##################################")
    
    groups$end <- groups$end[!groups_remove]
    groups$start <- groups$start[!groups_remove]
    do.call(encode_binary_vector, groups) %>% 
      plot(type = "s", main = "After")
    axis(1, groups_center_ind, group_names[!groups_remove])
    
    par(oldPar)
  }
  
  group[-groups_center_ind] <- ""
  stopifnot(all(group_set %in% unique(group)))
  group
}

skip_labels <- function(n) function(x) replace(x, which(!(1:length(x) %in% floor(seq(1, length(x), length.out = n)))), "")

label_prettify_scientific <- function(parse, digits) {
  Vectorize(function(x) {
    if (is.null(x) || is.na(x)) return(NA)
    if (is.numeric(x)) x <- as.character(signif(x, digits = digits))
    if (!is.character(x)) stop(paste0("Cannot parse type of \"", class(x), "\"."))
    if (!str_detect(x, "e\\+|e\\-")) return(x)
    out <- if (str_detect(x, "e\\+")) {
      parts <- str_split_fixed(x, "(?=e)", 2)
      parts[1] <- round(as.numeric(parts[1]), digits = digits)
      x <- paste0(parts, collapse = "")
      paste0("$", str_replace(x, "((?<=^)1)*e\\+", " \\\\times\\\\, 10^{"),"}$")
    }
    else if (str_detect(x, "e\\-")) {
      parts <- str_split_fixed(x, "(?=e)", 2)
      parts[1] <- round(as.numeric(parts[1]), digits = digits)
      x <- paste0(parts, collapse = "")
      paste0("$", str_replace(x, "((?<=^)1)*e\\-", " \\\\times\\\\, 10^{-"),"}$")
    }
    else {
      stop("Invalid format!")
    }
    out <- str_remove(out, "(?<=^\\$) \\\\times\\\\, ")
    if (parse) return(latex2exp::TeX(out)) else out
  })
}

prettify_scientific <- label_prettify_scientific(T, 2)

abslog <- function(base = 10, mult = 1) { 
  scales::trans_new(
    "abslog",
    transform = function(x) {
      sign(x) * log(abs(x * mult) + 1, base = base)
    },
    inverse =   inv <- function(x) {
      sign(x) * (base^abs(x) - 1) / mult
    }
    # breaks = scales::breaks_log(n, base)
  )
}

invlog <- function(base = 10, add = 0, mult = 1) { 
  scales::trans_new(
    "invlog",
    transform = function(x) {
      -log((1 - x + add), base) * mult
    },
    inverse =   inv <- function(x) {
      1 - base^(-x / mult) / mult + add
    }
    # breaks = scales::breaks_log(n, base)
  )
}

normlog <- function(base = 10, add = 0, mult = 1) { 
  scales::trans_new(
    "invlog",
    transform = function(x) {
      log((x + add), base) * mult
    },
    inverse =   inv <- function(x) {
      (base^(x / mult)) + add
    }
    # breaks = scales::breaks_log(n, base)
  )
}
