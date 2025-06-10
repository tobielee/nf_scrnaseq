args <- commandArgs(trailingOnly=TRUE)

# Define a function to extract parameter values
get_param_value <- function(param_name) {
  param_index <- match(param_name, args)
  if (is.na(param_index)) {
    stop(paste("Parameter", param_name, "not found."))
  }
  param_value <- args[param_index + 1]
  # parse bracketed list input as string vector 
  # (you will need to convert to numeric in Rscript if desired)
  if (grepl("\\[.*\\]", param_value) && grepl(",", param_value)) {
    # Remove square brackets and split by comma
    param_value <- strsplit(gsub("\\[|\\]", "", param_value), ",")[[1]]
  }
  return(param_value)
}
