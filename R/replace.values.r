#' Replace values
#'
#' Replace the values of a categorical variable with new values
#'
#' @param input_data A vector of values (N) to be replaced.
#' @param old_values A vector of length n containing the values to be replaced.
#' @param new_values A vector of length n containing the values to replace the old_values.
#' @note The length of old_values and new_values should be equal.
#' @export

replace.values <- function(input_data, old_values, new_values){

# Check if the length of replaced values equals with new values
if (length(old_values) != length(new_values)){
warning("The length of the old_values and new_values vectors does not match.")
}
output_data<-input_data
for (i in 1:length(old_values)){
output_data <- replace(output_data, output_data == old_values[i], new_values[i])
}
return(output_data)
}
