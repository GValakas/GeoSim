#' Replace values
#'
#' Replace the values of a categorical variable with new values
#'
#' @param input_data A vector of values (N) to be replaced
#' @param old_values A vector with n values to be replaced
#' @param new_values A vector with n values to replace the old_values
#' @note The length of old values and new values should be equal.
#' @export

replace.values <- function(input_data, old_values, new_values){

# Check if the length of replaced values equals with new values
if (length(old_values) != length(new_values)){
warning("The length of the old and new values does not equal" )
}
output_data<-input_data
for (i in 1:length(old_values)){
output_data <- replace(output_data, output_data == old_values[i], new_values[i])
}
return(output_data)
}
