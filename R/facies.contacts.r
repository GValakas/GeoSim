#' Calculate facies contacts
#'
#' This function calculates the contacts of facies in a borehole dataset and plots them.
#' @param data_Facies A data frame with columns 'id' and 'Facies', where 'id' determines the borehole ID and 'Facies' contains the facies categories. 
#' If no 'id' column is defined, insert the 'x' and 'y' columns instead. 
#' @param facies_color A vector with colors representing each facies category, arranged in alphabetical order. This argument is used when the contact probabilities 'prot' are produced (if not equal to NULL). 
#' @param display A value indicating whether to produce bar plots of contacts for each facies (1) or not (0).
#'
#' @return The frequencies of facies contacts and their respective plots.
#'
#' @examples
#' data_Facies <- SFM_data
#' facies_color <- c("orange", "black", "red", "green", "blue") 
#' # Colors in alphabetical order of facies: AL, CO, HD, MG, SG
#' facies.contacts(data_Facies, facies_color, display = 1)
#'
#' @export

facies.contacts <- function(data_Facies,facies_color,display){

facies_categories<-rownames(table(data_Facies$Facies))
nfacies<-length(facies_categories)

data_colnames <- colnames(data_Facies)
# Check the id column
if (length(which(data_colnames == "id")) == 0){

	if (length(which(colnames(data_Facies) == "x")) == 0 || length(which(colnames(data_Facies) == "y")) == 0){
	 stop("GeoSim Package: The data frame 'data_Facies' must include either the 'id' column for borehole ID or the 'x' and 'y' columns for coordinates. Please ensure that your data includes a column with the name 'id' or two columns containing the coordinates with names 'x' and 'y'.")
	}
	
index_col_x <- which(colnames(data_Facies) == "x")
index_col_y <- which(colnames(data_Facies) == "y")

ind <- duplicated(data_Facies[,c(index_col_x,index_col_y)])
unique_xy <- cbind(as.numeric(data_Facies$x[which(ind==FALSE)]), as.numeric(data_Facies$y[which(ind==FALSE)]))
create_id <- rep(NaN, nrow(data_Facies))
xy_data_facies <- cbind(as.numeric(data_Facies$x), as.numeric(data_Facies$y))
for (i in 1: nrow(unique_xy)){
create_id[which( xy_data_facies[ ,1] == unique_xy[i,1] & xy_data_facies[ ,2] == unique_xy[i,2])] <- paste("id_",i)
}

data_Facies <- cbind(create_id, data_Facies)
colnames(data_Facies)<- c("id",data_colnames)
}

# Create the possible combinations of facies contacts
comb_facies <- combn(facies_categories,2)
possible_contacts <- paste(comb_facies[1,],comb_facies[2,],sep=".")

# id and number of boreholes
index_unique_id <- which(!duplicated(data_Facies$id) == TRUE) #  unique id per borehole
n_boreholes <- length(index_unique_id) # number of boreholes

# Compute the frequency of possible facies contacts
data_contacts_boreholes <- matrix(0,nrow = n_boreholes, ncol = length(possible_contacts))
for (i_borehole in 1:n_boreholes){
data_Facies_borehole <- data_Facies[ which(data_Facies$id == data_Facies$id[index_unique_id[i_borehole]]), ] # data per borehole
contacts_borehole <- matrix(NaN, nrow = (nrow(data_Facies_borehole) - 1), 1)
for (i_contact in 1:(nrow(data_Facies_borehole) - 1)){
label_sort<- sort(c(data_Facies_borehole$Facies[i_contact],data_Facies_borehole$Facies[i_contact+1]))
contacts_borehole[i_contact] <- paste(label_sort[1],".",label_sort[2],sep="") # create the facies contacts per borehole
}

for (i_pos_contact in 1:length(possible_contacts)){
check_index <- which(rownames(table(contacts_borehole)) == possible_contacts[i_pos_contact])
if (length(check_index) == 1){
data_contacts_boreholes[i_borehole, i_pos_contact] <- table(contacts_borehole)[check_index]
} # end if
} # end loop i_pos_contact
} # end loop i_borehole
colnames(data_contacts_boreholes) <- possible_contacts
contacts_freq <- apply(data_contacts_boreholes, 2, sum)

# Convert the Facies contact into matrix
contacts_freq_matrix <- matrix(0, nrow = length(facies_categories), ncol = length(facies_categories))
colnames(contacts_freq_matrix) <- facies_categories
rownames(contacts_freq_matrix) <- facies_categories
sq<-c(0, seq(length(facies_categories)-1, 1, -1))

for (i in 1:(length(facies_categories) - 1)){
contacts_freq_matrix[i, (i + 1):length(facies_categories)] <- contacts_freq[(sum(sq[1:i]) +1 ):sum(sq[1:(i + 1)])]
}
contacts_matrix <- contacts_freq_matrix + t(contacts_freq_matrix)
diag(contacts_matrix) <- NaN

if (display==1){
# Check the number of colors
if (length(facies_color) != nfacies){
warning(paste("GeoSim Package: The function has identified", nfacies, "facies categories in the given data, but the number of colors provided does not match the number of facies categories.", call. = FALSE))
}
for (i in 1:nfacies){
dev.new(height = 4, width = 3.5)
barplot(((contacts_matrix[-i, i]/ sum(contacts_matrix[-i, i])) * 100),col = facies_color[-i],horiz = TRUE,
xlim = c(0, 100),
xlab=paste(facies_categories[i],"contacts in %"),
border = NaN, las=1)
box(col = facies_color[i])
} # end loop i
} 
return(contacts_freq)
}
