#' Co-simulation
#' Perform conditional or uncoditional co-simulation of a continuous (transformed multi-Gaussian) and categorical (plurigaussian) GRFs using the turning bands method.
#' The function simulates categorical (facies) and continuous variables using spatially correlated GRFs. The GRFs are jointly Gaussian,
#' with zero means and unit variances. The GRFs should be correlated.
#' Multiple realizations of the simulated facies can be produced in a regular grid.
#'
#' @param simu_grid A list of arguments to define the grid parameters. Use list(x0, y0, z0, nx, ny, nz, dx, dy, dz),
#' where:
#'   - x0, y0, z0: Single values indicating the minimum grid coordinates along the East-West, North-South, and vertical directions, respectively.
#'   - nx, ny, nz: Single values indicating the number of grid nodes along the East-West, North-South, and vertical directions, respectively.
#'   - dx, dy, dz: Single values indicating the grid mesh sizes along the East-West, North-South, and vertical directions, respectively.
#'
#' @param conditioning_data A matrix of conditioning data (categorical and continuous) with dimensions (number of data x 6). The columns are as follows:
#'   - Columns 1-3: Coordinates (x, y, z) of the conditioning data.
#'   - Column 4: Conditioning values of the categorical data (numerical codes of facies).
#'   - Column 5: Conditioning values of the continuous data.
#'   - Column 6: Normal scores of the continuous data.
#' Use NaN for the values in the 4th to 6th columns when categorical and continuous values are not at the same locations.
#' When no conditional data is imported, leave this argument empty to build an unconditional model.
#' 
#' @param truncation_rule A list of arguments to determine the truncation rule of Plurigaussian. Use list(nfield, nthres, thresholds, flag),
#' where:
#'   - nfield: Number of GRFs used for the truncation rule.
#'   - nthres: A vector with the number of thresholds for each GRF (1 x nfield).
#'   - thresholds: A vector of the thresholds for all GRFs (1 x sum(nthres)).
#'   - flag: A vector with category numbers codifying the truncation rule.
#'   - vpc_matrix: A matrix containing the proportions of facies for each z-level. It can be the matrix of vertical proportion curves obtained from the vpc function. Each row name in the matrix should define the corresponding z-level. Currently, the input is limited to a matrix that represents the vertical proportions across the entire examined area. If no available data is provided, leave it as NULL.
#'
#' @param variog_model A list of arguments to define the nested variogram theoretical models. Use list(model, cc, b, nugget),
#' where:
#'   - model: A matrix containing the covariance model for the GRFs (nested structures x 7 matrix). Each row corresponds to a nested structure and is codified as: type, scale factors, angles. Use the codes of the available types of variogram models (see details). There are three scale factors (along the rotated NS, EW, and vertical axes) and three angles to define the coordinate rotation (azimuth, dip, and plunge). See Deutsch & Journel (1997, p. 25) for more information.
#'   - cc: A matrix indicating the sills of nested structures (nested structures x (1 + nfield)^2).
#'   - b: A column vector with the additional parameters required for specific covariance types (nested structures x 1). See details for specific covariance types and their corresponding requirements.
#'   - nugget: A row vector with the nugget effect variance-covariance matrix of size (1 + nfield)^2.
#'
#' @param neigb_par A list of arguments to define the moving neighborhood to condition the data. Use list(radius, angles, octant, ndata),
#' where:
#'   - radius: A row vector with the maximum search radii along the rotated NS, EW, and vertical axes (1 x 3) for conditioning data.
#'   - angles: A row vector with the angles for anisotropic search (1 x 3), based on Deutsch and Journel (1997, p. 27).
#'   - octant: A value of 1 or 0, specifying if the neighborhood should be divided into octants (1) or not (0).
#'   - ndata: A single value indicating the number of conditioning data per octant or in total.
#'
#' @param simu_par A list of arguments to define the simulation parameters. Use list(nlines, nrealiz, seed, nnodes, itGibbs),
#' where:
#'   - nlines: A single value indicating the number of lines to use for simulating the nested structures by turning bands.
#'   - nrealiz: A single value indicating the number of realizations to produce.
#'   - seed: The seed number for generating random values.
#'   - nnodes: The number of nodes per line segment. 
#'   - itGibbs: The number of iterations in Gibbs' sampling. 
#'
#' @return A list of the simulated coordinates and the respective simulated categorical and continuous values.
#'   - coord: A matrix containing the center coordinates of gridded blocks (number of blocks x 3).
#'   - categoricalVar: A matrix containing the simulated values of the categorical variable (number of blocks x nrealiz).
#'     Each column represents a realization of the categorical variable.
#'   - continuousVar: A matrix containing the simulated values of the continuous variable (number of blocks x nrealiz).
#'
#' @references
#'   - Emery, X., & Silva, D. A. (2009). Conditional co-simulation of continuous and categorical variables for geostatistical applications. Computers & Geosciences, 35(6), 1234–1246. https://doi.org/10.1016/j.cageo.2008.07.005
#'
#' @seealso pluri.sim
#'
#' @export
					  
co.sim <- function(simu_grid, conditioning_data, truncation_rule, variog_model, neigb_par, simu_par){
simucoord <- c()
# Prepare and check the inputs for cosimulation
# =============================================
cc_unique <- -1
m0 <- 0
# Set data whose values are not in the trimming limits interval to NaN
if (!exists("limits")){
limits <- NULL
}
if (is.null(limits)){
limits <- c(-5, 5)
}
if (!all(is.na(conditioning_data[,1:4]))){ 
# Conditioning Data
m0 <- nrow(conditioning_data)
datacoord <- conditioning_data[ , 1:3]
idata <- matrix(conditioning_data[ , 4], ncol = 1)
ydata <- matrix(conditioning_data[ , 6], ncol = 1)
cc_unique <- 1

ydata[which((ydata < limits[1]) | (ydata > limits[2]))] <- NaN
# Anamorphosis Modeling (to transform back the continuous variable)
index_continuous_data <- which(!is.nan(ydata))
if (length(index_continuous_data) > 0){
tableYZ <- matrix(NaN, nrow = length(index_continuous_data), ncol = 2)
tableYZ[, 1] <- conditioning_data[index_continuous_data, 5]
tableYZ[, 2] <- conditioning_data[index_continuous_data, 6] + 1e-10 # avoid tied values in the conversion table
} else {
tableYZ <- matrix(NaN, nrow = m0, ncol = 2)
tableYZ[, 1] <- conditioning_data[ , 5]
tableYZ[, 2] <- conditioning_data[ , 6] + 1e-10
}
} else{ # unconditional simulation
index_YZ <- which(!is.nan(conditioning_data[ ,5]))
tableYZ <- matrix(NaN, nrow = length(index_YZ), ncol = 2)
tableYZ[, 1] <- conditioning_data[index_YZ, 5]
tableYZ[, 2] <- conditioning_data[index_YZ, 6] + 1e-10

index_YZ <- which(!is.nan(conditioning_data[ ,5]))
tableYZ <- matrix(NaN, nrow = length(index_YZ), ncol = 2)
tableYZ[, 1] <- conditioning_data[index_YZ, 5]
tableYZ[, 2] <- conditioning_data[index_YZ, 6] + 1e-10
}
tableYZ <- unique(tableYZ, fromLast = FALSE)
if (!exists("yzmin")){
yzmin <- NULL
}
if (is.null(yzmin)){
yzmin <- min(tableYZ[ , 1])
}
if (!exists("yzmax")){
yzmax <- NULL
}
if (is.null(yzmax)){
yzmax <- max(tableYZ[ , 1])
}
if (!exists("tail_anamor")){
tail_anamor <- NULL
}
if (is.null(tail_anamor)){
tail_anamor <- c(1, 1)
} 
# Truncation Rule
nfield <- truncation_rule$nfield
nthres <- truncation_rule$nthres
thresholds <- truncation_rule$thresholds
flag <- truncation_rule$flag
vpc_matrix <- truncation_rule$vpc_matrix
if (length(flag) != prod(nthres + 1)){stop('GeoSim Package: The truncation rule is inconsistent with the number of thresholds. Please make sure the length of the thresholds vector matches the total number of thresholds determined by the nthres vector.')}
# Coregionalization Model
model <- variog_model$model
nst <- nrow(model) # number of nested structures
cc <- variog_model$cc
nvar <- sqrt(ncol(cc)) #  number of variables (= 1 + nfield)
b <- matrix(variog_model$b)
nugget <- variog_model$nugget
sillnugget <- matrix(nugget, nvar, nvar, byrow = TRUE) 	
chk_variogram <- check_variogram(model, nfield, nst, cc, nugget, cosim = 1) # check and update
model <- chk_variogram$model
sill <- chk_variogram$sill
A1 <- chk_variogram$A1
A0 <- chk_variogram$A0
# Simulation Modeling
nlines <- matrix(simu_par$nlines)
if (length(nlines) != nst){nlines <- matrix(nlines[1] %*% rep(1, nst))} 
nrealiz <- simu_par$nrealiz
nrealiz <- nvar * nrealiz # number of realizations
seed <- simu_par$seed
if (is.null(simu_par$nnodes)){
nnodes <- 500
} else {
nnodes <- ceiling(simu_par$nnodes)
}
if (is.null(simu_par$itGibbs)){
itGibbs <- 1e4
} else {
itGibbs <- simu_par$itGibbs
}
# Moving Neighborhood Parameters 
if (!is.null(neigb_par)){
radius <- neigb_par$radius
angles <-  neigb_par$angles 
octant <- neigb_par$octant
ndata <- neigb_par$ndata
if (length(radius) != 3){
radius <- radius[1] * rep(1, 3)
} 
} else{
radius <- Inf  * rep(1, 3)
angles <- rep(0,3)
octant <- 0
ndata <- 0
}
# Conditional using moving neighborhood (unique = 0)  or conditional to all data data (unique = 1)
if (cc_unique != -1){
if (is.infinite(radius[1])){
cc_unique <- 1
} else {
cc_unique <- 0
} 
}
search_rotationmatrix <- setrot(t(c(1, radius, angles)), it = 1) # rotation-reduction matrix for data search 
  

# Locations for Co-simulation
dx <- simu_grid$dx; dy <- simu_grid$dy; dz <- simu_grid$dz
nx <- simu_grid$nx; ny <- simu_grid$ny; nz <- simu_grid$nz
if (!is.null(simucoord)){
x0 <- min(simucoord[ , 1]); y0 <- min(simucoord[ , 2]); z0 <- min(simucoord[ , 3])
xmax <- max(simucoord[ , 1]); ymax <- max(simucoord[ , 2]); zmax <- max(simucoord[ , 3])
nx <- length(seq(x0,xmax,dx)); ny <- length(seq(y0,ymax,dy)); nz <- length(seq(z0,zmax,dz))
} else {
x0 <- simu_grid$x0 + simu_grid$dx / 2 
y0 <- simu_grid$y0 + simu_grid$dy / 2 
z0 <- simu_grid$z0 + simu_grid$dz / 2
xmax <- (x0 + nx * dx) - dx / 2
ymax <- (y0 + ny * dy) - dy / 2
zmax <- (z0 + nz * dz) - dz / 2
simucoord <- data.matrix(expand.grid(seq(x0, xmax, dx),seq(y0, ymax, dy),seq(z0, zmax, dz),KEEP.OUT.ATTRS = FALSE))
}

# Remove data located too far from the locations to simulate
if (m0 > 0){
tmp <- datacoord %*% search_rotationmatrix
x <- matrix(c(x0, y0, z0, x0, y0, z0 + (nz - 1) * dz, x0, y0 + (ny - 1) * dy, z0, x0, y0 + (ny - 1) * dy, z0 + (nz - 1) * dz,
x0 + (nx - 1) * dx, y0, z0, x0 + (nx - 1) * dx, y0, z0 + (nz - 1) * dz, x0 + (nx - 1) * dx, y0 + (ny - 1) * dy, z0, 
x0 + (nx - 1) * dx, y0 + (ny - 1) * dy, z0 + (nz - 1) * dz),nrow = 8, ncol = 3, byrow = TRUE)
x <- x %*% search_rotationmatrix
minx = min(x[ ,1]); miny = min(x[ ,2]); minz = min(x[ ,3]);
maxx = max(x[ ,1]); maxy = max(x[ ,2]); maxz = max(x[ ,3]);
I <- which((tmp[ ,1] < (minx-1)) | (tmp[ ,1] > (maxx+1)) | (tmp[ ,2] < (miny-1)) | (tmp[ ,2] > (maxy+1)) | (tmp[ ,3] < (minz-1)) | (tmp[ ,3] > (maxz+1)) )
if (length(I) > 0){
datacoord <- datacoord[-I, ]
idata <- matrix(idata[-I, ],ncol=1) 
ydata <- matrix(ydata[-I, ],ncol=1) 
m0 <- nrow(datacoord)
}
} 

# Default values
eps <- 2.2204e-16
# Create seed numbers
set.seed(seed) # rand('state',seed);
seed_vdc <- ceiling(1e7 * runif(nst))
seed_line <- ceiling(1e7 * matrix(runif(nst * nrealiz * max(nlines)), nst * nrealiz, max(nlines)))

# Extremal coordinates to simulate (including data locations) in the original referential
if (m0 > 0){
minx <- min(c(datacoord[ ,1], x0)); miny <- min(c(datacoord[ ,2], y0)); minz <- min(c(datacoord[ ,3], z0)) 
maxx <- max(c(datacoord[ ,1], xmax-dx/2)); maxy <- max(c(datacoord[ ,2], ymax-dy/2)); maxz <- max(c(datacoord[ ,3], zmax-dz/2))
} else {
radius <- rep(0,3)  
minx <- x0; miny <- y0; minz <- z0
maxx <- x0 + (nx-1) * dx; maxy <- y0 + (ny - 1) * dy; maxz <- z0 + (nz - 1) * dz
} 
extreme_coord_vec <- c(minx, miny, minz, minx, miny, maxz, minx, maxy, minz, minx, maxy, maxz, maxx, miny, minz, maxx, miny, maxz, maxx, maxy, minz, maxx, maxy, maxz)
extreme_coord <- matrix(extreme_coord_vec, ncol = 3, byrow = TRUE); 
# Prepare the lines for tbmain
prepare_lines <- prepare_lines(model, nrealiz, nlines, nst, extreme_coord, seed_vdc)
max_nugget <- max(abs(nugget)) 

# NON CONDITIONAL SIMULATION AT DATA LOCATIONS
# ============================================
simudata <- matrix(0, nrow = m0, ncol = nrealiz)
cc_weights <- matrix(0, nrow = 1, ncol = nrealiz)
if (cc_unique >= 0){
nisb <- 0
nxsup <- 0
xmnsup <- 0
xsizsup <- 0
nysup <- 0
ymnsup <- 0
ysizsup <- 0
nzsup <- 0
zmnsup <- 0
zsizsup <- 0
ixsbtosr <- matrix(0, nrow=1)
iysbtosr <- matrix(0, nrow=1)
izsbtosr <- matrix(0, nrow=1)
if (cc_unique == 0){
# Prepare super-block search strategy
results_superblk <- superblk(datacoord,nx,ny,nz,x0,y0,z0,dx,dy,dz)
nisb <- results_superblk$nisb 
nxsup <- results_superblk$nxsup 
xmnsup <- results_superblk$xmnsup 
xsizsup <- results_superblk$xsizsup
nysup <- results_superblk$nysup
ymnsup <- results_superblk$ymnsup 
ysizsup <- results_superblk$ysizsup
nzsup <- results_superblk$nzsup
zmnsup <- results_superblk$zmnsup 
zsizsup <- results_superblk$zsizsup
# sort data according to ascending super-block number
datacoord <- datacoord[results_superblk$I,] 
idata <- matrix(idata[results_superblk$I,], ncol = 1) 
ydata <- matrix(ydata[results_superblk$I,], ncol = 1)  
# build template of super-blocks centered at the block containing the node to simulate
results_picksupr <- picksupr(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,search_rotationmatrix)
ixsbtosr <- matrix(results_picksupr$ixsbtosr, nrow=1) 
iysbtosr <- matrix(results_picksupr$iysbtosr, nrow=1)
izsbtosr <- matrix(results_picksupr$izsbtosr, nrow=1)
rm(results_picksupr)
}

# Prepare super-block search strategy 
t_start_neigb<-Sys.time() 
ydata2 <- Gibbs_cosim_cpp(datacoord, idata, ydata, nfield, flag, nthres, thresholds, model, sill, b, sillnugget, nrealiz / nvar, itGibbs, prepare_lines$model_rotationmatrix, search_rotationmatrix, cc_unique, octant, ndata, nxsup, nysup, nzsup, xmnsup, ymnsup, zmnsup, xsizsup, ysizsup, zsizsup, ixsbtosr, iysbtosr, izsbtosr,nisb)
t_finish_neigb<-Sys.time()
t_finish_neigb-t_start_neigb
# How many data locations can be simulated simultaneously?
m2 <- max(1 ,min(m0, nnodes)) 
cc_sequence <- c(seq(0,(m0 - 0.5), m2), m0) 
seed_nugget_data <- ceiling(1e7 * runif(1))   
# Loop over the sequences of data points
for (n in 1:(length(cc_sequence) - 1)){
index <- t(t((cc_sequence[n] + 1):cc_sequence[n + 1])) 
t_start_neigb<-Sys.time() 
simudata[index, ] <- tbmain_cosimu_cpp(datacoord[index, ],model, prepare_lines$cc_sigma, A1, nvar, nlines, nrealiz, seed_line, prepare_lines$all_lines, prepare_lines$all_offset, prepare_lines$all_r, prepare_lines$all_phi, prepare_lines$all_theta, prepare_lines$valid_lines)						
t_finish_neigb<-Sys.time()
t_finish_neigb-t_start_neigb
}
# Add nugget effect
if (max_nugget > eps){
set.seed(seed_nugget_data)
simunug <- matrix(rnorm(m0 * nrealiz), nrow = m0 * nrealiz / nvar, ncol = nvar)  
simunug <- matrix(t(simunug),nrow = nrealiz, ncol = m0) 
simudata <- simudata + t(simunug)
}
# Prepare conditioning co-kriging
simudata <- matrix(simudata, nrow = m0 * nvar, ncol = nrealiz / nvar)
cc_residuals <- rbind(ydata %*% matrix(1, nrow = 1, ncol = nrealiz / nvar), ydata2) - simudata
if (cc_unique == 1){# co-kriging in a unique neighborhood
index_missing <- which(is.nan(cc_residuals[ , 1]))
cc_residuals <- na.omit(cc_residuals)
cc_weights <- dual(datacoord, cc_residuals, index_missing, model, sill, b, sillnugget + 1e-7, prepare_lines$model_rotationmatrix)
} else{ # co-kriging in a moving neighborhood
cc_residuals <- matrix(cc_residuals, nrow = m0, ncol = nrealiz)
}
}
# Conditional Simulation 
# ======================
m1 <- nrow(simucoord) 
m2 <- max(1, min(m1, nnodes)) 
cc_sequence <- c(seq(0, (m1 - 0.5), m2), m1) 
nloops <- length(cc_sequence) - 1 
seed_nugget <- ceiling(1e7 * runif(nloops)) 
all_continuous_values <- matrix(NA, nrow = m1, ncol = (nrealiz / nvar))
all_categorical_values <- matrix(NA, nrow = m1, ncol = (nrealiz / nvar))
all_coord <- matrix(NA, nrow = m1, ncol = ncol(simucoord)) 

for (n in 1:nloops){ 
# Coordinates of the points to simulate
index <- t(t((cc_sequence[n] + 1):cc_sequence[n + 1])) 
m1 <- length(index) 
coord <- simucoord[index, ] 	
# Non-conditional simulation
t_start_neigb<-Sys.time() 
simu <- tbmain_cosimu_cpp(coord, model, prepare_lines$cc_sigma, A1, nvar, nlines, nrealiz, seed_line, prepare_lines$all_lines, prepare_lines$all_offset,
    prepare_lines$all_r, prepare_lines$all_phi, prepare_lines$all_theta, prepare_lines$valid_lines) 
t_finish_neigb<-Sys.time()
t_finish_neigb-t_start_neigb
# Add nugget effect
if (max_nugget > eps){ 
set.seed(seed_nugget[n])
simunug <- matrix(rnorm(m1 * nrealiz), nrow = (m1 * nrealiz / nvar), ncol = nvar) %*% A0 
simunug <- matrix(t(simunug), nrow = nrealiz, ncol = m1) 
simu <- simu + t(simunug)
}  
# Conditioning
if (cc_unique == 1){ # dual co-kriging
for (i in 1:m1){
# Substitution of residuals
k0 <- setdual(model, coord[i, ], sill, b, datacoord, index_missing, prepare_lines$model_rotationmatrix) 
simu[i, ] <- simu[i, ] + matrix(t(k0) %*% cc_weights, nrow = 1, ncol = nrealiz) 
}
} else if (cc_unique == 0) { 
t_start_neigb<-Sys.time()
simu <- cond_mooving_neigbor(m1, simu, nvar, nrealiz, datacoord, cc_residuals, coord, search_rotationmatrix,octant,ndata,nxsup, nysup, nzsup, xmnsup, ymnsup, zmnsup, xsizsup, ysizsup, zsizsup, ixsbtosr, iysbtosr, izsbtosr, nisb, model, sill, b, sillnugget, prepare_lines$model_rotationmatrix )
t_finish_neigb<-Sys.time()
t_finish_neigb-t_start_neigb
}
# Back transform to raw variable and rock types
# =============================================
if (!is.null(tableYZ)){ 
simu_seq <- simu[ , seq(1, nrealiz, nvar)]#matrix(simu[ , seq(1, nrealiz, nvar)], nrow = nrow(simu), ncol = length(seq(1, nrealiz, nvar)))
continuousvariable <- back_transform(simu_seq, tableYZ, yzmin, yzmax, tail_anamor)
continuousvariable <- matrix(continuousvariable,ncol = nrealiz/ nvar, byrow = FALSE)
}
simu <- simu[ ,-seq(1, nrealiz, nvar)] 
simu <- matrix(simu, nrow = nrow(coord), ncol = nrealiz - length(seq(1, nrealiz, nvar)))
if (is.null(vpc_matrix)){
categoricalvariable <- t(t(pluri_truncate(simu, nfield, flag, nthres, thresholds)))
categoricalvariable <- matrix(categoricalvariable, length(categoricalvariable) / (nrealiz / nvar), nrealiz / nvar, byrow = FALSE)
} else{
categoricalvariable <- t(t(vpc_truncate(coord,simu, nfield, flag, nthres, vpc_matrix)))
}
# Save results	
all_coord[index, ] <- coord
all_continuous_values[index, ] <- continuousvariable
all_categorical_values[index, ] <- categoricalvariable	     
}
colnames(all_coord) <- c("X", "Y", "Z")
cosimulation_results <- list(coord = all_coord, categoricalVar = all_categorical_values , continuousVar = all_continuous_values)	     
}
