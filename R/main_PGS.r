#' Back-transformation from Gaussian to raw scale
#' @noRd

back_transform <- function(z, tableYZ, ymin, ymax, cc_tail){
# Back-transformation from Gaussian to raw scale
p <- nrow(tableYZ)
m <- nrow(z)
n <- ncol(z)
y <- matrix(0, nrow = m, ncol = n)
# Do not perform back-transformation if the conversion table is empty
if (is.null(tableYZ)){
y <- z
stop()
}
# Values in the lower tail (exponential extrapolation)
y1 <- tableYZ[1, 1]
z1 <- tableYZ[1, 2]
I1 <- which(z < z1)

if (length(I1) > 0){
b0 <- (y1 - ymin) * exp(-cc_tail[1] * z1)
y[I1] <- ymin + b0 * exp(cc_tail[1] * z[I1])
}
# Values in the upper tail (exponential extrapolation)
yp <- tableYZ[p, 1]
zp <- tableYZ[p, 2]
I2 <- which(z > zp)
if (length(I2) > 0){
bp <- (yp - ymax) * exp(cc_tail[2] * zp)
y[I2] <- ymax + bp * exp(-cc_tail[2] * z[I2])
}
# Within-class values
I <- c(I1, I2)
I3 <- 1:(m * n)
I3 <- t(t(I3[-I]))
if (length(I3) > 0){
y[I3] <- approx(tableYZ[, 2], tableYZ[, 1], z[I3]) # Table lookup
}
return(y)
}

#' Check the user input of the variogram for inconsistency
#' @noRd

check_variogram <- function(model, nfield, nst, cc, nugget, b, cosim){
# Check The imput of variogam modeling
# This function use only R basic functions

nvar <- sqrt(ncol(cc))
sillnugget <- matrix(nugget, nvar, nvar, byrow = TRUE)
# Check the sill matrix
if (nvar > floor(nvar)){
stop('GeoSim Package: The number of columns in the sill matrix is inconsistent.')
}
if (!missing(cosim) && (nvar != 1 + nfield)){
stop('GeoSim Package: The number of columns in the sill matrix is inconsistent.')
}
sill<-array(0, dim = c(nvar, nvar, nst))
A1 <- array(0, dim = c(nvar,nvar,nst))
for (i in 1:nst){
sill[ , ,i] <- matrix(cc[i, ], nvar, nvar, byrow = TRUE)
eigenvalues <- eigen(sill[ , , i])$values
eigenvectors <- eigen(sill[ , , i])$vectors
A1[ , , i] <- sqrt(eigenvalues) * eigenvectors
if (min(eigenvalues) < 0){
stop('GeoSim Package: The sill matrix for structure ',i,' is not semi-positive definite.')
}
}
# Check the nugget
if (!isSymmetric(sillnugget)){
stop('GeoSim Package: The sill matrix for the nugget effect is not symmetric.')
}
eigenvalues <- eigen(sillnugget)$values
eigenvectors <- eigen(sillnugget)$vectors
A0 <- sqrt(eigenvalues) * eigenvectors
if (min(eigenvalues) < 0){
stop('GeoSim Package: The sill matrix for the nugget effect is not semi-positive definite.')
}
# Check the Variogram Models
for (i in 1:nst){
if (model[i,1] == 3){ #Gamma model
if (b[i] <= 0){
stop('GeoSim Package: The parameter of the gamma model must be positive.')
} # end
} else if (model[i,1] == 4){ # Stable model
if (b[i] > 2){
stop('GeoSim Package: The parameter of the stable model must be less than 2.')
} else if (b[i] <= 0){
stop('GeoSim Package: The parameter of the stable model must be positive.')
} else if(b[i] == 2){ # Gaussian model
model[i,1] <- 6
} else if (b[i] == 1){ # Exponential model
model[i,1] <- 2
} else if (b[i] > 1){ # Simulation via spectral method
model[i,1] <- 4.5
}
} else if(model[i,1] == 8){ # Bessel-J model
if (b[i] < 0.5){
stop('GeoSim Package: The parameter of the Bessel-J model must be greater than 0.5.')
} else if(b[i] == 0.5){ # Cardinal sine model
model[i,1] <- 7
}
} else if (model[i,1] == 9){ # Bessel-K model
if (b[i] <= 0){
stop('GeoSim Package: The parameter of the Bessel-K model must be positive.')
} else if (b[i] == 0.5){ # Exponential model
model[i,1] <- 2
} else if (b[i] > 0.5){ # Simulation via spectral method
model[i,1] <- 9.5
}
} else if (model[i,1] == 10){ # Generalized Cauchy model
if (b[i] <= 0){
stop('GeoSim Package: The parameter of the generalized Cauchy model must be positive.')
}
}
}
return(list("model" = model,"sill" = sill, "nugget" = nugget, "A1" = A1, "A0" = A0))

}

#' Compute covariance for reduced distance h
#' @noRd

cova <- function(it, h, b){
# This function use only R basic functions
epsilon <- 1e-12
if (it < 2){ # Spherical model
h[which(h > 1)] <- 1
C <- 1 - 1.5 * h + 0.5 * (h ^ 3)
} else if (it < 3){ # Exponential model
C <- exp(-h)
} else if (it < 4){ # Gamma model
C <- 1 / (1 + h) ^ b
} else if (it < 5){ # Stable model
C <- exp(-h ^ b)
} else if (it < 6){ # Cubic model
h1 <- h
h1[which(h > 1)] <- 1
C <- 1 - 7 * (h1 ^ 2) + 35 / 4 * (h1 ^ 3) - 7 / 2 * (h1 ^ 5) + 3 / 4 * (h1 ^ 7)
} else if (it < 7){ # Gaussian model
C <- exp(-h ^ 2)
} else if (it < 8){ # Cardinal sine model
C <- sin(h + epsilon) / (h + epsilon)
} else if (it < 9){ # J-Bessel model
C <- ((h + epsilon) / 2) ^ (-b) * gamma(b + 1) * besselJ(b, h + epsilon)
} else if (it < 10){ # K-Bessel model
if (b == 0){
stop('b cannot be 0')
}
C <- (((h + epsilon) / 2) ^ b) * 2 / gamma(b) * besselK(b, h + epsilon)
} else if (it < 11){ # Generalized Cauchy model
C <- 1/(1 + h ^ 2) ^ b
} else if (it < 12){ # Exponential sine model
C <- sin(pi / 2 * exp(-h))
} else{
stop('Unavailable covariance model')
}
return(C)
}

#' Compute dual co-kriging weights
#' @noRd

dual <- function(datacoord, ydata, index_missing, model, sill, b, sillnugget, model_rotationmatrix){
# Compute dual co-kriging weights
# This function uses the following GeoSim functions: cova

n <- length(datacoord[, 1]) # number of data
nst <- nrow(model) # number of nested structures
nvar <- nrow(sill) # number of variables
nrealiz <- ncol(ydata) # number of realizations

# Calculation of the left covariance matrix K
k <- matrix(0, nrow = nvar * n, ncol = nvar * n)
for (i in 1:nst){
# Calculation of matrix of reduced rotated distances
R <- model_rotationmatrix[ , ,i]
tt <- datacoord %*% R
tt <- tt %*% t(tt)
h <- -2 * tt + t(t(diag(tt))) %*% rep(1, n) + rep(1, n) %*% t(diag(tt))
h[which(h < 0)] <- 0
h <- sqrt(h)
# Evaluation of the current basic structure
C <- cova(model[i,1], h, b[i])
k <- k + kronecker(sill[ , ,i],C)
}
k <- k + kronecker(sillnugget, diag(n))
# Remove rows and columns corresponding to missing data
if (length(index_missing) > 0){
k <- k[-index_missing, ]
k <- k[ ,-index_missing]
}
# Co-kriging weights
cc_weights <- solve(k, ydata)
return(cc_weights)
}

#' Compute the Normalized Hermite polynomials
#' @noRd

hermite <- function(y, Pmax){
# Compute the Normalized Hermite polynomials
# This function use only R basic functions

y <- t(t(y))
h<-matrix(c(rep(0, length(y)), rep(1, length(y))), nrow = 2, ncol = length(y), byrow=TRUE)
for (p in 0:(Pmax - 1)){ # compute polynomial of degree p
h[, p + 3] <- (-y * h[ ,p + 2] - sqrt(p) * h[ ,p+1]) / sqrt(p + 1)
}
h <- h[ ,-1]
}

#' Converts Gaussian values into categorical values
#' @noRd

pluri_truncate <- function(y_simu, nfield, flag, nthres, thresholds){
# Converts Gaussian values into categorical values
# This function use only R basic functions

m1 <- nrow(y_simu) #[m1,m2] = size(simu);
m2 <- ncol(y_simu)
y_simu <- array(y_simu, dim = c(m1, nfield, m2 / nfield))
y_simu <- aperm(y_simu,c(1, 3, 2))
sthres <- c(0, cumsum(nthres))
pthres <- c(1, cumprod(nthres + 1))
i<-1
for (k in 1:sthres[nfield + 1]){
j <- sum(k > sthres)
I <- matrix(0, nrow = nrow(y_simu), ncol = ncol(y_simu))
I[which(y_simu[ , ,j] > thresholds[k])] <- 1
i <- i + pthres[j] * I
}
y_simu <- flag[i]
return(y_simu)
}

#' Converts Gaussian values into categorical values considering the vertical proportion matrix of facies
#' @noRd

vpc_truncate <- function(coord,y_simu, nfield, flag, nthres, vpc_matrix){
# Converts Gaussian values into categorical values
# This function use only R basic functions

vpc_thresholds<-matrix(NaN, nrow = nrow(vpc_matrix),ncol = sum(nthres))
for (i in 1:nrow(vpc_matrix)){
vpc_thresholds[i,] <- calculate.thresholds(flag,nthres,proportions = as.vector(vpc_matrix[i, ,1]))
}
rownames(vpc_thresholds)<-rownames(vpc_matrix)
m1 <- nrow(y_simu) 
m2 <- ncol(y_simu)
y_simu <- array(y_simu, dim = c(m1, nfield, m2 / nfield))
y_simu <- aperm(y_simu,c(1, 3, 2))
sthres <- c(0, cumsum(nthres))
pthres <- c(1, cumprod(nthres + 1))
y_simu_temp <- matrix(NaN, nrow = m1, ncol = m2 / nfield)
depth <- as.numeric(rownames(vpc_thresholds))
for (d in 1:(nrow(vpc_thresholds))){
if (d==1){
z_index <- which(coord[,3] <= depth[d])
} else if((d > 1) & (d <= nrow(vpc_thresholds))){
z_index <- which((depth[d-1]<coord[,3]) & (coord[,3]<= depth[d]))
} else if (d == nrow(vpc_thresholds)){
z_index <- which(coord[,3] > depth[d])
}
if (length(z_index)!=0){
i<-1
for (k in 1:sthres[nfield + 1]){
j <- sum(k > sthres)
I <- matrix(0, nrow = length(z_index), ncol = ncol(y_simu)) 
I[which(y_simu[z_index, ,j] > vpc_thresholds[d,k])] <- 1
i <- i + pthres[j] * I
}
temp_flag <- matrix(flag[i], nrow = length(z_index), ncol =  m2 / nfield, byrow = FALSE )
y_simu_temp[z_index,] <- temp_flag 
}
}
y_simu <- as.vector(y_simu_temp)
return(y_simu)
}

#' Calculate the right-hand side member of the dual co-kriging system
#' @noRd

setdual <- function(model, coord, sill, b, datacoord, index_missing, model_rotationmatrix){
# Calculate the right-hand side member of the dual co-kriging system
# This function uses the following subroutine: cova

if (is.null(datacoord)){
k0 <- 0
}
nst <- nrow(model) # number of nested structures
nvar <- nrow(sill) # number of variables
m0 <- nrow(datacoord) # number of data
k0 <- matrix(0, nrow = nvar * m0, ncol = nvar)
for (i in 1:nst){
# Calculation of matrix of reduced rotated distances
R <- model_rotationmatrix[,,i]
h <- (matrix(1, nrow = m0, ncol = 1) %*% coord - datacoord) %*% R
h <- h ^ 2
h <- sqrt(t(t(apply(h, 1, FUN = sum, na.rm = TRUE))))
# Evaluation of the current basic structure
C <- cova(model[i, 1], h, b[i])
k0 <- k0 + kronecker(sill[ , ,i], C)
}
# Remove rows corresponding to missing data
if (length(index_missing) > 0){
k0 <- k0[-index_missing, ]
}
return (k0)
}

#' Set up the matrix to transform the Cartesian coordinates to coordinates that account for angles and anisotropy
#' @noRd

setrot <- function(model, it){

deg2rad <- pi / 180
ranges <- model[it, 2:4]
angles <- model[it, 5:7]
eps <- 2.2204e-016
# Matrix of coordinate reduction
redmat <- diag(1 / (eps+ranges))
cosa <- cos((90 - angles[1]) * deg2rad)
sina <- sin((90 - angles[1]) * deg2rad)
cosb <- cos(-angles[2] * deg2rad)
sinb <- sin(-angles[2] * deg2rad)
cosc <- cos(angles[3] * deg2rad)
sinc <- sin(angles[3] * deg2rad)
rotmat <- matrix(0,3,3)
rotmat[1, 1] <- cosb * cosa
rotmat[1, 2] <- cosb * sina
rotmat[1, 3] <- -sinb
rotmat[2, 1] <- -cosc * sina + sinc * sinb * cosa
rotmat[2, 2] <- cosc * cosa + sinc * sinb * sina
rotmat[2, 3] <- sinc * cosb
rotmat[3, 1] <- sinc * sina + cosc * sinb * cosa
rotmat[3, 2] <- -sinc * cosa + cosc * sinb * sina
rotmat[3, 3] <- cosc * cosb
rotred_matrix <- solve(rotmat) %*% redmat
return(rotred_matrix)
}

#' Build template of super-blocks centered at the block containing the node to simulate
#' @noRd

picksupr <- function(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,search_rotationmatrix){
# Main Loop over all possible super blocks
ixsbtosr<-NaN
iysbtosr<-NaN
izsbtosr<-NaN
for (i in (-nxsup + 1):(nxsup - 1)){
for (j in (-nysup + 1):(nysup - 1)){
for (k in (-nzsup + 1):(nzsup - 1)){
xo <- i * xsizsup 
yo <- j * ysizsup
zo <- k * zsizsup
coord <- c(xo,yo,zo) 
# Calculate the closest distance between the corners of the super block and the block at the origin 
distance <- matrix(NaN, nrow = 27, ncol = 3)
distance[ ,1] <- coord[1] + kronecker(matrix(1,nrow = 9,ncol=1), c(-1:1)) * xsizsup 
distance[ ,2] <- coord[2] + kronecker(matrix(1,nrow = 3,ncol=1), kronecker(c(-1:1), matrix(1,nrow = 3,ncol = 1))) * ysizsup
distance[ ,3] <- coord[3] + kronecker(c(-1:1), matrix(1,nrow = 9,ncol = 1)) * zsizsup   
distance <- distance %*% search_rotationmatrix 
distance <- apply(distance^2, 1, sum) 
# Keep this super block if it is close enough:
if (min(distance) < 1){ 
ixsbtosr <- c(ixsbtosr,i)
iysbtosr <- c(iysbtosr,j) 
izsbtosr <- c(izsbtosr,k)
}
} 
} 
}
return(list('ixsbtosr' = ixsbtosr[2:length(ixsbtosr)],'iysbtosr' = iysbtosr[2:length(iysbtosr)],'izsbtosr' = izsbtosr[2:length(izsbtosr)]))
}

#' Set up super-block strategy
#' @noRd

superblk <- function(datacoord,nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz){

# Default parameters
MAXSBX <- 21 
MAXSBY <- 21 
MAXSBZ <- 21 
ndata <- length(datacoord[,1])

# Establish the number and size of the super blocks
nxsup <- min(MAXSBX,nx) 
nysup <- min(MAXSBY,ny) 
nzsup <- min(MAXSBZ,nz) 
xsizsup <- nx * xsiz / nxsup 
ysizsup <- ny * ysiz / nysup 
zsizsup <- nz * zsiz / nzsup 
xmnsup <- (xmn-0.5 * xsiz) + 0.5 * xsizsup 
ymnsup <- (ymn-0.5 * ysiz) + 0.5 * ysizsup 
zmnsup <- (zmn-0.5 * zsiz) + 0.5 * zsizsup 

# Assign the data to a super block

ix_min <- apply(matrix(c(rep(nxsup,ndata),floor((datacoord[,1] - xmnsup) / xsizsup + 1.5)), nrow = ndata, ncol = 2, byrow = FALSE), 1, FUN = min,na.rm = TRUE)
ix <- array(apply(matrix(c(rep(1,ndata),ix_min),nrow = ndata,ncol = 2,byrow = FALSE), 1, FUN = max, na.rm = TRUE), dim = c(ndata,1))

iy_min <- apply(matrix(c(rep(nysup,ndata),floor((datacoord[,2] - ymnsup) / ysizsup + 1.5)), nrow = ndata, ncol = 2, byrow = FALSE), 1, FUN = min,na.rm = TRUE)
iy <- array(apply(matrix(c(rep(1,ndata),iy_min), nrow = ndata, ncol = 2, byrow = FALSE), 1, FUN = max, na.rm = TRUE), dim = c(ndata,1)) 
 
iz_min <- apply(matrix(c(rep(nzsup,ndata),floor((datacoord[,3] - zmnsup) / zsizsup + 1.5)), nrow = ndata, ncol = 2,byrow = FALSE), 1, FUN = min ,na.rm = TRUE)
iz <- array(apply(matrix(c(rep(1,ndata),iz_min),nrow = ndata, ncol = 2, byrow = FALSE), 1, FUN = max, na.rm = TRUE), dim = c(ndata,1))
index <- ix + (iy - 1) * nxsup + (iz - 1) * nxsup * nysup 

# Accumulate how many data belong to each super block
nisb <- matrix(0, nrow = nxsup * nysup * nzsup, ncol =1 ) 
for (i in 1:ndata){
nisb[index[i]] <- nisb[index[i]] + 1 
}
nisb <- array(c(0,cumsum(nisb)),dim = c(length(nisb) + 1,1)) 

# Sort by ascending super block number
tmp <- array(sort(index),dim = c(length(index),1)) 
I <- array(order(index),dim = c(length(index),1))
return(list('I' = I,'nisb' = nisb,'nxsup' = nxsup,'xmnsup' = xmnsup,'xsizsup' = xsizsup,'nysup' = nysup,'ymnsup' = ymnsup,'ysizsup' = ysizsup,'nzsup' = nzsup,'zmnsup' = zmnsup,'zsizsup' = zsizsup))

}

#' Converts Gaussian values into categorical values
#' @noRd

cc_truncate<-function(y_simu,nfield,flag,nthres,thresholds){
m1 <- nrow(y_simu)
m2 <- ncol(y_simu)
y_simu <- array(y_simu,dim = c(m1,nfield,m2/nfield))
y_simu<-aperm(y_simu,c(1,3,2))
sthres <- c(0,cumsum(nthres))
pthres <- c(1,cumprod(nthres+1))
i <- 1
for (k in 1:sthres[nfield + 1]){
  j <- sum(k > sthres)
  I <- matrix(0,nrow = nrow(y_simu),ncol = ncol(y_simu))
  I[which(y_simu[,,j] > thresholds[k])] <- 1
  i <- i + pthres[j] * I 
}
y_simu <- flag[i]
return(y_simu)
}
