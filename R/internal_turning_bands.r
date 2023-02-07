#' Generation of equidistributed lines over the unit 3D sphere according to a van der Corput sequence
#' @noRd

vdc <- function(nlines, nrealiz, seed){
# This function uses the following functions: setrot
vdc_lines <- matrix(0, nrow = nlines * nrealiz, ncol = 3)
# rand('state',seed);
seed2 <- ceiling(1e7 * runif(1))
# randn('state',seed2);
i <- matrix(c(1:nlines), ncol = 1)
# binary decomposition of i
jj <- i
uu <- 0
pp <- 0
while (max(jj) > 0){
pp <- pp + 1
tt <- trunc(jj  / 2)
uu <- uu+2*(jj/2-tt)/(2^pp)
jj <- tt
}
# ternary decomposition of i
jj <- i
vv <- 0
pp <- 0
while (max(jj) > 0){
pp <- pp+1
tt <- trunc(jj / 3)
vv<- vv + 3 * (jj / 3 - tt) / 3 ^ pp
jj<-tt
}
# directing vector of the i-th line
x <- matrix(NaN,nrow=length(i),ncol=3)
x[, 1] <- cos(2 * pi * uu) * sqrt(1 - vv * vv)
x[, 2] <- sin(2 * pi * uu) * sqrt(1 - vv * vv)
x[, 3] <- vv
# random rotation
for (k in 1:nrealiz){
angles<-360 * runif(3)
R<-setrot(matrix(c(1, 1, 1, 1, angles), nrow = 1), 1)
vdc_lines[((k - 1) * nlines + 1):(k * nlines), ]<- x %*% R
}
return(vdc_lines)
}

# Prepare the lines for turning bands
#' @noRd

prepare_lines <- function(model, nrealiz, nlines, nst, extreme_coord, seed_vdc,b){
# This function uses the following functions: setrot, vdc

max_intervalnumber <- 1e4
# Initialization
model_rotationmatrix <- array(NaN, dim = c(3, 3, nst))
all_lines <- array(NaN, dim = c(nrealiz * nlines[1], 3, nst))
all_offset <- matrix(NaN, nrow = (nrealiz * nlines[1]), ncol = nst)
all_r <- array(NaN, dim = c(nlines[1], nrealiz, nst))
all_phi <- array(NaN, dim = c(nlines[1], nrealiz, nst))
all_theta <- array(NaN, dim = c(nlines[1], nrealiz, nst))
valid_lines <- matrix(NaN, nrow = (nrealiz * nlines[1]), ncol = nst)
cc_sigma <- rep(NaN, nst)
for (i in 1:nst){

# Line generation (Van der Corput sequence)
V <- vdc(nlines[i], nrealiz, seed_vdc[i])
R <- setrot(model, i)  # R = setrot(model,i);
model_rotationmatrix[ , , i] <- R
glines <- V %*% t(R)
all_lines[1:(nrealiz * nlines[i]), , i] <- glines
nottoosmall <- rep(NaN,nrow(glines))
# Spherical covariance model?
if (model[i, 1] == 1){
x <- extreme_coord %*% t(glines)
cc_sigma[i] <- 2 * sqrt(3 / nlines[i])
all_offset[1:(nrealiz * nlines[i]), i]<- t(t(apply(x, 2, min) - runif(nrealiz * nlines[i])))
# Identify the lines for which the scale factor is too small
# The simulation along these lines will be replaced by a nugget effect
interval_number <- t(t(apply(x, 2, max) - apply(x, 2, min)))
nottoosmall[which(interval_number < max_intervalnumber)] <- 1
nottoosmall[-which(interval_number < max_intervalnumber)] <- 0
valid_lines[1:(nrealiz * nlines[i]), i] <- nottoosmall
# Exponential covariance model?
} else if(model[i, 1] == 2){
v <- matrix(rnorm(6 * nlines[i] * nrealiz), nrow = 6,ncol = nlines[i] * nrealiz)
w <- 3 * matrix(runif(nlines[i] * nrealiz), nrow = 1, ncol = nlines[i] * nrealiz)
v[5, which(w < 1)] <- 0
v[6, which(w < 1)] <- 0
G <- 0.5 * t(t(apply(v ^ 2, 2, sum))) %*% matrix(1, nrow = 1, ncol = 3)
glines <- glines / G
x <- extreme_coord %*% t(glines)
cc_sigma[i] <- 2 * sqrt(3 / nlines[i])
all_lines[1:(nrealiz * nlines[i]), , i] <- glines
all_offset[1:(nrealiz * nlines[i]), i] <- t(t(apply(x, 2, min) - runif(nrealiz * nlines[i])))
# Identify the lines for which the scale factor is too small
# The simulation along these lines will be replaced by a nugget effect
interval_number <- t(t(apply(x, 2, max) - apply(x, 2, min)))
nottoosmall[which(interval_number < max_intervalnumber)] <- 1
nottoosmall[-which(interval_number < max_intervalnumber)] <- 0
valid_lines[1:(nrealiz * nlines[i]), i] <- nottoosmall
# Gamma covariance model?
} else if (model[i, 1] == 3){
v <- matrix(rnorm(6 * nlines[i] * nrealiz), nrow = 6, ncol = nlines[i] * nrealiz)
w <- 3 * matrix(runif(nlines[i] * nrealiz), nrow = 1, ncol = nlines[i] * nrealiz)
v[5, which(w < 1)] <- 0
v[6, which(w < 1)] <- 0
G <- 0.5 * (t(t(apply(v ^ 2, 2, sum))) / matrix(rgamma(nlines[i] * nrealiz, scale = b[i], shape = 1),nrow = nlines[i] * nrealiz, ncol=1)) %*% matrix(1, nrow = 1, ncol = 3)
glines <- glines / G
x <- extreme_coord %*% t(glines)
cc_sigma[i] <- 2 * sqrt(3 / nlines[i])
all_lines[1:(nrealiz * nlines[i]), , i] <- glines
all_offset[1:(nrealiz * nlines[i]), i] <- t(t(apply(x, 2, min) - runif(nrealiz * nlines[i])))
# Identify the lines for which the scale factor is too small
# The simulation along these lines will be replaced by a nugget effect
interval_number <- t(t(apply(x, 2, max) - apply(x, 2, min)))
nottoosmall[which(interval_number < max_intervalnumber)] <- 1
nottoosmall[-which(interval_number < max_intervalnumber)] <- 0
valid_lines[1:(nrealiz * nlines[i]), i] <- nottoosmall
# Stable covariance model with parameter <1?
} else if (model[i, 1] == 4){ # elseif (model(i,1) == 4)
v <- matrix(rnorm(6 * nlines[i] * nrealiz), nrow = 6, ncol = nlines[i] * nrealiz)
w <- 3 * matrix(runif(nlines[i] * nrealiz), nrow = 1, ncol = nlines[i] * nrealiz)
v[5, which(w < 1)] <- 0
v[6, which(w < 1)] <- 0
e <- -log(matrix(runif(nlines[i] * nrealiz), nrow = nlines[i] * nrealiz, ncol = 1))
f <- pi * matrix(runif(nlines[i] * nrealiz), nrow = nlines[i] * nrealiz, ncol = 1) - pi / 2
stable <- abs(sin(b[i] * (f - pi / 2)) / (cos(f) ^ (1 / b[i])) * (cos(f - b[i] * (f - pi / 2)) / e) ^ ((1 - b[i]) / b[i]))
G <- 0.5 * (t(t(apply(v ^ 2, 2, sum))) / stable) %*% matrix(1, nrow = 1, ncol = 3)
glines <- glines / G
x <- extreme_coord %*% t(glines)
cc_sigma[i] <- 2 * sqrt(3 / nlines[i])
all_lines[1:(nrealiz * nlines[i]), , i] <- glines
all_offset[1:(nrealiz * nlines[i]), i] <- t(t(apply(x, 2, min) - runif(nrealiz * nlines[i])))
# Identify the lines for which the scale factor is too small
# The simulation along these lines will be replaced by a nugget effect
interval_number <- t(t(apply(x, 2, max) - apply(x, 2, min)))
nottoosmall[which(interval_number < max_intervalnumber)] <- 1
nottoosmall[-which(interval_number < max_intervalnumber)] <- 0
valid_lines[1:(nrealiz * nlines[i]), i] <- nottoosmall
# Stable model with parameter >1?
} else if (model[i, 1] == 4.5){
cc_sigma[i] <- sqrt(2 / nlines[i])
tt <- matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2 + matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2 + matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2
all_r[1:nlines[i], 1:nrealiz, i] <- sqrt(6 * tt)
e <- -log(matrix(runif(nlines[i]*nrealiz),nrow=nlines[i]*nrealiz,ncol=1))
f <- pi * matrix(runif(nlines[i]*nrealiz),nrow=nlines[i]*nrealiz,ncol=1) - pi/2
G <- abs( sin(b[i] /2 * (f - pi / 2))/(cos(f) ^ (2 / b[i])) * (cos(f - b[i] / 2 * (f - pi / 2)) / e) ^ ((1 - b[i] / 2) / b[i] * 2) ) %*% matrix(1, nrow = 1, ncol = 3)
glines <- glines %*% sqrt(G / 3)
all_lines[1:(nrealiz * nlines[i]), , i] <- glines
all_phi[1:nlines[i], 1:nrealiz, i] <- 2 * pi * matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)
all_theta[1:nlines[i], 1:nrealiz, i] <- sqrt(-log(matrix(runif(nlines[i] * nrealiz),nrow = nlines[i], ncol = nrealiz)))
# Cubic covariance model?
} else if (model[i, 1] == 5){
x <- extreme_coord %*% t(glines)
cc_sigma[i] <- 2 * sqrt(210 / nlines[i])
all_lines[1:(nrealiz * nlines[i]), , i] <- glines
all_offset[1:(nrealiz * nlines[i]), i] <- t(t(apply(x, 2, min) - runif(nrealiz * nlines[i])))
# Identify the lines for which the scale factor is too small
# The simulation along these lines will be replaced by a nugget effect
interval_number <- t(t(apply(x, 2, max) - apply(x, 2, min)))
nottoosmall[which(interval_number < max_intervalnumber)] <- 1
nottoosmall[-which(interval_number < max_intervalnumber)] <- 0
valid_lines[1:(nrealiz * nlines[i]), i] <- nottoosmall
# Gaussian covariance model?
} else if (model[i, 1] == 6){
cc_sigma[i] <- sqrt(2 / nlines[i])
tt <- matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2 + matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2 + matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2
all_r[1:nlines[i], 1:nrealiz, i] <- sqrt(2 * tt)
all_phi[1:nlines[i], 1:nrealiz, i] <- 2 * pi * matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)
all_theta[1:nlines[i], 1:nrealiz, i] <- sqrt(-log(matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)))
# Cardinal sine covariance model?
} else if (model[i, 1] == 7){
cc_sigma[i] <- sqrt(2 / nlines[i])
tt <- 2 * matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) - 1
all_r[1:nlines[i], 1:nrealiz, i] <- sign(tt)
all_phi[1:nlines[i], 1:nrealiz, i] = 2 * pi * matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)
all_theta[1:nlines[i], 1:nrealiz, i] = sqrt(-log(matrix(runif(nlines[i] * nrealiz),nrow = nlines[i],ncol = nrealiz)))
# Bessel-J covariance model?
} else if (model[i, 1] == 8){
cc_sigma[i] <- sqrt(2 / nlines[i])
tt <- matrix(rbeta(nlines[i] * nrealiz, 1.5, b[i] - 0.5), nrow = nlines[i], ncol = nrealiz)
all_r[1:nlines[i], 1:nrealiz, i] <- sqrt(tt)
all_phi[1:nlines[i], 1:nrealiz, i] <- 2 * pi * matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)
all_theta[1:nlines[i],1:nrealiz, i] <- sqrt(-log(matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)))
# Bessel-K covariance model (b<0.5)?
} else if (model[i, 1] == 9){
v <- matrix(rnorm(6 * nlines[i] * nrealiz), nrow = 6, ncol = nlines[i] * nrealiz)
w <- 3 * matrix(runif(nlines[i] * nrealiz), nrow = 1, ncol = nlines[i] * nrealiz)
v[5, which(w < 1)] <- 0
v[6, which(w < 1)] <- 0
G <- 0.5 * ((t(t(apply(v ^ 2, 2, sum)))) * sqrt( matrix(rbeta(nlines[i] * nrealiz, b[i], 0.5 - b[i]), nrow = nlines[i] * nrealiz, ncol = 1))) %*% matrix(1, nrow = 1, ncol = 3)
glines <- glines / G
x <- extreme_coord %*% t(glines)
cc_sigma[i] <- 2 * sqrt(3 / nlines[i])
all_lines[1:(nrealiz * nlines[i]), , i] <- glines
all_offset[1:(nrealiz * nlines[i]), i] <- t(t(apply(x, 2, min) - runif(nrealiz * nlines[i])))
# Identify the lines for which the scale factor is too small
# The simulation along these lines will be replaced by a nugget effect
interval_number <- t(t(apply(x, 2, max) - apply(x, 2, min)))
nottoosmall[which(interval_number < max_intervalnumber)] <- 1
nottoosmall[-which(interval_number < max_intervalnumber)] <- 0
valid_lines[1:(nrealiz * nlines[i]), i] <- nottoosmall
# Bessel-K covariance model (b>0.5)?
} else if( model[i ,1] == 9.5){
cc_sigma[i] <- sqrt(2 / nlines[i])
tt <- matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2 + matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2 + matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2
all_r[1:nlines[i], 1:nrealiz, i] <- sqrt(6 * tt)
G <- matrix(rgamma(nlines[i] * nrealiz, scale = b[i], shape = 1), nrow = nlines[i] * nrealiz, ncol = 1) %*% matrix(1, nrow = 1, ncol = 3)
	glines <- glines / sqrt(G * 12)
all_lines[1:(nrealiz * nlines[i]), , i] <- glines
all_phi[1:nlines[i], 1:nrealiz, i] = 2 * pi * matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)
all_theta[1:nlines[i], 1:nrealiz, i] = sqrt(-log(matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)))
# Generalized Cauchy model?
} else if( model[i, 1] == 10){
cc_sigma[i] <- sqrt(2 / nlines[i])
tt <- matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2 + matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2 + matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2
all_r[1:nlines[i], 1:nrealiz, i] <- sqrt(6 * tt)
G <- matrix(rgamma(nlines[i] * nrealiz, scale = b[i], shape = 1), nrow = nlines[i] * nrealiz, ncol=1) %*% matrix(1, nrow = 1, ncol = 3)
all_lines[1:(nrealiz * nlines[i]), , i] <- glines
all_phi[1:nlines[i], 1:nrealiz, i] <- 2 * pi * matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)
all_theta[1:nlines[i], 1:nrealiz, i] <- sqrt(-log(matrix(runif(nlines[i] * nrealiz),nrow = nlines[i] ,ncol = nrealiz)))
# Exponential sine model?
} else if(model[i,1] == 11){
cc_sigma[i] <- sqrt(2 / nlines[i])
t2 <- matrix(NaN, nrow = nlines[i], ncol = nrealiz)
for (k in 1:nrealiz){
for (j in 1:nlines[i]){
iflag <- -1
while (iflag < 0){
t2[j, k] <- rgamma(1, scale = 0.5, shape = 1)
u <- runif(1)
c_expmodel <- 1
puis <- 1
n <- 0
iflag <- 0
while (iflag == 0){
n <- n + 1
puis <- puis * pi * pi / 4 / (2 * n * (2 * n - 1))
if (n / 2 == round(n / 2)){
c_expmodel <- c_expmodel + puis * exp(-4 * t2[j, k] * n * (n + 1))
if (u > c_expmodel){
iflag <- -1
} else{
c_expmodel <- c_expmodel - puis * exp(-4 * t2[j, k] * n * (n + 1))
if (u < c_expmodel){
iflag <- 1
}
}
}
}
}
}
}
t1 <- matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2 + matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2 + matrix(rnorm(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz) ^ 2
all_r[1:nlines[i], 1:nrealiz, i] <- sqrt(t1 / 2 / t2)
all_phi[1:nlines[i], 1:nrealiz, i] <- 2 * pi * matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)
all_theta[1:nlines[i], 1:nrealiz, i] <- sqrt(-log(matrix(runif(nlines[i] * nrealiz), nrow = nlines[i], ncol = nrealiz)))
}
}
return(list(model_rotationmatrix = model_rotationmatrix, cc_sigma = cc_sigma, all_r = all_r, all_phi = all_phi, all_theta = all_theta, all_lines = all_lines, all_offset = all_offset, valid_lines = valid_lines))
}
