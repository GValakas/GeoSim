#' Modeling Indicator Variograms
#'
#' The function helps to determine the coefficients of the corgionalization model. In the case of Plurigaussian moedling, the model indicator variograms and cross-variograms of facies are calculated to examine the fitting of the underlying GRFs variograms. 
#' In the case of co-simulation, the model cross variograms of facies and Normal scores are calculated. The calculation of simple and cross-variograms is performed by Monte Carlo simulation.
#'
#' @param truncation_rule A list of arguments to determine the truncation rule of Plurigaussian. Use list(nfield, nthres, thresholds, flag),
#' where:
#'   - nfield: Number of GRFs used for the truncation rule.
#'   - nthres: A vector with the number of thresholds for each GRF (1 x nfield).
#'   - thresholds: A vector of the thresholds for all GRFs (1 x sum(nthres)).
#'   - flag: A vector with category numbers codifying the truncation rule.
#'
#' @param variog_model A list of arguments to define the nested variogram theoretical models. Use list(model, cc, b, nugget),
#' where:
#'   - model: A matrix containing the covariance model for the GRFs (nested structures x 7 matrix). Each row corresponds to a nested structure and is codified as: type, scale factors, angles. Use the codes of the available types of variogram models (see details). There are three scale factors (along the rotated NS, EW, and vertical axes) and three angles to define the coordinate rotation (azimuth, dip, and plunge). See Deutsch & Journel (1997, p. 25) for more information.
#'   - cc: A matrix indicating the sills of nested structures (nested structures x (1 + nfield)^2).
#'   - b: A column vector with the additional parameters required for specific covariance types (nested structures x 1). See details for specific covariance types and their corresponding requirements.
#'   - nugget: A row vector with the nugget effect variance-covariance matrix of size (1 + nfield)^2.
#'
#' @param h_lag The lag distance (single value)
#'
#' @param nlag The number of lags (single value)
#'
#' @param azm The azimuth (single value)
#'
#' @param dip The dip (single value)
#'
#' @param nMC The number of Monte Carlo realizations (single value)
#'
#' @return In the case of Plurigaussian modeling the function returns an array of variograms (simple and cross) of facies indicators (ncategory * ncategory * nlag matrix). In the case of co-simulation modeling, the function returns cross-variograms between Gaussian field Z0 and facies indicators (ncategory * nlag matrix)
#'
#' @details To see the available types of variogram model see co.sim function. In the case of co-simulation, the number of rows of sills matrix differs from the number of GRFs of Plurigaussian modeling.  
#'
#' @examples
#' truncation_rule <- list(nfield = 2 , nthres = c(2,2),thresholds = 
#' c(-0.5769, 1.0762, -0.4895, 0.6592), flag = c(2,1,3,2,1,5,2,4,4))
#' # Gaussian Random Fields Modeling
#' nst_str <- 2
#' nugget <- c(0,0,0,0)
#' B_nst1 <- matrix(c(0.9986,-0.1,-0.1,0.8167),nst_str,truncation_rule$nfield)
#' B_nst2 <- matrix(c(0.00149,0,0,0.1833),nst_str,truncation_rule$nfield)
#' sills_matrix <- matrix(c(B_nst1,B_nst2),nrow=nst_str,ncol=truncation_rule$nfield^2,byrow=TRUE)
#' variogram_models <- matrix(c(1,400,400,140,0,0,0,2,300,300,90,0,0,0),nst_str,7,byrow=TRUE)			
#' variog_model <- list(model=variogram_models,cc=sills_matrix,b=matrix(0,nrow=nst_str,ncol=1),
#' nugget=nugget)
#' azm <- 0
#' dip <- 0
#' h_lag <- 100
#' nlag <- 10
#' nMC <- 1e4;

rvmodel <- function(truncation_rule,variog_model,h_lag,nlag,azm,dip,nMC){

# Define parameters
nfield <- truncation_rule$nfield
thresholds <- truncation_rule$thresholds
nthres <- truncation_rule$nthres
flag <- truncation_rule$flag
model <- variog_model$model
nst <- nrow(model) # number of nested structures
cc <-variog_model$cc
nvar <- nfield
b <- matrix(variog_model$b)
nugget <- variog_model$nugget
nfield <- length(nthres)

if (nfield < sqrt(ncol(cc))){
nvar <- 1 + nfield
} else {
nvar <- nfield
}

azm <- pi / 180 * azm;
dip <- pi / 180 * dip;
u <- c(sin(azm) * cos(dip), cos(azm) * cos(dip), sin(dip)) 
h <- matrix(c(0:nlag),ncol=1) %*% u * h_lag 
n <- length(flag)
ncategory <- max(flag)
nst <- nrow(model)
sthres <- c(0, cumsum(nthres))
pho<- matrix(0,nrow = (nlag + 1), ncol = (nvar * nvar))

# Check for model consistency
if (nfield < nvar){
chk_variogram <- check_variogram(model, nfield, nst, cc, nugget, cosim = 1) # check variogram
} else {
chk_variogram <- check_variogram(model, nfield, nst, cc, nugget) 
}

for (i in 1:nst){
R <- setrot(model,i)
ha<- h %*% R 
ha <- sqrt(colSums(t(ha)^2)) 
pho <- pho + cova(model[i,1], matrix(ha, ncol = 1),b[i]) %*% cc[i, ] 
}
pho[1, ] <- pho[1, ] + nugget

# Monte Carlo simulation 
variog_ind <- array(0, c(ncategory,ncategory,nlag))  
cross_variog <- matrix(0, nrow = ncategory, ncol = nlag)  
A1 = array(pho[1, ], c(nvar,nvar)) 
X <- matrix(rnorm(2 * nvar * nMC), nrow = 2 * nvar , ncol = nMC) 
for (k in 2:(nlag + 1)){  
Ak <- array(pho[k, ], c(nvar,nvar))  
A <- chol(rbind(cbind(A1,Ak),cbind(Ak,A1)))    
Y <- t(A) %*% X  
I1 <- cc_truncate(t(Y[2:nvar,]),nfield,flag,nthres,thresholds)   
Ik <- cc_truncate(t(Y[(nvar + 2):(2 * nvar), ]),nfield,flag,nthres,thresholds)  
Y1 <- Y[1, ]    
Yk <- Y[nvar + 1, ]    #
 for (i in 1:ncategory){   
I1_tmp_i <- rep(0,length(I1))
Ik_tmp_i <- rep(0,length(Ik))
I1_tmp_i[which(I1 == i) ] <- 1
Ik_tmp_i[which(Ik == i) ] <- 1
cross_variog[i,k-1] <- 0.5 * mean((I1_tmp_i - Ik_tmp_i) * t((Y1-Yk)))     
for (j in 1:ncategory){     
I1_tmp_j <- rep(0,length(I1))
Ik_tmp_j <- rep(0,length(Ik))
I1_tmp_j[which(I1 == j)] <- 1
Ik_tmp_j[which(Ik == j)] <- 1
variog_ind[i,j,k-1] <-  0.5 * mean((I1_tmp_i - Ik_tmp_i) * (I1_tmp_j - Ik_tmp_j)) 
} 
} 
}   
if (nvar > nfield){
vmodel_results <- list(pg_variograms = variog_ind , cross_variograms = cross_variog)
return(vmodel_results)
} else {
pg_variograms <- variog_ind
return(pg_variograms)
}
}





