#' Calculate the thresholds associated with a given truncation rule and given facies proportions when the number of GRFs are independent.
#'
#' Calculate Thresholds for Truncation Rule and Facies Proportions
#'
#' @param flag A vector with category numbers codifying the truncation rule.
#' @param nthres A vector with the number of thresholds for each GRF (1 x nfield).
#' @param proportions A vector with the proportions for each facies.
#' @param ymin A value representing negative infinity (-Inf). The 'ymin' value is updated recursively when the function is used iteratively.
#' @param ymax A value representing positive infinity (Inf). The 'ymax' value is updated recursively when the function is used iteratively.
#' @return A vector with the calculated thresholds.
#'
#' @details Use calculate_thresholds(flag, nthres, proportions) and do not define the arguments 'ymin' and 'ymax', which are used recursively.
#'
#' @references This function is a translation of MATLAB code written by Emery (2007).
#' - Emery, X. (2007). Simulation of geological domains using the plurigaussian model: New developments and computer programs. Computers & Geosciences, 33(9), 1189–1201. DOI: 10.1016/j.cageo.2007.01.006
#'
#' @examples
#' # Example usage
#' calculate.thresholds(flag = c(2,1,3,2,1,5,2,4,4), nthres = c(2,2), 
#' proportions = c(0.426, 0.253, 0.047, 0.205, 0.069))

calculate.thresholds <- function(flag,nthres,proportions,ymin,ymax){
nfield <- length(nthres) 
cthres <- c(0,cumsum(nthres)) 
# if (length(nthres) > 2){  
# flag <- array(flag,c(nthres+1))
# } else if (length(nthres) == 2){
# flag <- matrix(flag,nrow = nthres[1] + 1, ncol = nthres[2] + 1)
# } else if(length(nthres) == 1){
# flag <- matrix(flag,nrow = 1, ncol = nthres[1] + 1)
# }
  # Determine initial shape of flag array
  if (nfield > 1) {
    flag <- array(flag, dim = c(nthres + 1))
  } else {
    flag <- as.array(flag) # or flag <- as.vector(flag) if 1D is expected
  }

if (missing(ymin) && missing(ymax)){
ymin <- rep(-Inf, nfield) 
ymax <- rep(Inf, nfield) 
}
# Define an ordering of the Gaussian fields
ind <- matrix(c(1:nfield),nrow = length(1:nfield),ncol = 1) 
order_ind <- c();
for (i in 1:nfield){ 
order_ind <- cbind(order_ind ,ind)
ind <- matrix(c(ind,ind[1]),ncol = 1)
ind <- ind[-1,]
}
# Scan the thresholds corresponding to the k-th Gaussian field
for (k in 1:nfield){ 

# Order inputs so that the k-th Gaussian field comes first
if (nfield > 1){
ordered_flag <- aperm(flag,order_ind[k, ])
} else{
ordered_flag <- matrix(flag,nrow = length(flag),ncol = 1) 
} 
ordered_nthres <- nthres[order_ind[k, ]] 
ordered_ymin <- ymin[order_ind[k, ]] 
ordered_ymax <- ymax[order_ind[k, ]] 
prop <- proportions 

# Loop over the thresholds
if (ordered_nthres[1] > 0){
for (i in 1:ordered_nthres[1]){ 
# flag1 <- ordered_flag[1:i,,, drop = TRUE]
# flag2 <- ordered_flag[(i + 1):(ordered_nthres[1] + 1),,, drop = TRUE]
if (length(dim(ordered_flag)) >= 3) {
        flag1 <- ordered_flag[1:i, , , drop = FALSE]
        flag2 <- ordered_flag[(i + 1):(ordered_nthres[1] + 1), , , drop = FALSE]
      } else if (length(dim(ordered_flag)) == 2) {
        flag1 <- ordered_flag[1:i, , drop = FALSE]
        flag2 <- ordered_flag[(i + 1):(ordered_nthres[1] + 1), , drop = FALSE]
      } else {
        flag1 <- ordered_flag[1:i]
        flag2 <- ordered_flag[(i + 1):(ordered_nthres[1] + 1)]
      }

# Divide the flag into two sub-flags
# if (!is.na(dim(flag)[3])){
# flag1 <- array(NaN, c(1:i,dim(ordered_flag)[2],dim(ordered_flag)[3]))

# if(is.na(dim(ordered_flag)[3])){
# flag2 <- array(NaN, c(length(c((i + 1):(ordered_nthres[1] + 1 ))),dim(ordered_flag)[1]))
# } else{
# flag2 <- array(NaN, c(length(c((i + 1):(ordered_nthres[1] + 1 ))),dim(flag)[2],dim(ordered_flag)[3]))
# }
# for (j in 1:dim(ordered_flag)[3]){
# flag1[,,j] <- array(ordered_flag[i,,j])
# flag2[,,j] <- array(ordered_flag[c((i + 1):(ordered_nthres[1] +1 )),,j])
# }
# } else {
# flag1 <- ordered_flag[1:i, ];
# flag2 <- ordered_flag[(i+1):(ordered_nthres[1] + 1), ];  
# }
common <- 0
for (j in 1:length(flag1)){  
common <- common + length(which(flag1[j] == flag2)) 
}
if (common < 1){ #flag1 and flag2 have no common element (successful grouping)
# Determine the threshold that separates flag1 and flag2
n <- length(flag1)  
cumprop <- 0
for (j in 1:n){     
cumprop <- cumprop + prop[flag1[j]]; 
prop[flag1[j]] <- 0 
}

product <- prod(pnorm(ordered_ymax[2:nfield]) - pnorm(ordered_ymin[2:nfield]))
if(is.na(product)){
cumprop <- pnorm(ordered_ymin[1]) + cumprop 
}else{
cumprop <- pnorm(ordered_ymin[1]) + cumprop / product 
}
if (is.na(cumprop)){

yi <- -Inf
}else{
if (cumprop < 0) {
yi<--Inf
}else if(cumprop > 1){
yi<-Inf
}else{
yi <- qnorm(cumprop,mean = 0, sd = 1, lower.tail = TRUE)  
}
}
# Use the function recursively to determine the other thresholds
if (nfield > 1){
nthres1 <- c(i-1, ordered_nthres[2:nfield]) 
nthres2 <- c(ordered_nthres[1]-i, ordered_nthres[2:nfield]) 
} else {
nthres1 <- 0
nthres2 <- c(ordered_nthres[1]-i)
}    
sumthres1 <- sum(nthres1)     
sumthres2 <- sum(nthres2)     
ymin1 <- ordered_ymin      
ymin2 <- ordered_ymin     
ymin2[1] <- yi      
ymax1 <- ordered_ymax     
ymax1[1] <- yi      
ymax2 <- ordered_ymax     
thresholds1 <- calculate.thresholds(flag1,nthres1,proportions,ymin1,ymax1)	  
thresholds2 <- calculate.thresholds(flag2,nthres2,proportions,ymin2,ymax2)


t1 <- NA
t2 <- NA
if(sumthres1 >=  (nthres1[1] + 1)){
t1 <- thresholds1[(nthres1[1] + 1):sumthres1]
}
if(sumthres2 >=  (nthres2[1] + 1)){
t2 <- thresholds2[((nthres2[1] + 1):sumthres2)]
}
t_matrix<-rbind(t1,t2)
tmax<-(apply(t_matrix,2,max))


if (nthres1[1] >= 1){
temp_thresholds1 <- c(thresholds1[1:nthres1[1]])
} else {
temp_thresholds1 <- c()
}

if (nthres2[1] >= 1){
temp_thresholds2 <- c(thresholds2[1:nthres2[1]])
}else{
temp_thresholds2<-c()
}

thresholds <-rep(NaN,length(c(temp_thresholds1,temp_thresholds2,yi,tmax)))

indices1 <-c((cthres[k] + 1):cthres[nfield + 1])
if (k == 1){
indices2 <- c()
} else{
indices2 <- (2:cthres[k])
}
thresholds[unique(c(indices1, indices2))] <- c(temp_thresholds1, yi, temp_thresholds2,tmax)

return(thresholds)
}
}
}
}
 thresholds <- rep(-Inf,length(nthres))
}
