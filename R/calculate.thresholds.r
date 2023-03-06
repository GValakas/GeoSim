#' Calculate the thresholds associated with a given truncation rule and given facies proportions when the number of GRFs are independent.
#'
#' @param flag is a vector with category numbers codifying the truncation rule
#' @param nthres is a vector with the number of thresholds for each GRF (1  x nfield)
#' @param proportions is a vector with the proportions for each facies
#' @param ymin a value equals to -Inf. ymin is updated whren the function is used recursively
#' @param ymax a value equals to Inf. ymax is updated whren the function is used recursively
#' @return A vector with the thresholds 
#'
#' @details Use calculate_thresholds(flag,nthres,proportions) and do not define the arguments ymin and ymax that are used recursively

calculate.thresholds <- function(flag,nthres,proportions,ymin,ymax){
nfield <- length(nthres) 
cthres <- c(0,cumsum(nthres)) 
if (length(nthres) > 2){  
flag <- array(flag,c(nthres+1))
} else if (length(nthres) == 2){
flag <- matrix(flag,nrow = nthres[1] + 1, ncol = nthres[2] + 1)
} else if(length(nthres) == 1){
flag <- matrix(flag,nrow = 1, ncol = nthres[1] + 1)
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
# Divide the flag into two sub-flags
if (!is.na(dim(flag)[3])){
flag1 <- array(NaN, c(1:i,dim(flag)[2],dim(ordered_flag)[3]))

if(is.na(dim(ordered_flag)[3])){
flag2 <- array(NaN, c(length(c((i + 1):(ordered_nthres[1] + 1 ))),dim(flag)[1]))
} else{
flag2 <- array(NaN, c(length(c((i + 1):(ordered_nthres[1] + 1 ))),dim(flag)[2],dim(ordered_flag)[3]))
}
for (j in 1:dim(ordered_flag)[3]){
flag1[,,j] <- array(ordered_flag[i,,j])
flag2[,,j] <- array(ordered_flag[c((i + 1):(ordered_nthres[1] +1 )),,j])
}
} else {
flag1 <- ordered_flag[1:i, ];
flag2 <- ordered_flag[(i+1):(ordered_nthres[1] + 1), ];  
}
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
yi <- qnorm(cumprop)  
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

#thresholds2 <-c(1.0771,-0.4919,0.6591)
t1 <- NA
t2 <- NA
if(sumthres1 >=  (nthres1[1] + 1)){
t1 <- thresholds1[(nthres1[1] + 1):sumthres1]
}
if(sumthres2 >=  (nthres2[1] + 1)){
t2 <- thresholds2[((nthres2[1] + 1):sumthres2)]
}
t_matrix<-rbind(t1,t2)
if (nthres1[1] >= 1){
temp_thresholds <- c(thresholds1[1:nthres1[1]], yi)
} else {
temp_thresholds <- yi
}
if (nthres2[1] >= 1){
temp_thresholds <- c(temp_thresholds, thresholds2[1:nthres2[1]])
}

thresholds <- c(temp_thresholds,(sort(apply(t_matrix,2,max))))
if (all(flag2[1:length(flag2)]==flag2[length(flag2):1])){
index_thres <- c((cthres[k] + 1):cthres[nfield + 1],1:cthres[k])
thresholds <-thresholds[index_thres]
}

return(thresholds)
}
}
}
}
 thresholds <- rep(-Inf,length(nthres))
}


