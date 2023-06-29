#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]


struct add_multiple {
int incr;
int count;
add_multiple(int incr)
: incr(incr), count(0)
{}
inline int operator()(int d) {
return d + incr * count++;
}
};
/* using namespace Rcpp; */
using namespace arma;

//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector seq_rcpp(double from_, double to_, double by_ = 1.0) {
int adjust = std::pow(10, std::ceil(std::log10(10 / by_)) - 1);
int from = adjust * from_;
int to = adjust * to_;
int by = adjust * by_;
std::size_t n = ((to - from) / by) + 1;
Rcpp::IntegerVector res = Rcpp::rep(from, n);
add_multiple ftor(by);
std::transform(res.begin(), res.end(), res.begin(), ftor);
return Rcpp::NumericVector(res) / adjust;
}

//' @noRd
arma::mat cova_cpp(double it,arma::mat h,double b){

int epsilon = 1e-12;
arma::mat C(h.n_rows,h.n_cols);
if(it == 1){
arma::uvec idx = find(h > 1);
h.elem(idx).fill(1);
C = 1 - 1.5 * h + 0.5 * pow(h,3);
}else if (it == 2){
C = exp(-h);
}else if(it == 3){
C = 1/pow(1+h,b);
}else if(it == 4){
C = exp(-pow(h,b));
}else if(it == 5){
arma::mat h1(h.n_rows,h.n_cols);
h1 = h;
arma::uvec idx2 = find(h1 > 1);
h1.elem(idx2).fill(1);
C = 1 - 7 * pow(h1,2) + 8.75 * pow(h1,3) - 3.5 * pow(h1,5) + 0.75 * pow(h1,7);
}else if(it == 6){
C = exp(-pow(h,2));
}else if(it == 7){
C = sin(h + epsilon);
}
return(C);
}

//' @noRd
// [[Rcpp::export]]
arma::field<arma::mat> cokrige_cpp(arma::mat datacoord, arma::uvec index_missing, arma::mat coord, arma::mat model,
arma::cube sill, arma::mat b, arma::mat sillnugget, arma::cube model_rotationmatrix){

int n = datacoord.n_rows; 
int p = coord.n_rows; 
int nst = model.n_rows; 
int nvar = sill.n_rows; 
arma::mat D(n, n, fill::eye);
arma::mat k = kron(sillnugget,D);
arma::mat k0(nvar * n, nvar * p, fill::zeros);

arma::mat R;
arma::mat h;
arma::mat C;
arma::uvec idx;
/* 
arma::mat R(model_rotationmatrix.n_rows,model_rotationmatrix.n_cols);
arma::mat h(n,n);
arma::mat o1(1,n);
o1.ones(1,n);
arma::mat o2(n,1);
o2.ones(n,1);
arma::mat C(n,n); */
for (int i = 0; i < nst; ++i) {
R = model_rotationmatrix.slice(i); 
arma::mat tt = datacoord * R; 
tt = tt * trans(tt); 
h = -2 * tt + diagvec(tt) * ones<arma::rowvec>(n) + ones<arma::colvec>(n) * trans(diagvec(tt));
/* h = -2 * tt + tt.diag(0) * o1 + o2 * trans(tt.diag(0));  */
idx = find(h < 0);
/* arma::uvec idx = find(h < 0); */ 
h.elem(idx).fill(0);
h = sqrt(h); 
C = cova_cpp( model(i,0),h,b(i)); 
k += kron(sill.slice(i), C);
/* k = k + kron(sill.slice(i),C);  */

for (int jj = 0; jj < p; ++jj) {
h = (ones<arma::colvec>(n) * coord.row(jj) - datacoord) * R;
/* h = (o2 * coord.row(jj) - datacoord) * R; */ 
h = pow(h,2);
/* arma::mat h1(h.n_rows,1); */  
h = sum(h, 1);
h = sqrt(h);
/* h1 = sum(h,1); 
h1 = sqrt(h1); */ 
C = cova_cpp(model(i,0),h,b(i));
/* C = cova_cpp(model(i,0),h1,b(i)); */

/* arma::uvec colidx = regspace<arma::uvec>(jj, nvar * p - 1, p); */
Rcpp::NumericVector subsetCols = seq_rcpp(jj,nvar * p - 1,p);
arma::uvec colidx = Rcpp::as<arma::uvec>(subsetCols);
/* arma::uvec colidx = arma::shuffle(arma::linspace<arma::uvec>(jj, nvar * p - 1, p)); */ // Shuffle the indices
k0.cols(colidx) += kron(sill.slice(i), C);
/* k0.cols(colidx) = k0.cols(colidx) + kron(sill.slice(i),C); */ 
}
}
/* arma::uvec index_missing_shed = find((k.n_rows-1) < index_missing); */

if (!index_missing.is_empty() ){

k.shed_rows(index_missing);
k.shed_cols(index_missing);
k0.shed_rows(index_missing);
}

/* Co-kriging weights and covariances */
arma::mat cokr_weights = solve(k,k0,arma::solve_opts::fast);
arma::mat D2(p, p, fill::eye);
/* arma::mat D2(p,p);
D2.eye(p, p); */
arma::mat C0 = kron(sillnugget,D2);
arma::mat op1 = ones<arma::rowvec>(p);
/* arma::mat op1(1,p);
op1.ones(1,p); */
arma::mat op2 = ones<arma::colvec>(p);
/* arma::mat op2(p,1);
op2.ones(p,1); */
for (int i = 0; i < nst; ++i) {
R = model_rotationmatrix.slice(i);
arma::mat tt = coord * R ;
tt =  tt * trans(tt);
/* arma::mat tt2 =  tt*trans(tt); */
arma::mat h = -2 * tt + diagvec(tt) * op1 + op2 * trans(diagvec(tt));
idx = find(h < 0);
/* arma::uvec idx = find(h < 0); */
h.elem(idx).fill(0);
h = sqrt(h);
arma::mat C2 = cova_cpp(model(i,0),h,b(i));
C0 += kron(sill.slice(i), C2);
/* C0 = C0+kron(sill.slice(i),C2); */
}
arma::mat covariances = C0 - trans(cokr_weights) * k0;
arma::field<arma::mat> gb_cokrige(2,1);
gb_cokrige(0,0) = cokr_weights;
gb_cokrige(1,0) = covariances;
return(gb_cokrige);
}

//' @noRd
// [[Rcpp::export]]
int cc_truncate_cpp(arma::mat y_simu, int nfield, arma::vec flag, arma::vec nthres,arma::vec thresholds){

/* Converts Gaussian values into categorical values */
arma::rowvec sthres(nthres.n_elem + 1) ;
Rcpp::NumericVector idx_seq = seq_rcpp(1,(nthres.n_elem),1);
arma::uvec idx = Rcpp::as<arma::uvec>(idx_seq);
sthres(0) = -1;
sthres.elem(idx) = cumsum(nthres); 
arma::rowvec pthres(nthres.n_elem + 1);
pthres(0) = 1;
pthres.elem(idx)= cumprod(nthres + 1); 
int i = 0;
int n_sthres  = sthres.n_cols;
for (int k = 0; k < sthres(nfield); ++k) { 
int j = 0;
for (int kk = 0; kk < n_sthres; ++kk){
int v = k + 1;
if (v > sthres(kk)){
j = j + 1; 
}
}
int I = 0;
if (y_simu(0,j - 1) > thresholds(k)){

I = 1;
}
i = i + pthres(j-1) * I; 
}
int y_out = flag(i);
return(y_out);}

// [[Rcpp::export]]
arma::field<arma::mat>  cc_search_cpp(arma::mat datacoord, arma::mat datavalue, arma::rowvec coord, arma::mat search_rotationmatrix, 
int octant, int ndata,
 int nxsup, int nysup, int nzsup, int xmnsup, int ymnsup, int zmnsup, int xsizsup, int ysizsup, int zsizsup, 
 arma::rowvec ixsbtosr, arma::rowvec iysbtosr,  arma::rowvec izsbtosr, arma::colvec nisb ){

/* First selection of data according to super-block search */
/* Indices of the super-block containing the node to simulate */
arma::rowvec ix_vec(2);
ix_vec(0) = nxsup;
ix_vec(1) = floor((coord(0) - xmnsup) / xsizsup + 1.5);
int ix_min = min(ix_vec); 
ix_vec(0) = ix_min;
ix_vec(1) = 1;
int ix = max(ix_vec);  

arma::rowvec iy_vec(2);
iy_vec(0) = nysup;
iy_vec(1) = floor((coord(1) - ymnsup) / ysizsup + 1.5);
int iy_min = min(iy_vec); 
iy_vec(0) = iy_min;
iy_vec(1) = 1;
int iy = max(iy_vec); 
arma::rowvec iz_vec(2);
iz_vec(0) = nzsup;
iz_vec(1) = floor((coord(2) - zmnsup) / zsizsup + 1.5);
int iz_min = min(iz_vec);
iz_vec(0) = iz_min;
iz_vec(1) = 1;
int iz = max(iz_vec);

/* Indices of the acceptable super-blocks for searching the data */
arma::rowvec ix2 = ix + ixsbtosr; 
arma::rowvec iy2 = iy + iysbtosr; 
arma::rowvec iz2 = iz + izsbtosr; 

arma::vec index_ix = arma::conv_to<arma::vec>::from(find (ix2 >= 1));
arma::vec index_iy = arma::conv_to<arma::vec>::from(find (iy2 >= 1));
arma::vec index_iz = arma::conv_to<arma::vec>::from(find (iz2 >= 1));
arma::vec index2_ix = arma::conv_to<arma::vec>::from(find(ix2 <= nxsup));
arma::vec index2_iy = arma::conv_to<arma::vec>::from(find(iy2 <= nysup));
arma::vec index2_iz = arma::conv_to<arma::vec>::from(find(iz2 <= nzsup));

arma::vec index_ix_iy = intersect(index_ix,index_iy);
arma::vec index_ix_iy_iz = intersect(index_ix_iy,index_iz);

arma::vec index2_ix_iy = intersect(index2_ix,index2_iy);
arma::vec index2_ix_iy_iz = intersect(index2_ix_iy,index2_iz);

arma::vec I_vec = intersect(index_ix_iy_iz,index2_ix_iy_iz); 
arma::uvec I = arma::conv_to<arma::uvec>::from(I_vec);
arma::vec ind_block_vec = ix2(I) + (iy2(I)-1) * nxsup + (iz2(I)-1) * nxsup * nysup ;

/* Indices of acceptable data */
int n_block = ind_block_vec.n_elem;
arma::uvec ind_block = arma::conv_to<arma::uvec>::from(ind_block_vec);
arma::vec sr_start = nisb(ind_block-1) + 1; 
arma::vec sr_finish = nisb(ind_block ); 
arma::vec index(sr_finish.n_elem  + 1);
index(0) = 0;
index(arma::span(1,sr_finish.n_elem)) = cumsum(sr_finish - sr_start + 1); 

arma::vec I2;
I2.zeros(index(n_block));

for (int i = 0; i < n_block; ++i) { 
if (sr_start(i) <= sr_finish(i)){  
I2(arma::span(arma::span(index(i),index(i + 1) - 1))) =  Rcpp::as<arma::vec>(wrap(seq_rcpp(sr_start(i),sr_finish(i),1))); 
} 
}
/* Select the acceptable neighboring data */
arma::uvec I2_uvec = arma::conv_to<arma::uvec>::from(I2) - 1;
datacoord = datacoord.rows(I2_uvec); 
datavalue = datavalue.rows(I2_uvec);

/*  Second selection of data according to radius, angles and octant */
int n = datacoord.n_rows; 
arma::field<arma::mat> search_results(3,1); 
if (n == 0){ 
search_results(0,0) = R_NaN;
} else{
arma::colvec o1;
o1.ones(n);
arma::mat deltacoord = datacoord - o1 * coord ; 
deltacoord = trans(deltacoord * search_rotationmatrix); 
int nsector = 1;

arma::mat flag(8,n);
if ( octant==0 ){ 
flag.row(0) = Rcpp::as<arma::rowvec>(wrap(seq_rcpp(1,n,1)));  
} else{
nsector = 8; 
arma::vec ind0 = arma::conv_to<arma::vec>::from(find(deltacoord(0,arma::span(0,deltacoord.n_cols-1)) > 0));
arma::vec ind1 = arma::conv_to<arma::vec>::from(find(deltacoord(1,arma::span(0,deltacoord.n_cols-1)) > 0));
arma::vec ind2 = arma::conv_to<arma::vec>::from(find(deltacoord(2,arma::span(0,deltacoord.n_cols-1)) > 0));
arma::vec ind3 = arma::conv_to<arma::vec>::from(find(deltacoord(0,arma::span(0,deltacoord.n_cols-1)) <= 0));
arma::vec ind4 = arma::conv_to<arma::vec>::from(find(deltacoord(1,arma::span(0,deltacoord.n_cols-1)) <= 0));
arma::vec ind5 = arma::conv_to<arma::vec>::from(find(deltacoord(2,arma::span(0,deltacoord.n_cols-1)) <= 0));

arma::uvec ind0_ind1_ind2 = arma::conv_to<arma::uvec>::from(intersect(intersect(ind0,ind1),ind2));
arma::uvec ind0_ind1_ind5 = arma::conv_to<arma::uvec>::from(intersect(intersect(ind0,ind1),ind5));
arma::uvec ind0_ind4_ind2 = arma::conv_to<arma::uvec>::from(intersect(intersect(ind0,ind4),ind2));
arma::uvec ind0_ind4_ind5 = arma::conv_to<arma::uvec>::from(intersect(intersect(ind0,ind4),ind5));

arma::uvec ind3_ind1_ind2 = arma::conv_to<arma::uvec>::from(intersect(intersect(ind3,ind1),ind2));
arma::uvec ind3_ind1_ind5 = arma::conv_to<arma::uvec>::from(intersect(intersect(ind3,ind1),ind5));
arma::uvec ind3_ind4_ind2 = arma::conv_to<arma::uvec>::from(intersect(intersect(ind3,ind4),ind2));
arma::uvec ind3_ind4_ind5 = arma::conv_to<arma::uvec>::from(intersect(intersect(ind3,ind4),ind5));


arma::mat a;
arma::rowvec flag1(n, arma::fill::zeros);
arma::rowvec flag2(n, arma::fill::zeros);
arma::rowvec flag3(n, arma::fill::zeros);
arma::rowvec flag4(n, arma::fill::zeros);
arma::rowvec flag5(n, arma::fill::zeros);
arma::rowvec flag6(n, arma::fill::zeros);
arma::rowvec flag7(n, arma::fill::zeros);
arma::rowvec flag8(n, arma::fill::zeros);

flag1(ind0_ind1_ind2) = a.ones(size(ind0_ind1_ind2));
flag2(ind0_ind1_ind5) = a.ones(size(ind0_ind1_ind5));
flag3(ind0_ind4_ind2) = a.ones(size(ind0_ind4_ind2));
flag4(ind0_ind4_ind5) = a.ones(size(ind0_ind4_ind5));
flag5(ind3_ind1_ind2) = a.ones(size(ind3_ind1_ind2));
flag6(ind3_ind1_ind5) = a.ones(size(ind3_ind1_ind5));
flag7(ind3_ind4_ind2) = a.ones(size(ind3_ind4_ind2));
flag8(ind3_ind4_ind5) = a.ones(size(ind3_ind4_ind5));

flag.row(0) = flag1; 
flag.row(1) = flag2; 
flag.row(2) = flag3; 
flag.row(3) = flag4; 
flag.row(4) = flag5; 
flag.row(5) = flag6; 
flag.row(6) = flag7; 
flag.row(7) = flag8;

/* std::vector<int> flag_row = {flag1, flag2, flag3, flag4, flag5, flag6, flag7, flag8}; */
	
}	
	
/* Select the neighboring data */
arma::colvec I3(ndata * nsector,arma::fill::zeros);
int k =0;
for (int i = 0; i < nsector; ++i) {   
arma::uvec index_i = find(flag.row(i) > 0); 
if (!index_i.is_empty()){
arma::mat squareddist = pow(deltacoord.cols(index_i),2);
arma::rowvec squareddist_vec = sum(squareddist,0);
arma::rowvec sorteddist = sort(squareddist_vec);
arma::uvec J_sorted = sort_index(squareddist_vec);
index_i = index_i(J_sorted);
/* Discard the data located beyond the radius */
arma::uvec J = find(sorteddist < 1);   

if (!J.is_empty()){
arma::rowvec n_vec(2);
n_vec(0) = J.n_elem;
n_vec(1) = ndata;
int n_min = min(n_vec);
arma::uvec n_uvec = Rcpp::as<arma::uvec>(seq_rcpp(0,n_min-1,1));
index_i = index_i(n_uvec);
I3(arma::span(k,k + n_min-1)) = arma::conv_to<arma::colvec>::from(index_i); 
k = k + n_min ;	
}
}	
}
if (k > 0 ){
arma::vec I4_vec = I3(arma::span(0,k - 1));
arma::uvec I4 = arma::conv_to<arma::uvec>::from(I4_vec ); 
search_results(0,0) = I4_vec;
search_results(1,0) = datacoord.rows(I4);
search_results(2,0) = datavalue.rows(I4);
} else{
search_results(0,0) = R_NaN;
}
} 
return(search_results);
} 

//' @noRd
// [[Rcpp::export]]
arma::mat Gibbs_cosim_cpp(arma::mat datacoord, arma::mat idata, arma::mat ydata, int nfield, arma::vec flag, arma::vec nthres, arma::vec thresholds,
 arma::mat model, arma::cube sill, arma::mat b, arma::mat sillnugget, int nrealiz, int niterations,   
 arma::cube model_rotationmatrix, arma::mat search_rotationmatrix, 
 int cc_unique, int octant, int ndata,
 int nxsup, int nysup, int nzsup, int xmnsup, int ymnsup, int zmnsup, int xsizsup, int ysizsup, int zsizsup, 
 arma::rowvec ixsbtosr, arma::rowvec iysbtosr,  arma::rowvec izsbtosr, arma::colvec nisb){

int n = datacoord.n_rows; 
int nvar = 1 + nfield; 
int nsector = 1	; 
if (octant == 1){ 
nsector = 8; 
}
arma::mat index(n,n,arma::fill::ones); 
arma::mat gb_missing(nvar + n * nsector * nvar, n, arma::fill::zeros); 
arma::cube lambda(1 + n* nsector * nvar, nvar, n, arma::fill::zeros); 	

if (ndata > 0){
index.set_size(n,(1 + ndata * nsector));	
/* arma::mat index(n,(1 + ndata * nsector),arma::fill::ones); */
gb_missing.set_size(nvar + ndata * nsector * nvar, n); 
/* arma::mat gb_missing(nvar + ndata * nsector * nvar, n, arma::fill::zeros); */ 
lambda.set_size(1 + ndata * nsector * nvar, nvar, n); 
/* arma::cube lambda(1 + ndata * nsector * nvar, nvar, n, arma::fill::zeros); */ 

} 

arma::mat nb(n, 1, arma::fill::zeros); 
arma::mat nnb(n, 1, arma::fill::zeros); 
arma::cube sigma(nfield, nfield, n, arma::fill::ones); 
arma::mat idata2(n, 1); 

/* Rcpp::NumericVector I_seq = seq_rcpp(0, n-1, 1); */
/* arma::vec missing1; */
for (int i = 0; i < n; ++i) {
idata2 = idata; 
idata2(i) = R_NaN; 
arma::uvec I ; 
if (cc_unique == 0){ 
arma::field<arma::mat> search_results = cc_search_cpp(datacoord,idata2,datacoord.row(i),search_rotationmatrix,octant,ndata,nxsup,nysup,nzsup,xmnsup,ymnsup,zmnsup,xsizsup,ysizsup,zsizsup,ixsbtosr,iysbtosr,izsbtosr,nisb);
arma::vec I_search = search_results(0,0);
I = arma::conv_to<arma::uvec>::from(I_search);
nb(i) = I_search.n_elem;
} else{
I = arma::shuffle(arma::linspace<arma::uvec>(0, n - 1, n)); // Shuffle the indices
/* I = Rcpp::as<arma::uvec>(I_seq); */
nb(i) = n;
}

arma::vec missing1;
if (nb(i) > 0){
index(i,arma::span(0, (nb(i) - 1))) = arma::conv_to<arma::rowvec>::from(I);
missing1 = arma::conv_to<arma::vec>::from(arma::find_nonfinite(ydata(I)));
arma::vec missing2 = arma::conv_to<arma::vec>::from(arma::find_nonfinite(idata2(I)));
for (int j=0; j < nfield; ++j) {
missing2 +=  nb(i);
		
	
		if (!missing1.is_empty()){
		missing1 = arma::join_cols(missing1, missing2);
		}
		else{
		missing1 = missing2;
		}
}

if(missing1.is_empty()){
	nnb(i) = 1;
	  missing1 = arma::zeros(1);
	
	}
	else{
	nnb(i) = missing1.n_rows;
	}
if (nnb(i) == nvar * nb(i)){
nb(i) = 0;
} else{
gb_missing(arma::span(0, nnb(i) - 1), i) =  missing1;
arma::uvec missing1_uvec = arma::conv_to<arma::uvec>::from(missing1);
arma::field<arma::mat> gb_cokrige = cokrige_cpp(datacoord.rows(I), missing1_uvec, datacoord.row(i), model, sill, b, sillnugget, model_rotationmatrix);
lambda(arma::span(0, (nb(i) * nvar - nnb(i)) - 1),arma::span(0, nvar - 1),arma::span(i, i)) = gb_cokrige(0, 0);
sigma.slice(i) = trans(chol(gb_cokrige(1, 0)(arma::span(1, nvar-1),arma::span(1, nvar - 1))));	
}
}
if (nb(i) == 0){
arma::mat Q = sum(sill, 2);
arma::mat covariances =  sillnugget + Q;
sigma.slice(i) = trans(chol(covariances(arma::span(1, nvar - 1), arma::span(1, nvar - 1))));
}
}
/* LOOP OVER THE REALIZATIONS */
arma::mat yy_new1;
arma::mat yy_new2;
int accept = 0;
arma::cube simu(n, nfield, nrealiz);
for(int k = 0; k < nrealiz; ++k){ 
for(int ii = 0; ii < n; ++ii){ 
accept = 0; 
arma::mat y(1, nfield);
arma::mat y1(1, nfield);
while (accept < 1){ 
y.randn(1, nfield); 
y1 = y;
int gb_category = cc_truncate_cpp(y, nfield, flag, nthres, thresholds); 
if (arma::is_finite(idata(ii)) && (gb_category != idata(ii))){
accept=0;}
else{
accept = 1;
if (Rcpp::traits::is_nan<REALSXP>(idata(ii))){
y.fill(R_NaN);
}
}
}
for (int ifield = 0; ifield < nfield; ++ifield){
simu(ii,ifield,k) = y(ifield); 
}
}
/* Iterations */
for (int nit = 0; nit < niterations; ++nit){
arma::mat u;
u.randu(1, 1);
int ig = floor(n * u(0)); 
arma::mat unorm;
unorm.randn(nfield, 1);

if (arma::is_finite(idata(ig))){  
if (nb(ig) > 0){
 
arma::uvec idx1_gb = arma::conv_to<arma::uvec>::from(index(ig, arma::span(0, nb(ig) - 1)));
arma::mat QQ = simu.slice(k);
arma::mat yy = join_rows(ydata.rows(idx1_gb), QQ.rows(idx1_gb)); 
arma::mat yy_new = reshape(yy, 1, (yy.n_rows * yy.n_cols));  
arma::uvec idx2_gb = arma::conv_to<arma::uvec>::from(gb_missing(arma::span(0, nnb(ig) - 1), ig));  
yy_new.shed_cols(idx2_gb);
yy_new1 = trans(lambda.slice(ig)(arma::span(0, (nb(ig)) * (nvar) - nnb(ig) - 1),arma::span(1, nvar - 1))) * trans(yy_new) + sigma.slice(ig) * unorm;
}
else{
arma::mat yy_new1 = sigma.slice(ig) * unorm;  
}
arma::mat yy_new2 = trans(yy_new1); 
int gb_category2 = cc_truncate_cpp(yy_new2, nfield, flag, nthres, thresholds);
if (gb_category2 == idata(ig)){	
accept = 1;
}
if (accept == 1){
for (int ifield2 = 0; ifield2 < nfield; ++ifield2){
simu(ig, ifield2, k) = yy_new2(ifield2); 
}
}
}
}
}
simu.reshape(n*nfield,nrealiz,1);
return(simu);
}

//' @noRd
// [[Rcpp::export]]
arma::mat cond_mooving_neigbor(int m1,arma::mat simu, int nvar, int nrealiz, arma::mat datacoord, arma::mat cc_residuals, arma::mat coord, arma::mat search_rotationmatrix, 
int octant, int ndata,
 int nxsup, int nysup, int nzsup, int xmnsup, int ymnsup, int zmnsup, int xsizsup, int ysizsup, int zsizsup, 
 arma::rowvec ixsbtosr, arma::rowvec iysbtosr,  arma::rowvec izsbtosr, arma::colvec nisb,
arma::mat model,
arma::cube sill, arma::mat b, arma::mat sillnugget, arma::cube model_rotationmatrix ){

for (int i = 0; i < m1; ++i) {
arma::field<arma::mat> search_results = cc_search_cpp(datacoord,cc_residuals,coord.row(i),search_rotationmatrix,octant,ndata,nxsup,nysup,nzsup,xmnsup,ymnsup,zmnsup,xsizsup,ysizsup,zsizsup,ixsbtosr,iysbtosr,izsbtosr,nisb);
/* arma::mat check_empty = search_results(0,0); */
if (arma::is_finite(search_results(0, 0)(0))) {
/* if(arma::is_finite(check_empty(0))){ */
arma::mat datacoord_i = search_results(1,0);
arma::mat cc_residuals_i = search_results(2,0);
int n_i = datacoord_i.n_rows;

arma::uvec index_missing = arma::find_nonfinite(cc_residuals_i.cols(arma::span(0, nvar - 1)));
arma::field<arma::mat> gb_cokrige = cokrige_cpp(datacoord_i, index_missing, coord.row(i), model, sill, b, sillnugget+1e-7,model_rotationmatrix);
arma::mat cc_weights = gb_cokrige(0,0); 
cc_residuals_i = reshape(cc_residuals_i, n_i * nvar, nrealiz / nvar);
cc_residuals_i.shed_rows(index_missing);

Rcpp::NumericVector ii = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(i));
arma::uvec ii_uvec = Rcpp::as<arma::uvec>(ii);
for (int j = 0 ; j < nvar; j++){
Rcpp::NumericVector I_seq = seq_rcpp(j, nrealiz-1, nvar);
arma::uvec I = Rcpp::as<arma::uvec>(I_seq);
simu.submat(ii_uvec,I) = simu.submat(ii_uvec,I)+ trans(cc_weights.col(j)) * cc_residuals_i;
}
}
}
return(simu);
}

//' @noRd
// [[Rcpp::export]]
arma::mat Gibbs_plurisim_cpp(arma::mat datacoord, arma::mat idata, 
int nfield, arma::vec flag, arma::vec nthres, arma::vec thresholds, arma::mat model, arma::cube sill, arma::mat b, arma::mat sillnugget, 
int nrealiz, int niterations, arma::cube model_rotationmatrix, arma::mat search_rotationmatrix,
int cc_unique, int octant, int ndata,
int nxsup, int nysup, int nzsup, int xmnsup, int ymnsup, int zmnsup, int xsizsup, int ysizsup, int zsizsup, 
arma::rowvec ixsbtosr, arma::rowvec iysbtosr,  arma::rowvec izsbtosr, arma::colvec nisb){

int n = datacoord.n_rows; 
int nvar = nfield; 
int nsector = 1	;
if (octant == 1){
nsector = 8;
}
arma::mat index(n,n,arma::fill::ones);
arma::mat gb_missing(nvar + n * nsector * nvar, n, arma::fill::zeros);
arma::cube lambda(1 + n* nsector * nvar, nvar, n, arma::fill::zeros); 

if (ndata > 0){
arma::mat index(n,(1 + ndata * nsector),arma::fill::ones);
arma::mat gb_missing(nvar + ndata * nsector * nvar, n, arma::fill::zeros); 
arma::cube lambda(1 + ndata * nsector * nvar, nvar, n, arma::fill::zeros); 
}
arma::mat nb(n, 1, arma::fill::zeros);
arma::mat nnb(n, 1, arma::fill::zeros);
arma::cube sigma(nfield, nfield, n, arma::fill::ones);
arma::mat idata2(n, 1);
Rcpp::NumericVector I_seq = seq_rcpp(0, n-1, 1);

for (int i = 0; i < n; ++i) {
idata2 = idata;
idata2(i) = R_NaN;

arma::uvec I ;
if (cc_unique == 0){
arma::field<arma::mat> search_results = cc_search_cpp(datacoord,idata2,datacoord.row(i),search_rotationmatrix,octant,ndata,nxsup,nysup,nzsup,xmnsup,ymnsup,zmnsup,xsizsup,ysizsup,zsizsup,ixsbtosr,iysbtosr,izsbtosr,nisb);
arma::vec I_search = search_results(0,0);
I = arma::conv_to<arma::uvec>::from(I_search);
nb(i) = I_search.n_elem;
} else{
I = Rcpp::as<arma::uvec>(I_seq);
nb(i) = n;
}

if (nb(i) > 0){
	index(i,arma::span(0, (nb(i) - 1))) = arma::conv_to<arma::rowvec>::from(I);
arma::vec missing1 ;
nnb(i) = 0;
if (nnb(i) == nvar * nb(i)){
nb(i) = 0;
} else{

arma::uvec missing1_uvec = arma::conv_to<arma::uvec>::from(missing1);
arma::field<arma::mat> gb_cokrige = cokrige_cpp(datacoord.rows(I), missing1_uvec, datacoord.row(i), model, sill, b, sillnugget, model_rotationmatrix);
lambda(arma::span(0, (nb(i) * nvar - nnb(i)) - 1),arma::span(0, nvar - 1),arma::span(i, i)) = gb_cokrige(0, 0);
}
}

if (nb(i) == 0){
arma::mat Q = sum(sill, 2);
arma::mat covariances =  sillnugget + Q;
sigma.slice(i) = trans(chol(covariances(arma::span(1, nvar - 1), arma::span(1, nvar - 1))));
}
}
/* LOOP OVER THE REALIZATIONS */
arma::cube simu(n, nfield, nrealiz);
for(int k = 0; k < nrealiz; ++k){
for(int ii = 0; ii < n; ++ii){
int accept = 0;
arma::mat y(1, nfield);
arma::mat y1(1, nfield);
while (accept < 1){
y.randn(1, nfield);
y1 = y;
int gb_category = cc_truncate_cpp(y, nfield, flag, nthres, thresholds);
if (arma::is_finite(idata(ii)) && (gb_category != idata(ii))){
accept=0;}
else{
accept = 1;
if (Rcpp::traits::is_nan<REALSXP>(idata(ii))){
y.fill(R_NaN);
}
}
}
for (int ifield = 0; ifield < nfield; ++ifield){
simu(ii,ifield,k) = y(ifield);
}
}
/* Iterations */
for (int nit = 0; nit < niterations; ++nit){
arma::mat u;
u.randu(1, 1);
int ig = floor(n * u(0));
arma::mat unorm;
unorm.randn(nfield, 1);
arma::mat yy_new1;
if (arma::is_finite(idata(ig))){
if (nb(ig) > 0){
arma::uvec idx1_gb = arma::conv_to<arma::uvec>::from(index(ig, arma::span(0, nb(ig) - 1)));
arma::mat QQ = simu.slice(k);
arma::mat yy =  QQ.rows(idx1_gb);
arma::mat yy_new = reshape(yy, 1, (yy.n_rows * yy.n_cols));
yy_new1 = trans(lambda.slice(ig)(arma::span(0, (nb(ig)) * (nvar) - nnb(ig) - 1),arma::span(0, nvar - 1))) * trans(yy_new) + sigma.slice(ig) * unorm;
}
else{
arma::mat yy_new1 = sigma.slice(ig) * unorm;
}
arma::mat yy_new2 = trans(yy_new1);
int gb_category2 = cc_truncate_cpp(yy_new2, nfield, flag, nthres, thresholds);
if (gb_category2 == idata(ig)){
for (int ifield2 = 0; ifield2 < nfield; ++ifield2){
simu(ig, ifield2, k) = yy_new2(ifield2);
}
}
}
}
}
simu.reshape(n*nfield,nrealiz,1);
return(simu);
}

//' @noRd
// [[Rcpp::export]]
void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}

//' @noRd
// [[Rcpp::export]]
arma::mat tbmain_cosimu_cpp(arma::mat coord, arma::mat model, arma::vec cc_sigma, arma::cube A1, int nvar, arma::vec nlines, int nrealiz, arma::mat seed,arma::cube all_lines, arma::mat all_offset, arma::cube all_r, arma::cube all_phi, arma::cube all_theta, arma::mat valid_lines){
/* Non conditional simulation by the turning bands method (simulation along the lines) */
int nst = model.n_rows;
int m = coord.n_rows;
arma::mat simu(m,nrealiz,arma::fill::zeros);

/* Loop over each nested structure */
for (int i = 0; i < nst; ++i) {
arma::mat sim(m,nrealiz,arma::fill::zeros);
arma::mat cc_lines;
arma::uvec index;
if(model(i,0) < 4.5 || model(i,0) == 9){
/* Loop over the realizations */
arma::mat o1(m,1,arma::fill::ones);
for (int k = 0; k < nrealiz; ++k) {
Rcpp::NumericVector index_vec = seq_rcpp(((k) * (nlines(i) - 1)),((k + 1) * (nlines(i) - 1)),1);
index = Rcpp::as<arma::uvec>(index_vec);
cc_lines = all_lines.slice(i).rows(index);
arma::mat valid_lines_col = valid_lines.col(i);
arma::uvec valid = find(valid_lines_col.rows(index) > 0);
int nbvalid = valid.n_elem;
int nbinvalid = nlines(i) - nbvalid;
/* Project the points to simulate over the lines of the i-th nested structure */
arma::mat x = coord * trans(cc_lines);
arma::mat all_offset_col = all_offset.col(i);
arma::mat cc_offset = o1 * trans(all_offset_col.rows(index));
arma::mat interval_number = floor(x - cc_offset) + 1;
arma::mat position = x - cc_offset - interval_number + 0.5;
arma::mat state = seed((k * nst + i),arma::span(0,nlines(i) - 1));
/* Loop over the lines */
for (int j = 0; j < nbvalid; ++j) {
/* Simulate the values by a partition method */
set_seed(state(valid(j)));
int maxnum = max(interval_number(arma::span::all,valid(j)));
arma::mat norm_rand = arma::randn(1,maxnum);
norm_rand.elem(find(norm_rand > 0)).ones();
norm_rand.elem(find(norm_rand != 1)).zeros();
arma::mat slope = cc_sigma(i) * (2*norm_rand-1);
arma::mat A(m,1);
for (int jj=0; jj< m; ++jj){
A.row(jj) = slope.col(interval_number(jj,valid(j))-1);
}
sim.col(k) = sim.col(k) + A % position.col(valid(j));
}
sim.col(k) = sim.col(k) + sqrt(nbinvalid/nlines(i)) * arma::randn(m,1);
}
} else if(model(i,0) == 5){
/* Loop over the realizations */
arma::mat o1(m,1,arma::fill::ones);
for (int k = 0; k < nrealiz; ++k) {
Rcpp::NumericVector index_vec = seq_rcpp(((k)*(nlines(i)-1)),((k+1)*(nlines(i)-1)),1);
index = Rcpp::as<arma::uvec>(index_vec);
cc_lines = all_lines.slice(i).rows(index);
arma::mat valid_lines_col = valid_lines.col(i);
arma::uvec valid = find(valid_lines_col.rows(index) > 0);
int nbvalid = valid.n_elem;
int nbinvalid = nlines(i) - nbvalid;
/* Project the points to simulate over the lines of the i-th nested structure */
arma::mat x = coord * trans(cc_lines);
arma::mat all_offset_col = all_offset.col(i);
arma::mat cc_offset = o1*trans(all_offset_col.rows(index));
arma::mat interval_number = floor(x-cc_offset) + 1;
arma::mat position = x - cc_offset - interval_number + 0.5;
arma::mat state = seed((k * nst + i),arma::span(0,nlines(i)-1));
/* Loop over the lines */
for (int j = 0; j < nbvalid; ++j) {
/* Simulate the values by a partition method */
set_seed(state(valid(j)));
int maxnum = max(interval_number(arma::span::all,valid(j)));
arma::mat norm_rand = arma::randn(1,maxnum);
norm_rand.elem(find(norm_rand > 0)).ones();
norm_rand.elem(find(norm_rand != 1)).zeros();
arma::mat slope = cc_sigma(i) * (2 * norm_rand-1);
arma::mat A(m,1);
for (int jj = 0; jj< m; ++jj){
A.row(jj) = slope.col(interval_number(jj,valid(j))-1);
}
sim.col(k) = sim.col(k) + A % (0.25 * position.col(valid(j)) - pow(position.col(valid(j)),3));
}
sim.col(k) = sim.col(k) + sqrt(nbinvalid/nlines(i)) * arma::randn(m,1);
}
} else{
/* Loop over the realizations */
arma::mat o1(1,m,arma::fill::ones);
for (int k = 0; k < nrealiz; ++k) {
 /* Project the points to simulate over the lines of the i-th nested structure */
Rcpp::NumericVector index_vec = seq_rcpp(((k)*(nlines(i)-1)),((k+1)*(nlines(i)-1)),1);
index = Rcpp::as<arma::uvec>(index_vec);
cc_lines = all_lines.slice(i).rows(index);
arma::mat x = coord * trans(cc_lines);
 /* Simulate the values by a continuous spectral method */
arma::mat A_r = all_r(arma::span::all,arma::span(k),arma::span(i));
arma::mat r = A_r * o1 ;
arma::mat A_phi = all_phi(arma::span::all,arma::span(k),arma::span(i));
arma::mat phi =  A_phi  * o1;
arma::mat A_theta = all_theta(arma::span::all,arma::span(k),arma::span(i));
arma::mat theta = A_theta * o1;
arma::mat A = theta % cos((r % trans(x)) + phi);
sim.col(k) = sim.col(k) + cc_sigma(i) * trans(sum(A,0));
}
}
sim = reshape(trans(sim),nvar,(nrealiz / nvar * m));
sim = reshape(trans(A1.slice(i)) * sim,nrealiz,m);
simu = simu + trans(sim);
}
return(simu);
}

//' @noRd
// [[Rcpp::export]]
arma::mat tbmain_simu_cpp(arma::mat coord, arma::mat model, arma::vec cc_sigma, arma::cube A1, int nvar, arma::vec nlines, int nrealiz, arma::mat seed,arma::cube all_lines, arma::mat all_offset, arma::cube all_r, arma::cube all_phi, arma::cube all_theta, arma::mat valid_lines){
/* Non conditional simulation by the turning bands method (simulation along the lines) */
int nst = model.n_rows;
int m = coord.n_rows;
arma::mat simu(m,nrealiz,arma::fill::zeros);
/* Loop over each nested structure */
for (int i = 0; i < nst; ++i) {
arma::mat sim(m,nrealiz,arma::fill::zeros);
arma::mat cc_lines;
arma::uvec index;
if(model(i,0) < 4.5 || model(i,0) == 9){
/* Loop over the realizations */
arma::mat o1(m,1,arma::fill::ones);
for (int k = 0; k < nrealiz; ++k) {
Rcpp::NumericVector index_vec = seq_rcpp(((k) * (nlines(i) - 1)),((k + 1) * (nlines(i) - 1)),1);
index = Rcpp::as<arma::uvec>(index_vec);
cc_lines = all_lines.slice(i).rows(index);
arma::mat valid_lines_col = valid_lines.col(i);
arma::uvec valid = find(valid_lines_col.rows(index) > 0);
int nbvalid = valid.n_elem;
int nbinvalid = nlines(i) - nbvalid;
/* Project the points to simulate over the lines of the i-th nested structure */
arma::mat x = coord * trans(cc_lines);
arma::mat all_offset_col = all_offset.col(i);
arma::mat cc_offset = o1 * trans(all_offset_col.rows(index));
arma::mat interval_number = floor(x-cc_offset) + 1;
arma::mat position = x - cc_offset - interval_number + 0.5;
arma::mat state = seed((k * nst + i),arma::span(0,nlines(i) - 1));
/* Loop over the lines */
for (int j = 0; j < nbvalid; ++j) {
/* Simulate the values by a partition method */
set_seed(state(valid(j)));
int maxnum = max(interval_number(arma::span::all,valid(j)));
arma::mat norm_rand = arma::randn(1,maxnum);
norm_rand.elem(find(norm_rand > 0)).ones();
norm_rand.elem(find(norm_rand != 1)).zeros();
arma::mat slope = cc_sigma(i) * (2 * norm_rand - 1);
arma::mat A(m,1);
for (int jj = 0; jj< m; ++jj){
A.row(jj) = slope.col(interval_number(jj,valid(j)) - 1);
}
sim.col(k) = sim.col(k) + A % position.col(valid(j));
}
sim.col(k) = sim.col(k) + sqrt(nbinvalid/nlines(i)) * arma::randn(m,1);
}
} else if(model(i,0) == 5){
/* Loop over the realizations */
arma::mat o1(m,1,arma::fill::ones);
for (int k = 0; k < nrealiz; ++k) {
Rcpp::NumericVector index_vec = seq_rcpp(((k)*(nlines(i) - 1)),((k + 1) * (nlines(i) - 1)),1);
index = Rcpp::as<arma::uvec>(index_vec);
cc_lines = all_lines.slice(i).rows(index);
arma::mat valid_lines_col = valid_lines.col(i);
arma::uvec valid = find(valid_lines_col.rows(index) > 0);
int nbvalid = valid.n_elem;
int nbinvalid = nlines(i) - nbvalid;
/* Project the points to simulate over the lines of the i-th nested structure */
arma::mat x = coord * trans(cc_lines);
arma::mat all_offset_col = all_offset.col(i);
arma::mat cc_offset = o1*trans(all_offset_col.rows(index));
arma::mat interval_number = floor(x-cc_offset) + 1;
arma::mat position = x - cc_offset - interval_number + 0.5;
arma::mat state = seed((k * nst + i),arma::span(0,nlines(i) - 1));
// # % Loop over the lines
for (int j = 0; j < nbvalid; ++j) {
/* Simulate the values by a partition method */
set_seed(state(valid(j)));
int maxnum = max(interval_number(arma::span::all,valid(j)));
arma::mat norm_rand = arma::randn(1,maxnum);
norm_rand.elem(find(norm_rand > 0)).ones();
norm_rand.elem(find(norm_rand != 1)).zeros();
arma::mat slope = cc_sigma(i) * (2 * norm_rand - 1);
arma::mat A(m,1);
for (int jj = 0; jj< m; ++jj){
A.row(jj) = slope.col(interval_number(jj,valid(j)) - 1);
}
sim.col(k) = sim.col(k) + A % (0.25 * position.col(valid(j)) - pow(position.col(valid(j)),3));
}
sim.col(k) = sim.col(k) + sqrt(nbinvalid/nlines(i)) * arma::randn(m,1);
}
} else{
/* Loop over the realizations */
arma::mat o1(1,m,arma::fill::ones);
for (int k = 0; k < nrealiz; ++k) {
/* Project the points to simulate over the lines of the i-th nested structure */
Rcpp::NumericVector index_vec = seq_rcpp(((k) * (nlines(i) - 1)),((k + 1) * (nlines(i) - 1)),1);
index = Rcpp::as<arma::uvec>(index_vec);
cc_lines = all_lines.slice(i).rows(index);
arma::mat x = coord * trans(cc_lines);
/* Simulate the values by a continuous spectral method */
arma::mat A_r = all_r(arma::span::all,arma::span(k),arma::span(i));
arma::mat r = A_r * o1 ;
arma::mat A_phi = all_phi(arma::span::all,arma::span(k),arma::span(i));
arma::mat phi =  A_phi  * o1;
arma::mat A_theta = all_theta(arma::span::all,arma::span(k),arma::span(i));
arma::mat theta = A_theta * o1;
arma::mat A = theta % cos((r % trans(x)) + phi);//
sim.col(k) = sim.col(k) + cc_sigma(i) * trans(sum(A,0));
}
}
sim = reshape(trans(sim),nvar,(nrealiz / nvar * m));
sim = reshape(trans(A1.slice(i)) * sim,nrealiz,m);
simu = simu + trans(sim);
}
 return(simu);
}