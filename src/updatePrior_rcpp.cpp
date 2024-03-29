#include <RcppArmadillo.h>
#include <vector>
#include <cmath> // std::exp()
#ifdef _OPENMP
#include <omp.h> // openMP
#endif
//[[Rcpp::plugins(openmp)]]
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//' Computation of the sum sum_c(w_c*psi_c) parallelized version
//'@param  Hconfig list of  vector of 0 and 1, corresponding to the configurations
//'@param  NewPrior a double vector containing the prior w_c
//'@param  Logf0Mat a double matrix containing the log(f0(xi_q)) 
//'@param  Logf1Mat a double matrix containing the log(f1(xi_q)) 
//'@param  threads_nb an int the number of threads
//'@return a double vector containing sum_c(w_c*psi_c)
//'
// [[Rcpp::export]]
arma::vec fHconfig_sum_update_ptr_parallel(const List & Hconfig,
                                   const arma::vec & NewPrior,
                                   const arma::mat & Logf0Mat,
                                   const arma::mat & Logf1Mat,
                                   int threads_nb = 0 ) {

   if(threads_nb < 0) threads_nb = 1;
   // The number of thread is set to the number of core available by default.
   // and limited to the number of core available if a higher number is asked.
   #ifdef _OPENMP
   int max_threads_nb = omp_get_num_procs();
   #else
   int max_threads_nb = 1;
   #endif
   if(threads_nb > max_threads_nb || threads_nb == 0) threads_nb = max_threads_nb;

   // Get the dimensions
   arma::uword row_nb = Logf0Mat.n_rows;
   arma::uword col_nb = Logf0Mat.n_cols;
   arma::uword nb_config = Hconfig.size();
   // Check the dimensions
   if(NewPrior.size() != nb_config)
     stop("OldPrior and Hconfig lengths are different.");
   if(Logf1Mat.n_rows != row_nb || Logf1Mat.n_cols != col_nb)
     stop("Logf0Mat and Logf1Mat Dimensions are different.");

   std::vector<const int*> HconfigPtr(nb_config);
   // int** HconfigPtr = new int*[nb_config];
   for(size_t i=0; i<nb_config; ++i){
     const IntegerVector v = Hconfig[i];
     HconfigPtr[i] = &(v[0]);
   }

   // final fHconfig_sum, shared among all threads
   // => we cannot use this object in the parallel for loop,
   // See fHconfig_sum_local declaration below.
   arma::vec fHconfig_sum(row_nb, arma::fill::zeros);

   // beginning of the multi-thread section. Each thread executes this section.
   #pragma omp parallel num_threads(threads_nb)
   {
     // fHconfig_sum_local is a thread-local object, will accumulate the results for each thread,
     // each of these thread-local objects value will later be accumulated in the final fHconfig_sum.
     arma::vec fHconfig_sum_local(row_nb, arma::fill::zeros);

     // f_c is a thread-local object (each thread has its own f_c)
     arma::vec f_c(row_nb);

     //Loop on configurations: updating prior
     #pragma omp for
     for(arma::uword c=0; c<nb_config; ++c){

       //computation of the density f_c
       f_c.fill(0.0);

       const int* h = HconfigPtr[c];

       for(arma::uword i=0; i<row_nb; ++i) {
         double sum = 0;
         for (arma::uword j = 0; j < col_nb; ++j) {
           if ( h[(size_t)j] == 0) sum += Logf0Mat.at(i,j);
           else sum += Logf1Mat.at(i, j);
         }
         f_c.at(i) = sum;
       }

       // f_c = exp(f_c);
       for(arma::uword i=0; i<row_nb; ++i) { f_c.at(i) = exp(f_c.at(i)); }

       //update fHconfig_sum_local
       for(arma::uword i=0; i<row_nb; ++i){
         fHconfig_sum_local.at(i) += NewPrior.at(c) * f_c.at(i);
       }

     }
     // Sum the thread-local results in the final result within the critical section (one thread at a time).
     #pragma omp critical
     fHconfig_sum += fHconfig_sum_local;

   } // end of multi-threaded section

   return fHconfig_sum;
 }


 
//' Update of the prior estimate in EM algo parallelized version
//'@param  Hconfig list of  vector of 0 and 1, corresponding to the configurations
//'@param  fHconfig_sum a double vector containing sum_c(w_c*psi_c), obtained by fHconfig_sum_update_ptr_parallel()
//'@param  OldPrior a double vector containing the prior w_c 
//'@param  Logf0Mat a double matrix containing the log(f0(xi_q)) 
//'@param  Logf1Mat a double matrix containing the log(f1(xi_q)) 
//'@param  threads_nb an int the number of threads
//'@return a double vector containing the new estimate of prior w_c
//'
// [[Rcpp::export]]
arma::vec prior_update_arma_ptr_parallel(const List & Hconfig,
                                     const arma::vec & fHconfig_sum,
                                     const arma::vec & OldPrior,
                                     const arma::mat & Logf0Mat,
                                     const arma::mat & Logf1Mat,
                                     int threads_nb = 0){


  if(threads_nb < 0) threads_nb = 1;
  // The number of thread is set to the number of core available by default.
  // and limited to the number of core available if a higher number is asked.
  #ifdef _OPENMP
  int max_threads_nb = omp_get_num_procs();
  #else
  int max_threads_nb = 1;
  #endif
  if(threads_nb > max_threads_nb || threads_nb == 0) threads_nb = max_threads_nb;

  // Get the dimensions
  arma::uword row_nb = Logf0Mat.n_rows;
  arma::uword col_nb = Logf0Mat.n_cols;
  arma::uword nb_config = Hconfig.size();
  // Check the dimensions
  if(OldPrior.size() != nb_config)
    stop("OldPrior and Hconfig lengths are different.");
  if(fHconfig_sum.size() != row_nb)
    stop("fHconfig_sum length and Logf0Mat number of rows are different.");
  if(Logf1Mat.n_rows != row_nb || Logf1Mat.n_cols != col_nb)
    stop("Logf0Mat and Logf1Mat Dimensions are different.");

  // As it is not safe to involved Rcpp objects (here Rcpp::List and Rcpp::IntegerVector) in multi-threaded sections,
  // we first used arma::field<arma::uvec> for 1st parameter Hconfig, but this leads to deep copy of data when
  // the R list Hconfig is converted to arma::field.
  //
  // Therefore an alternative method is to define a standard C++ vector (HconfigPtr)
  // storing the memory addresses of the integer vectors included in Hconfig.
  std::vector<const int*> HconfigPtr(nb_config);
  for(size_t i=0; i<nb_config; ++i){
    const IntegerVector v = Hconfig[i];
    HconfigPtr[i] = &(v[0]);
  }
  // This allows to access the raw data located at these addresses in the multi-threaded section without involving
  // the Rcpp objects (Rcpp::IntegerVector) owning these data. See the schema below.

  //                      |   R list Hconfig            |
  //                      |                             |
  //                      |   [[1]] R integer vector    |
  //                      |   |   raw data           |  |
  //  HconfigPtr[0]-------|---|-->[0|0|0|0|0|0|0|0]  |  |
  //                      |   |                      |  |
  //                      |                             |
  //                      |   [[2]] R integer vector    |
  //                      |   |   raw data           |  |
  //  HconfigPtr[1]-------|---|-->[1|0|0|0|0|0|0|0]  |  |
  //                      |   |                      |  |
  //                      |                             |
  //

  // This should be safe as Hconfig and HconfigPtr are never modified in this function.
  // IMPORTANT: R side Integer vector uses a 32 bits (4 Bytes) integer type, even on 64 bits system.
  //            => We have to use a pointer type compatible with 32 bits integer !
  //              'long*' would be wrong as sizeof(long) = 8 on 64 bits system.
  //              'int*' is ok as sizeof(int) is always equal to 4 even on 64 bits system.
  //

  arma::vec NewPrior(nb_config);

  // beginning of the multi-thread section.
  // Every variable/object declared above this section is shared by all the threads.
  // Every variable/object declared inside this section is private (local) to each thread.
  #pragma omp parallel num_threads(threads_nb)
  {
    // f_c is a thread-local object (each thread has its own f_c)
    arma::vec f_c(row_nb);

    // Loop on configurations: updating prior
    // Iterations are equaly distributed among threads.
    // Example with threads_nb = 4 and nb_config = 256 (Q=8):
    // 1st thread executes iterations 0 to 63,
    // 2nd thread 64 to 127
    // 3rd thread 128 to 191
    // 4th thread 192 to 255
    #pragma omp for
    for(arma::uword c=0; c<nb_config; ++c){

      //computation of the density f_c
      f_c.fill(0.0);

      // arma::ivec advanced constructor allows to use the memory address of the vector data stored externally without copying it.
      // const arma::ivec h(HconfigPtr[c], col_nb, false, true);
      // However, we don't need to embed this memory block in a arma::ivec object since it is not implied in vector/matrix operations.
      //     => We can use it as is.
      // h can be used as a pure C array.
      const int* h = HconfigPtr[c];

      // New version:
      //   Swapping inner loop and outer loop allows to collect the sum in a temp variable (sum),
      //   and assign its contents to f_c.at(i) only after the inner loop is finished.
      //   => This divides the number of call to f_c.at(i) by col_nb, and also avoids repeated index lookups in multi-threaded context.
      for(arma::uword i=0; i<row_nb; ++i) {
        double sum = 0;
        for (arma::uword j = 0; j < col_nb; ++j) {
          if ( h[(size_t)j] == 0 ) sum += Logf0Mat.at(i, j);
          else sum += Logf1Mat.at(i, j);
        }
        f_c.at(i) = sum;
      }

      for(arma::uword i=0; i<row_nb; ++i) { f_c.at(i) = exp(f_c.at(i)); }

      //update the prior
      // NewPrior.at(c) = mean(OldPrior.at(c)*f_c/fHconfig_sum);
      double sum = 0;
      double oldPrior_c = OldPrior.at(c);
      for(arma::uword k = 0; k<row_nb;++k){
        sum += oldPrior_c*f_c.at(k)/fHconfig_sum.at(k);
      }
      NewPrior.at(c) = sum/row_nb;
    }
  } // end of multi-threaded section

  return(NewPrior);
}
 
 
//' Computation of the sum sum_c(w_c*psi_c) using Gaussian copula parallelized version
//'@param  Hconfig list of  vector of 0 and 1, corresponding to the configurations
//'@param  NewPrior a double vector containing the prior w_c
//'@param  Logf0Mat a double matrix containing the log(f0(xi_q)) 
//'@param  Logf1Mat a double matrix containing the log(f1(xi_q)) 
//'@param  zeta0 a double matrix containing the qnorm(F0(x_iq))
//'@param  zeta1 a double matrix containing the qnorm(F1(x_iq))
//'@param  R a double matrix corresponding to the copula parameter
//'@param  Rinv a double matrix corresponding to the inverse copula parameter
//'@param  threads_nb an int the number of threads
//'@return a double vector containing sum_c(w_c*psi_c)
//'
// [[Rcpp::export]]
arma::vec fHconfig_sum_update_gaussian_copula_ptr_parallel(const List & Hconfig,
                                                            const arma::vec & NewPrior,
                                                            const arma::mat & Logf0Mat,
                                                            const arma::mat & Logf1Mat,
                                                            const arma::mat & zeta0,
                                                            const arma::mat & zeta1,
                                                            const arma::mat & R,
                                                            const arma::mat & Rinv,
                                                            int threads_nb = 0 ) {
   
   if(threads_nb < 0) threads_nb = 1;
   // The number of thread is set to the number of core available by default.
   // and limited to the number of core available if a higher number is asked.
   #ifdef _OPENMP
   int max_threads_nb = omp_get_num_procs();
   #else
   int max_threads_nb = 1;
   #endif
   if(threads_nb > max_threads_nb || threads_nb == 0) threads_nb = max_threads_nb;
   
   // Get the dimensions
   arma::uword row_nb = Logf0Mat.n_rows;
   arma::uword col_nb = Logf0Mat.n_cols;
   arma::uword nb_config = Hconfig.size();
   // Check the dimensions
   if(NewPrior.size() != nb_config)
     stop("OldPrior and Hconfig lengths are different.");
   if(Logf1Mat.n_rows != row_nb || Logf1Mat.n_cols != col_nb)
     stop("Logf0Mat and Logf1Mat Dimensions are different.");
   
   std::vector<const int*> HconfigPtr(nb_config);
   
   for(size_t i=0; i<nb_config; ++i){
     const IntegerVector v = Hconfig[i];
     HconfigPtr[i] = &(v[0]);
   }
   
   
   // Initialization
   // final fHconfig_sum, shared among all threads
   // => we cannot use this object in the parallel for loop,
   // See fHconfig_sum_local declaration below.
   arma::vec fHconfig_sum(row_nb, arma::fill::zeros);
   arma::mat Rinv_I = Rinv - arma::eye(col_nb, col_nb);
   double sqrt_detR = sqrt(arma::det(R));
   
   // beginning of the multi-thread section. Each thread executes this section.
#pragma omp parallel num_threads(threads_nb)
{
  // fHconfig_sum_local is a thread-local object, will accumulate the results for each thread,
  // each of these thread-local objects value will later be accumulated in the final fHconfig_sum. 
  arma::vec fHconfig_sum_local(row_nb, arma::fill::zeros);
  
  // xxx_c is a thread-local object (each thread has its own xxx_c)
  arma::vec f_c(row_nb);
  arma::vec copula_c(row_nb);
  arma::mat zeta_c(row_nb,col_nb);
  
  //Loop on configurations: updating prior
#pragma omp for
  for(arma::uword c=0; c<nb_config; ++c){
    
    //computation of the density f_c and zeta_c
    f_c.fill(0.0);
    copula_c.fill(0.0);
    zeta_c.fill(0.0);
    
    
    const int* h = HconfigPtr[c];
    
    for(arma::uword i=0; i<row_nb; ++i) {
      double sum = 0;
      for (arma::uword j = 0; j < col_nb; ++j) {
        if ( h[(size_t)j] == 0) {sum += Logf0Mat.at(i,j);
          zeta_c.at(i,j) = zeta0.at(i,j);}
        else { // h[j] == 1
          sum += Logf1Mat.at(i, j);
          zeta_c.at(i,j) = zeta1.at(i,j);
        }
      }
      f_c.at(i) = sum;
    }
    
    
    // f_c = exp(f_c);
    for(arma::uword i=0; i<row_nb; ++i) { f_c.at(i) = exp(f_c.at(i)); }
    
    //computation of the gaussian copula density  
    for(arma::uword i=0; i<row_nb; ++i) {
      copula_c.at(i) = (1/sqrt_detR) * exp(-0.5*arma::dot(zeta_c.row(i), Rinv_I*zeta_c.row(i).t()));
    }
    
    //update fHconfig_sum_local
    for(arma::uword i=0; i<row_nb; ++i){
      fHconfig_sum_local.at(i) += NewPrior.at(c) * f_c.at(i) * copula_c.at(i);
    }
    
  }
  // Sum the thread-local results in the final result within the critical section (one thread at a time).
#pragma omp critical
  fHconfig_sum += fHconfig_sum_local;

} // end of multi-threaded section

return fHconfig_sum;
 }
 
 
 
//' Update of the prior estimate in EM algo using Gaussian copula, parallelized version
//'@param  Hconfig list of  vector of 0 and 1, corresponding to the configurations
//'@param  fHconfig_sum a double vector containing sum_c(w_c*psi_c), obtained by fHconfig_sum_update_ptr_parallel()
//'@param  OldPrior a double vector containing the prior w_c 
//'@param  Logf0Mat a double matrix containing the log(f0(xi_q)) 
//'@param  Logf1Mat a double matrix containing the log(f1(xi_q))
//'@param  zeta0 a double matrix containing the qnorm(F0(x_iq))
//'@param  zeta1 a double matrix containing the qnorm(F1(x_iq))
//'@param  R a double matrix corresponding to the copula parameter
//'@param  Rinv a double matrix corresponding to the inverse copula parameter 
//'@param  threads_nb an int the number of threads
//'@return a double vector containing the new estimate of prior w_c
//'
// [[Rcpp::export]]
arma::vec prior_update_gaussian_copula_ptr_parallel(const List & Hconfig,
                                                     const arma::vec & fHconfig_sum,
                                                     const arma::vec & OldPrior,
                                                     const arma::mat & Logf0Mat,
                                                     const arma::mat & Logf1Mat,
                                                     const arma::mat & zeta0,
                                                     const arma::mat & zeta1,
                                                     const arma::mat & R,
                                                     const arma::mat & Rinv,
                                                     int threads_nb = 0) {
   
   
   if(threads_nb < 0) threads_nb = 1;
   // The number of thread is set to the number of core available by default.
   // and limited to the number of core available if a higher number is asked.
   #ifdef _OPENMP
   int max_threads_nb = omp_get_num_procs();
   #else
   int max_threads_nb = 1;
   #endif
   if(threads_nb > max_threads_nb || threads_nb == 0) threads_nb = max_threads_nb;
   
   // Get the dimensions
   arma::uword row_nb = Logf0Mat.n_rows;
   arma::uword col_nb = Logf0Mat.n_cols;
   arma::uword nb_config = Hconfig.size();
   // Check the dimensions
   if(OldPrior.size() != nb_config)
     stop("OldPrior and Hconfig lengths are different.");
   if(fHconfig_sum.size() != row_nb)
     stop("fHconfig_sum length and Logf0Mat number of rows are different.");
   if(Logf1Mat.n_rows != row_nb || Logf1Mat.n_cols != col_nb)
     stop("Logf0Mat and Logf1Mat Dimensions are different.");
   
   std::vector<const int*> HconfigPtr(nb_config);
   for(size_t i=0; i<nb_config; ++i){
     const IntegerVector v = Hconfig[i];
     HconfigPtr[i] = &(v[0]);
   }
   
   // Initialization
   arma::vec NewPrior(nb_config);
   arma::mat Rinv_I = Rinv- arma::eye(col_nb, col_nb);
   double sqrt_detR = sqrt(arma::det(R));
   
   // beginning of the multi-thread section.
   // Every variable/object declared above this section is shared by all the threads.
   // Every variable/object declared inside this section is private (local) to each thread.
#pragma omp parallel num_threads(threads_nb)
{
  // xxx_c is a thread-local object (each thread has its own xxx_c)
  arma::vec f_c(row_nb);
  arma::vec copula_c(row_nb);
  arma::mat zeta_c(row_nb,col_nb);
  
  // Loop on configurations: updating prior
  // Iterations are equaly distributed among threads.
  // Example with threads_nb = 4 and nb_config = 256 (Q=8):
  // 1st thread executes iterations 0 to 63,
  // 2nd thread 64 to 127
  // 3rd thread 128 to 191
  // 4th thread 192 to 255
#pragma omp for
  for(arma::uword c=0; c<nb_config; ++c){
    
    //computation of the density f_c
    f_c.fill(0.0);
    copula_c.fill(0.0);
    zeta_c.fill(0.0);
    
    // arma::ivec advanced constructor allows to use the memory address of the vector data stored externally without copying it.
    // const arma::ivec h(HconfigPtr[c], col_nb, false, true);
    // However, we don't need to embed this memory block in a arma::ivec object since it is not implied in vector/matrix operations.
    //     => We can use it as is.
    // h can be used as a pure C array.
    const int* h = HconfigPtr[c];
    
    for(arma::uword i=0; i<row_nb; ++i) {
      double sum = 0;
      for (arma::uword j = 0; j < col_nb; ++j) {
        if ( h[(size_t)j] == 0){
          sum += Logf0Mat.at(i,j);
          zeta_c.at(i,j) = zeta0.at(i,j);
        }else { // h[j] == 1
          sum += Logf1Mat.at(i, j);
          zeta_c.at(i,j) = zeta1.at(i,j);
        }
        f_c.at(i) = sum;
      }
    }
    
    for(arma::uword i=0; i<row_nb; ++i) { f_c.at(i) = exp(f_c.at(i)); }
    
    //computation of the gaussian copula density  
    for(arma::uword i=0; i<row_nb; ++i){
      copula_c.at(i) = (1/sqrt_detR) * exp(-0.5*arma::dot(zeta_c.row(i), Rinv_I*zeta_c.row(i).t()));
    }
    
    //update the prior
    // NewPrior.at(c) = mean(OldPrior.at(c)*f_c/fHconfig_sum);
    double sum = 0;
    double oldPrior_c = OldPrior.at(c);
    for(arma::uword k = 0; k<row_nb;++k){
      sum += oldPrior_c*f_c.at(k)*copula_c.at(k)/fHconfig_sum.at(k);
    }
    NewPrior.at(c) = sum/row_nb;
  }
} // end of multi-threaded section
return NewPrior;
 }
 
//' Update the estimate of R correlation matrix of the gaussian copula, parallelized version 
//'@param  Hconfig list of vector of 0 and 1, corresponding to the configurations
//'@param  fHconfig_sum a double vector containing sum_c(w_c*psi_c), obtained by fHconfig_sum_update_ptr_parallel()
//'@param  OldPrior a double vector containing the prior w_c 
//'@param  Logf0Mat a double matrix containing the log(f0(xi_q)) 
//'@param  Logf1Mat a double matrix containing the log(f1(xi_q))
//'@param  zeta0 a double matrix containing the qnorm(F0(x_iq))
//'@param  zeta1 a double matrix containing the qnorm(F1(x_iq))
//'@param  OldR a double matrix corresponding to the copula parameter
//'@param  OldRinv a double matrix corresponding to the inverse copula parameter 
//'@param  RhoIndex a int matrix containing the index of lower triangular part of a matrix
//'@param  threads_nb an int the number of threads
//'@return a double vector containing the lower triangular part of the MLE of R
// [[Rcpp::export]]
arma::mat R_MLE_update_gaussian_copula_ptr_parallel(const List & Hconfig,
                                                     const arma::vec & fHconfig_sum,
                                                     const arma::vec & OldPrior,
                                                     const arma::mat & Logf0Mat,
                                                     const arma::mat & Logf1Mat,
                                                     const arma::mat & zeta0,
                                                     const arma::mat & zeta1,
                                                     const arma::mat & OldR,
                                                     const arma::mat & OldRinv,
                                                     const arma::mat & RhoIndex,
                                                     int threads_nb = 0) {
   
   if(threads_nb < 0) threads_nb = 1;
   // The number of thread is set to the number of core available by default.
   // and limited to the number of core available if a higher number is asked.
   #ifdef _OPENMP
   int max_threads_nb = omp_get_num_procs();
   #else
   int max_threads_nb = 1;
   #endif
   if(threads_nb > max_threads_nb || threads_nb == 0) threads_nb = max_threads_nb;
   
   // Get the dimensions
   arma::uword row_nb = Logf0Mat.n_rows;
   arma::uword col_nb = Logf0Mat.n_cols;
   arma::uword nb_config = Hconfig.size();
   // Check the dimensions
   if(OldPrior.size() != nb_config)
     stop("OldPrior and Hconfig lengths are different.");
   if(Logf1Mat.n_rows != row_nb || Logf1Mat.n_cols != col_nb)
     stop("Logf0Mat and Logf1Mat Dimensions are different.");
   
   std::vector<const int*> HconfigPtr(nb_config);
   // int** HconfigPtr = new int*[nb_config];
   for(size_t i=0; i<nb_config; ++i){
     const IntegerVector v = Hconfig[i];
     HconfigPtr[i] = &(v[0]);
   }
   
   // Initialization
   arma::vec NewRho(RhoIndex.n_rows, arma::fill::zeros);
   arma::mat OldRinv_I = OldRinv - arma::eye(col_nb, col_nb);
   double sqrt_detOldR = sqrt(arma::det(OldR));
   
   // beginning of the multi-thread section. Each thread executes this section.
#pragma omp parallel num_threads(threads_nb)
{
  // NewRho_local is a thread-local object, will accumulate the results for each thread,
  // each of these thread-local objects value will later be accumulated in the final NewRho. 
  arma::vec NewRho_local(RhoIndex.n_rows, arma::fill::zeros);
  
  // xxx_c is a thread-local object (each thread has its own xxx_c)
  arma::vec f_c(row_nb);
  arma::vec copula_c(row_nb);
  arma::mat zeta_c(row_nb,col_nb);
  
  //Loop on configurations: updating prior
#pragma omp for
  for(arma::uword c=0; c<nb_config; ++c){
    
    //computation of the density f_c and zeta_c
    f_c.fill(0.0);
    copula_c.fill(0.0);
    zeta_c.fill(0.0);
    
    const int* h = HconfigPtr[c];
    
    for(arma::uword i=0; i<row_nb; ++i) {
      double sum = 0;
      for (arma::uword j = 0; j < col_nb; ++j) {
        if ( h[(size_t)j] == 0) {sum += Logf0Mat.at(i,j);
          zeta_c.at(i,j) = zeta0.at(i,j);}
        else { // h[j] == 1
          sum += Logf1Mat.at(i, j);
          zeta_c.at(i,j) = zeta1.at(i,j);
        }
        f_c.at(i) = sum;
      }
    }
    
    // f_c = exp(f_c);
    for(arma::uword i=0; i<row_nb; ++i) { f_c.at(i) = exp(f_c.at(i)); }
    
    
    //computation of the gaussian copula density  
    for(arma::uword i=0; i<row_nb; ++i) {
      copula_c.at(i) = (1/sqrt_detOldR) * exp(-0.5*arma::dot(zeta_c.row(i), OldRinv_I*zeta_c.row(i).t()));
    }
    
    //Computation of NewRho_local
    for(arma::uword index=0; index<RhoIndex.n_rows; index++){
      for(arma::uword i = 0; i<row_nb; ++i){
        NewRho_local.at(index) +=  OldPrior.at(c)*f_c.at(i)*copula_c.at(i)/fHconfig_sum.at(i)*zeta_c.at(i,RhoIndex.at(index,0)-1)*zeta_c.at(i,RhoIndex.at(index,1)-1);
      }
    }
    
  }
  // Sum the thread-local results in the final result within the critical section (one thread at a time).
#pragma omp critical
  for(arma::uword index=0; index<RhoIndex.n_rows; index++){
    NewRho.at(index) +=  NewRho_local.at(index);
  }
} // end of multi-threaded section

// Reconstruction de la matrice
arma::mat NewR(col_nb,col_nb);
for(arma::uword index=0; index<RhoIndex.n_rows; index++){
  NewR.at(RhoIndex.at(index,0)-1,RhoIndex.at(index,1)-1) = NewRho.at(index)/row_nb ;
  NewR.at(RhoIndex.at(index,1)-1,RhoIndex.at(index,0)-1) = NewRho.at(index)/row_nb ;
}

return NewR;
 }
