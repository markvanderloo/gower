/*  gower - a C/R implementation of Gower's similarity (or distance) measure.
 *  Copyright (C) 2016  Mark van der Loo
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 *
 *  You can contact the author at: mark _dot_ vanderloo _at_ gmail _dot_ com
 *
 */


#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <R.h>
#include <Rdefines.h>

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

// presence or absence of a character. x and y are 0 (FALSE) or 1 (TRUE)
static inline void gower_logi(int x, int y, double *dijk, double *sijk){
   *dijk = (double) ((x | y) & !(x == NA_LOGICAL || y == NA_LOGICAL));
   *sijk = (double) (*dijk == 1) ? (x * y) : 0.0;    
}

// equality of categorical variables, encoded as x, y in {1,2,...,N}.
static inline void gower_cat(int x, int y, double *dijk, double *sijk){
   *dijk = (double) !(x == NA_INTEGER || y == NA_INTEGER);
   *sijk = (*dijk==1.0) ? (double) (x == y) : 0.0; 
}

// comparison of numerical variables, by absolute difference divided by range.
static inline void gower_num(double x, double y, double R, double *dijk, double *sijk){
  *dijk = (double) (isfinite(x) & isfinite(y));

  *sijk = (*dijk==1.0) ? (1.0-fabs(x-y)/R) : 0.0;
}



SEXP R_gower(SEXP x, SEXP y, SEXP logi_, SEXP cat_, SEXP num_, SEXP ranges_){

  int *logi = INTEGER(logi_)
    , *cat  =  INTEGER(cat_)
    , *num  =  INTEGER(num_);


  double *ranges = REAL(ranges_);



  int nrow_x = xlength(VECTOR_ELT(x, 0L))
    , nrow_y = xlength(VECTOR_ELT(y, 0L));

  SEXP out;
  out = PROTECT(allocVector(REALSXP, MAX(nrow_x, nrow_y)));
  double *dist = REAL(out);

  int n_logi = length(logi_)
    , n_cat  = length(cat_)
    , n_num  = length(num_)
    , i = 0L
    , j = 0L;

  for ( int i =0; i<n_logi; i++) logi[i]--; 
  for ( int i =0; i<n_cat; i++) cat[i]--; 
  for ( int i =0; i<n_num; i++) num[i]--; 

  double numerator, denominator, dijk, sijk;

  for ( int m=0; m<MAX(nrow_x, nrow_y); m++ ){
    numerator = denominator = 0.0;
    // compute gower similarity
    for ( int k=0; k < n_logi; k++){
      gower_logi( LOGICAL(VECTOR_ELT(x,logi[k]))[i]
                , LOGICAL(VECTOR_ELT(y,logi[k]))[j], &dijk, &sijk);
      numerator   += dijk * sijk;
      denominator += dijk;
    }

    for ( int k=0; k < n_cat; k++){
      gower_cat( INTEGER(VECTOR_ELT(x,cat[k]))[i]
             ,   INTEGER(VECTOR_ELT(y,cat[k]))[j], &dijk, &sijk);
      numerator += dijk * sijk;
      denominator += dijk;
    }
   
    for ( int k=0; k < n_num; k++){
      gower_num(REAL(VECTOR_ELT(x, num[k]))[i]
            , REAL(VECTOR_ELT(y, num[k]))[j]
            , ranges[k], &dijk, &sijk);
      numerator += dijk * sijk;
      denominator += dijk;
    }
  
    dist[m] = (denominator == 0) ? R_NaN : 1 - (numerator / denominator);
    // recycle shortest vector.
    i = (i+1 == nrow_x) ? 0 : i+1;
    j = (j+1 == nrow_y) ? 0 : j+1;

  }
  UNPROTECT(1);
  
  return out;

}



