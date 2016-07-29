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
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define RECYCLE(N, K) ((N) + 1L < (K) ? (N) + 1L : 0L )

// determine when something is numerically zero.
static double EPS = 1e-8;

// presence or absence of a character. x and y are 0 (FALSE) or 1 (TRUE)
static inline void gower_logi(int *x, int nx, int *y, int ny
   , double *num, double *den)
{
  int nt = MAX(nx,ny);
  double dijk, sijk;
  int i = 0, j = 0;
  double *inum = num,  *iden=den;

  for ( int k = 0; k < nt; k++, inum++, iden++){
    dijk = (double) ((x[i] | y[j]) & !((x[i] == NA_INTEGER) | (y[j] == NA_INTEGER)));
    sijk = (dijk == 1.0) ? (double) (x[i] * y[j]) : 0.0;
    *inum += dijk * sijk; 
    *iden += dijk;
    i = RECYCLE(i, nx);
    j = RECYCLE(j, ny);
  }
  
}

// equality of categorical variables, encoded as x, y in {1,2,...,N}.
static inline void gower_cat(int *x, int nx, int *y, int ny 
  , double *num, double *den)
{

  int nt = MAX(nx, ny);
  double dijk, sijk;
  int i=0, j=0;
  double *inum = num,  *iden=den;

  for ( int k=0; k<nt; k++, inum++, iden++ ){
    dijk = (double) !(x[i] == NA_INTEGER || y[j] == NA_INTEGER);
    sijk = (dijk==1.0) ? (double) (x[i] == y[j]) : 0.0; 
    *inum += dijk * sijk; 
    *iden += dijk;
    i = RECYCLE(i, nx);
    j = RECYCLE(j, ny);
  }

}

static inline void gower_str(SEXP x, int nx, SEXP y, int ny, double *num, double *den){
  int nt = MAX(nx, ny);
  double dijk, sijk;
  int i=0, j=0;
  double *inum = num,  *iden=den;
  SEXP xi, yj;

 
  for ( int k=0; k<nt; k++, inum++, iden++ ){
    xi = STRING_ELT(x,i);
    yj = STRING_ELT(y,j);
    dijk = (double) !(xi == NA_STRING || yj == NA_STRING);
    sijk = (dijk==1.0) ? (double) (CHAR(xi) == CHAR(yj)) : 0.0; 
    *inum += dijk * sijk; 
    *iden += dijk;
    i = RECYCLE(i, nx);
    j = RECYCLE(j, ny);
  }
}


// comparison of numerical variables, by absolute difference divided by range.
static inline void gower_num(double *x, int nx, double *y, int ny,double R
    , double *num, double *den)
{
  int nt = MAX(nx, ny);
  double dijk, sijk;
  int i=0, j=0;
  double *inum = num,  *iden=den;

  if ( !isfinite(R) || R < EPS ){
    warning("skipping variable with zero or non-finite range\n");
    return;
  }

  for ( int k=0; k<nt; k++, inum++, iden++ ){
    dijk = (double) (isfinite(x[i]) & isfinite(y[j]));
    sijk = (dijk==1.0) ? (1.0-fabs(x[i]-y[j])/R) : 0.0;
    (*inum) += dijk * sijk; 
    (*iden) += dijk;

    i = RECYCLE(i, nx);
    j = RECYCLE(j, ny);
  }
}


static inline void gower_dbl_int(double *x, int nx, int *y, int ny,double R
    , double *num, double *den)
{
  int nt = MAX(nx, ny);
  double dijk, sijk;
  int i=0, j=0;
  double *inum = num,  *iden=den;

  if ( !isfinite(R) || R < EPS ){
    warning("skipping variable with zero or non-finite range\n");
    return;
  }

  for ( int k=0; k<nt; k++, inum++, iden++ ){
    dijk = (double) (isfinite(x[i]) & (y[j] != NA_INTEGER));
    sijk = (dijk==1.0) ? (1.0-fabs(x[i] - ((double) y[j]) )/R) : 0.0;
    *inum += dijk * sijk; 
    *iden += dijk;
    i = RECYCLE(i, nx);
    j = RECYCLE(j, ny);
  }
}

static inline void gower_int(int *x, int nx, int *y, int ny, double R
    , double *num, double *den)
{
  int nt = MAX(nx, ny);
  double dijk, sijk;
  int i=0, j=0;
  double *inum = num,  *iden=den;

  if ( !isfinite(R) || R == 0 ){
    warning("skipping variable with zero or non-finite range\n");
    return;
  }

  for ( int k=0; k<nt; k++, inum++, iden++ ){
    dijk = (double) ( (x[i] !=NA_INTEGER) & (y[j] != NA_INTEGER));
    sijk = (dijk==1.0) ? (1.0-fabs( ((double)x[i]) - ((double)y[j]) )/R) : 0.0;
    *inum += dijk * sijk; 
    *iden += dijk;
    i = RECYCLE(i, nx);
    j = RECYCLE(j, ny);
  }
}

// range computations
static void get_dbl_range(double *x, int nx, double *min, double *max){

  double *ix = x;

  min[0] = x[0];
  max[0] = x[0];

  for ( int i=0; i<nx; i++, ix++ ){
    if (isfinite(*min)) break;
    *min = *ix; 
    *max = *ix;
  }
  
  // all non-finite, range not computable.
  if ( !isfinite(*min) ){
    return ;
  }

  ix = x;

  for ( int i=0; i<nx; i++, ix++){
    if (isfinite(*ix)){
      if (*ix > *max){
        *max = *ix;
      } else if ( *ix < *min ){
        *min = *ix;
      }
    }
  }

}


static void get_int_range(int *x, int nx, double *min, double *max){

  int *ix = x;

  int imin = x[0]
    , imax = x[0];

  for ( int i=0; i<nx; i++, ix++ ){
    if ( imin != NA_INTEGER ) break;
    imin = *ix; 
    imax = *ix;
  }
  
  // all missing, range not computable.
  if ( imin == NA_INTEGER ){
    return ;
  }

  ix = x;

  for ( int i=0; i<nx; i++, ix++){
    if ( *ix != NA_INTEGER ){
      if (*ix > imax){
        imax = *ix;
      } else if ( *ix < imin ){
        imin = *ix;
      }
    }
  }
  *min = (double) imin;
  *max = (double) imax;
}

static void get_range(SEXP x, double *min, double *max){

  switch( TYPEOF(x) ){
    case INTSXP : {
      get_int_range(INTEGER(x), length(x), min, max);
      break;
    }
    case REALSXP : {
      get_dbl_range(REAL(x), length(x), min, max);
      break;
    }
  }

}

static double get_xy_range(SEXP x, SEXP y){

  double x_min, x_max, y_min, y_max, min, max;

  get_range(x, &x_min, &x_max);
  get_range(y, &y_min, &y_max);

  if ( isfinite(x_min) & isfinite(y_min) ){
    min = MIN(x_min, y_min);
  } else if ( isfinite(x_min) & !(isfinite(y_min)) ){
    min = x_min;
  } else if ( (!isfinite(x_min)) & isfinite(y_min) ) {
    min = y_min;
  } else {
    min = NA_REAL;
  }

  if ( isfinite(x_max) & isfinite(y_min) ){
    max = MAX(x_max, y_max);
  } else if ( isfinite(x_max) & !isfinite(y_max) ){
    max = x_max;
  } else if ( (!isfinite(x_max)) & isfinite(y_max) ){
    max = y_max;
  } else {
    max = NA_REAL;
  }

  return max - min;

}


SEXP R_gower(SEXP x, SEXP y, SEXP pair_, SEXP factor_pair_, SEXP eps_){

  int *pair = INTEGER(pair_)
    , *factor_pair = INTEGER(factor_pair_);
  int npair = length(pair_);

  // set global epsilon
  EPS = REAL(eps_)[0];


  // from R [base-1] to C [base-0] index for columns.
  //for ( int j=0; j<npair; j++) pair[j]--;

  int nrow_x = length(VECTOR_ELT(x, 0L))
    , nrow_y = length(VECTOR_ELT(y, 0L));
  int nt = MAX(nrow_x, nrow_y);

  // output
  SEXP out;
  out = PROTECT(allocVector(REALSXP, nt));

  // numerator & denominator. 
  double *num = REAL(out)
       , *den = (double *) R_alloc(nt, sizeof(double));
 
  double *iden = den, *inum = num;
  for ( int j=0; j<nt; j++, iden++, inum++){
    *iden = 0.0;
    *inum = 0.0;
  }

  int type_y;
  double R;

  // loop over coluns of x, compare with paired columns in y.
  for ( int j = 0; j < npair; j++){
    if (pair[j] == -1L) continue; // no paired column.
    switch( TYPEOF(VECTOR_ELT(x,j)) ) {
      case LGLSXP : 
        gower_logi(INTEGER(VECTOR_ELT(x,j)), nrow_x
            , INTEGER(VECTOR_ELT(y,pair[j])), nrow_y
            ,num, den);
        break;
      case REALSXP : 
        R = get_xy_range(VECTOR_ELT(x,j), VECTOR_ELT(y,pair[j]));
        if (TYPEOF(VECTOR_ELT(y,pair[j])) == REALSXP){
          gower_num(REAL(VECTOR_ELT(x,j)), nrow_x
                , REAL(VECTOR_ELT(y,pair[j])), nrow_y
                , R, num, den);
        } else if (TYPEOF(VECTOR_ELT(y,pair[j])) == INTSXP) {
          gower_dbl_int(REAL(VECTOR_ELT(x,j)), nrow_x
                , INTEGER(VECTOR_ELT(y,pair[j])), nrow_y
                , R, num, den);
        }
        break;
      case INTSXP : 
        type_y = TYPEOF(VECTOR_ELT(y,pair[j]));
        if ( type_y == REALSXP ){ // treat as numeric
          R = get_xy_range(VECTOR_ELT(x,j), VECTOR_ELT(y,pair[j]));
          gower_dbl_int(REAL(VECTOR_ELT(y,pair[j])), nrow_y
                , INTEGER(VECTOR_ELT(x, j)), nrow_x
                , R, num, den);
        } else if ( type_y == INTSXP ){
          if ( factor_pair[j] ){ // factor variables
            gower_cat(INTEGER(VECTOR_ELT(x,j)), nrow_x
                    , INTEGER(VECTOR_ELT(y,pair[j])), nrow_y
                    , num, den);
          } else { // treat as integers
            R = get_xy_range(VECTOR_ELT(x,j), VECTOR_ELT(y,pair[j]));
            gower_int(INTEGER(VECTOR_ELT(x,j)), nrow_x
                    , INTEGER(VECTOR_ELT(y,pair[j])), nrow_y
                    , R, num, den);
          }
        } 
        break;
      case STRSXP :
        gower_str(VECTOR_ELT(x,j),nrow_x,VECTOR_ELT(y,pair[j]),nrow_y, num, den);
    } // end switch
  } // end for

  inum = num;
  iden = den;
  for (int i=0; i<nt; i++, inum++, iden++){
    (*inum) = (*iden == 0.0) ? R_NaN : (1.0 - (*inum)/(*iden));
  }

  UNPROTECT(1);
  
  return out;

}



/* Push a distance down the list while keeping things sorted.
 * 
 *
 */
static inline void push(double x, int ind, double *topn, int *index, int n){
  
  for ( int i=0; i<n; i++){
    if ( (x < topn[i]) | ((x == topn[i]) & (i < n-1))){
      // new entry encountered, bubble to keep sorted.
      if ( topn[i] == R_PosInf ){ // No entry yet, overwrite.
        topn[i] = x;
        index[i] = ind;
        break; // out of main loop
      } else { // Bubble up to insert new entry
        topn[n-1] = topn[n-2];
        index[n-1] = index[n-2];
        for ( int j=n-2; j>i; j-- ){
          topn[j] = topn[j-1];
          index[j] = index[j-1];
        }
        topn[i] = x;
        index[i] = ind;
        break; // out of main loop
      }
    }
  }

} 

/* For testing purposes
SEXP R_pushdown(SEXP entry_, SEXP index_, SEXP values_, SEXP indices_){
  double  entry   = REAL(entry_)[0];
  int     index   = INTEGER(index_)[0];
  double *values  = REAL(values_);
  int    *indices = INTEGER(indices_);
  int n = length(values_);
  
  push(entry, index, values, indices, n);
  return R_NilValue;
}
*/



static inline void copyrec(SEXP into, SEXP from, int i){
  int ncol = length(into);

  SEXP col_from, col_into;

  for ( int j = 0; j < ncol; j++){
    col_from = VECTOR_ELT(from,j);
    col_into = VECTOR_ELT(into,j);
    switch(TYPEOF(col_from)){
      case LGLSXP  : { INTEGER(col_into)[0] = INTEGER(col_from)[i]; break;}
      case INTSXP  : { INTEGER(col_into)[0] = INTEGER(col_from)[i]; break;}
      case REALSXP : { REAL(col_into)[0] = REAL(col_from)[i]; break;}
      case STRSXP  : { SET_STRING_ELT(col_from, 0, STRING_ELT(col_from,i)); break;}
    }
  }
}

/* for testing purposes only
void prvec(SEXP x){
  for (int i=0; i<length(x); i++){
    Rprintf("%8.4f",REAL(x)[i]);
  }
Rprintf("\n");
}
*/

SEXP R_gower_topn(SEXP x_, SEXP y_, SEXP pair_, SEXP factor_pair_, SEXP n_, SEXP eps_){

  int n = INTEGER(n_)[0];
  int ny = length(VECTOR_ELT(y_,0));
  int nrowx = length(VECTOR_ELT(x_,0)); 
  int nout = nrowx * n;

  // setup output list
  SEXP out = allocVector(VECSXP, 2L);
  PROTECT(out);
  SET_VECTOR_ELT(out, 0L, allocVector(INTSXP, nout));
  SET_VECTOR_ELT(out, 1L, allocVector(REALSXP, nout));

  // temporary record to pass to R_gower
  SEXP temprec_ = allocVector(VECSXP, length(x_));
  PROTECT(temprec_); 
 
  // initialize output distance.
  double *vv=REAL(VECTOR_ELT(out,1L));
  int *ii=INTEGER(VECTOR_ELT(out,0L));
  for(int i=0; i<nout; i++,vv++,ii++){
    *vv = R_PosInf;
    *ii = 0L;
  }


  // pointers to output
  int *index = INTEGER(VECTOR_ELT(out, 0L));
  double *value = REAL(VECTOR_ELT(out, 1L));

  // initialize spots in temporary record
  for ( int j=0; j<length(x_); j++){
    SET_VECTOR_ELT(temprec_, j, allocVector(TYPEOF(VECTOR_ELT(x_,j)), 1L));
  }

  SEXP d;
  double *dist; 
  // loop over records in x_
  for ( int i=0; i < nrowx; i++){
    // create a list to feed to R_gower
    copyrec(temprec_, x_, i);
    // compute distances
    d = R_gower(temprec_, y_, pair_, factor_pair_, eps_);
    // push down distances & indices.
    dist = REAL(d);
    for ( int k=0; k < ny; k++){
      push(dist[k], k+1, value, index, n);
    }
    value += n;
    index += n;
  }


  // return list with indices and distances.
  UNPROTECT(2);
  return(out);
}






