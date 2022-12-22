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

/* R-windows-oldrel (3.2.x) uses gcc 4.6.3  which we need to detect */
#ifdef __GNUC__
#if __GNUC__ <= 4 && __GNUC_MINOR__ <= 6
#else
#define HAS_REDUCTION
#endif
#endif




#include <math.h>
#include <R.h>
#include <Rdefines.h>

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define RECYCLE(N, K) ((N) + 1L < (K) ? (N) + 1L : 0L )

// determine when something is numerically zero.
static double EPS = 1e-8;
// set max nr of threads to use
static int NTHREAD = 1L;

// recycling in parallel region.
static inline int recycle(int i, int nthreads, int ni){
  i += nthreads;
  if ( i >= ni )
    i = (nthreads < ni) ? (i - ni) : (i % ni);
  return i;
}


SEXP R_get_max_threads(void){
  SEXP out = allocVector(INTSXP, 1L);
  PROTECT(out);
  INTEGER(out)[0] = 1L;
  #ifdef _OPENMP
  INTEGER(out)[0] = omp_get_max_threads();
  #endif
  UNPROTECT(1);
  return out;
}


SEXP R_get_num_procs(void){
  SEXP out = allocVector(INTSXP, 1L);
  PROTECT(out);
  INTEGER(out)[0] = 1L;
  #ifdef _OPENMP
  INTEGER(out)[0] = omp_get_num_procs();
  #endif
  UNPROTECT(1);
  return out;
}

SEXP R_get_thread_limit(void){
  SEXP out = allocVector(INTSXP, 1L);
  PROTECT(out);
  INTEGER(out)[0] = 1L;
  #ifdef _OPENMP
  INTEGER(out)[0] = omp_get_thread_limit();
  #endif
  UNPROTECT(1);
  return out;
}


// presence or absence of a character. x and y are 0 (FALSE) or 1 (TRUE)
static inline void gower_logi(int *x, int nx, int *y, int ny
   , double *num, double *den, double weight)
{

  #pragma omp parallel num_threads(NTHREAD) 
  {
    int nt = MAX(nx,ny);
    double dijk, sijk;
    int i = 0, j = 0;
    double *inum = num,  *iden=den;
    
    int ID = 0, num_threads=1;
    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, nx);
    j = recycle(ID-num_threads, num_threads, ny);
    inum += ID;
    iden += ID;
    #endif
    

    for ( int k = ID; k < nt; k += num_threads, inum += num_threads, iden += num_threads){
      dijk = (double) ((x[i] | y[j]) & !((x[i] == NA_INTEGER) | (y[j] == NA_INTEGER)));
      sijk = (dijk == 1.0) ? (double) (x[i] * y[j]) : 0.0;
      *inum += weight * dijk * sijk; 
      *iden += weight * dijk;
      i = recycle(i, num_threads, nx);
      j = recycle(j, num_threads, ny);
    }
  } // end parallel region.
}

// equality of categorical variables, encoded as x, y in {1,2,...,N}.
static inline void gower_cat(int *x, int nx, int *y, int ny 
  , double *num, double *den, double weight)
{

  #pragma omp parallel num_threads(NTHREAD)
  {
    int nt = MAX(nx,ny);
    double dijk, sijk;
    int i = 0, j = 0;
    double *inum = num,  *iden=den;
    
    int ID = 0, num_threads=1;
    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, nx);
    j = recycle(ID-num_threads, num_threads, ny);
    inum += ID;
    iden += ID;
    #endif

    for ( int k = ID; k < nt; k += num_threads, inum += num_threads, iden += num_threads){
      dijk = (double) !(x[i] == NA_INTEGER || y[j] == NA_INTEGER);
      sijk = (dijk==1.0) ? (double) (x[i] == y[j]) : 0.0; 
      *inum += weight * dijk * sijk; 
      *iden += weight * dijk;
      i = recycle(i, num_threads, nx);
      j = recycle(j, num_threads, ny);
    }
  } // end parallel region.

}

// strings. Treated as categories.
static inline void gower_str(SEXP x, int nx, SEXP y, int ny, double *num, double *den, double weight){
  #pragma omp parallel num_threads(NTHREAD)
  {
    int nt = MAX(nx, ny);
    double dijk, sijk;
    int i=0, j=0;
    double *inum = num,  *iden=den;
    SEXP xi, yj;

    int ID = 0, num_threads=1;
    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, nx);
    j = recycle(ID-num_threads, num_threads, ny);
    inum += ID;
    iden += ID;
    #endif
   
    for ( int k = ID; k < nt; k += num_threads, inum += num_threads, iden += num_threads){
      xi = STRING_ELT(x,i);
      yj = STRING_ELT(y,j);
      dijk = (double) !(xi == NA_STRING || yj == NA_STRING);
      sijk = (dijk==1.0) ? (double) (CHAR(xi) == CHAR(yj)) : 0.0; 
      *inum += weight * dijk * sijk; 
      *iden += weight * dijk;
      i = recycle(i, num_threads, nx);
      j = recycle(j, num_threads, ny);
    }
  } // end of parallel region 
}


// comparison of numerical variables, by absolute difference divided by range.
static inline void gower_num(double *x, int nx, double *y, int ny,double R
    , double *num, double *den, double weight)
{
  if ( !isfinite(R) || R < EPS ){
    warning("skipping variable with zero or non-finite range.");
    return;
  } 
  #pragma omp parallel num_threads(NTHREAD)
  {
    int nt = MAX(nx,ny);
    double dijk, sijk;
    int i = 0, j = 0;
    double *inum = num,  *iden=den;
    
    int ID = 0, num_threads=1;
    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, nx);
    j = recycle(ID-num_threads, num_threads, ny);
    inum += ID;
    iden += ID;
    #endif


    for ( int k = ID; k < nt; k += num_threads, inum += num_threads, iden += num_threads){
      dijk = (double) (isfinite(x[i]) & isfinite(y[j]));
      sijk = (dijk==1.0) ? (1.0-fabs(x[i]-y[j])/R) : 0.0;
      (*inum) += weight * dijk * sijk; 
      (*iden) += weight * dijk;
      i = recycle(i, num_threads, nx);
      j = recycle(j, num_threads, ny);
    }
    
  } // end of parallel region
}


static inline void gower_dbl_int(double *x, int nx, int *y, int ny,double R
    , double *num, double *den, double weight)
{

  if ( !isfinite(R) || R < EPS ){
    warning("skipping variable with zero or non-finite range\n");
    return;
  }

  #pragma omp parallel num_threads(NTHREAD)
  {
    int nt = MAX(nx, ny);
    double dijk, sijk;
    int i=0, j=0;
    double *inum = num,  *iden=den;

    int ID = 0, num_threads=1;
    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, nx);
    j = recycle(ID-num_threads, num_threads, ny);
    inum += ID;
    iden += ID;
    #endif
   
    for ( int k = ID; k < nt; k += num_threads, inum += num_threads, iden += num_threads){
      dijk = (double) (isfinite(x[i]) & (y[j] != NA_INTEGER));
      sijk = (dijk==1.0) ? (1.0-fabs(x[i] - ((double) y[j]) )/R) : 0.0;
      *inum += weight * dijk * sijk; 
      *iden += weight * dijk;
      i = recycle(i, num_threads, nx);
      j = recycle(j, num_threads, ny);
    }
  } // end of parallel region
}

static inline void gower_int(int *x, int nx, int *y, int ny, double R
    , double *num, double *den, double weight)
{
  if ( !isfinite(R) || R == 0 ){
    warning("skipping variable with zero or non-finite range\n");
    return;
  }
  #pragma omp parallel num_threads(NTHREAD)
  {
    int nt = MAX(nx, ny);
    double dijk, sijk;
    int i=0, j=0;
    double *inum = num,  *iden=den;

    int ID = 0, num_threads=1;
    #ifdef _OPENMP
    ID = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    i = recycle(ID-num_threads, num_threads, nx);
    j = recycle(ID-num_threads, num_threads, ny);
    inum += ID;
    iden += ID;
    #endif


    for ( int k = ID; k < nt; k += num_threads, inum += num_threads, iden += num_threads){
      dijk = (double) ( (x[i] !=NA_INTEGER) & (y[j] != NA_INTEGER));
      sijk = (dijk==1.0) ? (1.0-fabs( ((double)x[i]) - ((double)y[j]) )/R) : 0.0;
      *inum += weight * dijk * sijk; 
      *iden += weight * dijk;
      i = recycle(i, num_threads, nx);
      j = recycle(j, num_threads, ny);
    }
  } // end of parallel region
}

// range computations
static void get_dbl_range(double *x, int nx, double *min, double *max){

  double *ix = x;


  double imin=*ix, imax=*ix;

  for ( int i=0; i<nx; i++, ix++ ){
    imin = *ix; 
    imax = *ix;
    if (isfinite(*ix)) break;
  }
  
  // all non-finite, range not computable.
  if ( !isfinite(imin) ){
    return ;
  }

  #ifdef HAS_REDUCTION
  #pragma omp parallel num_threads(NTHREAD)
  #endif
  {
    #ifdef HAS_REDUCTION
    #pragma omp for reduction(min:imin), reduction(max:imax)
    #endif
    for ( int i=0; i<nx; i++){
      if (isfinite(x[i])){
        if ( x[i] > imax ) imax = x[i];
        if ( x[i] < imin ) imin = x[i];
      }
    }
  }// end parallel region
  *min = imin;
  *max = imax;
}


static void get_int_range(int *x, int nx, double *min, double *max){

  int *ix = x;

  int imin = *ix
    , imax = *ix;

  for ( int i=0; i<nx; i++, ix++ ){
    imin = *ix; 
    imax = *ix;
    if ( imin != NA_INTEGER ) break;
  }
  
  // all missing, range not computable.
  if ( imin == NA_INTEGER ){
    return ;
  }

   
  #ifdef HAS_REDUCTION
  #pragma omp parallel for reduction(min:imin), reduction(max:imax)
  #endif
  for ( int i=0; i<nx; i++){
    if ( x[i] != NA_INTEGER ){
      if ( x[i] > imax ) imax = x[i];
      if ( x[i] < imin ) imin = x[i];
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

  double x_min = R_NegInf
      , x_max = R_PosInf
      , y_min = R_NegInf
      , y_max = R_PosInf
      , min, max;

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

SEXP R_get_xy_range(SEXP x_, SEXP y_, SEXP nthread_){
  NTHREAD = INTEGER(nthread_)[0];
  SEXP out = allocVector(REALSXP,1L);
  PROTECT(out);
  REAL(out)[0] = get_xy_range(x_, y_);
  UNPROTECT(1);
  return out;
}


/*
static void print_vec(double *x, int n){
  for ( int i=0; i<n; i++){
    Rprintf("(%d) = %4.3f, ",i,x[i]);
  }
  Rprintf("\n");
}*/

static void do_gower(
  SEXP x               // a data.frame
  , SEXP y             // another data.frame
  , SEXP ranges_       // dbl; ranges[i] is range of variable (x[i],y[pair[i]])
  , SEXP pair_         // int; pair[i] is index in y
  , SEXP factor_pair_  // int; 0 if not factor
  , SEXP eps_          // dbl; numerical zero
  , SEXP weights_      // dbl; weights[i] is weight for variable i in distance
  , SEXP nthread_      // int; requested nr of threads
  , double *work       // dbl; of length max(nrow(x),nrow(y))
  , SEXP out_)         // dbl; output, length equals work.
{

  int *pair = INTEGER(pair_)
    , *factor_pair = INTEGER(factor_pair_);
  int npair = length(pair_);

  double *ranges = REAL(ranges_);
  double *weights = REAL(weights_);

  NTHREAD = INTEGER(nthread_)[0];

  // set global epsilon
  EPS = REAL(eps_)[0];


  int nrow_x = length(VECTOR_ELT(x, 0L))
    , nrow_y = length(VECTOR_ELT(y, 0L));
  int nt = MAX(nrow_x, nrow_y);

  // no point paralellizing over a small number of records,
  // so let's save some system calls.


  // numerator & denominator. 
  double *num = REAL(out_)
       , *den = work;

  // initialize work- and output space
  double *iden = den, *inum = num;
  for ( int j=0; j<nt; j++, iden++, inum++){
    *iden = 0.0;
    *inum = 0.0;
  }

  int type_y;
  double R;

  // loop over columns of x, compare with paired columns in y.
  for ( int j = 0; j < npair; j++){
    if (pair[j] == -1L) continue; // no paired column.
    
    // Get the weight value for this column
    double weight = weights[j];

    switch( TYPEOF(VECTOR_ELT(x,j)) ) {
      case LGLSXP : 
        gower_logi(INTEGER(VECTOR_ELT(x,j)), nrow_x
            , INTEGER(VECTOR_ELT(y,pair[j])), nrow_y
            ,num, den, weight);
        break;
      case REALSXP : 
        R = ranges[j];
        if (TYPEOF(VECTOR_ELT(y,pair[j])) == REALSXP){
          gower_num(REAL(VECTOR_ELT(x,j)), nrow_x
                , REAL(VECTOR_ELT(y,pair[j])), nrow_y
                , R, num, den, weight);
        } else if (TYPEOF(VECTOR_ELT(y,pair[j])) == INTSXP) {
          gower_dbl_int(REAL(VECTOR_ELT(x,j)), nrow_x
                , INTEGER(VECTOR_ELT(y,pair[j])), nrow_y
                , R, num, den, weight);
        }
        break;
      case INTSXP : 
        type_y = TYPEOF(VECTOR_ELT(y,pair[j]));
        if ( type_y == REALSXP ){ // treat as numeric
          R = ranges[j];
          gower_dbl_int(REAL(VECTOR_ELT(y,pair[j])), nrow_y
                , INTEGER(VECTOR_ELT(x, j)), nrow_x
                , R, num, den, weight);
        } else if ( type_y == INTSXP ){
          if ( factor_pair[j] ){ // factor variables
            gower_cat(INTEGER(VECTOR_ELT(x,j)), nrow_x
                    , INTEGER(VECTOR_ELT(y,pair[j])), nrow_y
                    , num, den, weight);
          } else { // treat as integers
            R = ranges[j];
            gower_int(INTEGER(VECTOR_ELT(x,j)), nrow_x
                    , INTEGER(VECTOR_ELT(y,pair[j])), nrow_y
                    , R, num, den, weight);
          }
        } 
        break;
      case STRSXP :
        gower_str(VECTOR_ELT(x,j),nrow_x,VECTOR_ELT(y,pair[j]),nrow_y, num, den, weight);
    } // end switch
  } // end for

  inum = num;
  iden = den;
  for (int i=0; i<nt; i++, inum++, iden++){
    (*inum) = (*iden == 0.0) ? R_NaN : (1.0 - (*inum)/(*iden));
  }

}



SEXP R_gower(SEXP x, SEXP y, SEXP ranges_, SEXP pair_
    , SEXP factor_pair_, SEXP eps_, SEXP _weights, SEXP nthread_){


  int nrow_x = length(VECTOR_ELT(x, 0L))
    , nrow_y = length(VECTOR_ELT(y, 0L));
  int nt = MAX(nrow_x, nrow_y);

  // Setup space for output.
  SEXP out_;
  out_ = PROTECT(allocVector(REALSXP, nt));
  // Make room for the workhorse
  double *work = (double *) R_alloc(nt, sizeof(double));
  // Make the horse work
  do_gower(x, y, ranges_, pair_, factor_pair_, eps_, _weights, nthread_, work, out_);

  // cleanup and return
  UNPROTECT(1);
  return(out_);

}



/* Push a value and an index down the list while keeping things sorted.
 * 
 * x    : value to push down
 * ind  : index value
 * topn : list, initialized with 1./0.
 * index: list, initialized with 0L
 * n    : nr of values in the list.
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

SEXP R_gower_topn(SEXP x_, SEXP y_, SEXP ranges_, SEXP pair_
    , SEXP factor_pair_, SEXP n_, SEXP eps_, SEXP weights_, SEXP nthread_){

  int n = INTEGER(n_)[0];
  int ny = length(VECTOR_ELT(y_,0));
  int nrowx = length(VECTOR_ELT(x_,0)); 
  int nout = nrowx * n;
  int nt = MAX(ny,nrowx);

  // setup output list
  SEXP out = allocVector(VECSXP, 2L);
  PROTECT(out);
  SET_VECTOR_ELT(out, 0L, allocVector(INTSXP, nout));
  SET_VECTOR_ELT(out, 1L, allocVector(REALSXP, nout));

  // temporary record to pass to R_gower
  SEXP temprec_ = allocVector(VECSXP, length(x_));
  PROTECT(temprec_);

  // Room for do_gower workhorse
  double *work = (double *) R_alloc(nt, sizeof(double));
 
  // initialize output distance.
  double *vv=REAL(VECTOR_ELT(out,1L));
  int *ii=INTEGER(VECTOR_ELT(out,0L));
  #pragma omp parallel num_threads(NTHREAD) 
  { 
    #pragma omp for schedule(static)
    for(int i=0; i<nout; i++){
      vv[i] = R_PosInf;
      ii[i] = 0L;
    }
  }


  // pointers to output
  int *index = INTEGER(VECTOR_ELT(out, 0L));
  double *value = REAL(VECTOR_ELT(out, 1L));

  // initialize spots in temporary record
  for ( int j=0; j<length(x_); j++){
    SET_VECTOR_ELT(temprec_, j, allocVector(TYPEOF(VECTOR_ELT(x_,j)), 1L));
  }

  // temporarily stores distance values.
  SEXP d_ = allocVector(REALSXP,nt);
  PROTECT(d_);
  double *dist; 
  // loop over records in x_
  for ( int i=0; i < nrowx; i++){
    // create a list to feed to R_gower
    copyrec(temprec_, x_, i);
    // compute distances
    do_gower(temprec_, y_, ranges_, pair_, factor_pair_, eps_, weights_, nthread_, work, d_);
    // push down distances & indices.
    dist = REAL(d_);
    for ( int k=0; k < ny; k++, dist++){
      push(*dist, k+1, value, index, n);
    }
    value += n;
    index += n;
  }


  // return list with indices and distances.
  UNPROTECT(3);
  return(out);
}






