#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP R_get_max_threads(void);
extern SEXP R_get_thread_limit(void);
extern SEXP R_get_xy_range(SEXP, SEXP, SEXP);
extern SEXP R_gower(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_gower_topn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"R_get_max_threads",  (DL_FUNC) &R_get_max_threads,  0},
    {"R_get_thread_limit", (DL_FUNC) &R_get_thread_limit, 0},
    {"R_get_xy_range",     (DL_FUNC) &R_get_xy_range,     3},
    {"R_gower",            (DL_FUNC) &R_gower,            8},
    {"R_gower_topn",       (DL_FUNC) &R_gower_topn,       9},
    {NULL, NULL, 0}
};

void R_init_gower(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
