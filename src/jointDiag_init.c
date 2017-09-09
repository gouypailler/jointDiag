#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void getW(void *, void *, void *, void *);
extern void jadiagw(void *, void *, void *, void *, void *, void *, void *, void *);
extern void sweepjedi(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"getW",      (DL_FUNC) &getW,      4},
  {"jadiagw",   (DL_FUNC) &jadiagw,   8},
  {"sweepjedi", (DL_FUNC) &sweepjedi, 6},
  {NULL, NULL, 0}
};

void R_init_jointDiag(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}