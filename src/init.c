#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */

/*
extern void comwr(void *CR, void *CI, void *LengthC, void *DR, void *DI, void *LengthD, void *HR, void *HI, void *GR, void *GI, void *LengthH, void *levels, void *firstC, void *lastC, void *offsetC, void *firstD, void *lastD, void *offsetD, void *type, void *bc, void *error);
*/
extern void wavereconsow(void *C, void *D, void *H, void *LengthH, void *levels, void *firstC, void *lastC, void *offsetC, void *firstD, void *lastD, void *offsetD, void *type, void *bc, void *error);

/* .Fortran calls */
extern void F77_NAME(ajv)(double *SNV, double *JVAL, int *ITYPE, double *GAMMA, double *DELTA, double *XLAM, double *XI, int *IFAULT);
extern void F77_NAME(jnsn)(double *XBAR, double *SD, double *RB1, double *BB2, int *ITYPE, double *GAMMA, double *DELTA, double *XLAM, double *XI, int *IFAULT);

static const R_CMethodDef CEntries[] = {
    {"wavereconsow", (DL_FUNC) &wavereconsow, 14},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"ajv",  (DL_FUNC) &F77_NAME(ajv),   8},
    {"jnsn", (DL_FUNC) &F77_NAME(jnsn), 10},
    {NULL, NULL, 0}
};

void R_init_waveband(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
