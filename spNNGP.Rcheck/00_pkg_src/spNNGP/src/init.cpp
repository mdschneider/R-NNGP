#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spNNGP.h"

static const R_CallMethodDef CallEntries[] = {
    {"rNNGP", (DL_FUNC) &rNNGP, 24},
    {"sNNGP", (DL_FUNC) &rNNGP, 24},
    {NULL, NULL, 0}
};

void R_init_spNNGP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
