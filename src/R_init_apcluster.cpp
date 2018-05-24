#include "apclusterCppHeaders.h"

extern "C"
{
#include "distanceL.h"
#include "aggExClusterC.h"
}

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static const R_CallMethodDef callMethods[] = {
    {"aggExClusterC", (DL_FUNC) &aggExClusterC, 13},
    {"apclusterC", (DL_FUNC) &apclusterC, 5},
    {"apclusterLeveragedC", (DL_FUNC) &apclusterLeveragedC, 5},
    {"apclusterSparseC", (DL_FUNC) &apclusterSparseC, 8},
    {"preferenceRangeC", (DL_FUNC) &preferenceRangeC, 2},
    {"preferenceRangeSparseC", (DL_FUNC) &preferenceRangeSparseC, 5},
    {"CdistR", (DL_FUNC) &CdistR, 4},
    {NULL, NULL, 0}
};

extern "C"
{
    void attribute_visible R_init_apcluster(DllInfo *info) {
		/* Register routines, allocate resources. */
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
    }

    void R_unload_apcluster(DllInfo *info) {
		/* Release resources. */
    }
}
