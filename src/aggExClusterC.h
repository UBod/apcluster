#ifndef __AGG_EX_CLUSTER_C_HEADERS__

#include <Rcpp.h>

#define __AGG_EX_CLUSTER_C_HEADERS__

RcppExport SEXP aggExClusterC(SEXP sR,
                              SEXP KR,
                              SEXP actClustR,
                              SEXP actExemR,
                              SEXP objMatR,
                              SEXP exeMatR,
                              SEXP actLabelsR,
                              SEXP selR,
                              SEXP clustersR,
                              SEXP exemplarsR,
                              SEXP mergeR,
                              SEXP heightR,
                              SEXP preserveNamesR);

#endif