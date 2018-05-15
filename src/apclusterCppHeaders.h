#ifndef __APCLUSTER_CPP_HEADERS__

#include <Rcpp.h>

#define __APCLUSTER_CPP_HEADERS__

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */

RcppExport SEXP apclusterC(SEXP sR,
			   SEXP maxitsR,
			   SEXP convitsR,
                           SEXP lamR,
			   SEXP detailsR);

RcppExport SEXP apclusterLeveragedC(SEXP sR,
				    SEXP selR,
				    SEXP maxitsR, 
                                    SEXP convitsR,
				    SEXP lamR);

RcppExport SEXP apclusterSparseC(SEXP iR,
				 SEXP jR,
				 SEXP valuesR,
				 SEXP maxitsR,
				 SEXP convitsR,
				 SEXP lamR,
				 SEXP nR,
				 SEXP detailsR);

RcppExport SEXP preferenceRangeC(SEXP sR,
				 SEXP exactR);

RcppExport SEXP preferenceRangeSparseC(SEXP iR,
				       SEXP jR,
				       SEXP valuesR,
				       SEXP nR,
				       SEXP exactR);


#endif
