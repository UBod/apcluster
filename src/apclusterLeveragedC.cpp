#include <limits.h>
#include <float.h>
#include <Rcpp.h>
#include <R_ext/Utils.h>

#include "apclusterCppHeaders.h"

using namespace Rcpp;


RcppExport SEXP apclusterLeveragedC(SEXP sR, SEXP selR, SEXP maxitsR, 
                                    SEXP convitsR, SEXP lamR)
{
    NumericMatrix s(sR);
    IntegerVector sel(selR);
    int maxits = as<int>(maxitsR);
    int convits = as<int>(convitsR);
    double lam = as<double>(lamR);
    int M = s.ncol(); int N = s.nrow();
    IntegerMatrix e(N, convits);
    IntegerVector I(N);
    IntegerVector se(N);
    NumericMatrix A(M, N);
    NumericMatrix R(M, N);
    NumericVector auxsum(M - 1);
    bool dn = false, unconverged = false;

    int i = 0, j, ii, K;

    while (!dn)
    {
        // first, compute responsibilities        
        for (ii = 0; ii < N; ii++)
        {
            double max1 = -DBL_MAX, max2 = -DBL_MAX;
            double avsim;
            int yMax;

            // determine largest and second largest element of A + S
            for (j = 0; j < M; j++)
            {
                if (j < M - 1 && sel[j] == ii)
                    continue;

                avsim = A(j, ii) + s(ii, j);

                if (avsim > max1)
                {
                    max2 = max1;
                    max1 = avsim;
                    yMax = j;
                }
                else if (avsim > max2)
                    max2 = avsim;
            }

            for (j = 0; j < M; j++) // R update including self responsibilities
            {
                if (j < M - 1 && sel[j] == ii)
                    continue;

                double newVal = (1 - lam) * (s(ii, j) -
                                (j == yMax ? max2 : max1)) + lam * R(j, ii);
                R(j, ii) = (newVal > DBL_MAX ? DBL_MAX : newVal);
                
                if (R(j, ii) > 0 && j < M - 1) 
                    auxsum[j] = auxsum[j] + R(j, ii);
            }
        }

        // correct auxsum with diag elements
        for (ii = 0; ii < M - 1; ii++)
                auxsum[ii] = auxsum[ii] + R(M - 1, sel[ii]);

        // secondly, compute availabilities        
        for (ii = 0; ii < M - 1; ii++)
        {
            for (j = 0; j < N; j++)
            {
                double newVal = auxsum[ii];

                if (R(ii, j) > 0)
                    newVal -= R(ii, j);

                if (sel[ii] == j)
                {
                    // update diagonal element back in last col
                    A(M - 1, j) = (1 - lam) * (newVal - R(M - 1, j)) + 
                                  lam * A(M-1, j);
                    newVal = 0; // set real diag elmenent to 0 - oldval is 0
                }
                else 
                {
                    if (newVal > 0)
                        newVal = 0;
                }

                A(ii, j) = (1 - lam) * newVal + lam * A(ii, j);
            }
            auxsum[ii] = 0;
        }

        // determine clusters and check for convergence
        unconverged = false;
        K = 0;

        for (j = 0; j < N; j++)
        {
            int ex = (A(M - 1, j) + R(M - 1, j) > 0 ? 1 : 0);
            se[j] = se[j] - e(j, i % convits) + ex;
            if (se[j] > 0 && se[j] < convits)
                unconverged = true;

            e(j, i % convits) = ex;
            if (ex)
                I[K] = j;
            K += ex;
        }

        if (i >= (convits - 1) || i >= (maxits - 1))
            dn = ((!unconverged && K > 0) || (i >= (maxits - 1)));

        i++;
    }

    List ret;

    ret["I"]      = I;
    ret["K"]      = K;
    ret["it"]     = IntegerVector::create(i - 1);
    ret["unconv"] = LogicalVector::create(unconverged);

    return(ret);
}
