#include <limits.h>
#include <float.h>
#include <Rcpp.h>

#include "apclusterCppHeaders.h"

using namespace Rcpp;

RcppExport SEXP apclusterC(SEXP sR, SEXP maxitsR, SEXP convitsR,
                           SEXP lamR, SEXP detailsR)
{
    NumericMatrix s(sR);
    int maxits = as<int>(maxitsR);
    int convits = as<int>(convitsR);
    double lam = as<double>(lamR);
    bool details = as<bool>(detailsR);
    int N = s.nrow();
    IntegerMatrix e(N, convits);
    IntegerVector E(N);
    IntegerVector I(N);
    IntegerVector se(N);
    NumericMatrix A(N, N);
    NumericMatrix R(N, N);
    NumericVector tmpidx(N);
    NumericVector netsimAll;
    NumericVector dpsimAll;
    NumericVector exprefAll;
    NumericMatrix idxAll;

    if (details)
    {
        netsimAll = NumericVector(maxits);
        dpsimAll  = NumericVector(maxits);
        exprefAll = NumericVector(maxits);
        idxAll    = NumericMatrix(N, maxits);
    }

    bool dn = false, unconverged = false;

    int i = 0, j, ii, K;

    while (!dn)
    {
        // first, compute responsibilities
        
        for (ii = 0; ii < N; ii++)
        {
            double max1 = -DBL_MAX, max2 = -DBL_MAX, avsim;
            int yMax;
            
            for (j = 0; j < N; j++) // determine second-largest element of AS
            {
		avsim = A(ii, j) + s(ii, j);

                if (avsim > max1)
                {
                    max2 = max1;
                    max1 = avsim;
                    yMax = j;
                }
                else if (avsim > max2)
                    max2 = avsim;
            }
            
            for (j = 0; j < N; j++) // perform update
            {
                double oldVal = R(ii, j);
                double newVal = (1 - lam) * (s(ii, j) -
                                (j == yMax ? max2 : max1)) + lam * oldVal;
                R(ii, j) = (newVal > DBL_MAX ? DBL_MAX : newVal);
            }
        }
        
        // secondly, compute availabilities
        
        for (ii = 0; ii < N; ii++)
        {
            NumericVector Rp(N);
            double auxsum = 0;
            
            for (j = 0; j < N; j++)
            {
                if (R(j, ii) < 0 && j != ii)
                    Rp[j] = 0;
                else
                    Rp[j] = R(j, ii);
                
                auxsum += Rp[j];
            }
            
            for (j = 0; j < N; j++)
            {
                double oldVal = A(j, ii);
                double newVal = auxsum - Rp[j];
                
                if (newVal > 0 && j != ii)
                    newVal = 0;
                
                A(j, ii) = (1 - lam) * newVal + lam * oldVal;
            }
        }
        
        // determine clusters and check for convergence
        
        unconverged = false;
        
        K = 0;
        
        for (ii = 0; ii < N; ii++)
        {
            int ex = (A(ii, ii) + R(ii, ii) > 0 ? 1 : 0);
            se[ii] = se[ii] - e(ii, i % convits) + ex;
            if (se[ii] > 0 && se[ii] < convits)
                unconverged = true;
            E[ii] = ex;
            e(ii, i % convits) = ex;
            K += ex;
        }
        
        if (i >= (convits - 1) || i >= (maxits - 1))
            dn = ((!unconverged && K > 0) || (i >= (maxits - 1)));
        
        if (K == 0)
        {
            if (details)
            {
                netsimAll[i] = R_NaN;
                dpsimAll[i]  = R_NaN;
                exprefAll[i] = R_NaN;
                
                for (ii = 0; ii < N; ii++)
                    idxAll(ii, i) = R_NaN;
            }
        }
        else
        {
            int cluster = 0;
            
            for (ii = 0; ii < N; ii++)
            {
                if (E[ii])
                {
                    I[cluster] = ii;
                    cluster++;
                }
            }
            
            for (ii = 0; ii < N; ii++)
            {
                if (E[ii])
                    tmpidx[ii] = (double)ii;
                else
                {
                    double maxSim = s(ii, I[0]);
                    tmpidx[ii] = (double)I[0];
                    
                    for (j = 1; j < K; j++)
                    {
                        if (s(ii, I[j]) > maxSim)
                        {
                            maxSim = s(ii, I[j]);
                            tmpidx[ii] = (double)I[j];
                        }
                    }
                }
            }
            
            if (details)
            {
                double sumPref = 0;
                
                for (j = 0; j < K; j++)
                    sumPref += s(I[j], I[j]);
                
                double sumSim = 0;
                
                for (ii = 0; ii < N; ii++)
                {
                    if (!E[ii])
                        sumSim += s(ii, (int)tmpidx[ii]);
                }
                
                netsimAll[i] = sumSim + sumPref;
                dpsimAll[i]  = sumSim;
                exprefAll[i] = sumPref;
                
                NumericMatrix::Column idxLocal = idxAll(_, i);
                idxLocal = tmpidx;
            }
        }
        
        i++;
    }

    List ret;

    ret["I"]      = I;
    ret["K"]      = K;
    ret["it"]     = IntegerVector::create(i - 1);
    ret["unconv"] = LogicalVector::create(unconverged);

    if (details)
    {
        ret["netsimAll"]  = netsimAll;
        ret["dpsimAll"]   = dpsimAll;
        ret["exprefAll"]  = exprefAll;
        ret["idxAll"]     = idxAll;
    }

    return(ret);
}
