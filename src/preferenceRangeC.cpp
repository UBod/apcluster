#include <limits.h>
#include <float.h>
#include <Rcpp.h>

#include "apclusterCppHeaders.h"

using namespace Rcpp;


RcppExport SEXP preferenceRangeC(SEXP sR, SEXP exactR)
{
    NumericMatrix s(sR);
    bool exact = as<bool>(exactR);
    int N = s.nrow();
    double dpsim1 = R_NegInf, pmin = R_NegInf, pmax = R_NegInf;

    NumericVector colS(N);

    for (int j = 0; j < N; j++)
    {
	double sumOfCol = R_NegInf;

	for (int i = 0; i < N; i++)
	{
	    if (s(i, j) > R_NegInf)
	    {
		if (sumOfCol == R_NegInf)
		    sumOfCol = s(i, j);
		else
		    sumOfCol += s(i, j);

		if (s(i, j) > pmax && i != j)
		    pmax = s(i, j);
	    }
	}

	if (sumOfCol > dpsim1)
	    dpsim1 = sumOfCol;
    }

    if (dpsim1 == R_NegInf)
        pmin = R_NaN;
    else if (exact)
    {
        double dpsim2 = R_NegInf;

	for (int j21 = 0; j21 < N - 1; j21++)
	{
	    for (int j22 = j21 + 1; j22 < N; j22++)
	    {
		double tmpSum = R_NegInf;

		for (int k = 0; k < N; k++)
		{
		    double maxi = R_NegInf;

		    if (s(k, j21) > s(k, j22))
			maxi = s(k, j21);
		    else if (s(k, j22) > R_NegInf)
			maxi = s(k, j22);

		    if (maxi > R_NegInf)
		    {
			if (tmpSum == R_NegInf)
			    tmpSum = maxi;
			else
			    tmpSum += maxi;
		    }
		}

		if (tmpSum > dpsim2)
		    dpsim2 = tmpSum;
	    }
	}

        pmin = dpsim1 - dpsim2;
    }
    else
    {
	double sumM = R_NegInf, sm1 = R_PosInf, sm2 = R_PosInf;

	for (int i = 0; i < N; i++)
	{
	    colS[i] = R_NegInf;

	    for (int j = 0; j < i; j++)
		if (s(i, j) > colS[i])
		    colS[i] = s(i, j);

	    for (int j = i + 1; j < N; j++)
		if (s(i, j) > colS[i])
		    colS[i] = s(i, j);

	    if (colS[i] > R_NegInf)
	    {
		if (sumM == R_NegInf)
		    sumM = colS[i];
		else
		    sumM += colS[i];

		if (colS[i] < sm1)
		{
		    sm2 = sm1;
		    sm1 = colS[i];
		}
		else if (colS[i] < sm2)
		    sm2 = colS[i];
	    }
	}

	if (sm2 == R_PosInf || sumM == R_NegInf)
	    pmin = R_NegInf;
	else
	    pmin = dpsim1 - sumM + sm1 + sm2;
    }

    return NumericVector::create(pmin, pmax);
}
