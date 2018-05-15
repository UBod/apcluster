#include <limits.h>
#include <float.h>
#include <Rcpp.h>

#include "apclusterCppHeaders.h"


using namespace Rcpp;


RcppExport SEXP preferenceRangeSparseC(SEXP iR, SEXP jR, SEXP valuesR, SEXP nR,
				       SEXP exactR)
{
    IntegerVector s_i(iR), s_j(jR);
    NumericVector s_values(valuesR);
    int N = as<int>(nR), M = s_i.length();
    bool exact = as<bool>(exactR);
    IntegerVector ind1(M), ind1s(N), ind1e(N);
    IntegerVector ind2(M), ind2s(N), ind2e(N);
    int i, j, ii, K,  temp1, temp2, length;

    // build ind1e
    for (i = 0; i < M; i++)
    {
	ind1e[s_i[i]]++; // count ind. occurance.
	ind2e[s_j[i]]++;
    }

    temp1 = 0;
    temp2 = 0;
  
    for (i = 0; i < N; i++) // cumsum
    {
	temp1 += ind1e[i];
	ind1e[i] = temp1 - 1;
      
	temp2 += ind2e[i];
	ind2e[i] = temp2 - 1;
    }
     
    //build ind1s
    ind1s[0] = 0;
    ind2s[0] = 0;

    for (i = 1; i < N; i++)
    {
	ind1s[i] = ind1e[i - 1] + 1;
	ind2s[i] = ind2e[i - 1] + 1;
    }

    temp1 = 0;
    temp2 = 0;
    
    //build ind1
    for(i = 0; i < M; i++)
    {
	temp1 = s_i[i];
	ind1[ind1s[temp1]] = i;
	ind1s[temp1] += 1;
       
	temp2 = s_j[i];
	ind2[ind2s[temp2]] = i;
	ind2s[temp2] += 1;
    }

    //rebuild ind1s changed when build ind1
    ind1s[0] = 0;
    ind2s[0] = 0;

    for(i = 1; i < N; i++)
    {
	ind1s[i] = ind1e[i - 1] + 1;
	ind2s[i] = ind2e[i - 1] + 1;
    }

    double dpsim1 = R_NegInf, pmin = R_NegInf, pmax = R_NegInf;

    NumericVector colS(N);

    for (int j = 0; j < N; j++)
    {
	double sumOfCol = R_NegInf;

	for (int i = ind2s[j]; i <= ind2e[j]; i++)
	{
	    if (sumOfCol == R_NegInf)
		sumOfCol = s_values[ind2[i]];
	    else
		sumOfCol += s_values[ind2[i]];

	    if (s_values[ind2[i]] > pmax)
		pmax = s_values[ind2[i]];
	}

	if (sumOfCol > dpsim1)
	    dpsim1 = sumOfCol;
    }

    if (dpsim1 == R_NegInf)
        pmin = R_NaN;
    else if (exact)
    {
        double dpsim2 = R_NegInf;
	IntegerVector Index(N, -1);

	for (int j21 = 0; j21 < N - 1; j21++)
	{
	    double j21sum = R_NegInf;

	    for (int k = ind2s[j21]; k <= ind2e[j21]; k++)
	    {
		Index[s_i[ind2[k]]] = ind2[k];
		
		if (j21sum == R_NegInf)
		    j21sum = s_values[ind2[k]];
		else
		    j21sum += s_values[ind2[k]];
	    }

	    for (int j22 = j21 + 1; j22 < N; j22++)
	    {
		double tmpSum = j21sum;

		for (int k22 = ind2s[j22]; k22 <= ind2e[j22]; k22++)
		{
		    if (Index[s_i[ind2[k22]]] >= 0)
		    {
			if (s_values[ind2[k22]] >
			    s_values[Index[s_i[ind2[k22]]]])
			    tmpSum += (s_values[ind2[k22]] -
				       s_values[Index[s_i[ind2[k22]]]]);
		    }
		    else
		    {
			double maxi;

			if (s_i[ind2[k22]] == j21)
			{
			    if (s_values[ind2[k22]] > 0)
				maxi = s_values[ind2[k22]];
			    else
				maxi = 0;
			}
			else
			    maxi = s_values[ind2[k22]];

			if (tmpSum == R_NegInf)
			    tmpSum = maxi;
			else
			    tmpSum += maxi;
		    }
		}

		if (Index[j22] >= 0 && tmpSum > R_NegInf &&
		    s_values[Index[j22]] < 0)
		    tmpSum -= s_values[Index[j22]];

		if (tmpSum > dpsim2)
		    dpsim2 = tmpSum;
	    }

	    for (int k = ind2s[j21]; k <= ind2e[j21]; k++)
		    Index[s_i[ind2[k]]] = -1;
	}

        pmin = dpsim1 - dpsim2;
    }
    else
    {
	double sumM = R_NegInf, sm1 = R_PosInf, sm2 = R_PosInf;

	for (int i = 0; i < N; i++)
	{
	    colS[i] = R_NegInf;

	    for (int j = ind1s[i]; j <= ind1e[i]; j++)
		if (s_values[ind1[j]] > colS[i])
		    colS[i] = s_values[ind1[j]];

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
