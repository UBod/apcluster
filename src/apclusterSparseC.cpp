#include <limits.h>
#include <float.h>
#include <Rcpp.h> 

#include "apclusterCppHeaders.h"

using namespace Rcpp;



RcppExport SEXP apclusterSparseC(SEXP iR, SEXP jR, SEXP valuesR, SEXP maxitsR,
				 SEXP convitsR, SEXP lamR, SEXP nR,
				 SEXP detailsR)
{
    IntegerVector s_i(iR), s_j(jR);
    NumericVector s_values(valuesR);
    int maxits = as<int>(maxitsR);
    int convits = as<int>(convitsR);
    double lam = as<double>(lamR);
    int N = as<int>(nR), M = s_i.length();
    bool details = as<bool>(detailsR);
    IntegerVector ind1(M), ind1s(N), ind1e(N);
    IntegerVector ind2(M), ind2s(N), ind2e(N);
    NumericVector A(M);
    NumericVector R(M);
    IntegerVector se(N);
    IntegerMatrix e(N, convits);
    IntegerVector E(N);
    IntegerVector I(N);
    NumericVector netsimAll;
    NumericVector dpsimAll;
    NumericVector exprefAll;
    IntegerMatrix idxAll;
    double tmpnetsim, tmpdpsim, tmpexpref; 
    IntegerVector tmpidx(N);
    
    bool dn = false, unconverged = false;
    int i, j, ii, K,  temp1, temp2, length;
    
    if (details)
    {
        netsimAll = NumericVector(maxits + 1);
        dpsimAll  = NumericVector(maxits + 1);
        exprefAll = NumericVector(maxits + 1);
        idxAll    = IntegerMatrix(N, maxits + 1);
    }

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

    temp1 = 0;
    temp2 = 0;

    int count_loop = 0;

    while (!dn)
    {
        // first, compute responsibilities

        for (ii = 0; ii < N; ii++)
        {
            double max1 = -DBL_MAX, max2 = -DBL_MAX, avsim;
            int yMax = 0;
            
            for (j = ind1s[ii]; j <= ind1e[ii]; j++)
            {
		temp1 = ind1[j];

		avsim = A[temp1] + s_values[temp1];        
		if (avsim > max1) // determine second-largest element of AS 
                {
                    max2 = max1;
                    max1 = avsim;
                    yMax = j;
                }
                else if (avsim > max2)
                    max2 = avsim;
                    
            }

            for (j = ind1s[ii]; j <= ind1e[ii]; j++)
            {
		temp1 = ind1[j];
		double oldVal = R(temp1);
		double newVal = (1 - lam) *
		    (s_values[temp1] - (j == yMax ? max2 : max1)) +
		    lam * oldVal;
                
                R(temp1) = newVal;
	    }
        }
        
        // secondly, compute availabilities
        NumericVector rp(M);

        for (ii = 0; ii < N; ii++)
        {
	    double auxsum = 0;
	    temp2 = ind2e[ii];

	    for(j = ind2s[ii]; j <= ind2e[ii]; j++)
            {
		temp1 = ind2[j];

		if (R(temp1) < 0 && j != temp2)
		    rp[j] = 0;
		else
		    rp[j] = R[temp1];
              
		auxsum += rp[j];
            }

	    for (j = ind2s[ii]; j <= ind2e[ii]; j++)
            {
                temp1 = ind2[j];

                double oldVal = A(temp1);
                double newVal = auxsum - rp[j];
                
                if (newVal > 0 && j != temp2)
                    newVal = 0;
                
                A(temp1) = (1 - lam) * newVal + lam * oldVal;
            }
        }

        // find exemplars
        temp2 = 0;
        unconverged = false;
        K = 0;
        
        for (j = (M - N); j < M; j++) //loop through the diagonal
        {
            int ex = (A(j) + R(j) > 0 ? 1 : 0);

            se[temp2] = se[temp2] - e(temp2, count_loop % convits) + ex;

            if (se[temp2] > 0 && se[temp2] < convits)
		unconverged = true;

            E[temp2] = ex;
            e(temp2, count_loop % convits) = ex;
            K += ex;
            temp2++;
        }
        
        if (count_loop >= (convits - 1) || count_loop >= (maxits - 1))
	    dn = ((!unconverged && K > 0) || (count_loop >= (maxits - 1)));
    
	// ******storage of details**********

	if (details)
	{         
	    if (K == 0)
	    {
		tmpnetsim = R_NaN;
		tmpdpsim  = R_NaN;
		tmpexpref = R_NaN;
		
		for (ii = 0; ii < N; ii++)
		    tmpidx(ii) = R_NaInt;
	    }
	    else
	    {
		length = 0;
		tmpnetsim = 0;
		tmpdpsim  = 0;
		tmpexpref = 0;
		double maxSim = 0;
		int discon = 0;

		IntegerVector I_temp(N); 

		for (ii = 0; ii < N; ii++)
		{
		    if (E[ii]) // if it is exemplar
		    {            
			tmpidx[ii] = ii;
			I_temp[length] = ii; // I = find(E)  
			length++;
		    }
		    else // non-exemplar points
		    { 
			tmpidx[ii] = R_NaInt;  
			temp1 = 0;
			temp2 = 0;
			length = 0;
			IntegerVector ee(N);

			NumericVector temp_ss(N);
			NumericVector temp_j(N);

			for(j = ind1s[ii]; j <= ind1e[ii]; j++)
			{
			    temp1 = ind1[j];
			    temp_ss[temp2]= s_values[temp1];
			    temp_j[temp2]= s_j[temp1];

			    if (E[temp_j[temp2]])
                            {  
                                ee[length] = temp2; // I = find(E(temp_j))
                                length++;
                            }

			    temp2++;
			}
                     
			if (length == 0)
			    discon = 1;
			else
			{
			    maxSim = temp_ss[ee[0]];
			    tmpidx[ii] = temp_j(ee[0]);
                        
			    for (int jj = 1; jj < length; jj++)
			    {
				temp1 = ee[jj];

				if (temp_ss(temp1) > maxSim)
                                {
                                    maxSim = temp_ss[temp1];
                                    tmpidx[ii] = temp_j(temp1);
                                }
			    }
                        
			    tmpdpsim = tmpdpsim+maxSim;
			}
		    }
		}

		I=I_temp;

		// preference sum             
		if (discon == 1)
		{
		    tmpnetsim = R_NaN;
		    tmpdpsim  = R_NaN;
		    tmpexpref = R_NaN;
		    
		    for (int jj = 0; jj < N; jj++)
			tmpidx[jj] = R_NaInt;
		}
		else
		{
		    temp1 = 0;
		    
		    for (ii = (M - N); ii < M; ii++)
		    {
			for(j = 0; j < K; j++)
			    if (temp1 == I_temp[j])
				tmpexpref += s_values(ii);

			temp1++;
		    }

                  
		    tmpnetsim = tmpdpsim + tmpexpref;
		}


	    }
	}
    
	if (details)
	{
	    exprefAll[count_loop] = tmpexpref;
	    dpsimAll[count_loop]  = tmpdpsim;
	    netsimAll[count_loop] = tmpnetsim;

	    for (ii = 0; ii < N; ii++)
		idxAll(ii, count_loop) = tmpidx[ii];
	}

	count_loop++;
    }
    // end of AP main loop
 
    // final refinement
    temp2 = 0;
    K = 0;
    
    for (j = (M - N); j < M; j++) // I think is not needed (exists in matlab)
    {
        int ex = (A(j) + R(j) > 0 ? 1 : 0);
        E[temp2] = ex;
        K += ex;
        temp2++;
    }
    
    if (K > 0)
    {
        tmpnetsim = 0;
	tmpdpsim  = 0;
	tmpexpref = 0;

        double maxSim = 0;

	// first loop finds the tmpidx if the user ask for details we can
	// skip this step and take the last tmpidx which already is computed
        for (ii = 0; ii < N; ii++)
	{
	    if (E[ii])
	    {            
		tmpidx[ii] = ii;
	    }
	    else
	    {
                tmpidx[ii] = R_NaInt;  
                NumericVector temp_ss(N);
                NumericVector temp_j(N);
                IntegerVector ee(N);
                temp1 = 0; // store the idx through loop
                temp2 = 0; // just counter
                length = 0;
               
                for(j = ind1s[ii]; j <= ind1e[ii]; j++)
		{
		    temp1 = ind1[j];
		    temp_ss[temp2] = s_values[temp1];
		    temp_j[temp2] = s_j[temp1];
                      
		    if (E[temp_j[temp2]])
		    {  
			ee[length] = temp2; // I = find(E(temp_j))
			length++;
		    }
        
		    temp2++;
		}

		maxSim = temp_ss[ee[0]];
		tmpidx[ii]=temp_j(ee[0]);

		for (int jj = 1; jj < length; jj++)
		{
		    temp1 = ee[jj];
                              
		    if (temp_ss(temp1) > maxSim)
		    {
			maxSim = temp_ss[temp1];
			tmpidx[ii] = temp_j(temp1);
		    }
		}
	    }
	}

          
	IntegerVector E_new(N);
          
	//*********************
	for (ii = 0; ii < N; ii++)
	{  
	    if (E[ii])
	    {
		IntegerVector temp_II(N);
		length = 0;

		for (int jj = 0; jj < N; jj++)
		{
                    if (tmpidx[jj] == ii)
		    {
                        temp_II[length] = jj; // I = find(E)
                        length++;
                    }
		}
                   
		NumericVector ns(N);
		NumericVector msk(N);
                  
		for (int jj = 0; jj < length; jj++) // loop only over exemplars
		{
                    temp1 = 0;
                    temp2 = 0;
                    NumericVector temp_j(N);
                    NumericVector temp_ss(N);

                    for (j = ind1s[temp_II[jj]]; j <= ind1e[temp_II[jj]]; j++)
                    {
			temp1 = ind1[j];
			temp_ss[temp2] = s_values[temp1];
			temp_j[temp2] = s_j[temp1];
			msk[temp_j[temp2]] += 1;
			ns[temp_j[temp2]] += temp_ss[temp2];
			temp2++;
		    }
		}

		IntegerVector II(length);
		IntegerVector III(length);
		int newcounter = 0, minuslength = 0;

		for (int jj = 0; jj < length; jj++)
		{
                    if (msk[temp_II[jj]] == length)
		    {
			II[newcounter] = jj;
			III[newcounter] = temp_II[II[newcounter]];
			newcounter++;
                    }
                    else
			minuslength++;
		}
                  
		maxSim = ns[III[0]];
		int index_max = 0;
		
		for (int jj = 1; jj < length - minuslength; jj++)
		{
		    temp1 = III[jj];

		    if (ns(temp1) > maxSim)
		    {
			maxSim = ns[temp1];
			index_max = jj;
		    }
		}

		E_new[III[index_max]] = 1;
            }
	}
         
	// **************************************
         
	E = E_new;
	length = 0;
	int lengthI = 0;
	tmpnetsim = 0;
	tmpdpsim  = 0;
	tmpexpref = 0;
	maxSim = 0;

	IntegerVector I_tempfinal(N);

	for (ii = 0; ii < N; ii++)
        {
            if (E[ii])
            {            
                tmpidx[ii] = ii;
                I_tempfinal[lengthI] = ii; // I = find(E) final
                lengthI++;
            }
	    else
	    {
		tmpidx[ii] = R_NaInt;  
		NumericVector temp_ss(N);
		NumericVector temp_j(N);
		IntegerVector ee(N);
		temp1 = 0;
		temp2 = 0;
		length = 0;
              
		for(j = ind1s[ii]; j <= ind1e[ii]; j++)
		{
		    temp1 = ind1[j];
		    temp_ss[temp2] = s_values[temp1];
		    temp_j[temp2] = s_j[temp1];
	
                    
		    if (E[temp_j[temp2]])
		    {  
			ee[length] = temp2;// I = find(E(temp_j))
			length++;
		    }
                      
		    temp2++;
		}
		//find max and update idx         
		maxSim = temp_ss[ee[0]];
		tmpidx[ii] = temp_j(ee[0]);

		for (int jj = 1; jj < length; jj++)
		{
		    temp1 = ee[jj];
		  
		    if (temp_ss(temp1) > maxSim)
		    {
			maxSim = temp_ss[temp1];
			tmpidx[ii] = temp_j(temp1);
		    }
		}

		tmpdpsim = tmpdpsim+maxSim;
	    }
	}

	// preference sum
        temp1 = 0;

        for (ii = (M - N); ii < M; ii++)
	{
            for(j = 0; j < K; j++)
		if(temp1 == I_tempfinal[j])
		    tmpexpref += s_values(ii);
            
            temp1++;
	}

	I = I_tempfinal;
	tmpnetsim = tmpdpsim + tmpexpref;
    }
    else
    {
	tmpnetsim = R_NaN;
	tmpdpsim  = R_NaN;
	tmpexpref = R_NaN;
	
	for (ii = 0; ii < N; ii++)
            tmpidx(ii) = R_NaInt;
    }

    if (details)
    {
	exprefAll[count_loop] = tmpexpref;
	dpsimAll[count_loop] = tmpdpsim;
	netsimAll[count_loop] = tmpnetsim;

	for (ii = 0; ii < N; ii++)
	    idxAll(ii, count_loop) = tmpidx[ii];
    }
    
    List ret; 
 
    ret["I"]         = I;
    ret["E"]         = E;
    ret["tmpidx"]    = tmpidx;
    ret["tmpnetsim"] = tmpnetsim;
    ret["tmpdpsim"]  = tmpdpsim;
    ret["tmpexpref"] = tmpexpref;
    ret["K"]         = K;
    ret["it"]        = IntegerVector::create(count_loop - 1);
    ret["unconv"]    = unconverged;
    
    if (details)
    {
        ret["netsimAll"] = netsimAll;
        ret["dpsimAll"]  = dpsimAll;
        ret["exprefAll"] = exprefAll;
        ret["idxAll"]    = idxAll;
    }
  
    return(ret);
}
