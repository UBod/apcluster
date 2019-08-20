#include <Rcpp.h>
#include "aggExClusterC.h"

using namespace Rcpp;

IntegerVector concat(IntegerVector x, IntegerVector y)
{
    IntegerVector res(x.size() + y.size());
  
    std::copy(x.begin(), x.end(), res.begin());
    std::copy(y.begin(), y.end(), res.begin() + x.size());
  
    return res;
}

NumericMatrix subsetMatrix(NumericMatrix x, IntegerVector col,
                            IntegerVector row)
{
    NumericMatrix res(col.length(), row.length());
  
    for (int i = 0; i < col.length(); i++)
    {
        for (int j = 0; j < row.length(); j++)
        {
            res(i,j) = x(col[i] - 1, row[j] - 1);
        }
    }
    
    return res;
}

NumericVector subsetMatrixToVec(NumericMatrix x, int col, IntegerVector row)
{
    NumericVector res(row.length());
    
    for (int i = 0; i < row.length(); i++)
    {
        res[i] = x(col - 1, row[i] - 1);
    }
    
    return res;
}

NumericVector colMeans(NumericMatrix x)
{
    NumericVector res(x.ncol());
    
    for (int i = 0; i < x.ncol(); i++)
    {
        res[i] = mean(x(_, i));
    }
    return res;
}

int which_max_NoNA(NumericVector x)
{
    int index;
    double max = R_NegInf;
    
    for (int i = 0; i < x.size(); i++)
    {
        double value = x[i];

        if(R_IsNA(value))
        {
            continue;
        }
        if(value > max)
        {
            max = value;
            index = i;
        }
    }
    return index;
}


RcppExport SEXP aggExClusterC(SEXP sR, SEXP KR, SEXP actClustR, SEXP actExemR,
    SEXP objMatR, SEXP exeMatR, SEXP actLabelsR, SEXP selR, SEXP clustersR,
    SEXP exemplarsR, SEXP mergeR, SEXP heightR, SEXP preserveNamesR)
{
    NumericMatrix s(sR);
    int K = as<int>(KR);
    List actClust(actClustR);
    IntegerVector actExem(actExemR);
    NumericMatrix objMat(objMatR);
    IntegerMatrix exeMat(exeMatR);
    IntegerVector actLabels(actLabelsR);
    IntegerVector sel(selR);
    List clusters(clustersR);
    List exemplars(exemplarsR);
    IntegerMatrix merge(mergeR);
    NumericVector height(heightR);
    bool preserveNames(preserveNamesR);
    
    IntegerVector colInd(s.nrow());
    
    if (sel.length() > 0)
    {
        for (int i = 0; i < sel.length(); i++)
        {
            colInd[sel[i] - 1] = i + 1;
        }
    }
    
    List ret;
    
    // compute complete matrices before starting joining
    for (int i = 0; i < K - 1; i++)
    {
        for (int j = i + 1; j < K; j++)
        {
            IntegerVector joint = concat(actClust[i], actClust[j]);
            
            if(sel.length() > 0)
            {
                IntegerVector inters = intersect(sel, joint);
                IntegerVector ci = colInd[inters - 1];

                if(ci.length() > 0)
                {
                    NumericVector cM = colMeans(subsetMatrix(s, joint, ci));
                    
                    int ex = inters[which_max(cM)];
                    exeMat(i, j) = ex;
                    objMat(i, j) = (mean(subsetMatrixToVec(s, ex, colInd[
                        intersect(sel, IntegerVector(actClust[i])) - 1])) +
                        mean(subsetMatrixToVec(s, ex, colInd[intersect(sel,
                        IntegerVector(actClust[j])) - 1]))) / 2;
                }
                else
                {
                    // joining not possible - no similarities available
                    ret["error"] = 1;
                    return(ret);
                }
            }
            else
            {
                NumericVector cM = colMeans(subsetMatrix(s, joint, joint));
                
                int ex = joint[which_max(cM)];
                
                exeMat(i, j) = ex;
                objMat(i, j) = (mean(subsetMatrixToVec(s, ex, IntegerVector(
                actClust[i]))) + mean(subsetMatrixToVec(s, ex, IntegerVector(
                actClust[j])))) / 2;
            }
        }
    }
    
    // agglomeration loop
    for (int k = K - 1; k > 0; k--)
    {
        int tojoin = which_max_NoNA(objMat);

        int I = tojoin % K;
        int J = std::floor(tojoin / K);
        
	IntegerVector newClust = concat(actClust[I], actClust[J]);
	IntegerVector newClustNM =
	    MAYBE_REFERENCED(newClust) ? clone(newClust) : newClust;
	newClust.names() = CharacterVector(newClustNM);
	
        LogicalVector rem(actClust.length(), true);
        rem[I] = false;
	rem[J] = false;
        actClust = actClust[rem];
        
        if (actClust.length() < (k - 1))
        {
            actClust[k - 1] = newClust;
        }
        else
        {
            actClust.insert(k - 1, newClust);
        }
        
        actExem = actExem[(actExem != actExem[I]) & (actExem != actExem[J])];
        actExem.push_back(exeMat(I, J));
        
        clusters[k - 1] = actClust;
        
        merge((K - k - 1), 0) = actLabels[I];
        merge((K - k - 1), 1) = actLabels[J];
        
        actLabels = actLabels[(actLabels != actLabels[I]) &
            (actLabels != actLabels[J])];
        actLabels.push_back(K - k);
        height[K - k - 1] = objMat(I, J);
        exemplars[k - 1] = actExem;

        if (preserveNames && !Rf_isNull(colnames(s)) && 
            (Rf_length(colnames(s)) > 0))
        {
            IntegerVector(exemplars[k-1]).names() = ifelse(
                actExem <= as<CharacterVector>(colnames(s)).length(),
                CharacterVector(actExem), NA_STRING);
        }
        
        if (k == 1)
        {
            break;
        }
        
        // rearrange matrices objMat and exeMat
        // put values for unchanged clusters in the first k-1 rows/columns
        IntegerVector indexVec = seq_len(k + 1);
        indexVec = indexVec[(indexVec != indexVec[I]) &
            (indexVec != indexVec[J])];
        
        for (int i = 0; i < k - 1; i++)
        {
            for (int j = 0; j < k - 1; j++)
            {
                exeMat(i,j) = exeMat(indexVec[i] - 1, indexVec[j] - 1);
                objMat(i,j) = objMat(indexVec[i] - 1, indexVec[j] - 1);
            }
        }
        
        // wipe out k+1-st column
        for (int i = 0; i < exeMat.nrow(); i++)
        {
            exeMat(i, k) = NA_INTEGER;
            objMat(i, k) = NA_REAL;
        }
        
        // update k-th column with objective values and joint exemplars of
        // unchanged clusters and the newly joined cluster
        for (int i = 1; i < k; i++)
        {
            IntegerVector joint = concat(actClust[i-1], actClust[k-1]);
            
            if(sel.length() > 0)
            {
                IntegerVector inters = intersect(sel, joint);
                IntegerVector ci = colInd[inters - 1];
                
                if(ci.length() > 0)
                {
                    NumericVector cM = colMeans(subsetMatrix(s, joint, ci));
                    
                    int ex = inters[which_max(cM)];
                    
                    exeMat(i - 1, k - 1) = ex;
                    objMat(i - 1, k - 1) = (mean(subsetMatrixToVec(s, ex, colInd[
                        intersect(sel, IntegerVector(actClust[i - 1])) - 1])) +
                        mean(subsetMatrixToVec(s, ex, colInd[intersect(sel,
                        IntegerVector(actClust[k - 1])) - 1]))) / 2;
                }
                else
                {
                    // joining not possible - no similarities available
                    ret["error"] = 2;
                    return(ret);
                }
            }
            else
            {
                NumericVector cM = colMeans(subsetMatrix(s, joint, joint));
                
                int ex = joint[which_max(cM)];
                
                exeMat(i - 1, k - 1) = ex;
                objMat(i - 1, k - 1) = (mean(subsetMatrixToVec(s, ex, 
                    IntegerVector(actClust[i - 1]))) + mean(subsetMatrixToVec(s,
                    ex, IntegerVector(actClust[k - 1])))) / 2;
            }
        }
    }
    
    ret["exeMat"] = exeMat;
    ret["objMat"] = objMat;
    ret["merge"] = merge;
    ret["height"] = height;
    ret["clusters"] = clusters;
    
    if(sel.length() > 0)
        ret["colInd"] = colInd;
    
    return(ret);
}
