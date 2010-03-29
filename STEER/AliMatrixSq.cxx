/**********************************************************************************************/
/* Abstract class for matrix used for millepede2 operation.                                   */
/* Author: ruben.shahoyan@cern.ch                                                             */
/*                                                                                            */ 
/**********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
//
#include "TClass.h"
#include "TMath.h"
#include "AliMatrixSq.h"
//

using namespace std;

ClassImp(AliMatrixSq)



//___________________________________________________________
void AliMatrixSq::MultiplyByVec(const Double_t *vecIn,Double_t *vecOut) const
{
  // fill vecOut by matrix*vecIn
  // vector should be of the same size as the matrix
  for (int i=GetSize();i--;) {
    vecOut[i] = 0.0;
    for (int j=GetSize();j--;) vecOut[i] += vecIn[j]*(*this)(i,j);
  }
  //
}

//___________________________________________________________
void AliMatrixSq::PrintCOO() const
{
  // print matrix in COO sparse format
  //
  // get number of non-zero elements
  int nnz = 0;
  int sz = GetSize();
  for (int ir=0;ir<sz;ir++) for (int ic=0;ic<sz;ic++) if (Query(ir,ic)!=0) nnz++;
  //
  printf("%d %d %d\n",sz,sz,nnz);
  double vl;
  for (int ir=0;ir<sz;ir++) for (int ic=0;ic<sz;ic++) if ((vl=Query(ir,ic))!=0) printf("%d %d %f\n",ir,ic,vl);
  //
}
