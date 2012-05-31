
/*********************************************************************************/
/* Symmetric Band Diagonal matrix with half band width W (+1 for diagonal)       */
/* Only lower triangle is stored in the "profile" format                         */ 
/*                                                                               */ 
/*                                                                               */
/* Author: ruben.shahoyan@cern.ch                                                */
/*                                                                               */ 
/*********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <float.h>
//
#include "TClass.h"
#include "TMath.h"
#include "AliSymBDMatrix.h"
//

using namespace std;

ClassImp(AliSymBDMatrix)


//___________________________________________________________
AliSymBDMatrix::AliSymBDMatrix() 
: fElems(0)
{
  // def. c-tor
  fSymmetric = kTRUE;
}

//___________________________________________________________
AliSymBDMatrix::AliSymBDMatrix(Int_t size, Int_t w)
  : AliMatrixSq(),fElems(0)
{
  // c-tor for given size
  //
  fNcols = size;    // number of rows
  if (w<0) w = 0;
  if (w>=size) w = size-1;
  fNrows = w;
  fRowLwb = w+1;
  fSymmetric = kTRUE;
  //
  // total number of stored elements
  fNelems = size*(w+1) - w*(w+1)/2;
  //
  fElems = new Double_t[fNcols*fRowLwb];
  memset(fElems,0,fNcols*fRowLwb*sizeof(Double_t));
  //
}

//___________________________________________________________
AliSymBDMatrix::AliSymBDMatrix(const AliSymBDMatrix &src) 
  : AliMatrixSq(src),fElems(0)
{
  // copy c-tor
  if (src.GetSize()<1) return;
  fNcols = src.GetSize();
  fNrows = src.GetBandHWidth();
  fRowLwb = fNrows+1;
  fNelems = src.GetNElemsStored();
  fElems = new Double_t[fNcols*fRowLwb];
  memcpy(fElems,src.fElems,fNcols*fRowLwb*sizeof(Double_t));
  //
}

//___________________________________________________________
AliSymBDMatrix::~AliSymBDMatrix() 
{
  // d-tor
  Clear();
}

//___________________________________________________________
AliSymBDMatrix&  AliSymBDMatrix::operator=(const AliSymBDMatrix& src)
{
  // assignment
  //
  if (this != &src) {
    TObject::operator=(src);
    if (fNcols!=src.fNcols) {
      // recreate the matrix
      if (fElems) delete[] fElems;
      fNcols = src.GetSize();
      fNrows = src.GetBandHWidth();
      fNelems = src.GetNElemsStored();
      fRowLwb = src.fRowLwb;
      fElems = new Double_t[fNcols*fRowLwb];
    }
    memcpy(fElems,src.fElems,fNcols*fRowLwb*sizeof(Double_t));     
    fSymmetric = kTRUE;
  }
  //
  return *this;
  //
}

//___________________________________________________________
void AliSymBDMatrix::Clear(Option_t*)
{
  // clear dynamic part
  if (fElems) {delete[] fElems; fElems = 0;}
  //  
  fNelems = fNcols = fNrows = fRowLwb = 0;
  //
}

//___________________________________________________________
Float_t AliSymBDMatrix::GetDensity() const
{
  // get fraction of non-zero elements
  if (!fNelems) return 0;
  Int_t nel = 0;
  for (int i=fNelems;i--;) if (!IsZero(fElems[i])) nel++;
  return nel/fNelems;
}

//___________________________________________________________
void AliSymBDMatrix::Print(Option_t* option) const
{
  // print data
  printf("Symmetric Band-Diagonal Matrix : Size = %d, half bandwidth = %d\n",
	 GetSize(),GetBandHWidth());
  TString opt = option; opt.ToLower();
  if (opt.IsNull()) return;
  opt = "%"; opt += 1+int(TMath::Log10(double(GetSize()))); opt+="d|";
  for (Int_t i=0;i<GetSize();i++) {
    printf(opt,i);
    for (Int_t j=TMath::Max(0,i-GetBandHWidth());j<=i;j++) printf("%+.3e|",GetEl(i,j));
    printf("\n");
  }
}

//___________________________________________________________
void AliSymBDMatrix::MultiplyByVec(const Double_t *vecIn,Double_t *vecOut) const
{
  // fill vecOut by matrix*vecIn
  // vector should be of the same size as the matrix
  if (IsDecomposed()) {
    for (int i=0;i<GetSize();i++) {
      double sm = 0;
      int jmax = TMath::Min(GetSize(),i+fRowLwb);
      for (int j=i+1;j<jmax;j++) sm += vecIn[j]*Query(j,i);
      vecOut[i] = QueryDiag(i)*(vecIn[i]+sm);
    }
    for (int i=GetSize();i--;) {
      double sm = 0;
      int jmin = TMath::Max(0,i - GetBandHWidth());
      int jmax = i-1;
      for (int j=jmin;j<jmax;j++) sm += vecOut[j]*Query(i,j);
      vecOut[i] += sm;
    }
  }
  else { // not decomposed
    for (int i=GetSize();i--;) {
      vecOut[i] = 0.0;
      int jmin = TMath::Max(0,i - GetBandHWidth());
      int jmax = TMath::Min(GetSize(),i + fRowLwb);
      for (int j=jmin;j<jmax;j++) vecOut[i] += vecIn[j]*Query(i,j);
    }
  }
  //
}

//___________________________________________________________
void AliSymBDMatrix::Reset()
{
  // set all elems to 0
  if (fElems) memset(fElems,0,fNcols*fRowLwb*sizeof(Double_t));
  SetDecomposed(kFALSE);
  //
}

//___________________________________________________________
void AliSymBDMatrix::AddToRow(Int_t r, Double_t *valc,Int_t *indc,Int_t n)
{
  // add list of elements to row r
  for (int i=0;i<n;i++) (*this)(r,indc[i]) = valc[i];
}

//___________________________________________________________
void AliSymBDMatrix::DecomposeLDLT()
{
  // decomposition to L Diag L^T
  if (IsDecomposed()) return;
  //
  Double_t eps = std::numeric_limits<double>::epsilon()*std::numeric_limits<double>::epsilon();
  //
  Double_t dtmp,gamma=0.0,xi=0.0;
  int iDiag;
  //
  // find max diag and number of non-0 diag.elements
  for (dtmp=0.0,iDiag=0;iDiag<GetSize();iDiag++) {
    if ( (dtmp=QueryDiag(iDiag)) <=0.0) break;
    if (gamma < dtmp) gamma = dtmp; 
  }
  //
  // find max. off-diag element
  for (int ir=1;ir<iDiag;ir++) {
    for (int ic=ir-GetBandHWidth();ic<ir;ic++) {
      if (ic<0) continue;
      dtmp = TMath::Abs(Query(ir,ic));
      if (xi<dtmp) xi = dtmp;
    }
  }
  double delta = eps*TMath::Max(1.0,xi+gamma);
  //
  double sn = GetSize()>1 ? 1.0/TMath::Sqrt( GetSize()*GetSize() - 1.0) : 1.0;
  double beta = TMath::Sqrt(TMath::Max(eps,TMath::Max(gamma,xi*sn))); 
  //
  for (int kr=1;kr<GetSize();kr++) {
    int colKmin = TMath::Max(0,kr - GetBandHWidth());
    double theta = 0.0;
    //
    for (int jr=colKmin;jr<=kr;jr++) {
      int colJmin = TMath::Max(0,jr - GetBandHWidth());
      //
      dtmp = 0.0;
      for (int i=TMath::Max(colKmin,colJmin);i<jr;i++) 
	dtmp += Query(kr,i)*QueryDiag(i)*Query(jr,i);
      dtmp = (*this)(kr,jr) -= dtmp;
      //
      theta = TMath::Max(theta, TMath::Abs(dtmp));
      //
      if (jr!=kr) {
	if (!IsZero(QueryDiag(jr))) (*this)(kr,jr) /= QueryDiag(jr);
	else                        (*this)(kr,jr) = 0.0;
      }
      else if (kr<iDiag) {
	dtmp = theta/beta;
	dtmp *= dtmp;
	dtmp = TMath::Max(dtmp, delta);
	(*this)(kr,jr) = TMath::Max( TMath::Abs(Query(kr,jr)), dtmp); 
      }
    } // jr
  } // kr
  //
  for (int i=0;i<GetSize();i++) {
    dtmp = QueryDiag(i);
    if (!IsZero(dtmp)) DiagElem(i) = 1./dtmp;
  }
  //
  SetDecomposed();
}

//___________________________________________________________
void AliSymBDMatrix::Solve(Double_t *rhs)
{
  // solve matrix equation
  //
  if (!IsDecomposed()) DecomposeLDLT();
  //
  for (int kr=0;kr<GetSize();kr++) 
    for (int jr=TMath::Max(0,kr-GetBandHWidth());jr<kr;jr++) rhs[kr] -= Query(kr,jr)*rhs[jr];
  //
  for (int kr=GetSize();kr--;) rhs[kr] *= QueryDiag(kr);
  //
  for (int kr=GetSize();kr--;)
    for (int jr=TMath::Max(0,kr - GetBandHWidth());jr<kr;jr++) rhs[jr] -= Query(kr,jr)*rhs[kr];
  //
}

//___________________________________________________________
void AliSymBDMatrix::Solve(const Double_t *rhs,Double_t *sol)
{
  // solve matrix equation
  memcpy(sol,rhs,GetSize()*sizeof(Double_t));
  Solve(sol);
  //
}
