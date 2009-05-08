#include <stdlib.h>
#include <stdio.h>
#include <iostream>
//
#include "TClass.h"
#include "TMath.h"
#include "AliSymMatrix.h"
//

using namespace std;

ClassImp(AliSymMatrix)


AliSymMatrix* AliSymMatrix::fgBuffer = 0; 
Int_t         AliSymMatrix::fgCopyCnt = 0; 
//___________________________________________________________
AliSymMatrix::AliSymMatrix() 
: fElems(0),fElemsAdd(0)
{
  fSymmetric = kTRUE;
  fgCopyCnt++;
}

//___________________________________________________________
AliSymMatrix::AliSymMatrix(Int_t size)
  : AliMatrixSq(),fElems(0),fElemsAdd(0)
{
  //
  fNrows = 0;
  fNrowIndex = fNcols = size;
  fElems     = new Double_t[fNcols*(fNcols+1)/2];
  fSymmetric = kTRUE;
  Reset();
  fgCopyCnt++;
  //
}

//___________________________________________________________
AliSymMatrix::AliSymMatrix(const AliSymMatrix &src) 
  : AliMatrixSq(src),fElems(0),fElemsAdd(0)
{
  fNrowIndex = fNcols = src.GetSize();
  fNrows = 0;
  if (fNcols) {
    int nmainel = fNcols*(fNcols+1)/2;
    fElems     = new Double_t[nmainel];
    nmainel = src.fNcols*(src.fNcols+1)/2;
    memcpy(fElems,src.fElems,nmainel*sizeof(Double_t));
    if (src.fNrows) { // transfer extra rows to main matrix
      Double_t *pnt = fElems + nmainel;
      int ncl = src.fNcols + 1;
      for (int ir=0;ir<src.fNrows;ir++) {
	memcpy(pnt,src.fElemsAdd[ir],ncl*sizeof(Double_t));
	pnt += ncl;
	ncl++; 
      }
    }
  }
  else fElems = 0;
  fElemsAdd = 0;
  fgCopyCnt++;
  //
}

//___________________________________________________________
AliSymMatrix::~AliSymMatrix() 
{
  Clear();
  if (--fgCopyCnt < 1 && fgBuffer) {delete fgBuffer; fgBuffer = 0;}
}

//___________________________________________________________
AliSymMatrix&  AliSymMatrix::operator=(const AliSymMatrix& src)
{
  //
  if (this != &src) {
    TObject::operator=(src);
    if (fNcols!=src.fNcols && fNrows!=src.fNrows) {
      // recreate the matrix
      if (fElems) delete[] fElems;
      for (int i=0;i<fNrows;i++) delete[] fElemsAdd[i]; 
      delete[] fElemsAdd;
      //
      fNrowIndex = src.GetSize(); 
      fNcols = src.GetSize();
      fNrows = 0;
      fElems     = new Double_t[fNcols*(fNcols+1)/2];
      int nmainel = src.fNcols*(src.fNcols+1);
      memcpy(fElems,src.fElems,nmainel*sizeof(Double_t));
      if (src.fNrows) { // transfer extra rows to main matrix
	Double_t *pnt = fElems + nmainel*sizeof(Double_t);
	int ncl = src.fNcols + 1;
	for (int ir=0;ir<src.fNrows;ir++) {
	  ncl += ir; 
	  memcpy(pnt,src.fElemsAdd[ir],ncl*sizeof(Double_t));
	  pnt += ncl*sizeof(Double_t);
	}
      }
      //
    }
    else {
      memcpy(fElems,src.fElems,fNcols*(fNcols+1)/2*sizeof(Double_t));
      int ncl = fNcols + 1;
      for (int ir=0;ir<fNrows;ir++) { // dynamic rows
	ncl += ir; 
	memcpy(fElemsAdd[ir],src.fElemsAdd[ir],ncl*sizeof(Double_t));
      }
    }
  }
  //
  return *this;
}

//___________________________________________________________
void AliSymMatrix::Clear(Option_t*)
{
  if (fElems) {delete[] fElems; fElems = 0;}
  //  
  if (fElemsAdd) {
    for (int i=0;i<fNrows;i++) delete[] fElemsAdd[i]; 
    delete[] fElemsAdd;
    fElemsAdd = 0;
  }
  fNrowIndex = 0;
  fNcols = 0;
  fNrows = 0;
  //
}

//___________________________________________________________
Float_t AliSymMatrix::GetDensity() const
{
  // get fraction of non-zero elements
  Int_t nel = 0;
  for (int i=GetSize();i--;) for (int j=i+1;j--;) if (GetEl(i,j)!=0) nel++;
  return 2.*nel/( (GetSize()+1)*GetSize() );
}

//___________________________________________________________
void AliSymMatrix::Print(Option_t* option) const
{
  printf("Symmetric Matrix: Size = %d (%d rows added dynamically)\n",GetSize(),fNrows);
  TString opt = option; opt.ToLower();
  if (opt.IsNull()) return;
  opt = "%"; opt += 1+int(TMath::Log10(double(GetSize()))); opt+="d|";
  for (Int_t i=0;i<fNrowIndex;i++) {
    printf(opt,i);
    for (Int_t j=0;j<=i;j++) printf("%+.3e|",GetEl(i,j));
    printf("\n");
  }
}

//___________________________________________________________
void AliSymMatrix::MultiplyByVec(Double_t *vecIn,Double_t *vecOut) const
{
  // fill vecOut by matrix*vecIn
  // vector should be of the same size as the matrix
  for (int i=fNrowIndex;i--;) {
    vecOut[i] = 0.0;
    for (int j=fNrowIndex;j--;) vecOut[i] += vecIn[j]*GetEl(i,j);
  }
  //
}

//___________________________________________________________
AliSymMatrix* AliSymMatrix::DecomposeChol() 
{
  // Return a matrix with Choleski decomposition
  // Adopted from Numerical Recipes in C, ch.2-9, http://www.nr.com
  // consturcts Cholesky decomposition of SYMMETRIC and
  // POSITIVELY-DEFINED matrix a (a=L*Lt)
  // Only upper triangle of the matrix has to be filled.
  // In opposite to function from the book, the matrix is modified:
  // lower triangle and diagonal are refilled.
  //
  if (!fgBuffer || fgBuffer->GetSize()!=GetSize()) {
    delete fgBuffer; 
    try {
      fgBuffer = new AliSymMatrix(*this);
    }
    catch(bad_alloc&) {
      printf("Failed to allocate memory for Choleski decompostions\n");
      return 0;
    }
  }
  else (*fgBuffer) = *this;
  //
  AliSymMatrix& mchol = *fgBuffer;
  //
  for (int i=0;i<fNrowIndex;i++) {
    Double_t *rowi = mchol.GetRow(i);
    for (int j=i;j<fNrowIndex;j++) {
      Double_t *rowj = mchol.GetRow(j);
      double sum = rowj[i];
      for (int k=i-1;k>=0;k--) if (rowi[k]&&rowj[k]) sum -= rowi[k]*rowj[k];
      if (i == j) {
	if (sum <= 0.0) { // not positive-definite
	  printf("The matrix is not positive definite [%e]\n"
		 "Choleski decomposition is not possible\n",sum);
	  return 0;
	}
	rowi[i] = TMath::Sqrt(sum);
	//
      } else rowj[i] = sum/rowi[i];
    }
  }
  return fgBuffer;
}

//___________________________________________________________
Bool_t AliSymMatrix::InvertChol() 
{
  // Invert matrix using Choleski decomposition
  //
  AliSymMatrix* mchol = DecomposeChol();
  if (!mchol) {
    printf("Failed to invert the matrix\n");
    return kFALSE;
  }
  //
  InvertChol(mchol);
  return kTRUE;
  //
}

//___________________________________________________________
void AliSymMatrix::InvertChol(AliSymMatrix* pmchol) 
{
  // Invert matrix using Choleski decomposition, provided the Cholseki's L matrix
  //
  Double_t sum;
  AliSymMatrix& mchol = *pmchol;
  //
  // Invert decomposed triangular L matrix (Lower triangle is filled)
  for (int i=0;i<fNrowIndex;i++) { 
    mchol(i,i) =  1.0/mchol(i,i);
    for (int j=i+1;j<fNrowIndex;j++) { 
      Double_t *rowj = mchol.GetRow(j);
      sum = 0.0; 
      for (int k=i;k<j;k++) if (rowj[k]) { 
	double &mki = mchol(k,i); if (mki) sum -= rowj[k]*mki;
      }
      rowj[i] = sum/rowj[j];
    } 
  }
  //
  // take product of the inverted Choleski L matrix with its transposed
  for (int i=fNrowIndex;i--;) {
    for (int j=i+1;j--;) {
      sum = 0;
      for (int k=i;k<fNrowIndex;k++) {
	double &mik = mchol(i,k); 
	if (mik) {
	  double &mjk = mchol(j,k);
	  if (mjk) sum += mik*mjk;
	}
      }
      (*this)(j,i) = sum;
    }
  }
  //
}


//___________________________________________________________
Bool_t AliSymMatrix::SolveChol(Double_t *b, Bool_t invert) 
{
  // Adopted from Numerical Recipes in C, ch.2-9, http://www.nr.com
  // Solves the set of n linear equations A x = b, 
  // where a is a positive-definite symmetric matrix. 
  // a[1..n][1..n] is the output of the routine CholDecomposw. 
  // Only the lower triangle of a is accessed. b[1..n] is input as the 
  // right-hand side vector. The solution vector is returned in b[1..n].
  //
  Int_t i,k;
  Double_t sum;
  //
  AliSymMatrix *pmchol = DecomposeChol();
  if (!pmchol) {
    printf("SolveChol failed\n");
    //    Print("l");
    return kFALSE;
  }
  AliSymMatrix& mchol = *pmchol;
  //
  for (i=0;i<fNrowIndex;i++) {
    Double_t *rowi = mchol.GetRow(i);
    for (sum=b[i],k=i-1;k>=0;k--) if (rowi[k]&&b[k]) sum -= rowi[k]*b[k];
    b[i]=sum/rowi[i];
  }
  //
  for (i=fNrowIndex-1;i>=0;i--) {
    for (sum=b[i],k=i+1;k<fNrowIndex;k++) if (b[k]) {
      double &mki=mchol(k,i); if (mki) sum -= mki*b[k];
    }
    b[i]=sum/mchol(i,i);
  }
  //
  if (invert) InvertChol(pmchol);
  return kTRUE;
  //
}

//___________________________________________________________
Bool_t AliSymMatrix::SolveChol(TVectorD &b, Bool_t invert) 
{
  return SolveChol((Double_t*)b.GetMatrixArray(),invert);
}


//___________________________________________________________
Bool_t AliSymMatrix::SolveChol(Double_t *brhs, Double_t *bsol,Bool_t invert) 
{
  memcpy(bsol,brhs,GetSize()*sizeof(Double_t));
  return SolveChol(bsol,invert);  
}

//___________________________________________________________
Bool_t AliSymMatrix::SolveChol(TVectorD &brhs, TVectorD &bsol,Bool_t invert) 
{
  bsol = brhs;
  return SolveChol(bsol,invert);
}

//___________________________________________________________
void AliSymMatrix::AddRows(int nrows)
{
  if (nrows<1) return;
  Double_t **pnew = new Double_t*[nrows+fNrows];
  for (int ir=0;ir<fNrows;ir++) pnew[ir] = fElemsAdd[ir]; // copy old extra rows
  for (int ir=0;ir<nrows;ir++) {
    int ncl = GetSize()+1;
    pnew[fNrows] = new Double_t[ncl];
    memset(pnew[fNrows],0,ncl*sizeof(Double_t));
    fNrows++;
    fNrowIndex++;
  }
  delete[] fElemsAdd;
  fElemsAdd = pnew;
  //
}

//___________________________________________________________
void AliSymMatrix::Reset()
{
  // if additional rows exist, regularize it
  if (fElemsAdd) {
    delete[] fElems;
    for (int i=0;i<fNrows;i++) delete[] fElemsAdd[i]; 
    delete[] fElemsAdd; fElemsAdd = 0;
    fNcols = fNrowIndex;
    fElems = new Double_t[fNcols*(fNcols+1)/2];
    fNrows = 0;
  }
  if (fElems) memset(fElems,0,fNcols*(fNcols+1)/2*sizeof(Double_t));
  //
}

//___________________________________________________________
/*
void AliSymMatrix::AddToRow(Int_t r, Double_t *valc,Int_t *indc,Int_t n)
{
  //   for (int i=n;i--;) {
  //     (*this)(indc[i],r) += valc[i];
  //   }
  //   return;

  double *row;
  if (r>=fNrowIndex) {
    AddRows(r-fNrowIndex+1); 
    row = &((fElemsAdd[r-fNcols])[0]);
  }
  else row = &fElems[GetIndex(r,0)];
  //
  int nadd = 0;
  for (int i=n;i--;) {
    if (indc[i]>r) continue;
    row[indc[i]] += valc[i];
    nadd++;
  }
  if (nadd == n) return;
  //
  // add to col>row
  for (int i=n;i--;) {
    if (indc[i]>r) (*this)(indc[i],r) += valc[i];
  }
  //
}
*/

//___________________________________________________________
Double_t* AliSymMatrix::GetRow(Int_t r)
{
  if (r>=fNrowIndex) {
    int nn = fNrowIndex;
    AddRows(r-fNrowIndex+1); 
    printf("create %d of %d\n",r, nn);
    return &((fElemsAdd[r-fNcols])[0]);
  }
  else return &fElems[GetIndex(r,0)];
}


//___________________________________________________________
int AliSymMatrix::SolveSpmInv(double *vecB, Bool_t stabilize)
{
  //   Solution a la MP1: gaussian eliminations
  ///  Obtain solution of a system of linear equations with symmetric matrix 
  ///  and the inverse (using 'singular-value friendly' GAUSS pivot)
  //

  Int_t nRank = 0;
  int iPivot;
  double vPivot = 0.;
  double eps = 0.00000000000001;
  int nGlo = GetSize();
  bool   *bUnUsed = new bool[nGlo];
  double *rowMax,*colMax=0;
  rowMax  = new double[nGlo];
  //
  if (stabilize) {
    colMax   = new double[nGlo];
    for (Int_t i=nGlo; i--;) rowMax[i] = colMax[i] = 0.0;
    for (Int_t i=nGlo; i--;) for (Int_t j=i+1;j--;) { 
	double vl = TMath::Abs(Query(i,j));
	if (vl==0) continue;
	if (vl > rowMax[i]) rowMax[i] = vl; // Max elemt of row i
	if (vl > colMax[j]) colMax[j] = vl; // Max elemt of column j
	if (i==j) continue;
	if (vl > rowMax[j]) rowMax[j] = vl; // Max elemt of row j
	if (vl > colMax[i]) colMax[i] = vl; // Max elemt of column i
      }
    //
    for (Int_t i=nGlo; i--;) {
      if (0.0 != rowMax[i]) rowMax[i] = 1./rowMax[i]; // Max elemt of row i
      if (0.0 != colMax[i]) colMax[i] = 1./colMax[i]; // Max elemt of column i
    }
    //
  }
  //
  for (Int_t i=nGlo; i--;) bUnUsed[i] = true;
  //  
  if (!fgBuffer || fgBuffer->GetSize()!=GetSize()) {
    delete fgBuffer; 
    try {
      fgBuffer = new AliSymMatrix(*this);
    }
    catch(bad_alloc&) {
      printf("Failed to allocate memory for matrix inversion buffer\n");
      return 0;
    }
  }
  else (*fgBuffer) = *this;
  //
  if (stabilize) for (int i=0;i<nGlo; i++) { // Small loop for matrix equilibration (gives a better conditioning) 
      for (int j=0;j<=i; j++) {
	double vl = Query(i,j);
	if (vl!=0) SetEl(i,j, TMath::Sqrt(rowMax[i])*vl*TMath::Sqrt(colMax[j]) ); // Equilibrate the V matrix
      }
      for (int j=i+1;j<nGlo;j++) {
	double vl = Query(j,i);
	if (vl!=0) fgBuffer->SetEl(j,i,TMath::Sqrt(rowMax[i])*vl*TMath::Sqrt(colMax[j]) ); // Equilibrate the V matrix
      }
    }
  //
  for (Int_t j=nGlo; j--;) fgBuffer->DiagElem(j) = TMath::Abs(QueryDiag(j)); // save diagonal elem absolute values 
  //
  for (Int_t i=0; i<nGlo; i++) {
    vPivot = 0.0;
    iPivot = -1;
    //
    for (Int_t j=0; j<nGlo; j++) { // First look for the pivot, ie max unused diagonal element       
      double vl;
      if (bUnUsed[j] && (TMath::Abs(vl=QueryDiag(j))>TMath::Max(TMath::Abs(vPivot),eps*fgBuffer->QueryDiag(j)))) {    
	vPivot = vl;
	iPivot = j;
      }
    }
    //
    if (iPivot >= 0) {   // pivot found          
      nRank++;
      bUnUsed[iPivot] = false; // This value is used
      vPivot = 1.0/vPivot;
      DiagElem(iPivot) = -vPivot; // Replace pivot by its inverse
      //
      for (Int_t j=0; j<nGlo; j++) {      
	for (Int_t jj=0; jj<nGlo; jj++) {  
	  if (j != iPivot && jj != iPivot) {// Other elements (!!! do them first as you use old matV[k][j]'s !!!)	  
	    double &r = j>=jj ? (*this)(j,jj) : (*fgBuffer)(jj,j);
	    r -= vPivot* ( j>iPivot  ? Query(j,iPivot)  : fgBuffer->Query(iPivot,j) )
	      *          ( iPivot>jj ? Query(iPivot,jj) : fgBuffer->Query(jj,iPivot));
	  }
	}
      }
      //
      for (Int_t j=0; j<nGlo; j++) if (j != iPivot) { // Pivot row or column elements 
	  (*this)(j,iPivot)     *= vPivot;
	  (*fgBuffer)(iPivot,j) *= vPivot;
	}
      //
    }
    else {  // No more pivot value (clear those elements)
      for (Int_t j=0; j<nGlo; j++) {
	if (bUnUsed[j]) {
	  vecB[j] = 0.0;
	  for (Int_t k=0; k<nGlo; k++) {
	    (*this)(j,k) = 0.;
	    if (j!=k) (*fgBuffer)(j,k) = 0;
	  }
	}
      }
      break;  // No more pivots anyway, stop here
    }
  }
  //
  for (Int_t i=0; i<nGlo; i++) for (Int_t j=0; j<nGlo; j++) {
      double vl = TMath::Sqrt(colMax[i])*TMath::Sqrt(rowMax[j]); // Correct matrix V
      if (i>=j) (*this)(i,j) *= vl;
      else      (*fgBuffer)(j,i) *= vl;
    }
  //
  for (Int_t j=0; j<nGlo; j++) {
    rowMax[j] = 0.0;
    for (Int_t jj=0; jj<nGlo; jj++) { // Reverse matrix elements
      double vl;
      if (j>=jj) vl = (*this)(j,jj)     = -Query(j,jj);
      else       vl = (*fgBuffer)(j,jj) = -fgBuffer->Query(j,jj);
      rowMax[j] += vl*vecB[jj];
    }		
  }
  
  for (Int_t j=0; j<nGlo; j++) {
    vecB[j] = rowMax[j]; // The final result
  }
  //
  delete [] bUnUsed;
  delete [] rowMax;
  if (stabilize) delete [] colMax;

  return nRank;
}


