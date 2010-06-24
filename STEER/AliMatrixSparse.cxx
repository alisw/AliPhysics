#include "AliMatrixSparse.h"
#include "TStopwatch.h"

/**********************************************************************************************/
/* Sparse matrix class, used as a global matrix for AliMillePede2                             */
/*                                                                                            */ 
/* Author: ruben.shahoyan@cern.ch                                                             */
/*                                                                                            */ 
/**********************************************************************************************/

//___________________________________________________________
ClassImp(AliMatrixSparse)

//___________________________________________________________
AliMatrixSparse::AliMatrixSparse(Int_t sz)
: AliMatrixSq(),fVecs(0)
{
  fNcols=fNrows=sz;
  //
  fVecs = new AliVectorSparse*[sz];
  for (int i=GetSize();i--;) fVecs[i] = new AliVectorSparse();
}

//___________________________________________________________
AliMatrixSparse::AliMatrixSparse(const AliMatrixSparse& src)
  : AliMatrixSq(src),fVecs(0)
{
  fVecs = new AliVectorSparse*[src.GetSize()];
  for (int i=GetSize();i--;) fVecs[i] = new AliVectorSparse( *src.GetRow(i));
}

//___________________________________________________________
AliVectorSparse* AliMatrixSparse::GetRowAdd(Int_t ir)
{
  if (ir>=fNrows) {
    AliVectorSparse** arrv = new AliVectorSparse*[ir+1];
    for (int i=GetSize();i--;) arrv[i] = fVecs[i];
    delete fVecs;
    fVecs = arrv;    
    for (int i=GetSize();i<=ir;i++) fVecs[i] = new AliVectorSparse();
    fNrows = ir+1;
    if (IsSymmetric() && fNcols<fNrows) fNcols = fNrows;
  }
  return fVecs[ir];
}

//___________________________________________________________
AliMatrixSparse& AliMatrixSparse::operator=(const AliMatrixSparse& src)
{
  if (*this == src) return *this;
  Clear();
  fNcols = src.GetNCols();
  fNrows = src.GetNRows();
  SetSymmetric(src.IsSymmetric());
  fVecs = new AliVectorSparse*[fNrows];
  for (int i=fNrows;i--;) fVecs[i] = new AliVectorSparse( *src.GetRow(i));
  return *this;
}

//___________________________________________________________
void AliMatrixSparse::Clear(Option_t*) 
{
  for (int i=fNrows;i--;) delete GetRow(i);
  delete[] fVecs;
  fNcols = fNrows = 0;
}

//___________________________________________________________
void AliMatrixSparse::Print(Option_t* opt)  const
{
  printf("Sparse Matrix of size %d x %d %s\n",fNrows,fNcols,IsSymmetric() ? " (Symmetric)":"");
  for (int i=0;i<fNrows;i++) {
    AliVectorSparse* row = GetRow(i);
    if (!row->GetNElems()) continue;
    printf("%3d: ",i); 
    row->Print(opt);
  }
}

//___________________________________________________________
void AliMatrixSparse::MultiplyByVec(const Double_t* vecIn, Double_t* vecOut) const
{
  // fill vecOut by matrix*vecIn
  // vector should be of the same size as the matrix
  //
  memset(vecOut,0,GetSize()*sizeof(Double_t));
  //
  for (int rw=GetSize();rw--;) {  // loop over rows >>>
    const AliVectorSparse* rowV = GetRow(rw);
    Int_t nel  = rowV->GetNElems();
    if (!nel) continue;
    //
    UShort_t* indV = rowV->GetIndices();
    Double_t* elmV = rowV->GetElems();
    //
    if (IsSymmetric()) {
      // treat diagonal term separately. If filled, it should be the last one
      if (indV[--nel]==rw) vecOut[rw] += vecIn[rw]*elmV[nel];
      else nel = rowV->GetNElems(); // diag elem was not filled
      //
      for (int iel=nel;iel--;) {          // less element retrieval for symmetric case
	if (elmV[iel]) {        
	  vecOut[rw]        += vecIn[indV[iel]]*elmV[iel];
	  vecOut[indV[iel]] += vecIn[rw]*elmV[iel];
	}
      }
    }
    else for (int iel=nel;iel--;) if (elmV[iel]) vecOut[rw] += vecIn[indV[iel]]*elmV[iel];
    //
  } // loop over rows <<<
  //
}

//___________________________________________________________
void AliMatrixSparse::SortIndices(Bool_t valuesToo)
{
  TStopwatch sw; 
  sw.Start();
  printf("AliMatrixSparse:sort>>\n");
  // sort columns in increasing order. Used to fix the matrix after ILUk decompostion
  for (int i=GetSize();i--;) GetRow(i)->SortIndices(valuesToo);
  sw.Stop();
  sw.Print();
  printf("AliMatrixSparse:sort<<\n");
}

//___________________________________________________________
void AliMatrixSparse::AddToRow(Int_t r, Double_t *valc,Int_t *indc,Int_t n)
{
  // for sym. matrix count how many elems to add have row>=col and assign excplicitly 
  // those which have row<col
  //   
  // range in increasing order of indices
  for (int i=n;i--;) {
    for (int j=i;j>=0;j--) {
      if (indc[j]>indc[i]) { // swap
	int    ti = indc[i]; indc[i] = indc[j]; indc[j] = ti;
	double tv = valc[i]; valc[i] = valc[j]; valc[j] = tv;
      }
    }
  }
  //
  int ni=n;
  if (IsSymmetric()) {
    while(ni--) {
      if (indc[ni]>r) (*this)(indc[ni],r) += valc[ni]; 
      else break;  // use the fact that the indices are ranged in increasing order
    }
  }
  //
  if (ni<0) return;
  AliVectorSparse* row = GetRowAdd(r);
  row->Add(valc,indc,ni+1);
}

//___________________________________________________________
Float_t AliMatrixSparse::GetDensity() const
{
  // get fraction of non-zero elements
  Int_t nel = 0;
  for (int i=GetSize();i--;) nel += GetRow(i)->GetNElems();
  int den = IsSymmetric() ? (GetSize()+1)*GetSize()/2 : GetSize()*GetSize();
  return float(nel)/den;
}

