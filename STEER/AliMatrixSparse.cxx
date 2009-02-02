#include "AliMatrixSparse.h"

//___________________________________________________________
AliVectorSparse::AliVectorSparse()
  : fNElems(0),fIndex(0),fElems(0) {}

//___________________________________________________________
AliVectorSparse::AliVectorSparse(const AliVectorSparse& src)
  : fNElems(src.fNElems),fIndex(0),fElems(0)
{
  fIndex = new UShort_t[fNElems];
  fElems = new Double_t[fNElems];
  memcpy(fIndex,src.fIndex,fNElems*sizeof(UShort_t));
  memcpy(fElems,src.fElems,fNElems*sizeof(Double_t));
}

//___________________________________________________________
void AliVectorSparse::Clear()
{
  delete[] fIndex; fIndex = 0;
  delete[] fElems; fElems = 0;
  fNElems = 0;
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

//___________________________________________________________
AliVectorSparse& AliVectorSparse::operator=(const AliVectorSparse& src)
{
  if (&src==this) return *this;
  Clear();
  fNElems = src.fNElems;
  fIndex = new UShort_t[fNElems];
  fElems = new Double_t[fNElems];
  memcpy(fIndex,src.fIndex,fNElems*sizeof(UShort_t));
  memcpy(fElems,src.fElems,fNElems*sizeof(Double_t));
  //
  return *this;
}

//___________________________________________________________
Double_t AliVectorSparse::FindIndex(Int_t ind) const
{
  // return an element with given index
  //printf("V: findindex\n");
  int first = 0;
  int last = fNElems-1;
  while (first<=last) {
    int mid = (first+last)>>1;
    if (ind>fIndex[mid]) first = mid+1;
    else if (ind<fIndex[mid]) last = mid-1;
    else return fElems[mid];
  }
  return 0.0;
}

//___________________________________________________________
void AliVectorSparse::SetToZero(Int_t ind)
{
  // set element to 0 if it was already defined
  int first = 0;
  int last = fNElems-1;
  while (first<=last) {
    int mid = (first+last)>>1;
    if (ind>fIndex[mid]) first = mid+1;
    else if (ind<fIndex[mid]) last = mid-1;
    else {fElems[mid] = 0.; return;}
  }
}

//___________________________________________________________
Double_t& AliVectorSparse::FindIndexAdd(Int_t ind)
{
  // increment an element with given index
  //printf("V: findindexAdd\n");
  int first = 0;
  int last = fNElems-1;
  while (first<=last) {
    int mid = (first+last)>>1;
    if (ind>fIndex[mid]) first = mid+1;
    else if (ind<fIndex[mid]) last = mid-1;
    else return fElems[mid];
  }
  // need to insert a new element
  UShort_t *arrI = new UShort_t[fNElems+1];
  memcpy(arrI,fIndex,first*sizeof(UShort_t));
  arrI[first] = ind;
  memcpy(arrI+first+1,fIndex+first,(fNElems-first)*sizeof(UShort_t));
  delete[] fIndex;
  fIndex = arrI;
  //
  Double_t   *arrE = new Double_t[fNElems+1];
  memcpy(arrE,fElems,first*sizeof(Double_t));
  arrE[first] = 0;
  memcpy(arrE+first+1,fElems+first,(fNElems-first)*sizeof(Double_t));
  delete[] fElems;
  fElems = arrE;
  //
  fNElems++;
  return fElems[first];
  //
}

//__________________________________________________________
void AliVectorSparse::ReSize(Int_t sz,Bool_t copy)
{
  if (sz<1) {Clear(); return;}
    // need to insert a new element
  UShort_t *arrI = new UShort_t[sz];
  Double_t *arrE = new Double_t[sz];
  memset(arrI,0,sz*sizeof(UShort_t));
  memset(arrE,0,sz*sizeof(Double_t));
  //
  if (copy && fIndex) {
    int cpsz = TMath::Min(fNElems,sz);
    memcpy(arrI,fIndex,cpsz*sizeof(UShort_t));
    memcpy(arrE,fElems,cpsz*sizeof(Double_t));
  }
  delete[] fIndex;
  delete[] fElems;
  fIndex = arrI;
  fElems = arrE;
  fNElems = sz;
  //
}

//__________________________________________________________
void AliVectorSparse::SortIndices(Bool_t valuesToo)
{
  // sort indices in increasing order. Used to fix the row after ILUk decomposition
  for (int i=fNElems;i--;) for (int j=i;j--;) if (fIndex[i]<fIndex[j]) { //swap
	UShort_t tmpI = fIndex[i]; fIndex[i] = fIndex[j]; fIndex[j]=tmpI;
	if (valuesToo) {Double_t tmpV = fElems[i];fElems[i]=fElems[j];fElems[j]=tmpV;}
      }
}

//__________________________________________________________
void AliVectorSparse::Print(Option_t* )  const
{
  printf("|");
  for (int i=0;i<fNElems;i++) printf("%2d:%+.2e|",fIndex[i],fElems[i]);
  printf("|\n");
}


///////////////////////////////////////////////////////////////////////////////////////////
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
  fVecs = new AliVectorSparse*[fNcols];
  for (int i=GetSize();i--;) fVecs[i] = new AliVectorSparse( *src.GetRow(i));
}

//___________________________________________________________
AliVectorSparse* AliMatrixSparse::GetRowAdd(Int_t ir)
{
  if (ir>=fNcols) {
    AliVectorSparse** arrv = new AliVectorSparse*[ir+1];
    for (int i=GetSize();i--;) arrv[i] = fVecs[i];
    delete fVecs;
    fVecs = arrv;    
    for (int i=GetSize();i<=ir;i++) fVecs[i] = new AliVectorSparse();
    fNcols = fNrows = ir+1;
  }
  return fVecs[ir];
}

//___________________________________________________________
AliMatrixSparse& AliMatrixSparse::operator=(const AliMatrixSparse& src)
{
  if (*this == src) return *this;
  Clear();
  fNcols = fNrows = src.GetSize();
  SetSymmetric(src.IsSymmetric());
  fVecs = new AliVectorSparse*[fNcols];
  for (int i=fNcols;i--;) fVecs[i] = new AliVectorSparse( *src.GetRow(i));
  return *this;
}

//___________________________________________________________
void AliMatrixSparse::Clear(Option_t*) 
{
  for (int i=fNcols;i--;) delete GetRow(i);
  delete[] fVecs;
  fNcols = fNrows = 0;
}

//___________________________________________________________
void AliMatrixSparse::Print(Option_t*)  const
{
  printf("Sparse Matrix of size %d %s\n",fNcols,IsSymmetric() ? " (Symmetric)":"");
  for (int i=0;i<fNcols;i++) {
    AliVectorSparse* row = GetRow(i);
    if (!row->GetNElems()) continue;
    printf("%3d: ",i); 
    row->Print();
  }
}

//___________________________________________________________
void AliMatrixSparse::MultiplyByVec(Double_t* vecIn, Double_t* vecOut) const
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
  // sort columns in increasing order. Used to fix the matrix after ILUk decompostion
  for (int i=GetSize();i--;) GetRow(i)->SortIndices(valuesToo);
}
