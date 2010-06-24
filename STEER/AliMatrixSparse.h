#ifndef ALIMATRIXSPARSE_H
#define ALIMATRIXSPARSE_H

#include "AliMatrixSq.h"
#include "AliVectorSparse.h"


/**********************************************************************************************/
/* Sparse matrix class, used as a global matrix for AliMillePede2                             */
/*                                                                                            */ 
/* Author: ruben.shahoyan@cern.ch                                                             */
/*                                                                                            */ 
/**********************************************************************************************/


//
class AliMatrixSparse : public AliMatrixSq 
{
 public:
  AliMatrixSparse() : fVecs(0) {}
  AliMatrixSparse(Int_t size);
  AliMatrixSparse(const AliMatrixSparse &mat);
  virtual ~AliMatrixSparse() {Clear();}
  //
  AliVectorSparse* GetRow(Int_t ir)        const {return (ir<fNcols) ? fVecs[ir] : 0;}
  AliVectorSparse* GetRowAdd(Int_t ir);
  //
  virtual Int_t   GetSize()                            const {return fNrows;}
  virtual Int_t   GetNRows()                           const {return fNrows;}
  virtual Int_t   GetNCols()                           const {return fNcols;}
  //
  void Clear(Option_t* option="");
  void Reset()                            {for (int i=fNcols;i--;) GetRow(i)->Reset();}
  void Print(Option_t* option="")                      const;
  AliMatrixSparse& operator=(const AliMatrixSparse& src);
  Double_t&        operator()(Int_t row,Int_t col);
  Double_t         operator()(Int_t row,Int_t col)     const;
  void             SetToZero(Int_t row,Int_t col);
  Float_t          GetDensity()                        const;
  //
  Double_t         DiagElem(Int_t r)                   const;
  Double_t&        DiagElem(Int_t r);
  void             SortIndices(Bool_t valuesToo=kFALSE);
  //
  void MultiplyByVec(const TVectorD &vecIn, TVectorD &vecOut) const; 
  void MultiplyByVec(const Double_t* vecIn, Double_t* vecOut) const;
  //
  void AddToRow(Int_t r, Double_t *valc,Int_t *indc,Int_t n);
  //
 protected:
  //
  AliVectorSparse** fVecs;
  //
  ClassDef(AliMatrixSparse,0)
};

//___________________________________________________
inline void AliMatrixSparse::MultiplyByVec(const TVectorD &vecIn, TVectorD &vecOut) const 
{
  MultiplyByVec((Double_t*)vecIn.GetMatrixArray(),(Double_t*)vecOut.GetMatrixArray());
}

//___________________________________________________
inline void AliMatrixSparse::SetToZero(Int_t row,Int_t col)
{
  //  set existing element to 0
  if (IsSymmetric() && col>row) Swap(row,col); 
  AliVectorSparse* rowv = GetRow(row);
  if (rowv) rowv->SetToZero(col);
}

//___________________________________________________
inline Double_t AliMatrixSparse::operator()(Int_t row,Int_t col) const
{
  //  printf("M: find\n");
  if (IsSymmetric() && col>row) Swap(row,col); 
  AliVectorSparse* rowv = GetRow(row);
  if (!rowv) return 0;
  return rowv->FindIndex(col);
}

//___________________________________________________
inline Double_t& AliMatrixSparse::operator()(Int_t row,Int_t col)
{
  //  printf("M: findindexAdd\n");
  if (IsSymmetric() && col>row) Swap(row,col); 
  AliVectorSparse* rowv = GetRowAdd(row);
  if (col>=fNcols) fNcols = col+1;
  return rowv->FindIndexAdd(col);
}

//___________________________________________________
inline Double_t AliMatrixSparse::DiagElem(Int_t row) const
{
  AliVectorSparse* rowv = GetRow(row);
  if (!rowv) return 0;
  if (IsSymmetric()) return (rowv->GetNElems()>0 && rowv->GetLastIndex()==row) ? rowv->GetLastElem() : 0.;
  else return rowv->FindIndex(row);
  //
}

//___________________________________________________
inline Double_t &AliMatrixSparse::DiagElem(Int_t row)
{
  AliVectorSparse* rowv = GetRowAdd(row);
  if (row>=fNcols) fNcols = row+1;
  if (IsSymmetric()) {
    return (rowv->GetNElems()>0 && rowv->GetLastIndex()==row) ? 
		       rowv->GetLastElem() : rowv->FindIndexAdd(row);
  }
  else return rowv->FindIndexAdd(row);
  //
}


#endif

