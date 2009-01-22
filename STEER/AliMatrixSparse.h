#ifndef ALIMATRIXSPARSE_H
#define ALIMATRIXSPARSE_H

#include "AliMatrixSq.h"


///////////////////////////////////////////////////////////////////////////////////////
class AliVectorSparse {
 public:
  AliVectorSparse();
  AliVectorSparse(const AliVectorSparse& src);
  virtual ~AliVectorSparse() {Clear();}
  virtual void Print(Option_t* option="")                 const;
  //
  Int_t     GetNElems()                                   const {return fNElems;}
  UShort_t *GetIndices()                                  const {return fIndex;}
  Double_t *GetElems()                                    const {return fElems;}
  UShort_t& GetIndex(Int_t i)                                   {return fIndex[i];}
  Double_t& GetElem(Int_t i)                              const {return fElems[i];}
  void      Clear();
  void      Reset()                                             {memset(fElems,0,fNElems*sizeof(Double_t));}
  void      ReSize(Int_t sz,Bool_t copy=kFALSE);
  void      SortIndices(Bool_t valuesToo=kFALSE);
  //
  AliVectorSparse& operator=(const AliVectorSparse& src);
  //
  virtual Double_t         operator()(Int_t ind)         const;
  virtual Double_t&        operator()(Int_t ind);
  virtual void             Zero(Int_t ind);
  Double_t                 FindIndex(Int_t ind)          const;
  Double_t&                FindIndexAdd(Int_t ind);
  //
  Int_t     GetLastIndex()                               const {return fIndex[fNElems-1];}
  Double_t  GetLastElem()                                const {return fElems[fNElems-1];}
  Double_t &GetLastElem()                                      {return fElems[fNElems-1];}
  //
 protected:
  Int_t            fNElems;   // 
  UShort_t*        fIndex;    // Index of stored elems
  Double_t*        fElems;    // pointer on elements
};


//___________________________________________________
inline Double_t AliVectorSparse::operator()(Int_t ind) const
{
  return FindIndex(ind);
}

//___________________________________________________
inline Double_t& AliVectorSparse::operator()(Int_t ind)
{
  return FindIndexAdd(ind);
}

//////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////
// Sparse matrix class
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

  void Clear(Option_t* option="");
  void Reset()                            {for (int i=fNcols;i--;) GetRow(i)->Reset();}
  void Print(Option_t* option="")                      const;
  AliMatrixSparse& operator=(const AliMatrixSparse& src);
  Double_t&        operator()(Int_t row,Int_t col);
  Double_t         operator()(Int_t row,Int_t col)     const;
  void             Zero(Int_t row,Int_t col);
  Float_t          GetDensity()                        const;
  //
  Double_t         DiagElem(Int_t r)                   const;
  Double_t&        DiagElem(Int_t r);
  void             SortIndices(Bool_t valuesToo=kFALSE);
  //
  void MultiplyByVec(Double_t* vecIn, Double_t* vecOut) const;
  void MultiplyByVec(TVectorD &vecIn, TVectorD &vecOut) const {MultiplyByVec((Double_t*)vecIn.GetMatrixArray(),(Double_t*)vecOut.GetMatrixArray());}
  //
 protected:
  //
  AliVectorSparse** fVecs;
  //
  ClassDef(AliMatrixSparse,0)
};

//___________________________________________________
inline void AliMatrixSparse::Zero(Int_t row,Int_t col)
{
  //  set existing element to 0
  if (IsSymmetric() && col>row) Swap(row,col); 
  AliVectorSparse* rowv = GetRow(row);
  if (rowv) rowv->Zero(col);
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
  if (IsSymmetric()) return (rowv->GetNElems()>0 && rowv->GetLastIndex()==row) ? 
		       rowv->GetLastElem() : rowv->FindIndexAdd(row);
  else return rowv->FindIndexAdd(row);
  //
}

#endif

