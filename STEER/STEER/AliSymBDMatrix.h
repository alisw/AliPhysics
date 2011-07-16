#ifndef ALISYMBDMATRIX_H
#define ALISYMBDMATRIX_H

/*********************************************************************************/
/* Symmetric Band Diagonal matrix with half band width W (+1 for diagonal)       */
/* Only lower triangle is stored in the "profile" format                         */ 
/*                                                                               */ 
/*                                                                               */
/* Author: ruben.shahoyan@cern.ch                                                */
/*                                                                               */ 
/*********************************************************************************/

//#include <string.h>
#include <TObject.h>
#include <TVectorD.h>
#include "AliMatrixSq.h"


class AliSymBDMatrix : public AliMatrixSq {
  //
 public:
  enum {kDecomposedBit = 0x1};
  //
  AliSymBDMatrix();
  AliSymBDMatrix(Int_t size, Int_t w = 0);
  AliSymBDMatrix(const AliSymBDMatrix &mat);
  virtual ~AliSymBDMatrix();
  //
  Int_t         GetBandHWidth()                                  const {return fNrows;}
  Int_t         GetNElemsStored()                                const {return fNelems;}
  void          Clear(Option_t* option="");
  void          Reset();
  //
  Float_t       GetDensity()                                     const;
  AliSymBDMatrix& operator=(const AliSymBDMatrix& src);
  Double_t      operator()(Int_t rown, Int_t coln)               const;
  Double_t&     operator()(Int_t rown, Int_t coln);
  Double_t      operator()(Int_t rown)                           const;
  Double_t&     operator()(Int_t rown);
  //
  Double_t      DiagElem(Int_t r)                                const {return (*(const AliSymBDMatrix*)this)(r,r);}
  Double_t&     DiagElem(Int_t r)                                      {return (*this)(r,r);}
  void          DecomposeLDLT();
  void          Solve(Double_t *rhs);
  void          Solve(const Double_t *rhs,Double_t *sol);
  void          Solve(TVectorD &rhs)                                   {Solve(rhs.GetMatrixArray());}
  void          Solve(const TVectorD &rhs,TVectorD &sol)               {Solve(rhs.GetMatrixArray(),sol.GetMatrixArray());}
  //
  void          Print(Option_t* option="")                       const;
  void          SetDecomposed(Bool_t v=kTRUE)                          {SetBit(kDecomposedBit,v);}
  Bool_t        IsDecomposed()                                   const {return TestBit(kDecomposedBit);}
  //
  void          MultiplyByVec(const Double_t* vecIn, Double_t* vecOut) const;
  void          MultiplyByVec(const TVectorD &vecIn, TVectorD &vecOut) const;
  void          AddToRow(Int_t r, Double_t *valc,Int_t *indc,Int_t n);
  //
  // protected:
  virtual             Int_t GetIndex(Int_t row,Int_t col)      const;
  virtual             Int_t GetIndex(Int_t diagID)             const;
  Double_t            GetEl(Int_t row, Int_t col)              const {return operator()(row,col);}
  void                SetEl(Int_t row, Int_t col,Double_t val)       {operator()(row,col) = val;}
  //
 protected:
  Double_t*  fElems;     //   Elements booked by constructor
  //
  ClassDef(AliSymBDMatrix,0) //Symmetric Matrix Class
};


//___________________________________________________________
inline Int_t AliSymBDMatrix::GetIndex(Int_t row,Int_t col) const
{
  // lower triangle band is actually filled
  if (row<col) Swap(row,col);
  col -= row;
  if (col < -GetBandHWidth()) return -1;
  return GetIndex(row) + col;
}

//___________________________________________________________
inline Int_t AliSymBDMatrix::GetIndex(Int_t diagID) const
{
  // Get index of the diagonal element on row diagID
  return (diagID+1)*fRowLwb-1;
}

//___________________________________________________________
inline Double_t AliSymBDMatrix::operator()(Int_t row, Int_t col) const
{
  // query element
  int idx = GetIndex(row,col);
  return (const Double_t&) idx<0 ? 0.0 : fElems[idx];
}

//___________________________________________________________
inline Double_t& AliSymBDMatrix::operator()(Int_t row, Int_t col)
{
  // get element for assingment; assignment outside of the stored range has no effect
  int idx = GetIndex(row,col);  
  if (idx>=0) return fElems[idx];
  fTol = 0; 
  return fTol;
}

//___________________________________________________________
inline Double_t AliSymBDMatrix::operator()(Int_t row) const
{
  // query diagonal 
  return (const Double_t&) fElems[GetIndex(row)];
}

//___________________________________________________________
inline Double_t& AliSymBDMatrix::operator()(Int_t row)
{
  // get diagonal for assingment; assignment outside of the stored range has no effect
  return fElems[GetIndex(row)];
}

//___________________________________________________________
inline void AliSymBDMatrix::MultiplyByVec(const TVectorD &vecIn, TVectorD &vecOut) const
{
  MultiplyByVec(vecIn.GetMatrixArray(), vecOut.GetMatrixArray());
}

#endif
