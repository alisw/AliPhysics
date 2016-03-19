#ifndef ALISYMMATRIX_H
#define ALISYMMATRIX_H
/**********************************************************************************************/
/* Fast symmetric matrix with dynamically expandable size.                                    */
/* Only part can be used for matrix operations. It is defined as:                             */ 
/* fNCols: rows built by constructor (GetSizeBooked)                                          */ 
/* fNRows: number of rows added dynamically (automatically added on assignment to row)        */ 
/*         GetNRowAdded                                                                       */ 
/* fNRowIndex: total size (fNCols+fNRows), GetSize                                            */ 
/* fRowLwb   : actual size to used for given operation, by default = total size, GetSizeUsed  */ 
/*                                                                                            */ 
/* Author: ruben.shahoyan@cern.ch                                                             */
/*                                                                                            */ 
/**********************************************************************************************/

#include <TVectorD.h>
#include "AliMatrixSq.h"



class AliSymMatrix : public AliMatrixSq {
  //
 public:
  AliSymMatrix();
  AliSymMatrix(Int_t size);
  AliSymMatrix(const AliSymMatrix &mat);
  virtual ~AliSymMatrix();
  //
  void          Clear(Option_t* option="");
  void          Reset();
  //
  Int_t         GetSize()                                        const {return fNrowIndex;}
  Int_t         GetSizeUsed()                                    const {return fRowLwb;}
  Int_t         GetSizeBooked()                                  const {return fNcols;}
  Int_t         GetSizeAdded()                                   const {return fNrows;}
  Float_t       GetDensity()                                     const;
  AliSymMatrix& operator=(const AliSymMatrix& src);
  AliSymMatrix& operator+=(const AliSymMatrix& src);
  AliSymMatrix& operator-=(const AliSymMatrix& src);
  Double_t      operator()(Int_t rown, Int_t coln)               const;
  Double_t&     operator()(Int_t rown, Int_t coln);
  //
  Double_t      DiagElem(Int_t r)                                const {return (*(const AliSymMatrix*)this)(r,r);}
  Double_t&     DiagElem(Int_t r)                                      {return (*this)(r,r);}
  //
  Double_t*     GetRow(Int_t r);
  //
  void          Print(const Option_t* option="")                 const;
  void          AddRows(int nrows=1);
  void          SetSizeUsed(Int_t sz)                                  {fRowLwb = sz;}
  //
  void          Scale(Double_t coeff);
  Bool_t        Multiply(const AliSymMatrix& right);
  void          MultiplyByVec(const Double_t* vecIn, Double_t* vecOut) const;
  void          MultiplyByVec(const TVectorD &vecIn, TVectorD &vecOut) const;
  void          AddToRow(Int_t r, Double_t *valc, Int_t *indc,Int_t n);
  //
  // ---------------------------------- Dummy methods of MatrixBase
  virtual       const Double_t   *GetMatrixArray  () const {return fElems;};
  virtual             Double_t   *GetMatrixArray  ()       {return (Double_t*)fElems;}
  virtual       const Int_t      *GetRowIndexArray() const {Error("GetRowIndexArray","Dummy"); return 0;};
  virtual             Int_t      *GetRowIndexArray()       {Error("GetRowIndexArray","Dummy"); return 0;};
  virtual       const Int_t      *GetColIndexArray() const {Error("GetColIndexArray","Dummy"); return 0;};
  virtual             Int_t      *GetColIndexArray()       {Error("GetColIndexArray","Dummy"); return 0;};
  virtual             TMatrixDBase &SetRowIndexArray(Int_t *) {Error("SetRowIndexArray","Dummy"); return *this;}
  virtual             TMatrixDBase &SetColIndexArray(Int_t *) {Error("SetColIndexArray","Dummy"); return *this;}
  virtual             TMatrixDBase &GetSub(Int_t,Int_t,Int_t,Int_t,TMatrixDBase &,Option_t *) const {Error("GetSub","Dummy"); return *((TMatrixDBase*)this);}
  virtual             TMatrixDBase &SetSub(Int_t,Int_t,const TMatrixDBase &) {Error("GetSub","Dummy"); return *this;}
  virtual             TMatrixDBase &ResizeTo (Int_t,Int_t,Int_t) {Error("ResizeTo","Dummy"); return *this;}
  virtual             TMatrixDBase &ResizeTo (Int_t,Int_t,Int_t,Int_t,Int_t) {Error("ResizeTo","Dummy"); return *this;}
  //
  // ----------------------------- Choleski methods ----------------------------------------
  AliSymMatrix*       DecomposeChol();                  //Obtain Cholesky decomposition L matrix 
  void                InvertChol(AliSymMatrix* mchol);  //Invert using provided Choleski decomposition
  Bool_t              InvertChol();                     //Invert 
  Bool_t              SolveChol(Double_t *brhs, Bool_t invert=kFALSE);
  Bool_t              SolveChol(Double_t *brhs, Double_t *bsol,Bool_t invert=kFALSE);
  Bool_t              SolveChol(TVectorD &brhs, Bool_t invert=kFALSE);
  Bool_t              SolveChol(const TVectorD &brhs, TVectorD& bsol,Bool_t invert=kFALSE);
  Bool_t              SolveCholN(Double_t *bn, int nRHS, Bool_t invert=kFALSE);
  //
  int                 SolveSpmInv(double *vecB, Bool_t stabilize=kTRUE);

 protected:
  virtual             Int_t GetIndex(Int_t row,Int_t col)      const;
  Double_t            GetEl(Int_t row, Int_t col)              const {return operator()(row,col);}
  void                SetEl(Int_t row, Int_t col,Double_t val)       {operator()(row,col) = val;}
  //
 protected:
  Double_t*  fElems;     //   Elements booked by constructor
  Double_t** fElemsAdd;  //   Elements (rows) added dynamicaly
  //
  static AliSymMatrix* fgBuffer;  // buffer for fast solution
  static Int_t         fgCopyCnt;  // matrix copy counter
  ClassDef(AliSymMatrix,0) //Symmetric Matrix Class
};


//___________________________________________________________
inline Int_t AliSymMatrix::GetIndex(Int_t row,Int_t col) const
{
  // lower triangle is actually filled
  return ((row*(row+1))>>1) + col;
}

//___________________________________________________________
inline Double_t AliSymMatrix::operator()(Int_t row, Int_t col) const
{
  //
  if (row<col) Swap(row,col);
  if (row>=fNrowIndex) return 0;
  return (const Double_t&) (row<fNcols ? fElems[GetIndex(row,col)] : (fElemsAdd[row-fNcols])[col]);
}

//___________________________________________________________
inline Double_t& AliSymMatrix::operator()(Int_t row, Int_t col)
{
  if (row<col) Swap(row,col);
  if (row>=fNrowIndex) AddRows(row-fNrowIndex+1);
  return (row<fNcols ? fElems[GetIndex(row,col)] : (fElemsAdd[row-fNcols])[col]);
}

//___________________________________________________________
inline void AliSymMatrix::MultiplyByVec(const TVectorD &vecIn, TVectorD &vecOut) const
{
  MultiplyByVec(vecIn.GetMatrixArray(), vecOut.GetMatrixArray());
}

//___________________________________________________________
inline void AliSymMatrix::Scale(Double_t coeff)
{
  for (int i=fNrowIndex;i--;) for (int j=i;j--;) { double& el = operator()(i,j); if (el) el *= coeff;}
}

//___________________________________________________________
inline void AliSymMatrix::AddToRow(Int_t r, Double_t *valc, Int_t *indc,Int_t n)
{
  for (int i=n;i--;) (*this)(indc[i],r) += valc[i];
}

#endif
