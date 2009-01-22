#ifndef ALIMATRIXSQ_H
#define ALIMATRIXSQ_H

#include <TMatrixDBase.h>
#include <TVectorD.h>

class AliMatrixSq : public TMatrixDBase {
  //
 public:
  AliMatrixSq(): fSymmetric(kFALSE) {}
  AliMatrixSq(const AliMatrixSq &src) : TMatrixDBase(src), fSymmetric(src.fSymmetric) {}
  virtual ~AliMatrixSq() {}
  virtual Int_t   GetSize()                            const {return fNcols;}
  virtual Float_t GetDensity()                         const     = 0;
  //
  virtual void  Clear(Option_t* option="")                       = 0;//{Error("Clear","Dummy");}
  //
  virtual       Double_t      Querry(Int_t rown, Int_t coln)     const {return operator()(rown,coln);}
  virtual       Double_t      operator()(Int_t rown, Int_t coln) const = 0;//{Error("(i,j)","Dummy");return 0;}
  virtual       Double_t&     operator()(Int_t rown, Int_t coln) = 0;//{Error("(i,j)","Dummy");return 0;}
  //
  virtual       Double_t      QuerryDiag(Int_t rc)               const {return DiagElem(rc);}
  virtual       Double_t      DiagElem(Int_t r)                  const = 0;
  virtual       Double_t&     DiagElem(Int_t r)                  = 0;
  //
  virtual void  Print(Option_t* option="")           const       = 0;//{Error("Print","Dummy");}
  virtual void  Reset()                                          = 0;
  virtual void  PrintCOO()                           const;          // print in COO format
  //
  virtual void  MultiplyByVec(Double_t* vecIn, Double_t* vecOut) const;
  virtual void  MultiplyByVec(TVectorD &vecIn, TVectorD &vecOut) const;
  //
  Bool_t        IsSymmetric()                       const {return fSymmetric;}
  void          SetSymmetric(Bool_t v=kTRUE)              {fSymmetric = v;}
  //
  // ---------------------------------- Dummy methods of MatrixBase
  virtual       const Double_t   *GetMatrixArray  () const {Error("GetMatrixArray","Dummy"); return 0;};
  virtual             Double_t   *GetMatrixArray  ()       {Error("GetMatrixArray","Dummy"); return 0;};
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
  virtual void Allocate      (Int_t ,Int_t ,Int_t , Int_t ,Int_t ,Int_t ) 
  {Error("Allocate","Dummy"); return;}
  //
 protected:
  //
  void    Swap(int &r,int &c) const {int t=r;r=c;c=t;}
  //
 protected:
  //
  Bool_t        fSymmetric;     // is the matrix symmetric? Only lower triangle is filled
  //
  ClassDef(AliMatrixSq,1) //Square Matrix Class
};


//___________________________________________________________
inline void AliMatrixSq::MultiplyByVec(TVectorD &vecIn, TVectorD &vecOut) const
{
  MultiplyByVec(vecIn.GetMatrixArray(), vecOut.GetMatrixArray());
}


#endif
