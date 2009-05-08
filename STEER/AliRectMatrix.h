#ifndef ALIRECTMATRIX_H
#define ALIRECTMATRIX_H

#include "TObject.h"
class TString;

class AliRectMatrix : public TObject {
  //
 public:
  AliRectMatrix();
  AliRectMatrix(Int_t nrow,Int_t ncol);
  AliRectMatrix(const AliRectMatrix &src);
  virtual ~AliRectMatrix();
  //
  Int_t         GetNRows()                            const {return fNRows;}
  Int_t         GetNCols()                            const {return fNCols;}
  //
  Double_t      Query(Int_t rown, Int_t coln)         const {return operator()(rown,coln);}

  AliRectMatrix& operator=(const AliRectMatrix& src);
  Double_t      operator()(Int_t rown, Int_t coln)    const;
  Double_t&     operator()(Int_t rown, Int_t coln);
  Double_t*     operator()(Int_t row)                 const {return GetRow(row);}
  Double_t*     GetRow(Int_t row)                     const {return fRows[row];}
  //
  void          Reset();
  //
  virtual void  Print(Option_t* option="")           const;
  //
 protected:
  //
  Int_t   fNRows;       // Number of rows
  Int_t   fNCols;       // Number of columns
  Double_t **fRows;     // pointers on rows
  //
  ClassDef(AliRectMatrix,0) //Rectangular Matrix Class
};

//___________________________________________________________
inline Double_t AliRectMatrix::operator()(Int_t row, Int_t col) const
{
  return (const Double_t&) GetRow(row)[col];
}

//___________________________________________________________
inline Double_t& AliRectMatrix::operator()(Int_t row, Int_t col)
{
  return (Double_t&) fRows[row][col];
}


#endif
