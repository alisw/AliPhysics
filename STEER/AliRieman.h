#ifndef ALIRIEMANN_H
#define ALIRIEMANN_H

#include "TMatrixDSym.h"

class AliRieman : public TObject{
 public:
  AliRieman();
  AliRieman(Int_t capacity);
  AliRieman(const AliRieman &rieman);
  ~AliRieman();
  void Reset();
  void AddPoint(Float_t x, Float_t y, Float_t z, Float_t sy, Float_t sz);
  Int_t GetN() const {return fN;}
  Int_t GetCapacity() const {return fCapacity;}
  Float_t * GetX(){return fX;}
  Float_t * GetY(){return fY;}
  Float_t * GetZ(){return fZ;}
  Float_t * GetSy(){return fSy;}
  Float_t * GetSz(){return fSz;}
  void Update();
  void UpdatePol();
  Double_t*  GetParam(){return fParams;}
  const TMatrixDSym &  GetCovariance(){return *fCovar;}
  Double_t GetC(); 
  Double_t GetYat(Double_t x);
  Double_t GetZat(Double_t x);
  Double_t GetDYat(Double_t x);
  Double_t GetDZat(Double_t x);
 protected:
  // public:
  Int_t         fCapacity;  // capacity
  Int_t         fN;         // numebr of points
  Float_t      *fX;         //[fN] x coordinate
  Float_t      *fY;         //[fN] y coordinate
  Float_t      *fZ;         //[fN] z coordinate
  Float_t      *fSy;        //[fN] sigma y coordinate
  Float_t      *fSz;        //[fN] sigma z coordinate
  Double_t      fParams[6]; //Parameters
  TMatrixDSym  *fCovar;     //Covariance
  Double_t      fSumXY[9];  //sums for XY part
  Double_t      fSumXZ[9];  //sums for XZ part
  Bool_t        fConv;      // indicates convergation
 protected:
  
 private:
  ClassDef(AliRieman,1)  // Fast fit of helices on ITS RecPoints
};



#endif
