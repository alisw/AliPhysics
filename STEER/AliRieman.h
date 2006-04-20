#ifndef ALIRIEMAN_H
#define ALIRIEMAN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
// Class for global helix fit of a track
// Author: M.Ivanov
// This class uses decomposition of the chi2 based on the fact that
// one can rotate the coordinate system and provide xi >> yi for each
// space point


#include <TMatrixDSymfwd.h>

class AliRieman : public TObject{
 public:
  AliRieman();
  AliRieman(Int_t capacity);
  AliRieman(const AliRieman &rieman);
  ~AliRieman();
  void Reset();
  void AddPoint(Double_t x, Double_t y, Double_t z, Double_t sy, Double_t sz);
  Int_t GetN() const {return fN;}
  Int_t GetCapacity() const {return fCapacity;}
  Double_t * GetX(){return fX;}
  Double_t * GetY(){return fY;}
  Double_t * GetZ(){return fZ;}
  Double_t * GetSy(){return fSy;}
  Double_t * GetSz(){return fSz;}
  void Update();
  void UpdatePol();
  Double_t*  GetParam(){return fParams;}
  const TMatrixDSym &  GetCovariance() const {return *fCovar;}
  Double_t GetC() const;
  Double_t GetYat(Double_t x) const;
  Double_t GetZat(Double_t x) const; 
  Double_t GetDYat(Double_t x) const;
  Double_t GetDZat(Double_t x) const;
  //
  Double_t GetChi2Y() const { return fChi2Y;}
  Double_t GetChi2Z() const { return fChi2Z;}
  Double_t GetChi2() const  { return fChi2; }

  Double_t CalcChi2Y() const;  
  Double_t CalcChi2Z() const;
  Double_t CalcChi2() const;
  AliRieman * MakeResiduals() const;
  //
 protected:
  // public:
  Int_t         fCapacity;  // capacity
  Int_t         fN;         // numebr of points
  Double_t      *fX;         //[fN] x coordinate
  Double_t      *fY;         //[fN] y coordinate
  Double_t      *fZ;         //[fN] z coordinate
  Double_t      *fSy;        //[fN] sigma y coordinate
  Double_t      *fSz;        //[fN] sigma z coordinate
  Double_t      fParams[6]; //Parameters
  TMatrixDSym  *fCovar;     //Covariance
  Double_t      fSumXY[9];  //sums for XY part
  Double_t      fSumXZ[9];  //sums for XZ part
  Double_t      fChi2;      //sums of chi2
  Double_t      fChi2Y;     //sums of chi2 for y coord
  Double_t      fChi2Z;     //sums of chi2 foz z coord 
  Bool_t        fConv;      // indicates convergation
 protected:  
 private:
  AliRieman& operator=(const AliRieman &rieman){return *this;}
  ClassDef(AliRieman,1)  // Fast fit of helices on ITS RecPoints
};



#endif
