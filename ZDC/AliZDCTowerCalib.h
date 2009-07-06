#ifndef ALIZDCTOWERCALIB_H
#define ALIZDCTOWERCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for ZDC calibration -> CALIBCOEFF   //
////////////////////////////////////////////////

#include "TNamed.h"
#include "AliCDBEntry.h"

class AliZDC;

class AliZDCTowerCalib: public TNamed {

 public:
  AliZDCTowerCalib();
  AliZDCTowerCalib(const char* name);
  AliZDCTowerCalib(const AliZDCTowerCalib &calibda);
  AliZDCTowerCalib& operator= (const AliZDCTowerCalib &calibda);
  virtual ~AliZDCTowerCalib();
  void Reset();
  virtual void  Print(Option_t *) const; 
  //
  Float_t  GetZN1EqualCoeff(Int_t channel) const {return fZN1EqualCoeff[channel];}
  Float_t* GetZN1EqualCoeffs()		   const {return (float*)fZN1EqualCoeff;}
  Float_t  GetZP1EqualCoeff(Int_t channel) const {return fZP1EqualCoeff[channel];}
  Float_t* GetZP1EqualCoeffs()		   const {return (float*)fZP1EqualCoeff;}
  Float_t  GetZN2EqualCoeff(Int_t channel) const {return fZN2EqualCoeff[channel];}
  Float_t* GetZN2EqualCoeffs()		   const {return (float*)fZN2EqualCoeff;}
  Float_t  GetZP2EqualCoeff(Int_t channel) const {return fZP2EqualCoeff[channel];}
  Float_t* GetZP2EqualCoeffs()		   const {return (float*)fZP2EqualCoeff;}

  void 	SetZN1EqualCoeff(Int_t channel, Float_t val) {fZN1EqualCoeff[channel]=val;}
  void 	SetZN1EqualCoeff(Float_t* EqualCoeff);
  void 	SetZP1EqualCoeff(Int_t channel, Float_t val) {fZP1EqualCoeff[channel]=val;}
  void 	SetZP1EqualCoeff(Float_t* EqualCoeff);
  void 	SetZN2EqualCoeff(Int_t channel, Float_t val) {fZN2EqualCoeff[channel]=val;}
  void 	SetZN2EqualCoeff(Float_t* EqualCoeff);
  void 	SetZP2EqualCoeff(Int_t channel, Float_t val) {fZP2EqualCoeff[channel]=val;}
  void 	SetZP2EqualCoeff(Float_t* EqualCoeff);
  
 protected:
  // --- Coefficients for tower calibration
  Float_t  fZN1EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  Float_t  fZP1EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  Float_t  fZN2EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  Float_t  fZP2EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  //
  ClassDef(AliZDCTowerCalib,3)    // ZDC calibration calibration data
};

#endif
