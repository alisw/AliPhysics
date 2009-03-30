#ifndef ALIZDCCALIB_H
#define ALIZDCCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for ZDC calibration -> CALIBCOEFF   //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TH1.h"
#include "AliCDBEntry.h"

class AliZDC;

class AliZDCCalib: public TNamed {

 public:
  AliZDCCalib();
  AliZDCCalib(const char* name);
  AliZDCCalib(const AliZDCCalib &calibda);
  AliZDCCalib& operator= (const AliZDCCalib &calibda);
  virtual ~AliZDCCalib();
  void Reset();
  virtual void  Print(Option_t *) const; 
  //
  Float_t  GetEnCalib(Int_t channel)	const {return fEnCalibration[channel];}
  Float_t* GetEnCalib()   		const {return (float*)fEnCalibration;}
  //
  Float_t  GetZN1EqualCoeff(Int_t channel) const {return fZN1EqualCoeff[channel];}
  Float_t* GetZN1EqualCoeffs()		   const {return (float*)fZN1EqualCoeff;}
  Float_t  GetZP1EqualCoeff(Int_t channel) const {return fZP1EqualCoeff[channel];}
  Float_t* GetZP1EqualCoeffs()		   const {return (float*)fZP1EqualCoeff;}
  Float_t  GetZN2EqualCoeff(Int_t channel) const {return fZN2EqualCoeff[channel];}
  Float_t* GetZN2EqualCoeffs()		   const {return (float*)fZN2EqualCoeff;}
  Float_t  GetZP2EqualCoeff(Int_t channel) const {return fZP2EqualCoeff[channel];}
  Float_t* GetZP2EqualCoeffs()		   const {return (float*)fZP2EqualCoeff;}

  void 	SetEnCalib(Int_t channel, Float_t val) {fEnCalibration[channel]=val;}
  void 	SetEnCalib(Float_t* EnCalib);
  //
  void 	SetZN1EqualCoeff(Int_t channel, Float_t val) {fZN1EqualCoeff[channel]=val;}
  void 	SetZN1EqualCoeff(Float_t* EqualCoeff);
  void 	SetZP1EqualCoeff(Int_t channel, Float_t val) {fZP1EqualCoeff[channel]=val;}
  void 	SetZP1EqualCoeff(Float_t* EqualCoeff);
  void 	SetZN2EqualCoeff(Int_t channel, Float_t val) {fZN2EqualCoeff[channel]=val;}
  void 	SetZN2EqualCoeff(Float_t* EqualCoeff);
  void 	SetZP2EqualCoeff(Int_t channel, Float_t val) {fZP2EqualCoeff[channel]=val;}
  void 	SetZP2EqualCoeff(Float_t* EqualCoeff);
  
 protected:
  // --- E calibration
  Float_t  fEnCalibration[6];	 // Coeff. for energy calibration
  // --- Coefficients for tower calibration
  Float_t  fZN1EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  Float_t  fZP1EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  Float_t  fZN2EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  Float_t  fZP2EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  //
  ClassDef(AliZDCCalib,5)    // ZDC calibration calibration data
};

#endif
