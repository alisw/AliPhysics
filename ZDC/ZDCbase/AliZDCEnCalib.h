#ifndef ALIZDCENCALIB_H
#define ALIZDCENCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for ZDC calibration -> CALIBCOEFF   //
////////////////////////////////////////////////

#include "TNamed.h"
#include "AliCDBEntry.h"

class AliZDC;

class AliZDCEnCalib: public TNamed {

 public:
  AliZDCEnCalib();
  AliZDCEnCalib(const char* name);
  AliZDCEnCalib(const AliZDCEnCalib &calibda);
  AliZDCEnCalib& operator= (const AliZDCEnCalib &calibda);
  virtual ~AliZDCEnCalib();
  void Reset();
  virtual void  Print(Option_t *) const; 
  //
  Float_t  GetEnCalib(Int_t channel)	const {return fEnCalibration[channel];}
  Float_t* GetEnCalib()   		const {return (float*)fEnCalibration;}

  void 	SetEnCalib(Int_t channel, Float_t val) {fEnCalibration[channel]=val;}
  void 	SetEnCalib(Float_t* EnCalib);
  
 protected:
  // --- E calibration
  Float_t  fEnCalibration[6];	 // Coeff. for energy calibration
  //
  ClassDef(AliZDCEnCalib,3)    // ZDC calibration calibration data
};

#endif
