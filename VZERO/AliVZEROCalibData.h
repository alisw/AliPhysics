#ifndef ALIVZEROCALIBDATA_H
#define ALIVZEROCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//                                            // 
//  class for VZERO calibration               //
//                                            //
////////////////////////////////////////////////

#include "TNamed.h"
#include "AliVZERO.h"

class AliVZEROCalibData: public TNamed {

 public:
  AliVZEROCalibData();
  AliVZEROCalibData(const char* name);
  AliVZEROCalibData(const AliVZEROCalibData &calibda);
  AliVZEROCalibData& operator= (const AliVZEROCalibData &calibda);
  virtual ~AliVZEROCalibData();

  Float_t  GetPedestal(Int_t channel)   const {return fPedestal[channel];}
  Float_t* GetPedestal()   const {return (float*)fPedestal;}
  Float_t  GetGain(Int_t channel)	const {return fGain[channel];}
  Float_t* GetGain()   const {return (float*)fGain;}
  //
  void     SetPedestal(Float_t val, Int_t channel) {fPedestal[channel]=val;}
  void     SetPedestal(Float_t* Pedestal);
  void 	   SetGain(Float_t val, Int_t channel) {fGain[channel]=val;}
  void 	   SetGain(Float_t* Gain);

 protected:
  Float_t  fPedestal[80];     // Mean pedestal values
  Float_t  fGain[80];	      // PM gains

  ClassDef(AliVZEROCalibData,1)    // VZERO Calibration data
};

#endif
