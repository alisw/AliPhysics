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
  void Reset();

  Float_t  GetPedestal(Int_t channel)   const {return fPedestal[channel];}
  Float_t* GetPedestal()   const {return (float*)fPedestal;}
  Float_t  GetSigma(Int_t channel)   const {return fSigma[channel];}
  Float_t* GetSigma()   const {return (float*)fSigma;}
  Float_t  GetGain(Int_t channel)	const {return fGain[channel];}
  Float_t* GetGain()   const {return (float*)fGain;}
  
  Float_t  GetTimeOffset(Int_t channel)	const {return fTimeOffset[channel];}
  Float_t* GetTimeOffset()   const {return (float*)fTimeOffset;}
  Float_t  GetTimeGain(Int_t channel)	const {return fTimeGain[channel];}
  Float_t* GetTimeGain()   const {return (float*)fTimeGain;}
  
  Float_t  GetMeanHV(Int_t channel)	const {return fMeanHV[channel];}
  Float_t* GetMeanHV()   const {return (float*)fMeanHV;}
 
  Float_t  GetWidthHV(Int_t channel)	const {return fWidthHV[channel];}
  Float_t* GetWidthHV()   const {return (float*)fWidthHV;}
 
  
  void     SetPedestal(Float_t val, Int_t channel) {fPedestal[channel]=val;}
  void     SetPedestal(Float_t* Pedestal);
  void     SetSigma(Float_t val, Int_t channel) {fSigma[channel]=val;}
  void     SetSigma(Float_t* Sigma);
  void 	   SetGain(Float_t val, Int_t channel) {fGain[channel]=val;}
  void 	   SetGain(Float_t* Gain);

  void     SetTimeOffset(Float_t val, Int_t channel) {fTimeOffset[channel]=val;}
  void     SetTimeOffset(Float_t* TimeOffset);
  void     SetTimeGain(Float_t val, Int_t channel) {fTimeGain[channel]=val;}
  void     SetTimeGain(Float_t* TimeGain);
  
  void     SetMeanHV(Float_t val, Int_t channel) {fMeanHV[channel]=val;}
  void     SetMeanHV(Float_t* MeanHVs);
  
  void     SetWidthHV(Float_t val, Int_t channel) {fWidthHV[channel]=val;}
  void     SetWidthHV(Float_t* WidthHVs);

  

 protected:
  Float_t  fPedestal[128];     // Mean pedestal values
  Float_t  fSigma[128];        // Sigmas of pedestal peaks
  Float_t  fGain[128];	       // PM gains
  
  Float_t  fTimeOffset[64];
  Float_t  fTimeGain[64];

  Float_t  fMeanHV[64];    		// Mean PMT HV
  Float_t  fWidthHV[64];		// Width of the PMT HV

  ClassDef(AliVZEROCalibData,2)    // VZERO Calibration data
};

#endif
