#ifndef AliFITCalibTimeEq_H
#define AliFITCalibTimeEq_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for FIT calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TH1F.h"

class AliFITCalibTimeEq: public TNamed {

 public:
  AliFITCalibTimeEq();
  AliFITCalibTimeEq(const char* name);
  AliFITCalibTimeEq(const AliFITCalibTimeEq &calibda);
  AliFITCalibTimeEq& operator= (const AliFITCalibTimeEq &calibda);
  virtual ~AliFITCalibTimeEq();
  void Reset();
  
  virtual void  Print(Option_t* option= "") const; 
  
  Float_t  GetCFDvalue(Int_t channel)  const {return fCFDvalue[channel];}
  Float_t* GetCFDvalue()          const {return (Float_t*)fCFDvalue;}
  Float_t  GetTimeEq(Int_t channel)        const {return fTimeEq[channel];}
  Float_t* GetTimeEq()          const {return (float*)fTimeEq;}
 
  void SetCFDvalue(Int_t channel, Float_t val) {fCFDvalue[channel]=val;}
  void SetTimeEq(Int_t channel, Float_t val) {fTimeEq[channel]=val;}
 
  

 protected:

  void GetMeanAndSigma(TH1F* hist, Float_t &mean, Float_t &sigma); 
  Float_t     *fCFDvalue;       // CFD values
  Float_t     *fTimeEq;	        // Time Equalized for OCDB	 
  
 //
  ClassDef(AliFITCalibTimeEq,1)    // T0 Sensor Calibration data
};

#endif

