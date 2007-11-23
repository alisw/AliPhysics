#ifndef AliT0CalibTimeEq_H
#define AliT0CalibTimeEq_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for T0 calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"

class AliT0CalibTimeEq: public TNamed {

 public:
  AliT0CalibTimeEq();
  AliT0CalibTimeEq(const char* name);
  AliT0CalibTimeEq(const AliT0CalibTimeEq &calibda);
  AliT0CalibTimeEq& operator= (const AliT0CalibTimeEq &calibda);
  virtual ~AliT0CalibTimeEq();
  void Reset();
  
  virtual void  Print(Option_t* option= "") const; 
  void     SetMeanT0(Int_t mean=500) { fMeanT0 = mean; };
  Int_t    GetMeanT0 () {return fMeanT0;};
  void     SetTimeDelayTVD(Int_t r=150)   { fTimeDelayTVD = r; };
  Float_t  GetTimeDelayTVD()   { return fTimeDelayTVD; }
  
  void ComputeOnlineParams(char* name1, char* name2, char* canv, Int_t npeaks, Double_t sigma);
  Float_t  GetCFDvalue(Int_t channel,Int_t number)        const {return fCFDvalue[channel][number];}
  Float_t* GetCFDvalue()          const {return (float*)fCFDvalue;}
  Float_t  GetTimeEq(Int_t channel)        const {return fTimeEq[channel];}
  Float_t* GetTimeEq()          const {return (float*)fTimeEq;}
  void SetCFDvalue(Int_t channel, Int_t number, Float_t val) {fCFDvalue[channel][number]=val;}
  void SetTimeEq(Int_t channel, Float_t val) {fTimeEq[channel]=val;}
   

 protected:

  Float_t     fTimeDelayTVD;          //time delay for TVD (vertex trigger channel)
  Int_t       fMeanT0;                //mean of T0distribution with vertex=0;
  Float_t     fCFDvalue[24][5];       // CFD values
  Float_t     fTimeEq[24];	      // Time Equalized for OCDB	 

  //
  ClassDef(AliT0CalibTimeEq,2)    // T0 Sensor Calibration data
};

typedef AliT0CalibTimeEq AliSTARTCalibTimeEq; // for backward compatibility

#endif

