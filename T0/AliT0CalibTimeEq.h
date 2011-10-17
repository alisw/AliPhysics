#ifndef AliT0CalibTimeEq_H
#define AliT0CalibTimeEq_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for T0 calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TH1F.h"

class AliT0CalibTimeEq: public TNamed {

 public:
  AliT0CalibTimeEq();
  AliT0CalibTimeEq(const char* name);
  AliT0CalibTimeEq(const AliT0CalibTimeEq &calibda);
  AliT0CalibTimeEq& operator= (const AliT0CalibTimeEq &calibda);
  virtual ~AliT0CalibTimeEq();
  void Reset();
  
  virtual void  Print(Option_t* option= "") const; 
  
  Bool_t ComputeOnlineParams(const char* filePhys);
  Int_t ComputeOfflineParams(const char* filePhys, Float_t *timecdb,Float_t *cfdcdb, Int_t badpmt);
  Float_t  GetCFDvalue(Int_t channel,Int_t number)  const {return fCFDvalue[channel][number];}
  Float_t* GetCFDvalue()          const {return (float*)fCFDvalue;}
  Float_t  GetTimeEq(Int_t channel)        const {return fTimeEq[channel];}
  Float_t* GetTimeEq()          const {return (float*)fTimeEq;}
 
  Float_t  GetTimeEqRms(Int_t channel)        const {return fTimeEqRms[channel];}

  Float_t  GetMeanT0() const {return 1.;} // WARNING: USED IN AliT0Parameters!!!!
   void SetCFDvalue(Int_t channel, Int_t number, Float_t val) {fCFDvalue[channel][number]=val;}
  void SetTimeEq(Int_t channel, Float_t val) {fTimeEq[channel]=val;}
  void  SetTimeEqRms(Int_t channel, Float_t rms)  { fTimeEqRms[channel]=rms;}
  
  void SetMeanVertex(Float_t mean=0) { fMeanVertex = mean; };
  Float_t GetMeanVertex () {return fMeanVertex;};

  void SetRmsVertex(Float_t rms=0) { fRmsVertex = rms; };
  Float_t GetRmsVertex () {return fRmsVertex;};


 protected:

  void GetMeanAndSigma(TH1F* hist, Float_t &mean, Float_t &sigma); 
  Float_t     fCFDvalue[24][5];       // CFD values
  Float_t     fTimeEq[24];	      // Time Equalized for OCDB	 
  Float_t     fTimeEqRms[24];	      // RMS of Time Equalized for OCDB	 
  Float_t     fMeanVertex;            // mean of vertex distribution   
  Float_t     fRmsVertex;            // RMS of vertex distribution   
 
 //
  ClassDef(AliT0CalibTimeEq,4)    // T0 Sensor Calibration data
};

typedef AliT0CalibTimeEq AliSTARTCalibTimeEq; // for backward compatibility

#endif

