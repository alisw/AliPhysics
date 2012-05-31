#ifndef AliT0TriggerParameters_H
#define AliT0TriggerParameters_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for T0 calibration                 //
////////////////////////////////////////////////

#include "TObject.h"

class AliT0TriggerParameters: public TObject {

 public:
  AliT0TriggerParameters();
  AliT0TriggerParameters& operator= (const AliT0TriggerParameters &);
  AliT0TriggerParameters(const AliT0TriggerParameters &calibda);
  virtual ~AliT0TriggerParameters();
  
  //Switched on/off
  void     SetPMTstatus(Int_t i, Int_t val);
  Int_t    GetPMTstatus(Int_t i) const;
 //Thresholds
  void     SetThreshold(Int_t i, Int_t val) {fThreshold[i]=val;}
  Int_t    GetThreshold(Int_t i) {return fThreshold[i];}
  Int_t*   GetThreshold()  {return  fThreshold;}
  //mult threshold
  Float_t GetTimeWindowLow()  const  {return fTimeWindowLow;}
  void    SetTimeWindowLow(Float_t low)   { fTimeWindowLow = low;}
  Float_t GetTimeWindowHigh()  const  {return fTimeWindowHigh;}
  void    SetTimeWindowHigh(Float_t high)   { fTimeWindowHigh = high;}
  Int_t   GetAmpCentr() const {return fAmpCentr;}
  void    SetAmpCentr(Int_t ref) {fAmpCentr = ref;}
  Int_t   GetAmpSemiCentr() const {return fAmpSemiCentr;}
  void    SetAmpSemiCentr(Int_t ref) {fAmpSemiCentr = ref;}
  

  void Reset();
  virtual void  Print(Option_t* option= "") const; 

 protected:

  Int_t fSwtPmt;    // PMT on/off
  Int_t fThreshold[24]; // thresholds
  Int_t fAmpCentr;         // threshold for central event
  Int_t fAmpSemiCentr;        // threshold for semi-central event
  Float_t fTimeWindowLow;      //low border for TVDC
  Float_t fTimeWindowHigh;     //low border for TVDC


  // AliT0TriggerParameters& operator= (const AliT0TriggerParameters &);
  // AliT0TriggerParameters(const AliT0TriggerParameters &calibda);
 //
  ClassDef(AliT0TriggerParameters,2)    // T0 Sensor Calibration data
};


#endif

