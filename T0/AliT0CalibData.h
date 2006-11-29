#ifndef AliT0CalibData_H
#define AliT0CalibData_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for T0 calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TF1.h"
#include "AliT0CalibData.h"
#include "TMap.h"
#include "TGraph.h"
#include "TString.h"
#include "TObjArray.h"
#include "AliT0.h"

class AliT0CalibData: public TNamed {

 public:
  AliT0CalibData();
  AliT0CalibData(const char* name);
  AliT0CalibData(const AliT0CalibData &calibda);
  AliT0CalibData& operator= (const AliT0CalibData &calibda);
  virtual ~AliT0CalibData();
  void Reset();
  
  virtual void  Print(Option_t* option= "") const; 
  Float_t  GetTimeDelayCFD(Int_t channel) const {return fTimeDelayCFD[channel];}
  Float_t* GetTimeDelayCFD()  const  {return(float*) fTimeDelayCFD;}
  Float_t  GetTimeDelayLED(Int_t channel) const {return fTimeDelayLED[channel];}
  Float_t* GetTimeDelayLED()  const  {return(float*) fTimeDelayLED;}

  Float_t   GetGain(Int_t channel) const {return fGain[channel];}
  Float_t*  GetGain()  const {return (float*)fGain;}
  void     SetGain(Float_t val, Int_t channel)  {fGain[channel]=val;}
  void     SetGain(Float_t* Gain);
  
  Float_t  GetWalk(Int_t ipmt, Float_t mv )  const {return ((TF1*)fWalk.At(ipmt))->Eval(mv);}
  void SetWalk(Int_t ipmt, const Char_t *filename="calibr/re.root") ;

   TGraph *  GetSlew(Int_t ipmt) const   {return (TGraph*)fSlewingLED.At(ipmt);}
  Float_t  GetSlewingLED(Int_t ipmt, Float_t mv)  const 
      {return((TGraph*)fSlewingLED.At(ipmt))->Eval(mv);}
   TGraph *  GetSlewRec(Int_t ipmt) const   {return (TGraph*)fSlewingRec.At(ipmt);}
  Float_t  GetSlewingRec(Int_t ipmt, Float_t mv)  const 
      {return((TGraph*)fSlewingRec.At(ipmt))->Eval(mv);}

  void SetSlewingLED(Int_t ipmt, const Char_t *filename) ;
  void SetSlewingRec(Int_t ipmt, const Char_t *filename) ;

  void     SetTimeDelayCFD(Float_t val, Int_t channel) {fTimeDelayCFD[channel]=val;}
  void     SetTimeDelayCFD(Float_t* TimeDelay);
  void     SetTimeDelayLED(Float_t val, Int_t channel) {fTimeDelayLED[channel]=val;}
  void     SetTimeDelayLED(Float_t* TimeDelay);

  void SetTimeDelayTVD(Int_t r=150)   { fTimeDelayTVD = r; };
  Float_t GetTimeDelayTVD()   { return fTimeDelayTVD; }

 protected:

  Float_t  fTimeDelayCFD[24]; // Coeff. for time delay (24 different cables & CFD )
  Float_t  fTimeDelayLED[24]; // Coeff. for time delay (24 different cables & CFD )
  Float_t  fGain[24]; // Coeff. for gain (24 different cables & CFD )
  Float_t fTimeDelayTVD; //time delay for TVD (vertex trigger channel)
  TObjArray fWalk;  //time - amp. walk
  TObjArray fSlewingLED;  //time - amp.LED slew
  TObjArray fSlewingRec;  //time - amp. LED slew for reconstruction
  //
  ClassDef(AliT0CalibData,1)    // T0 Sensor Calibration data
};

typedef AliT0CalibData AliSTARTCalibData; // for backward compatibility

#endif

