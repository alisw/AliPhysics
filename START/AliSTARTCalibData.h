#ifndef AliSTARTCalibData_H
#define AliSTARTCalibData_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for START calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TF1.h"
#include "AliSTARTCalibData.h"
#include "TMap.h"
#include "TGraph.h"
#include "TString.h"
#include "TObjArray.h"
#include "AliSTART.h"

class AliSTARTCalibData: public TNamed {

 public:
  AliSTARTCalibData();
  AliSTARTCalibData(const char* name);
  AliSTARTCalibData(const AliSTARTCalibData &calibda);
  AliSTARTCalibData& operator= (const AliSTARTCalibData &calibda);
  virtual ~AliSTARTCalibData();
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

  void SetSlewingLED(Int_t ipmt, const Char_t *filename) ;

  void     SetTimeDelayCFD(Float_t val, Int_t channel) {fTimeDelayCFD[channel]=val;}
  void     SetTimeDelayCFD(Float_t* TimeDelay);
  void     SetTimeDelayLED(Float_t val, Int_t channel) {fTimeDelayLED[channel]=val;}
  void     SetTimeDelayLED(Float_t* TimeDelay);


 protected:

  Float_t  fTimeDelayCFD[24]; // Coeff. for time delay (24 different cables & CFD )
  Float_t  fTimeDelayLED[24]; // Coeff. for time delay (24 different cables & CFD )
  Float_t  fGain[24]; // Coeff. for gain (24 different cables & CFD )
  TObjArray fWalk;  //time - amp. walk
  TObjArray fSlewingLED;  //time - amp. walk
  //
  ClassDef(AliSTARTCalibData,1)    // START Sensor Calibration data
};

#endif

