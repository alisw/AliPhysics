#ifndef AliT0CalibWalk_H
#define AliT0CalibWalk_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for T0 amplitude calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TGraph.h"
#include "TObjArray.h"

class AliT0CalibWalk: public TNamed {

 public:
  AliT0CalibWalk();
  AliT0CalibWalk(const char* name);
  AliT0CalibWalk(const AliT0CalibWalk &calibda);
  AliT0CalibWalk& operator= (const AliT0CalibWalk &calibda);
  virtual ~AliT0CalibWalk();
  
  
  TGraph *GetWalk(Int_t ipmt )  const {return (TGraph*)fWalk.At(ipmt);}
  void SetWalk(Int_t ipmt) ;
  TObjArray* GetfWalk() {return &fWalk;}    

  TGraph *GetQTC(Int_t ipmt )  const {return (TGraph*)fQTC.At(ipmt);}
  Float_t  GetQTCpar(Int_t channel, Int_t ipar)    const 
  {return fQTCpar[channel][ ipar ];}
  void SetQTCpar(Int_t channel,Int_t ipar, Float_t val) 
  {fQTCpar[channel][ipar]=val;}
 
  TGraph *GetAmpLED(Int_t ipmt )  const {return (TGraph*)fAmpLED.At(ipmt);}
  Float_t  GetLEDpar(Int_t channel, Int_t ipar)    const 
  {return fAmpLEDpar[channel][ ipar ];}
  void SetAmpLEDpar(Int_t channel,Int_t ipar, Float_t val)
  {fAmpLEDpar[channel][ipar]=val;}
 
  void MakeWalkCorrGraph(const char *laserFile);
  
  

  TGraph *  GetAmpLEDRec(Int_t ipmt) const   {return (TGraph*)fAmpLEDRec.At(ipmt);}
  void      SetAmpLEDRec(Int_t ipmt) ;
   
 protected:
   
   TObjArray   fWalk;  //time - amp. walk
   TObjArray   fAmpLEDRec;  //time - amp. LED-CFD for reconstruction
   TObjArray   fQTC;  //time - amp. walk
   TObjArray   fAmpLED;  //time - amp. LED-CFD for reconstruction
   Float_t   fQTCpar[24][2]; // fitted parameters QTC amplitude
   Float_t   fAmpLEDpar[24][2]; // fitted parameters fAmpLEDpar amplitude
  
   //
   ClassDef(AliT0CalibWalk,4)    // T0 Amplitude Calibration data
     };

     typedef AliT0CalibWalk AliSTARTCalibWalk; // for backward compatibility

#endif

