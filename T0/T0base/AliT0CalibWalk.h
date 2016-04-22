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
  TObjArray* GetfWalk() {return &fWalk;}    

  TGraph *GetQTC(Int_t ipmt )  const {return (TGraph*)fQTC.At(ipmt);} 
  TGraph *GetAmpLED(Int_t ipmt )  const {return (TGraph*)fAmpLED.At(ipmt);} 
  Bool_t MakeWalkCorrGraph(const char *laserFile);
  TGraph *  GetAmpLEDRec(Int_t ipmt) const   {return (TGraph*)fAmpLEDRec.At(ipmt);}
  void    GetMeanAndSigma(TH1F* hist, Float_t &mean, Float_t &sigma);
  void SetCalibByData(Bool_t calib) {fCalibByData=calib;}
  void SetWalk2015(TString filename);
  void SetWalkDima(TString filename);
  void SetWalkZero();
 protected:
   
   TObjArray   fWalk;  //time - amp. walk
   TObjArray   fAmpLEDRec;  //time - amp. LED-CFD for reconstruction
   TObjArray   fQTC;  //time - amp. walk
   TObjArray   fAmpLED;  //time - amp. LED-CFD for reconstruction
   Bool_t fCalibByData;
   //
   ClassDef(AliT0CalibWalk,6)    // T0 Amplitude Calibration data
     };

     typedef AliT0CalibWalk AliSTARTCalibWalk; // for backward compatibility

#endif

