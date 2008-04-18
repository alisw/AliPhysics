#ifndef AliT0CalibWalk_H
#define AliT0CalibWalk_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for T0 calibration                 //
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
  // void SetQTC(Int_t ipmt) ;
  
  TGraph *GetAmpLED(Int_t ipmt )  const {return (TGraph*)fAmpLED.At(ipmt);}
  //  void SetAmpLED(Int_t ipmt) ;
 
  void MakeWalkCorrGraph(const char *laserFile);
  
  

  TGraph *  GetAmpLEDRec(Int_t ipmt) const   {return (TGraph*)fAmpLEDRec.At(ipmt);}
  void     SetAmpLEDRec(Int_t ipmt) ;
   
 protected:
   
   TObjArray   fWalk;  //time - amp. walk
   TObjArray fAmpLEDRec;  //time - amp. LED-CFD for reconstruction
   TObjArray   fQTC;  //time - amp. walk
   TObjArray fAmpLED;  //time - amp. LED-CFD for reconstruction
  
   //
   ClassDef(AliT0CalibWalk,2)    // T0 Sensor Calibration data
     };

     typedef AliT0CalibWalk AliSTARTCalibWalk; // for backward compatibility

#endif

