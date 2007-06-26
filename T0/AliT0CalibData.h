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
#include "AliT0LookUpValue.h"

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
  Float_t  GetTimeDelayDA(Int_t channel) const {return fTimeDelayDA[channel];}
  Float_t* GetTimeDelayDA()  const  {return(float*) fTimeDelayDA;}

  
  TGraph *GetWalk(Int_t ipmt )  const {return ((TGraph*)fWalk.At(ipmt));}
  Float_t  GetWalkVal(Int_t ipmt, Float_t mv )  const {return ((TGraph*)fWalk.At(ipmt))->Eval(mv);}
  void SetWalk(Int_t ipmt) ;

  //   TGraph *  GetAmpLED(Int_t ipmt) const   {return (TGraph*)fAmpLED.At(ipmt);}
  //  Float_t  GetAmpLEDVal(Int_t ipmt, Float_t mv)  const 
  //  {return((TGraph*)fAmpLED.At(ipmt))->Eval(mv);}
   TGraph *  GetAmpLEDRec(Int_t ipmt) const   {return (TGraph*)fAmpLEDRec.At(ipmt);}
  Float_t  GetAmpLEDRecVal(Int_t ipmt, Float_t mv)  const 
      {return((TGraph*)fAmpLEDRec.At(ipmt))->Eval(mv);}

  //  void SetAmpLED(Int_t ipmt) ;
  void     SetAmpLEDRec(Int_t ipmt) ;

  void     SetTimeDelayCFD(Float_t val, Int_t channel) {fTimeDelayCFD[channel]=val;}
  void     SetTimeDelayCFD(Float_t* TimeDelay);
  void     SetTimeDelayDA(Float_t val, Int_t channel) {fTimeDelayDA[channel]=val;}
  void     SetTimeDelayDA(Float_t* TimeDelay);

  void     SetTimeDelayTVD(Int_t r=150)   { fTimeDelayTVD = r; };
  Float_t  GetTimeDelayTVD()   { return fTimeDelayTVD; }

  void     ReadAsciiLookup(const Char_t *filename);
  Int_t    GetChannel(Int_t trm,  Int_t tdc, Int_t chain, Int_t channel);
  void     PrintLookup(Option_t* option= "", Int_t iTRM=0, Int_t iTDC=0, Int_t iChannel=0) const;
  TMap    *GetMapLookup(void) {return &fLookup;}
  Int_t    GetNumberOfTRMs() const {return fNumberOfTRMs;}
  void     SetNumberOfTRMs(Int_t ntrms=2) {fNumberOfTRMs = ntrms;}

  void     SetMeanT0(Int_t mean=500) { fMeanT0 = mean; };
  Int_t    GetMeanT0 () {return fMeanT0;};

 protected:

  Float_t     fTimeDelayCFD[24]; // Coeff. for time delay (24 different cables & CFD )
  Float_t     fTimeDelayDA[24]; // number of channel with mean time+delay if vertex=0 )
  Float_t     fTimeDelayTVD; //time delay for TVD (vertex trigger channel)
  Int_t       fMeanT0; //mean of T0distribution with vertex=0;
  TObjArray   fWalk;  //time - amp. walk
  TObjArray fAmpLEDRec;  //time - amp. LED-CFD for reconstruction
  TMap fLookup;           //lookup table
  Int_t fNumberOfTRMs;    // number of TRMs in setup

  //
  ClassDef(AliT0CalibData,7)    // T0 Sensor Calibration data
};

typedef AliT0CalibData AliSTARTCalibData; // for backward compatibility

#endif

