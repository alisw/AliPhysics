#ifndef ALIT0RECOPARAM_H
#define ALIT0RECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with T0 reconstruction parameters                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliDetectorRecoParam.h"

class AliT0RecoParam : public AliDetectorRecoParam
{
 public: 
  AliT0RecoParam();
  AliT0RecoParam(const AliT0RecoParam &p); //copy constructor
  AliT0RecoParam& operator=(const AliT0RecoParam &p);
  virtual ~AliT0RecoParam();
 
  
 static   AliT0RecoParam *GetLowFluxParam();        // make reco parameters for low  flux env
  static   AliT0RecoParam *GetHighFluxParam();       // make reco parameters for high flux env 
  static   AliT0RecoParam *GetLaserTestParam();  // special setting for laser SetLaserTestParam 
  //for monitoring
  //  static   AliT0RecoParam *GetHistRange();  //  limit of monitoring histograms

  Float_t GetRefAmp()  const  {return fRefAmp;}
  void    SetRefAmp(Float_t amp)   { fRefAmp = amp;}
  Int_t   GetRefPoint() const {return fRefPoint;}
  void    SetRefPoint(Int_t ref) {fRefPoint = ref;}
 
  Float_t   GetLow(Int_t numhist) const {return fLow[numhist];}
  //  Float_t   GetLow() {return *fLow;}
  void      SetLow(Int_t numhist, Float_t low) {fLow[numhist] = low;}

  Float_t   GetHigh(Int_t numhist) const  {return fHigh[numhist];}
  //  Float_t   GetHigh()  {return *fHigh;}
  void      SetHigh(Int_t numhist, Float_t high) {fHigh[numhist] = high;}
  
  void PrintParameters() const;
  
 protected:
  Float_t fRefAmp;
  Int_t   fRefPoint;   
  Float_t   fLow[500];
  Float_t   fHigh[500];



  ClassDef(AliT0RecoParam, 1);

};
#endif
