#ifndef ALIZDCPEDESTALS_H
#define ALIZDCPEDESTALS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for ZDC calibration -> PEDESTALS    //
////////////////////////////////////////////////

#include "TNamed.h"

class AliZDC;

class AliZDCPedestals: public TNamed {

 public:
  AliZDCPedestals();
  AliZDCPedestals(const char* name);
  AliZDCPedestals(const AliZDCPedestals &calibda);
  AliZDCPedestals& operator= (const AliZDCPedestals &calibda);
  virtual ~AliZDCPedestals();
  void Reset();
  enum ESubPedModeBit{
     kPedSubModeFromOCDB = 14
  };

  virtual void  Print(Option_t *) const; 
  //
  Float_t  GetMeanPed(Int_t channel)   	   const {return fMeanPedestal[channel];}
  Float_t* GetMeanPed()   		   const {return (float*)fMeanPedestal;}
  Float_t  GetMeanPedWidth(Int_t channel)  const {return fMeanPedWidth[channel];}
  Float_t* GetMeanPedWidth()   		   const {return (float*)fMeanPedWidth;}
  Float_t  GetOOTPed(Int_t channel)   	   const {return fOOTPedestal[channel];}
  Float_t* GetOOTPed()   		   const {return (float*)fOOTPedestal;}
  Float_t  GetOOTPedWidth(Int_t channel)   const {return fOOTPedWidth[channel];}
  Float_t* GetOOTPedWidth()   	           const {return (float*)fOOTPedWidth;}
  Float_t  GetPedCorrCoeff0(Int_t channel) const {return fPedCorrCoeff[0][channel];}
  Float_t  GetPedCorrCoeff1(Int_t channel) const {return fPedCorrCoeff[1][channel];}
  Float_t* GetPedCorrCoeff()		   const {return (float*)fPedCorrCoeff;}
  
  UInt_t   GetPedSubModefromOCDB()         const {return fPedSubModefromOCDB;}
  Bool_t   TestPedModeBit()	   	   const {return TESTBIT(fPedSubModefromOCDB, kPedSubModeFromOCDB);}
  Bool_t   GetUseCorrFit(int ich)	   const 
  	   {if(AliZDCPedestals::TestPedModeBit()) return fUseCorrFit[ich]; else return kFALSE;}

  void  SetMeanPed(Int_t channel, Float_t val) {fMeanPedestal[channel]=val;}
  void  SetMeanPed(Float_t* MeanPed);
  void  SetMeanPedWidth(Int_t channel, Float_t val) {fMeanPedWidth[channel]=val;}
  void  SetMeanPedWidth(Float_t* MeanPedWidth);
  void  SetOOTPed(Int_t channel, Float_t val) {fOOTPedestal[channel]=val;}
  void  SetOOTPed(Float_t* OOTPed);
  void  SetOOTPedWidth(Int_t channel, Float_t val) {fOOTPedWidth[channel]=val;}
  void  SetOOTPedWidth(Float_t* OOTPedWidth);
  void  SetPedCorrCoeff(Int_t channel, Float_t valCoeff0, Float_t valCoeff1)
  	{fPedCorrCoeff[0][channel]=valCoeff0; fPedCorrCoeff[1][channel]=valCoeff1;}
  void  SetPedCorrCoeff(Float_t* PedCorrCoeff);
  void  SetPedCorrCoeff(Float_t* PedCorrCoeff0, Float_t* PedCorrCoeff1);
  
  void  SetPedModeBit(Bool_t on=kTRUE) 
  	{on ? SETBIT(fPedSubModefromOCDB, kPedSubModeFromOCDB) : CLRBIT(fPedSubModefromOCDB, kPedSubModeFromOCDB);}
  void  SetSubFromCorr(int ich) {fUseCorrFit[ich] = kTRUE;}
  
 protected:
  // --- Pedestals
  Float_t  fMeanPedestal[48];	 // Mean pedestal values 
  Float_t  fMeanPedWidth[48];	 // Mean pedestal widths 
  Float_t  fOOTPedestal[48];	 // "Out of Time" pedestal values
  Float_t  fOOTPedWidth[48];	 // "Out of Time" pedestal widths
  Float_t  fPedCorrCoeff[2][48]; // Fit of correlation in-time vs. out-of-time
  UInt_t   fPedSubModefromOCDB;  // test whether the OCDB provides the ped sub mode ch.bych. (from RUN2)
  Bool_t   fUseCorrFit[24];      // if pedestal subtraction mode is from OCDB decide the mode!
  //
  ClassDef(AliZDCPedestals,4)    // ZDC pedestal calibration data
};

#endif
