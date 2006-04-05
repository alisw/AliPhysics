#ifndef ALIZDCCALIBDATA_H
#define ALIZDCCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for ZDC calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TH1.h"
#include "AliCDBEntry.h"

class AliZDC;

class AliZDCCalibData: public TNamed {

 public:
  AliZDCCalibData();
  AliZDCCalibData(const char* name);
  AliZDCCalibData(const AliZDCCalibData &calibda);
  AliZDCCalibData& operator= (const AliZDCCalibData &calibda);
  virtual ~AliZDCCalibData();
  void Reset();
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
  //
  Float_t  GetEnCalib(Int_t channel)	const {return fEnCalibration[channel];}
  Float_t* GetEnCalib()   		const {return (float*)fEnCalibration;}
  //
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
  void  SetPedCorrCoeff(Float_t* PedCorrCoeff0, Float_t* PedCorrCoeff1);
  //
  void 	SetEnCalib(Int_t channel, Float_t val) {fEnCalibration[channel]=val;}
  void 	SetEnCalib(Float_t* EnCalib);
  //
//  void     PrepHistos();
//  TH1F*    GetHistMeanPed() const {return fHistMeanPed;}
//  void     CleanHistos();

 protected:
  // --- Pedestals
  Float_t  fMeanPedestal[47];	 // Mean pedestal values 
  Float_t  fMeanPedWidth[47];	 // Mean pedestal values 
  Float_t  fOOTPedestal[44];	 // "Out of Time" pedestal values
  Float_t  fOOTPedWidth[44];	 // "Out of Time" pedestal values
  Float_t  fPedCorrCoeff[2][44]; // Fit of correlation in-time vs. out-of-time
  // --- E calibration
  Float_t  fEnCalibration[6];	 // Coeff. for energy calibration
  //TH1F*    fHistMeanPed;        //! histos for drawing mean pedestals
  //
  ClassDef(AliZDCCalibData,3)    // ZDC  Calibration data
};

#endif
