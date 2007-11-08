#ifndef ALIZDCRECPARAM_H
#define ALIZDCRECPARAM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for ZDC calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TH1.h"
#include "AliCDBEntry.h"

class AliZDC;

class AliZDCRecParam: public TNamed {

 public:
  AliZDCRecParam();
  AliZDCRecParam(const char* name);
  AliZDCRecParam(const AliZDCRecParam &calibda);
  AliZDCRecParam& operator= (const AliZDCRecParam &calibda);
  virtual ~AliZDCRecParam();
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
  Float_t  GetZN1EqualCoeff(Int_t channel) const {return fZN1EqualCoeff[channel];}
  Float_t* GetZN1EqualCoeffs()		   const {return (float*)fZN1EqualCoeff;}
  Float_t  GetZP1EqualCoeff(Int_t channel) const {return fZP1EqualCoeff[channel];}
  Float_t* GetZP1EqualCoeffs()		   const {return (float*)fZP1EqualCoeff;}
  Float_t  GetZN2EqualCoeff(Int_t channel) const {return fZN2EqualCoeff[channel];}
  Float_t* GetZN2EqualCoeffs()		   const {return (float*)fZN2EqualCoeff;}
  Float_t  GetZP2EqualCoeff(Int_t channel) const {return fZP2EqualCoeff[channel];}
  Float_t* GetZP2EqualCoeffs()		   const {return (float*)fZP2EqualCoeff;}
  //
  Float_t GetZEMEndValue()     const {return fZEMEndValue;}
  Float_t GetZEMCutFraction()  const {return fZEMCutFraction;}
  Float_t GetDZEMSup()	       const {return fDZEMSup;}
  Float_t GetDZEMInf()	       const {return fDZEMInf;}
  //
  Float_t GetEZN1MaxValue()  const {return fEZN1MaxValue;}
  Float_t GetEZP1MaxValue()  const {return fEZP1MaxValue;}
  Float_t GetEZDC1MaxValue() const {return fEZDC1MaxValue;}
  Float_t GetEZN2MaxValue()  const {return fEZN2MaxValue;}
  Float_t GetEZP2MaxValue()  const {return fEZP2MaxValue;}
  Float_t GetEZDC2MaxValue() const {return fEZDC2MaxValue;}

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
  //
  void 	SetEnCalib(Int_t channel, Float_t val) {fEnCalibration[channel]=val;}
  void 	SetEnCalib(Float_t* EnCalib);
  //
  void 	SetZN1EqualCoeff(Int_t channel, Float_t val) {fZN1EqualCoeff[channel]=val;}
  void 	SetZN1EqualCoeff(Float_t* EqualCoeff);
  void 	SetZP1EqualCoeff(Int_t channel, Float_t val) {fZP1EqualCoeff[channel]=val;}
  void 	SetZP1EqualCoeff(Float_t* EqualCoeff);
  void 	SetZN2EqualCoeff(Int_t channel, Float_t val) {fZN2EqualCoeff[channel]=val;}
  void 	SetZN2EqualCoeff(Float_t* EqualCoeff);
  void 	SetZP2EqualCoeff(Int_t channel, Float_t val) {fZP2EqualCoeff[channel]=val;}
  void 	SetZP2EqualCoeff(Float_t* EqualCoeff);
  //  
  void  SetZEMEndValue(Float_t ZEMEndValue) {fZEMEndValue = ZEMEndValue;}
  void  SetZEMCutFraction(Float_t ZEMCutFraction) {fZEMCutFraction = ZEMCutFraction;}
  void  SetDZEMSup(Float_t DZEMSup) {fDZEMSup = DZEMSup;}
  void  SetDZEMInf(Float_t DZEMInf) {fDZEMInf = DZEMInf;}
  //
  void	SetEZN1MaxValue(Float_t value)  {fEZN1MaxValue = value;}
  void	SetEZP1MaxValue(Float_t value)  {fEZP1MaxValue = value;}
  void	SetEZDC1MaxValue(Float_t value) {fEZDC1MaxValue = value;}
  void	SetEZN2MaxValue(Float_t value)  {fEZN2MaxValue = value;}
  void	SetEZP2MaxValue(Float_t value)  {fEZP2MaxValue = value;}
  void	SetEZDC2MaxValue(Float_t value) {fEZDC2MaxValue = value;}
  
 protected:
  // --- Pedestals
  Float_t  fMeanPedestal[48];	 // Mean pedestal values 
  Float_t  fMeanPedWidth[48];	 // Mean pedestal widths 
  Float_t  fOOTPedestal[48];	 // "Out of Time" pedestal values
  Float_t  fOOTPedWidth[48];	 // "Out of Time" pedestal widths
  Float_t  fPedCorrCoeff[2][48]; // Fit of correlation in-time vs. out-of-time
  // --- E calibration
  Float_t  fEnCalibration[6];	 // Coeff. for energy calibration
  // --- Coefficients for tower calibration
  Float_t  fZN1EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  Float_t  fZP1EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  Float_t  fZN2EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  Float_t  fZP2EqualCoeff[5];	 // Equalization coefficients for ZN1 PTMs
  // --- Coefficients for centrality selection from ZEM signal
  Float_t  fZEMEndValue;    	 // End point value of ZEM energy spectrum
  Float_t  fZEMCutFraction; 	 // Fraction of ZEM energy spectrum used to cut
  Float_t  fDZEMSup;// Upper value of EZDCvs.ZEM correlation where ZEM signal is used
  Float_t  fDZEMInf;// Lower value of EZDCvs.ZEM correlation where ZEM signal is used
  // --- Parameters from EZDC vs. Nspec correlation
  Float_t  fEZN1MaxValue;	 // Max value of ZN1 vs. Nspec n correlation
  Float_t  fEZP1MaxValue;	 // Max value of ZP1 vs. Nspec p correlation
  Float_t  fEZDC1MaxValue;	 // Max value of ZDC1 vs. Nspec n+p correlation
  Float_t  fEZN2MaxValue;	 // Max value of ZN2 vs. Nspec n correlation
  Float_t  fEZP2MaxValue;	 // Max value of ZP2 vs. Nspec p correlation
  Float_t  fEZDC2MaxValue;	 // Max value of ZDC2 vs. Nspec n+p correlation
  //
  ClassDef(AliZDCRecParam,11)    // ZDC  Calibration data
};

#endif
