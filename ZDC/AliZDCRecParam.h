#ifndef ALIZDCRECPARAM_H
#define ALIZDCRECPARAM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for ZDC calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"

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
  ClassDef(AliZDCRecParam,3)    // ZDC  Calibration data
};

#endif
