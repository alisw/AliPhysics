#ifndef ALIPHOSV3_H
#define ALIPHOSV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Implementation version v1 of PHOS Manager class 
// The main goal of this version of AliPHOS is to calculte the 
//  induced charged in the PIN diode, taking into account light
//  tracking in the PbWO4 crystal, induced signal in the 
//  PIN due to MIPS particle and electronic noise.
// This is done in the StepManager 
//                  
//*-- Author:  Odd Harald Oddland & Gines Martinez (SUBATECH)

// --- ROOT system ---


// --- AliRoot header files ---
#include "AliPHOSv1.h"


class AliPHOSv3 : public AliPHOSv1 {
  
 public:
  
  AliPHOSv3(void) ; 
  AliPHOSv3(const char *name, const char *title="") ;
  //  AliPHOSv3(AliPHOSReconstructioner * Reconstructioner, const char *name, const char *title="") ;
  virtual ~AliPHOSv3(void) {
    // dtor
  } 
  
  virtual Int_t   IsVersion(void) const { 
    // Gives the version number 
    return 3 ; 
  }
  virtual const TString Version(void)const { 
    // returns the version number 
    return TString("v3") ; 
  }   
  virtual void   StepManager(void) ;
  
  Float_t GetLightYieldMean()         const { return  fLightYieldMean ;}
  Float_t GetLightYieldAttenuation()  const { return  fLightYieldAttenuation ;}
  Float_t GetRecalibrationFactor()    const { return  fRecalibrationFactor ;}
  Float_t GetAPDGain()                const { return  fAPDGain ;}
  Float_t GetIntrinsicPINEfficiency() const { return  fIntrinsicPINEfficiency ;}
  Float_t GetElectronsPerGeV()        const { return  fElectronsPerGeV ;}

  void    SetLightYieldMean(Float_t LightYieldMean) 
    {fLightYieldMean = LightYieldMean;}
  void    SetLightYieldAttenuation(Float_t LightYieldAttenuation)
    {fLightYieldAttenuation = LightYieldAttenuation;}
  void    SetIntrinsicPINEfficiency(Float_t IntrinsicPINEfficiency) 
    {fIntrinsicPINEfficiency = IntrinsicPINEfficiency;}
  void    SetRecalibrationFactor(Float_t RecalibrationFactor) 
    {fRecalibrationFactor = RecalibrationFactor;}
  void    SetElectronsPerGeV(Float_t ElectronsPerGeV) 
    {fElectronsPerGeV = ElectronsPerGeV;}
  void    SetAPDGain(Float_t APDGain) 
    {fAPDGain = APDGain;}
 
 private:
  
  Float_t fLightYieldMean ;         // Mean lightyield in the PbOW4 xtal per GeV (Poisson distribution)
  Float_t fIntrinsicPINEfficiency ; // Photo efficiency of the PIN diode   
  Float_t fLightYieldAttenuation ;  // Attenuation of the light through the crystal
  Float_t fRecalibrationFactor ;    // Recalibration factor
  Float_t fElectronsPerGeV ;        // Number of electrons per GeV created in the PIN by a ionizing particle
  Float_t fAPDGain ;                // APD Gain
  
  ClassDef(AliPHOSv3,1)  // Implementation of PHOS manager class for layout EMC+PPSD with light transport, MIPS in PIN and electronic noise

};

#endif // AliPHOSV3_H
