#ifndef ALIPHOSV1_H
#define ALIPHOSV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Implementation version v0 of PHOS Manager class 
// Layout EMC + PPSD has name GPS2  
// The main goal of this version of AliPHOS is to calculte the 
//  induced charged in the PIN diode, taking into account light
//  tracking in the PbWO4 crystal, induced signal in the 
//  PIN due to MIPS particle and electronic noise.
// This is done in the StepManager 
//                  
//*-- Author:  Odd Harald Oddland & Gines Martinez (SUBATECH)

// --- ROOT system ---


// --- AliRoot header files ---
#include "AliPHOSv0.h"


class AliPHOSv1 : public AliPHOSv0 {

public:

  AliPHOSv1(void) : AliPHOSv0() {} // ctor 
  AliPHOSv1(const char *name, const char *title="") ;
  AliPHOSv1(AliPHOSReconstructioner * Reconstructioner, const char *name, const char *title="") ;
  virtual ~AliPHOSv1(void) {} 
                            
  virtual void   StepManager(void) ;                                // does the tracking through PHOS and a preliminary digitalization
 
  
private:
  
  Float_t fLightYieldMean ;         // Mean lightyield in the PbOW4 xtal per GeV (Poisson distribution)
  Float_t fIntrinsicPINEfficiency ; // Photo efficiency of the PIN diode   
  Float_t fLightYieldAttenuation ;  // Attenuation of the light through the crystal
  Float_t fRecalibrationFactor ;    // Recalibration factor
  Float_t fElectronsPerGeV ;        // Number of electrons per GeV created in the PIN by a ionizing particle

  ClassDef(AliPHOSv1,1)  // Implementation of PHOS manager class for layout EMC+PPSD with light transport, MIPS in PIN and electronic noise

};

#endif // AliPHOSV1_H
