#ifndef ALIPHOSV1_H
#define ALIPHOSV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager class  for PHOS                   //
//  Version SUBATECH
//  Author : Odd Harald Oddland & 
//           Gines Martinez        Feb-2000 
// The main goal of this version of AliPHOS is to calculted the 
// induced charged in the PIN diode, taking into account light
// tracking in the PbWO4 crystal, induced signal in the 
// PIN due to MIPS particle and electronic noise.
// In this respect, this class derived from AliPHOSv0 and 
// only the StepManager function has been "surcharged"
////////////////////////////////////////////////////////////////////            

// --- ROOT system ---


// --- AliRoot header files ---
#include "AliPHOSv0.h"


class AliPHOSv1 : public AliPHOSv0 {

public:

  AliPHOSv1(void) ;
  AliPHOSv1(const char *name, const char *title="") ;
  AliPHOSv1(AliPHOSReconstructioner * Reconstructioner, const char *name, const char *title="") ;
  virtual ~AliPHOSv1(void) ;
                            
  virtual void   StepManager(void) ;                                // does the tracking through PHOS and a preliminary digitalization
 
  
private:
  //
  // Number of electrons created in the PIN due to light collected in the PbWo4 crystal is calculated using 
  // following formula
  // NumberOfElectrons = EnergyLost * LightYield * PINEfficiency * 
  //                     exp (-LightYieldAttenuation * DistanceToPINdiodeFromTheHit) *
  //                     RecalibrationFactor ;
  // LightYield is obtained as a Poissonian distribution with a mean at 700000 photons per GeV fromValery Antonenko
  // PINEfficiency is 0.1875 from Odd Harald Odland work
  // k_0 is 0.0045 from Valery Antonenko 
  //
  Float_t fLightYieldMean ;   // Mean of the Poisson distribution which is the mean lightyield in the PbOW4 xtal per GeV
  Float_t fIntrinsicPINEfficiency ;    
  Float_t fLightYieldAttenuation ; 
  Float_t fRecalibrationFactor ;
  Float_t fElectronsPerGeV ;   //Number of electrons per GeV created in the PIN by a ionizing particle

  ClassDef(AliPHOSv1,1)  // PHOS v1 main class , version subatech with light transportation, MIPS in PIN and electronic noise

};

#endif // AliPHOSV1_H
