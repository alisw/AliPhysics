#ifndef ALIEMCALV1_H
#define ALIEMCALV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                      
         */
/* $Id$ */

//_________________________________________________________________________
// Implementation version v1 of EMCAL Manager class 
//*--                  
//*-- Author: Sahal Yacoob (LBL / UCT) 
//*--  and  : Jennifer Klay (LBL)
//#include <assert.h>

// --- ROOT system ---
class TClonesArray;
class TLorentzVector;
class TFile;

// --- AliRoot header files ---
#include "AliEMCALv0.h"

class AliEMCALv1 : public AliEMCALv0 {
  
public:

  AliEMCALv1(void) ; 
  AliEMCALv1(const char *name, const char *title="") ;
  // cpy ctor: no implementation yet
  // requested by the Coding Convention
  AliEMCALv1(const AliEMCALv0 & emcal) {abort();}
  virtual ~AliEMCALv1(void) ;
  virtual void  AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t iparent, Float_t ienergy,
			Int_t id, Float_t *hits, Float_t *p);
  // Gives the version number 
  virtual Int_t  IsVersion(void) const {return 1;}
  virtual void StepManager(void) ;
  virtual const TString Version(void)const {return TString("v0");}
  // assignement operator requested by coding convention but not needed  
  AliEMCALv1 & operator = (const AliEMCALv0 & rvalue){abort();return *this;}
 
    
private:

  Float_t fLightYieldMean ;         // Mean lightyield in a plastic layer per GeV (Poisson distribution)
  Float_t fIntrinsicAPDEfficiency ; // Photo efficiency of the APD diode   
  Float_t fLightYieldAttenuation ;  // Attenuation of the light through the light fiber
  Float_t fRecalibrationFactor ;    // Recalibration factor
  Float_t fAPDGain ;                // APD Gain
  Float_t fLightFactor ;            //! a calculated factor
  Float_t fAPDFactor ;              //! a calculated factor

  ClassDef(AliEMCALv1,4)//Implementation of EMCAL manager class to produce hits in a Central Calorimeter 

};

#endif // AliEMCALV1_H
