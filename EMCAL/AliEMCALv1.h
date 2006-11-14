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
  virtual ~AliEMCALv1(void) ;

  using AliEMCALv0::AddHit;
  virtual void  AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t iparent, Float_t ienergy,
			Int_t id, Float_t *hits, Float_t *p);
  // Gives the version number 
  virtual Int_t  IsVersion(void) const {return 1;}
  virtual void StepManager(void) ;
  virtual void RemapTrackHitIDs(Int_t *map);
  virtual void FinishPrimary();
  virtual const TString Version(void)const {return TString("v0");}
  virtual void SetTimeCut(Float_t tc){ fTimeCut = tc;}
  virtual Float_t GetTimeCut() const {return fTimeCut;} 
    
protected:
  Int_t fCurPrimary;  // Current primary track
  Int_t fCurParent;   // Current parent 
  Int_t fCurTrack;    // Current track
  Float_t fTimeCut;   // Cut to remove the background from the ALICE system

 private:
  AliEMCALv1(const AliEMCALv1 & emcal);
  AliEMCALv1 & operator = (const AliEMCALv1  & /*rvalue*/);

  ClassDef(AliEMCALv1,9) // Implementation of EMCAL manager class to produce hits in a Central Calorimeter 
    
};

#endif // AliEMCALV1_H
