#ifndef ALIEMCALJETFINDERALGOUA1UNIT_H
#define ALIEMCALJETFINDERALGOUA1UNIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  *  * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Unit used by UA1 algorithm
//
//*-- Author: Sarah Blyth (LBL/UCT)


#include <TObject.h>
#include "AliEMCALJetFinderTypes.h"

class AliEMCALJetFinderAlgoUA1Unit : public TObject
{
 public:
  AliEMCALJetFinderAlgoUA1Unit();
  ~AliEMCALJetFinderAlgoUA1Unit();
  void SetUnitEnergy(Float_t energy)  {fUnitEnergy = energy;}
  const Float_t GetUnitEnergy()             { return fUnitEnergy;}
  void SetUnitEta(Float_t eta)        {fUnitEta = eta;} 
  const Float_t GetUnitEta()                {return fUnitEta;}
  void SetUnitPhi(Float_t phi)        {fUnitPhi = phi;}
  const Float_t GetUnitPhi()                {return fUnitPhi;}         
  void SetUnitID(Int_t id)            {fUnitID = id;}
  const Int_t GetUnitID()                   {return fUnitID;}
  
  void SetUnitFlag(AliEMCALJetFinderAlgoUA1UnitFlagType_t flag)   
  {
	  fUnitFlag = flag;
  }
  
  const AliEMCALJetFinderAlgoUA1UnitFlagType_t GetUnitFlag()     
  {
	  return fUnitFlag;
  }
  
  const Bool_t operator> ( AliEMCALJetFinderAlgoUA1Unit unit1);
  const Bool_t operator< ( AliEMCALJetFinderAlgoUA1Unit unit1);
  const Bool_t operator== ( AliEMCALJetFinderAlgoUA1Unit unit1);

 protected:
  Float_t         fUnitEnergy;        // Energy of the unit 
  Float_t         fUnitEta;           // Eta of the unit
  Float_t         fUnitPhi;           // Phi of the unit
  Int_t           fUnitID;            // ID of the unit
  AliEMCALJetFinderAlgoUA1UnitFlagType_t      fUnitFlag; //Flag of the unit

  ClassDef(AliEMCALJetFinderAlgoUA1Unit,2)
};

#endif

