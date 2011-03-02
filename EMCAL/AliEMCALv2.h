#ifndef ALIEMCALV2_H
#define ALIEMCALV2_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                      
         */
/* $Id$ */

//_________________________________________________________________________
// Implementation version v2 of EMCAL Manager class for Shish-Kebab case 
//*--                  
//*-- Author:  Aleksei Pavlinov

// --- ROOT system ---
class TBrowser;

// --- AliRoot header files ---
#include "AliEMCALv1.h"

class AliEMCALv2 : public AliEMCALv1 {
  
public:

  AliEMCALv2(void) ; 
  AliEMCALv2(const char *name, const char *title="") ;
  virtual ~AliEMCALv2(void) ;

  using AliEMCALv1::AddHit;
  virtual void  AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t iparent, Float_t ienergy,
			Int_t id, Float_t *hits, Float_t *p);

  virtual void StepManager(void) ;

  // Gives the version number 
  virtual Int_t  IsVersion(void) const {return 2;}
  virtual const TString Version(void)const {return TString("v2");}

 protected:

 private:
  AliEMCALv2(const AliEMCALv2 & emcal);
  AliEMCALv2 & operator = (const AliEMCALv2  & /*rvalue*/);
 
  ClassDef(AliEMCALv2,3)    //Implementation of EMCAL manager class to produce hits in a Shish-Kebab
    
};

#endif // AliEMCALV2_H
