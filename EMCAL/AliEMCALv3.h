#ifndef ALIEMCALV3_H
#define ALIEMCALV3_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                      
         */
/* $Id$ */

//_________________________________________________________________________
// Implementation version v3 of EMCAL Manager class for Shish-Kebab case 
// Save all hits inside of Sc - Nov 25, 05
//*--                  
//*-- Author:  Aleksei Pavlinov

class AliEMCALGeometry;

// --- AliRoot header files ---
#include "AliEMCALv1.h"

class AliEMCALv3 : public AliEMCALv1 {
  
public:

  AliEMCALv3(void) ; 
  AliEMCALv3(const char *name, const char *title="") ;
  AliEMCALv3(const AliEMCALv3 & emcal); //copy ctor

  virtual ~AliEMCALv3(void) ;

  using AliEMCALv1::AddHit;
  virtual void  AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t iparent, Float_t ienergy,
			Int_t id, Float_t *hits, Float_t *p);

  virtual void StepManager(void) ;

  // Gives the version number 
  virtual Int_t  IsVersion(void) const {return 3;}
  virtual const TString Version(void)const {return TString("v3");}
  AliEMCALv3 & operator = (const AliEMCALv3 & /*rvalue*/){
    Fatal("operator =", "not implemented") ;  
    return *this;}

  virtual void Browse(TBrowser* b);

  AliEMCALGeometry* fGeometry; //!

  ClassDef(AliEMCALv3,1)    //Implementation of EMCAL manager class to produce hits in a Shish-Kebab
    
};

#endif // AliEMCALV3_H
