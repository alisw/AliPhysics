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

  // 30-aug-04
  virtual void Browse(TBrowser* b);
  // drawing
  void DrawCalorimeterCut(const char *name="SMOD", int axis=3, double dcut=1.); // *MENU*
  void DrawSuperModuleCut(const char *name="EMOD", int axis=2, double dcut=0.03, int fill = 6);//  *MENU*
  void DrawTowerCut(const char *name="SCMY", int axis=2, double dcut=0., int fill=1, const char *optShad="on");   //  *MENU*
  void DrawAlicWithHits(int mode=1);                            // *MENU*
  void SetVolumeAttributes(const char *name="SCM0", int seen=1, int color=1, int fill=1); // *MENU*
  void TestIndexTransition(int pri=0, int idmax=0); // *MENU*

 protected:

 private:
  AliEMCALv2(const AliEMCALv2 & emcal);
  AliEMCALv2 & operator = (const AliEMCALv2  & /*rvalue*/);
 
  ClassDef(AliEMCALv2,2)    //Implementation of EMCAL manager class to produce hits in a Shish-Kebab
    
};

#endif // AliEMCALV2_H
