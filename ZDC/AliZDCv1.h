#ifndef ALIZDCV1_H
#define ALIZDCV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:ZDC      //
////////////////////////////////////////////////

#include "AliZDC.h"

//____________________________________________________________________________ 
class AliZDCv1 : public AliZDC {

public:
  AliZDCv1();
  AliZDCv1(const char *name, const char *title);
  virtual      ~AliZDCv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateBeamLine();
  virtual void  CreateZDC();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  DrawModule();
  virtual void  Init();
  virtual void  InitTables();
  virtual void  StepManager();
 
protected:
  //Sensitive media
  Int_t   fMedSensF1;  // Sensitive medium F1
  Int_t   fMedSensF2;  // Sensitive medium F2
  Int_t   fMedSensZP;  // Sensitive medium for Protons
  Int_t   fMedSensZN;  // Sensitive medium for Neutrons
  Int_t   fMedSensGR;  // Other sensitive medium
  //Parameter for light tables
  Int_t   fNalfan;             // Number of Alfa neutrons
  Int_t   fNalfap;             // Number of Alfa protons
  Int_t   fNben;               // Number of beta neutrons
  Int_t   fNbep;               // Number of beta protons
  Float_t fTablen[4][90][18];  // Table neutrons
  Float_t fTablep[4][90][28];  // Table protons

   ClassDef(AliZDCv1,1)  // Zero Degree Calorimeter version 1
}; 
 
#endif
