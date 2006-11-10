#ifndef STARTV1_H
#define STARTV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager and hits classes for set:START     //
////////////////////////////////////////////////
 
#include "AliSTART.h"
 
class AliSTARTv1 : public AliSTART {
  
public:

  enum constants {kAir=1, kSc=2, kVac=3, kCer=4, kGlass=6, kSteel=8, kRibber=9, kBrass=11, kLucite=12, kC=13, kPP=14, kAl=15, kOpGlass=16, kOpAir=17, kOpAirNext=18, kOpGlassCathode=19};

 
  AliSTARTv1() {};
  AliSTARTv1(const char *name, const char *title);
  virtual       ~AliSTARTv1();
  virtual void   CreateGeometry();
  virtual void   AddAlignableVolumes() const;
  virtual void   CreateMaterials();
  virtual void   DrawDetector();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  Bool_t RegisterPhotoE(Int_t impt, Double_t energy);
  virtual void   StepManager();


protected:
  Int_t fIdSens1; // Sensetive volume  in START
  TObjArray fEffPMT; //pmt registration effeicincy
 
  ClassDef(AliSTARTv1,2)  //Class for START version 1
};

#endif


