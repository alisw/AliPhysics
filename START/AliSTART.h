#ifndef START_H
#define START_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:START     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"
 
 
class AliSTART : public AliDetector {
 
public:
  AliSTART();
  AliSTART(const char *name, const char *title);
  virtual       ~AliSTART() {}
  virtual void   AddHit(Int_t, Int_t*, Float_t*);
  virtual void   BuildGeometry();
  virtual void   CreateGeometry()=0;
  virtual void   CreateMaterials()=0; 
  virtual Int_t  DistanceToPrimitive(Int_t px, Int_t py);
  virtual Int_t  IsVersion() const =0;
  virtual void   Init();
  virtual void   MakeBranch(Option_t *opt=" ");
  virtual void   DrawModule()=0;
  virtual void   StepManager()=0;
  
 protected:
  Int_t fIdSens1;
  ClassDef(AliSTART,1)  //Class for the START detector
};

//_____________________________________________________________________________
 
class AliSTARThit : public AliHit {
public:
  Int_t      fVolume;
  Int_t      fPmt;
  Int_t      fParticle;     //Particle identificator
  Float_t    fEdep;    //Energy deposition
  Float_t    fEtot;    //Energy of particle 
  Float_t    fTime;    //Time of flight 
 
public:
  AliSTARThit() {}
  AliSTARThit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliSTARThit() {}
  
  ClassDef(AliSTARThit,1)  //Hits for detector START
};

#endif
