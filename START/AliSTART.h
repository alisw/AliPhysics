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
#include "TNamed.h"
#include "TTree.h"
class TDirectory;
R__EXTERN TDirectory *  gDirectory;
 
 
 
class AliSTART : public AliDetector {
 
public:
  AliSTART();
  AliSTART(const char *name, const char *title);
  virtual       ~AliSTART() {}
  virtual void   AddHit(Int_t, Int_t*, Float_t*);
  virtual void   AddDigit( Int_t*, Int_t*);
  virtual void   BuildGeometry();
  virtual void   CreateGeometry() = 0;
  virtual void   CreateMaterials() = 0;
  virtual Int_t  DistanceToPrimitive(Int_t px, Int_t py);
  virtual Int_t  IsVersion() const = 0;
  virtual void   Init();
  virtual void   MakeBranch(Option_t *opt=" ");
  virtual void   DrawModule() = 0;
  virtual void   StepManager() = 0;

  void   Hit2digit(Int_t iEventNum);
  void   Hit2digit(){return;}
public:
  TTree   *fTreeD;        //tree
  TTree * GetTree() { return fTreeD;}
  //return refeence to actual tree 
  Bool_t  SetTree(Int_t nevent=0, TDirectory *dir = gDirectory);
  //map tree from given directory
  Bool_t  MakeTree(Int_t nevent=0);
  //map tree from given directory
protected:
  Int_t fIdSens1;
  ClassDef(AliSTART,1)  //Class for the START detector
};

//____________________________________________________________

#endif
