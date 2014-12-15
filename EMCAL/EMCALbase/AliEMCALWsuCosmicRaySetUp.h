#ifndef ALIBODY_H
#define ALIBODY_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////
//  Manager class for detector:  AliEMCALWsuCosmicRaySetUp        //
//   This is the envelop for Alice                                //
///////////////////////////////////////////////////////////////////
 
#include "AliModule.h"

class TList;

class AliEMCALWsuCosmicRaySetUp : public AliModule {
 
public:
  AliEMCALWsuCosmicRaySetUp();
  AliEMCALWsuCosmicRaySetUp(const char *name, const char *title);
  virtual     ~AliEMCALWsuCosmicRaySetUp() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  void DefineCuts(const Int_t idtmed=1);
  virtual Int_t IsVersion() const {return 0;}
  // GetMethod
  Float_t* GetMasterVolume() {return fMasterVolume;}
  TList*   GetLhists() {return fLHists;}
  TList*   GetLhists(Int_t ind) {return ind<0?fLHists:dynamic_cast<TList *>(fLHists->At(ind));}
  // Dec 1,2010
  virtual void StepManager(void) ;
  virtual void FinishEvent();

  TList*  BookKineHists(const Double_t p=1., const Char_t *opt="kine");
  //
  virtual Bool_t  IsFolder() const {return kTRUE;}
  virtual void Browse(TBrowser* b);

  protected:
  TList *fLHists;           // list of hists
  Float_t fMasterVolume[3]; // size of MASTER volume
private:
  // Keep for convention only
  AliEMCALWsuCosmicRaySetUp(const AliEMCALWsuCosmicRaySetUp &var);
  AliEMCALWsuCosmicRaySetUp& operator = (const AliEMCALWsuCosmicRaySetUp & /*rvalue*/);

  ClassDef(AliEMCALWsuCosmicRaySetUp,1)  // Class manager for the Wsu Cosmic Ray or TB CERN SetUp
};

#endif
