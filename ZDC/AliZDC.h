#ifndef ALIZDC_H
#define ALIZDC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and classes for set ZDC           //
////////////////////////////////////////////////

class AliZDCMerger;
 
#include "AliDetector.h"
 
class AliZDC : public AliDetector {

public:
  AliZDC();
  AliZDC(const char *name, const char *title);
  virtual       ~AliZDC();
  virtual void  AddHit(Int_t track, Int_t *vol, Float_t *hits);
  virtual void  AddDigit(Int_t *sector, Int_t digit);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials() {}
  Int_t         DistancetoPrimitive(Int_t px, Int_t py);
  virtual Int_t IsVersion() const =0;
  virtual Float_t ZMin() const;	// Minimum overall dimension of the ZDC
  virtual Float_t ZMax() const;	// Maximum overall dimension of the ZDC
  virtual void  SetTreeAddress();
  virtual void  MakeBranch(Option_t* opt, const char *file=0);
  virtual void  MakeBranchInTreeS(TTree *treeS, const char *file=0);
  virtual void  MakeBranchInTreeD(TTree *treeD, const char *file=0);
  virtual void  MakeBranchInTreeR(TTree *treeD, const char *file=0);
  virtual void  Hits2SDigits();
  virtual void  SDigits2Digits();
  virtual void  Hits2Digits();
  virtual void  Digits2Reco();
  TClonesArray  *Reconstructed()   const {return fRecPoints;}
  virtual void  SetMerger(AliZDCMerger* merger) {fMerger = merger;}
  virtual AliZDCMerger* Merger() {return fMerger;}
  virtual void  StepManager() {}
    
  // Switching off the shower development in ZDCs
  void  NoShower(){fNoShower=1;}
  void  Shower()  {fNoShower=0;}

protected:

  Int_t        fNoShower;	// Flag to switch off the shower	
  AliZDCMerger *fMerger;   	// ! pointer to merger
  
  Int_t        fNMergedhits;    // Number of Merged hits for background
  TClonesArray *fMergedHits;    //! TCA for "merged" hits  
  
  Int_t        fNRecPoints;	// Number of RecPoints
  TClonesArray *fRecPoints;	//! List of RecPoints

  ClassDef(AliZDC,3)  	// Zero Degree Calorimeter base class
};
 
#endif
