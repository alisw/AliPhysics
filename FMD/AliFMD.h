#ifndef ALIFMD_H
#define ALIFMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager and hits classes for set:Si-FMD     //
////////////////////////////////////////////////
 
#include <AliDetector.h>
#include <TBranch.h>
#include <AliLoader.h>

class TClonesArray;
 class AliFMD : public AliDetector {
 
public:
  AliFMD();
  AliFMD(const char *name, const char *title);
  virtual       ~AliFMD(); 
  virtual void   AddHit(Int_t, Int_t*, Float_t*);
  virtual void   AddDigit(Int_t*);
   virtual void   BuildGeometry();
  virtual void   CreateGeometry() {}
  virtual void   CreateMaterials()=0; 
  virtual Int_t  DistanceToPrimitive(Int_t px, Int_t py);
  virtual Int_t  IsVersion() const =0;
  virtual void   Init();
  virtual void   MakeBranch(Option_t *opt=" ");
  virtual void   MakeBranchInTreeD(TTree *treeD, const char *file=0);
  virtual void   SetTreeAddress();
  virtual void   ResetHits();
  virtual void   ResetDigits();
  virtual void   DrawDetector()=0;
  virtual void   StepManager() {}
   
  void SetEventNumber(Int_t i)     {fEvNrSig = i;}
  void  Eta2Radius(Float_t, Float_t, Float_t*);

 
   // Digitisation
  TClonesArray *ReconParticles() const {return fReconParticles;}   
  virtual void SetHitsAddressBranch(TBranch *b){b->SetAddress(&fHits);}
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;

 protected:
//Granularity
  Int_t fRingsSi1;       // Number of rings
  Int_t fSectorsSi1;    // Number of sectors
  Int_t fRingsSi2;       // Number of rings
  Int_t fSectorsSi2;    // Number of sectors

  Int_t   fNevents ;        // Number of events to digitize
  Int_t fEvNrSig;                 // signal     event number


  TClonesArray *fReconParticles;

 ClassDef(AliFMD,6)  //Class for the FMD detector
};
#endif // AliFMD_H


