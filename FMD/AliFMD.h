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
  AliFMD(const AliFMD& FMD):AliDetector(FMD) {;}  //copy ctor
  virtual       ~AliFMD(); 
 AliFMD&  operator=(const AliFMD&)                 {return *this;}
  virtual void   AddHit(Int_t track, Int_t * vol, Float_t * hits);
  virtual void   AddDigit(Int_t* digits);
   virtual void   BuildGeometry();
  virtual void   CreateGeometry() {}
  virtual void   CreateMaterials()=0; 
  virtual const Int_t  DistanceToPrimitive(Int_t px, Int_t py);
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


 
   // Digitisation
  TClonesArray *ReconParticles() const {return fReconParticles;}   
  virtual const void SetHitsAddressBranch(TBranch *b){b->SetAddress(&fHits);}
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;

 protected:
//Granularity
  Int_t fRingsSi1;       // Number of rings
  Int_t fSectorsSi1;    // Number of sectors
  Int_t fRingsSi2;       // Number of rings
  Int_t fSectorsSi2;    // Number of sectors

  Int_t   fNevents ;        // Number of events to digitize
  Int_t fEvNrSig;                 // signal     event number


  TClonesArray *fReconParticles; //array of reconstructed multiplicity in 0.1eta

 ClassDef(AliFMD,6)  //Class for the FMD detector
};
#endif // AliFMD_H


