#ifndef ALISTART_H
#define ALISTART_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager and hits classes for set:START     //
////////////////////////////////////////////////
 
#include <AliDetector.h>
#include <TTree.h>

class TDirectory;
class TFile;
class AliESD;
R__EXTERN TDirectory *  gDirectory;
 
 
 
class AliSTART : public AliDetector {



public:
   AliSTART();
   AliSTART(const char *name, const char *title);
   virtual       ~AliSTART();
   virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits);
   virtual void   AddDigit(Int_t *tracks, Int_t *digits);
   virtual void   BuildGeometry();
   virtual void   CreateGeometry(){}
   virtual void   CreateMaterials(){} 
   virtual Int_t  DistanceToPrimitive(Int_t px, Int_t py);
   virtual void   DrawDetector(){}
   virtual Int_t  IsVersion()const {return 0;}
   virtual void   Init();
   virtual void SetHitsAddressBranch(TBranch *b1)
     {b1->SetAddress(&fHits);}
   virtual void   MakeBranch(Option_t *opt=" ");
   virtual void   StepManager(){}
   virtual void   ResetHits();
   virtual void   SetTreeAddress();
   virtual AliLoader* MakeLoader(const char* topfoldername);
   virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
   void  Digits2Raw ();

protected:
   Int_t fIdSens;    // Sensetive Cherenkov radiator
 
private:
  ClassDef(AliSTART,4)  //Base class for the T0 aka START detector
};

//_____________________________________________________________________________
 
#endif































