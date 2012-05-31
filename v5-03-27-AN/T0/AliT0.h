#ifndef ALIT0_H
#define ALIT0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager and hits classes for set:T0     //
////////////////////////////////////////////////
 
#include <AliDetector.h>
#include <TTree.h>
#include "AliT0RecPoint.h"
#include "AliT0digit.h"
#include "AliT0Trigger.h"
#include "AliT0RawReader.h"

class TDirectory;
class TFile;
class AliESD;


class AliT0 : public AliDetector {



public:
   AliT0();
   AliT0(const char *name, const char *title);
   virtual       ~AliT0();
   virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits);
   virtual void AddDigit(Int_t *, Int_t *) {};
   virtual void   AddDigit(Int_t besttimeright, Int_t besttimeleft, Int_t meantime, 
			   Int_t timediff, Int_t sumMult, Int_t refpoint,
			   TArrayI *time, TArrayI *adc, TArrayI *timeAmp, TArrayI *adcAmp);
   virtual void   CreateGeometry(){}
   virtual void   CreateMaterials(){} 
   virtual Int_t  IsVersion()const {return 0;}
   virtual void   Init();
   virtual void   DefineOpticalProperties(){}
   virtual void SetHitsAddressBranch(TBranch *b1)
     {b1->SetAddress(&fHits);}
   virtual void   MakeBranch(Option_t *opt=" ");
   virtual void   StepManager(){}
   virtual void   ResetHits();
   virtual void   ResetDigits();
    virtual void   SetTreeAddress();
   virtual void   MakeBranchInTreeD(TTree *treeD, const char *file=0);
   // virtual AliLoader* MakeLoader(const char* topfoldername);
   virtual AliDigitizer* CreateDigitizer(AliDigitizationInput* digInput) const;
   void  Digits2Raw ();
   void  Raw2Digits (AliRawReader *reader,TTree* digitsTree);
   virtual AliTriggerDetector* CreateTriggerDetector() const 
     { return new  AliT0Trigger(); }
   
   virtual void  Raw2Digits() {}

protected:
   Int_t           fIdSens;    // Sensetive Cherenkov photocathode
   AliT0digit     *fDigits;    // pointer to T0digits
   AliT0RecPoint  *fRecPoints; // pointer to T0 RecPoints

 private:
   AliT0(const AliT0&);
   AliT0& operator=(const AliT0&);
 
  ClassDef(AliT0,4)  //Base class for the T0 aka T0 detector
};

typedef AliT0 AliSTART; // for backward compatibility

//_____________________________________________________________________________
 
#endif































