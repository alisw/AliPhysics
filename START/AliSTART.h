#ifndef START_H
#define START_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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
  Int_t fZposit;


public:
   AliSTART();
   AliSTART(const char *name, const char *title);
   virtual       ~AliSTART() {}
   virtual void   AddHit(Int_t, Int_t*, Float_t*);
   virtual void   AddDigit(Int_t*, Int_t*);
   virtual void   BuildGeometry();
   virtual void   CreateGeometry(){}
   virtual void   CreateMaterials(){} 
   virtual Int_t  DistanceToPrimitive(Int_t px, Int_t py);
   virtual void   DrawDetector(){}
   virtual Int_t  IsVersion()const {return 0;}
   virtual void   Init();
   void Hit2digit(Int_t iEventNum);
   void Hit2digit(){return;}
   virtual void   MakeBranch(Option_t *opt=" ");
   virtual void   StepManager(){}
   void PrintMedium(Int_t iMediumId=0);//Prints "iMediumId" TMED properties
   /*
   TTree   *fTreeD;        //tree
   TTree * GetTree() { return fTreeD;}//return reference to actual tree 
   Bool_t  SetTree(Int_t nevent=0, TDirectory *dir = gDirectory);//map tree from given directory
   Bool_t  MakeTree(Int_t nevent=0);//map tree from given directory
   */   
   
   //  void Fill();
   //   void Write();
   
protected:
   Int_t fIdSens;    // Sensetive Cherenkov radiator
  ClassDef(AliSTART,1)  //Base class for the T0 aka START detector
};

//_____________________________________________________________________________
 
#endif































