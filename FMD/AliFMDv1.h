#ifndef ALIFMDV1_H
#define ALIFMDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////
//  Manager and hits classes for set:FMD     //
////////////////////////////////////////////////
 
#include "AliFMD.h"
#include "AliFMDSDigitizer.h"
 
class AliFMDv1 : public AliFMD {
  
public:
  AliFMDv1() {};
  AliFMDv1(const char *name, const char *title);
  virtual       ~AliFMDv1() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawDetector();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();

protected:
   Int_t fIdSens1; // Sensetive volume  in FMD
   
// Background event for event mixing
   Text_t *fFileName;           // ! File with background hits
   TTree *fTrH1;                // Hits Tree for background event
   TClonesArray *fHits2;        // List of hits for one track only
  
   ClassDef(AliFMDv1,2)  //Class for FMD version 0
};

#endif


