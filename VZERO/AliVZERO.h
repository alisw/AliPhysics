#ifndef ALIVZERO_H
#define ALIVZERO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//////////////////////////////////////////////////
//  Manager and hits classes for set : VZERO    //
//////////////////////////////////////////////////

/*
#include "AliRun.h"
#include "AliMC.h"
#include "AliDetector.h"
#include "AliVZEROLoader.h"
 
#include <TNamed.h>
#include <TTree.h>
*/
#include "AliDetector.h"
#include "AliVZEROTrigger.h"

class TNamed;
class TTree;

class AliVZEROLoader;
class AliVZEROhit; 
class AliVZEROdigit;
  
class AliVZERO : public AliDetector {
 
public:

  AliVZERO() {}
  AliVZERO(const char *name, const char *title);
  virtual       ~AliVZERO();
//  virtual void   AddHit(Int_t track , Int_t *vol, Float_t *hits)=0; 
//  virtual void   AddDigits(Int_t* tracks, Int_t* digits)=0;
  virtual void   BuildGeometry();
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   AddAlignableVolumes() const {}
  virtual Int_t  DistanceToPrimitive(Int_t px, Int_t py);
  virtual Int_t  IsVersion() const = 0;
  virtual void   Init();
  virtual AliLoader* MakeLoader(const char* topfoldername);
  virtual void   Hits2Digits();
  virtual void   Digits2Raw();
  virtual void   SetTreeAddress();  
  virtual void   MakeBranch(Option_t *option) =0;
  virtual void   DrawModule() const {};
  virtual void   StepManager() {};
// Trigger 
  virtual AliTriggerDetector* CreateTriggerDetector() const 
  { return new AliVZEROTrigger(); }
  
  virtual void   SetThickness(Float_t thick)  {fThickness = thick;};
  virtual void   SetThickness1(Float_t thick) {fThickness1 = thick;};
// Set Stepping Parameters
  virtual void   SetMaxStepQua(Float_t p1);
  virtual void   SetMaxStepAlu(Float_t p1);
  virtual void   SetMaxDestepQua(Float_t p1);
  virtual void   SetMaxDestepAlu(Float_t p1);

  AliDigitizer*  CreateDigitizer(AliRunDigitizer* manager) const;

protected:

   Int_t   fIdSens1;      // Sensitive volume  in VZERO
   Float_t fThickness;    // Total thickness of box holding Right detector V0R i.e. 4.1 cm
   Float_t fThickness1;   // Thickness of elementary cells i.e. 0.7 cm
  
// Stepping Parameters
   Float_t fMaxStepQua;   // Maximum step size inside the quartz volumes
   Float_t fMaxStepAlu;   // Maximum step size inside the  aluminum volumes
   Float_t fMaxDestepQua; // Maximum relative energy loss in quartz
   Float_t fMaxDestepAlu; // Maximum relative energy loss in aluminum
  
  ClassDef(AliVZERO,1)  //Class for the VZERO detector
};

//____________________________________________________________

#endif
