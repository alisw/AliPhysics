#ifndef ALIMUONV1_H
#define ALIMUONV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
 
#include "AliMUON.h"
#include <TF1.h>

class AliMUONv1 : public AliMUON {
public:
   AliMUONv1();
   AliMUONv1(const char *name, const char *title);
   virtual  ~AliMUONv1() {}
   virtual void   CreateGeometry();
   virtual void   CreateMaterials();
   virtual void   Init();
   virtual Int_t  IsVersion() const {return 1;}
   virtual void   StepManager();
   void StepManagerOld();
   void StepManagerNew();
   void StepManagerTest();



   void SetStepManagerVersionOld(Bool_t Opt) 
     { fStepManagerVersionOld = Opt; }
   void SetStepManagerVersionNew(Bool_t Opt) 
     { fStepManagerVersionNew = Opt; }
   void SetStepManagerVersionTest(Bool_t Opt) 
     { fStepManagerVersionTest = Opt; }
   void SetStepMaxInActiveGas(Float_t StepMax)
     {fStepMaxInActiveGas = StepMax; }
protected:
   Int_t*  fStations; //! allow to externally set which station to create
   Bool_t  fStepManagerVersionOld; // Version of StepManager, Default is false
   Bool_t  fStepManagerVersionNew; // Version of StepManager, Default is false
   Bool_t  fStepManagerVersionTest; // Version of StepManager, Default is false
   Float_t fStepMaxInActiveGas; // Step mas in active gas default 0.6cm
   virtual Int_t  GetChamberId(Int_t volId) const;
   

private:
   ClassDef(AliMUONv1,1)  // MUON Detector class Version 1


};
#endif







