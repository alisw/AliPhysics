#ifndef ALIMUONV3_H
#define ALIMUONV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 3    //
/////////////////////////////////////////////////////////
//
// Old MUONv1 class (to be removed later)
// - now replaced with a new one where geometry and materials
// are created using new geometry builders
// (See ALIMUON*GeometryBuilder classes)
 
#include <TLorentzVector.h>

#include "AliMUON.h"

class TF1;

class AliMUONv3 : public AliMUON 
{
 public:
   AliMUONv3();
   AliMUONv3(const char *name, const char *title);
   virtual  ~AliMUONv3() {}
   virtual void   CreateGeometry();
   virtual void   CreateMaterials();
   virtual void   Init();
   virtual Int_t  IsVersion() const {return 3;}
   virtual void   StepManager();
   void StepManagerOld();
   void SetStepManagerVersionOld(Bool_t Opt) 
     { fStepManagerVersionOld = Opt; }
   void SetAngleEffect(Bool_t Opt) 
     { fAngleEffect = Opt; }
   void SetStepMaxInActiveGas(Float_t StepMax)
     {fStepMaxInActiveGas = StepMax; }

 protected:
   AliMUONv3(const AliMUONv3& right);
   AliMUONv3&  operator = (const AliMUONv3& right);
   Int_t*  fStations;              //! allow to externally set which station to create
   Bool_t  fStepManagerVersionOld; // Version of StepManager, Default is false
   Bool_t  fAngleEffect; // Angle Effect along wires, Default is true
   Float_t fStepMaxInActiveGas;    // Step max in active gas default 0.6cm
   virtual Int_t  GetChamberId(Int_t volId) const;

   // StepManager 
   Float_t *  fStepSum; //!
   Float_t *  fDestepSum; //!
  
   TLorentzVector fTrackMomentum; // Momentum of the particle entering in the active gas of chamber
   TLorentzVector fTrackPosition; // Position of the particle exiting the active gas of chamber
   TF1 *          fElossRatio;    // Ratio of particle mean eloss with respect MIP's 
   TF1 *          fAngleEffect10; // Angle effect in tracking chambers at theta =10 degres as a function of ElossRatio (Khalil BOUDJEMLINE sep 2003 Ph.D Thesis) (in micrometers)  
   TF1 *          fAngleEffectNorma;// Angle effect: Normalisation form theta=10 degres to theta between 0 and 10 (Khalil BOUDJEMLINE sep 2003 Ph.D Thesis)

   ClassDef(AliMUONv3,1)  // MUON Detector class Version 1
};
#endif







