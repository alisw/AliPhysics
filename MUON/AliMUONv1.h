#ifndef ALIMUONV1_H
#define ALIMUONV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 1    //
/////////////////////////////////////////////////////////
 
#include "TLorentzVector.h"

#include "AliMUON.h"

class TF1;
class TGeoCombiTrans;

class TString;
class TGeoHMatrix;

class AliMUONv1 : public AliMUON {
public:
   AliMUONv1();
   AliMUONv1(const char *name, const char *title);
   virtual  ~AliMUONv1();
   virtual void   CreateGeometry();
   virtual void   CreateMaterials();
   virtual void   Init();
   virtual Int_t  IsVersion() const {return 1;}
   virtual void   StepManager();
   void StepManagerOld();
   void SetStepManagerVersionOld(Bool_t Opt) 
     { fStepManagerVersionOld = Opt; }
   void SetAngleEffect(Bool_t Opt) 
     { fAngleEffect = Opt; }
   void SetStepMaxInActiveGas(Float_t StepMax)
     {fStepMaxInActiveGas = StepMax; }
protected:
   AliMUONv1(const AliMUONv1& right);
   AliMUONv1&  operator = (const AliMUONv1& right);

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
   TGeoCombiTrans* fGlobalTransformation; // global transformation 
                                  // applied to the whole geometry 
private:
   // method
   void PlaceVolume(const TString& name, const TString& mName, Int_t copyNo, 
             const TGeoHMatrix& matrix, Int_t npar, Double_t* param,
	     const char* only) const;

   ClassDef(AliMUONv1,2)  // MUON Detector class Version 1


};
#endif







