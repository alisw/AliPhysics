#ifndef AliPHOSRecCpvManager_H
#define AliPHOSRecCpvManager_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Class for the management of the CPV reconstruction.
// Author  : Boris Polichtchouk (IHEP, Protvino)
// 6 March 2001

#include "AliPHOSRecManager.h"
//class AliPHOSGeometry ; 

class AliPHOSRecCpvManager : public AliPHOSRecManager {

 public: 

  AliPHOSRecCpvManager();
  ~AliPHOSRecCpvManager(void);


  void AG(Float_t e, Float_t dx, Float_t dy, Float_t& a, Float_t& gradx, Float_t& grady );
  Float_t Dispersion(Float_t etot, Float_t ai) const;

  Float_t OneGamChi2(Float_t ai, Float_t ei, Float_t etot, Float_t& gi) const ;
  Float_t TwoGamChi2(Float_t ai, Float_t ei, Float_t etot, Float_t& gi) const ;

  Float_t OneGamChisqCut() const { return fOneGamChisqCut; } 
  Float_t OneGamInitialStep() const { return fOneGamInitialStep; }
  Float_t OneGamChisqMin() const { return fOneGamChisqMin; }
  Float_t OneGamStepMin() const { return fOneGamStepMin; }
  Int_t OneGamNumOfIterations() const { return fOneGamNumOfIterations; }

  Float_t TwoGamInitialStep() const { return fTwoGamInitialStep; }
  Float_t TwoGamChisqMin() const { return fTwoGamChisqMin; }
  Float_t TwoGamEmin() const { return fTwoGamEmin; }
  Float_t TwoGamStepMin() const { return fTwoGamStepMin; } 
  Int_t TwoGamNumOfIterations() const { return fTwoGamNumOfIterations; }

  Float_t KillGamMinEnergy() const { return fThr0; } 
  Float_t MergeGammasMinDistanceCut() const { return fSqdCut; } 

  void SetTwoPointsMinDistance(Float_t dist) { fSqdCut=dist; }
  void SetPointMinEnergy(Float_t emin) { fThr0=emin; }

 private:

  Float_t Fcml(Float_t x, Float_t y); // what is it ?
  Float_t GradX(Float_t x, Float_t y); // what is it ?
  Float_t GradY(Float_t x, Float_t y); // what is it ?

  Float_t fOneGamChisqCut; // what is it ?

  Float_t fOneGamInitialStep; // what is it ?
  Float_t fOneGamChisqMin; // what is it ?
  Float_t fOneGamStepMin; // what is it ?
  Int_t fOneGamNumOfIterations; // what is it ?

  Float_t fTwoGamInitialStep; // what is it ?
  Float_t fTwoGamChisqMin; // what is it ?
  Float_t fTwoGamEmin; // what is it ?
  Float_t fTwoGamStepMin; // what is it ?
  Int_t fTwoGamNumOfIterations; // what is it ?

  Float_t fThr0; // what is it ?
  Float_t fSqdCut; // what is it ?

  ClassDef(AliPHOSRecCpvManager,1)        // CPV reconstruction management class 

} ;

#endif // AliPHOSRecCpvManager_H



