#ifndef AliPHOSRecCpvManager_H
#define AliPHOSRecCpvManager_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Class for the management of the CPV reconstruction.
// Author  : Boris Polichtchouk (IHEP, Protvino)
// 6 March 2001

#include "AliPHOSRecManager.h"
#include "AliPHOSGeometry.h"

class AliPHOSRecCpvManager : public AliPHOSRecManager {

 public: 

  AliPHOSRecCpvManager();
  ~AliPHOSRecCpvManager(void);


  void AG(Float_t E, Float_t dx, Float_t dy, Float_t& A, Float_t& grad_x, Float_t& grad_y );
  Float_t Dispersion(Float_t Etot, Float_t Ai, Float_t Ei);

  Float_t OneGamChi2(Float_t Ai, Float_t Ei, Float_t Etot, Float_t& Gi);
  Float_t TwoGamChi2(Float_t Ai, Float_t Ei, Float_t Etot, Float_t& Gi);

  Float_t OneGamChisqCut() { return fOneGamChisqCut; }
  Float_t OneGamInitialStep() { return fOneGamInitialStep; }
  Float_t OneGamChisqMin() { return fOneGamChisqMin; }
  Float_t OneGamStepMin() { return fOneGamStepMin; }
  Int_t OneGamNumOfIterations() { return fOneGamNumOfIterations; }

  Float_t TwoGamInitialStep() { return fTwoGamInitialStep; }
  Float_t TwoGamChisqMin() { return fTwoGamChisqMin; }
  Float_t TwoGamEmin() { return fTwoGamEmin; }
  Float_t TwoGamStepMin() { return fTwoGamStepMin; } 
  Int_t TwoGamNumOfIterations() { return fTwoGamNumOfIterations; }

  Float_t KillGamMinEnergy() { return fThr0; } 
  Float_t MergeGammasMinDistanceCut() { return fSqdCut; } 

  void SetTwoPointsMinDistance(Float_t dist) { fSqdCut=dist; }
  void SetPointMinEnergy(Float_t emin) { fThr0=emin; }

 private:

  Float_t Fcml(Float_t x, Float_t y);
  Float_t GradX(Float_t x, Float_t y);
  Float_t GradY(Float_t x, Float_t y);

 private:

  Float_t fOneGamChisqCut;

  Float_t fOneGamInitialStep;
  Float_t fOneGamChisqMin;
  Float_t fOneGamStepMin;
  Int_t fOneGamNumOfIterations;

  Float_t fTwoGamInitialStep;
  Float_t fTwoGamChisqMin;
  Float_t fTwoGamEmin;
  Float_t fTwoGamStepMin;
  Int_t fTwoGamNumOfIterations;

  Float_t fThr0;
  Float_t fSqdCut;

 public:

  ClassDef(AliPHOSRecCpvManager,1)        // CPV reconstruction management class 

} ;

#endif // AliPHOSRecCpvManager_H



