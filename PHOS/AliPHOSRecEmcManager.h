#ifndef AliPHOSRecEmcManager_H
#define AliPHOSRecEmcManager_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Class for the management by the Emc reconstruction.
// Author  : Boris Polichtchouk (IHEP, Protvino)
// 6 March 2001

#include "AliPHOSRecManager.h"
#include "AliPHOSGeometry.h"

class AliPHOSRecEmcManager : public AliPHOSRecManager {

 public: 

  AliPHOSRecEmcManager();
  ~AliPHOSRecEmcManager(void);


  void AG(Float_t E, Float_t dx, Float_t dy, Float_t& A, Float_t& grad_x, Float_t& grad_y );
  Float_t Dispersion(Float_t Etot, Float_t Ai, Float_t Ei) const;

  Float_t OneGamChi2(Float_t Ai, Float_t Ei, Float_t Etot, Float_t& Gi);
  Float_t TwoGamChi2(Float_t Ai, Float_t Ei, Float_t Etot, Float_t& Gi)const;

  Float_t OneGamChisqCut() const{ return fOneGamChisqCut; }
  Float_t OneGamInitialStep() const{ return fOneGamInitialStep; }
  Float_t OneGamChisqMin() const{ return fOneGamChisqMin; }
  Float_t OneGamStepMin() const{ return fOneGamStepMin; }
  Int_t OneGamNumOfIterations() const{ return fOneGamNumOfIterations; }

  Float_t TwoGamInitialStep() const{ return fTwoGamInitialStep; }
  Float_t TwoGamChisqMin() const{ return fTwoGamChisqMin; }
  Float_t TwoGamEmin() const{ return fTwoGamEmin; }
  Float_t TwoGamStepMin() const{ return fTwoGamStepMin; } 
  Int_t TwoGamNumOfIterations() const{ return fTwoGamNumOfIterations; }

  Float_t KillGamMinEnergy() const{ return fThr0; } 
  Float_t MergeGammasMinDistanceCut() const{ return fSqdCut; } 

  void SetTwoPointsMinDistance(Float_t dist) { fSqdCut=dist; }
  void SetPointMinEnergy(Float_t emin) { fThr0=emin; }

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

  ClassDef(AliPHOSRecEmcManager,1)        // Emc reconstruction management class 

} ;

#endif // AliPHOSRecEmcManager_H



