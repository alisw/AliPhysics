#ifndef AliPHOSRecManager_H
#define AliPHOSRecManager_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Base class for the management by the PHOS reconstruction.
// It contains only virtual member functions
// which will be implemented for the Emc and CPV reconstruction
// in the appropriate derived classes.
// Author  : Boris Polichtchouk (IHEP, Protvino)
// 6 March 2001

#include "TNamed.h"

class AliPHOSRecManager : public TNamed {

 public:

  AliPHOSRecManager();
  virtual ~AliPHOSRecManager(void) {}

  virtual void AG(Float_t E, Float_t dx, Float_t dy, Float_t& A, Float_t& grad_x, Float_t& grad_y ) = 0;
  virtual Float_t Dispersion(Float_t Etot, Float_t Ai, Float_t Ei) = 0;

  virtual Float_t OneGamChi2(Float_t Ai, Float_t Ei, Float_t Etot, Float_t& Gi) = 0;
  virtual Float_t TwoGamChi2(Float_t Ai, Float_t Ei, Float_t Etot, Float_t& Gi) = 0;

  virtual Float_t OneGamChisqCut() = 0 ;
  virtual Float_t OneGamInitialStep() = 0;
  virtual Float_t OneGamChisqMin() = 0;
  virtual Float_t OneGamStepMin() = 0;
  virtual Int_t OneGamNumOfIterations() = 0;

  virtual Float_t TwoGamInitialStep() = 0;
  virtual Float_t TwoGamChisqMin() = 0;
  virtual Float_t TwoGamEmin() = 0;
  virtual Float_t TwoGamStepMin() = 0;
  virtual Int_t TwoGamNumOfIterations() = 0;

  virtual Float_t KillGamMinEnergy() = 0;
  virtual Float_t MergeGammasMinDistanceCut() = 0;

  virtual void SetTwoPointsMinDistance(Float_t dist) = 0;
  virtual void SetPointMinEnergy(Float_t emin) = 0;

  ClassDef(AliPHOSRecManager,1)

} ;

#endif // AliPHOSRecManager_H













