//
// Class AliRsnCutRange
//
// General implementation of cuts which check a value inside a range.
// This range can be defined by two integers or two doubles.
// A user-friendly enumeration allows to define what is checked.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTBETHEBLOCK_H
#define ALIRSNCUTBETHEBLOCK_H

#include "AliPID.h"
#include "AliRsnCut.h"

class AliRsnCutBetheBloch : public AliRsnCut
{
 public:

  AliRsnCutBetheBloch();
  AliRsnCutBetheBloch(const char *name, Double_t fractionRange, Double_t mass, Double_t mip = 50.0, Bool_t correct = kTRUE);
  AliRsnCutBetheBloch(const char *name, Double_t fractionRange, AliPID::EParticleType type, Double_t mip = 50.0, Bool_t correct = kTRUE);

  void           SetMass(Double_t mass) {fMass = mass;}
  void           SetMass(AliPID::EParticleType type) {AliPID pid; fMass = pid.ParticleMass(type);}
  void           SetMIP(Double_t mip) {fMIP = mip;}
  void           SetCalibConstant(Int_t i, Double_t value) {if (i>=0&&i<5) fConst[i] = value;}
  Double_t       BetheBloch(AliRsnDaughter *track);

  virtual Bool_t IsSelected(ETarget tgt, AliRsnDaughter *daughter);
  virtual Bool_t IsSelected(ETarget tgt, AliRsnPairParticle *pair);
  virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *event);
  virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *ev1, AliRsnEvent *ev2);

protected:

  Bool_t   fCorrect;   // apply or not the saturation corrections
  Double_t fMass;      // mass hypothesis
  Double_t fMIP;       // MIP normalization
  Double_t fConst[5];  // calibration constants

  ClassDef(AliRsnCutBetheBloch, 1)
};

#endif
