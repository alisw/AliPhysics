#ifndef ALIANALYSISTASKPARTICLERANDOMIZER_H
#define ALIANALYSISTASKPARTICLERANDOMIZER_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TClonesArray;
class TString;
class TRandom3;

#include "AliAnalysisTaskSE.h"

// Analysis task to copy (and randomize) a given TClonesArray of AliVParticle-derived objects

class AliAnalysisTaskParticleRandomizer : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskParticleRandomizer();
  virtual ~AliAnalysisTaskParticleRandomizer();

  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *);
  virtual void  Terminate(Option_t *) {}
  void          ExecOnce();

  // ###### Configuration setters
  void          SetRandomizeInPhi(Bool_t val)               {fRandomizeInPhi = val;}
  void          SetRandomizeInEta(Bool_t val)               {fRandomizeInEta = val;}
  void          SetRandomizeInPt(Bool_t val)                {fRandomizeInPt = val;}

  void          SetPhiMin(Double_t val)                     {fMinPhi = val;}
  void          SetPhiMax(Double_t val)                     {fMaxPhi = val;}
  void          SetEtaMin(Double_t val)                     {fMinEta = val;}
  void          SetEtaMax(Double_t val)                     {fMaxEta = val;}
  void          SetPtMin(Double_t val)                      {fMinPt = val;}
  void          SetPtMax(Double_t val)                      {fMaxPt = val;}

  void          SetInputArrayName(const char* name)         {fInputArrayName = name;}
  void          SetOutputArrayName(const char* name)        {fOutputArrayName = name;}


private:
  Bool_t              fInitialized;               // internal state when ExecOnce has been executed
  Bool_t              fRandomizeInPhi;            // randomize the particle's position in azimuth
  Bool_t              fRandomizeInEta;            // randomize the particle's position in pseudorap
  Bool_t              fRandomizeInPt;             // randomize the particle's position in Pt

  Double_t            fMinPhi;                    // range for phi for randomization
  Double_t            fMaxPhi;                    // range for phi for randomization
  Double_t            fMinEta;                    // range for eta for randomization
  Double_t            fMaxEta;                    // range for eta for randomization
  Double_t            fMinPt;                     // range for Pt for randomization
  Double_t            fMaxPt;                     // range for Pt for randomization

  TString             fInputArrayName;            // Name of the TClonesArray that will be loaded
  TString             fOutputArrayName;           // Name of the destination TClonesArray

  TClonesArray*       fInputArray;                //! TClonesArray that will be loaded
  TClonesArray*       fOutputArray;               //! Destination TClonesArray

  TRandom3*           fRandom;                    //! random number generator

  ClassDef(AliAnalysisTaskParticleRandomizer, 1);
};

#endif
