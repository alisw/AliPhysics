#ifndef ALIANALYSISTASKPARTICLERANDOMIZER_H
#define ALIANALYSISTASKPARTICLERANDOMIZER_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
 * \class AliAnalysisTaskParticleRandomizer
 * \brief Particle randomization task
 *
 * This task clones the tracks and randomize them in a given acceptance
 * Use ActivateJetRemoval() to remove the given jets from the event before randomization
 * 
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Apr 21, 2016
 */
// 
class TClonesArray;
class TString;
class TRandom3;
class TH2D;

class AliAnalysisTaskParticleRandomizer : public AliAnalysisTaskEmcal {
public:
  AliAnalysisTaskParticleRandomizer();
  virtual ~AliAnalysisTaskParticleRandomizer();
  virtual void  UserCreateOutputObjects();
  virtual void  Terminate(Option_t *) {}

  // ###### Configuration setters
  void          SetRandomizeInPhi(Bool_t val)               {fRandomizeInPhi = val;}
  void          SetRandomizeInEta(Bool_t val)               {fRandomizeInEta = val;}
  void          SetRandomizeInTheta(Bool_t val)             {fRandomizeInTheta = val;}
  void          SetRandomizeInPt(Bool_t val)                {fRandomizeInPt = val;}

  void          SetPhiMin(Double_t val)                     {fMinPhi = val;}
  void          SetPhiMax(Double_t val)                     {fMaxPhi = val;}
  void          SetEtaMin(Double_t val)                     {fMinEta = val;}
  void          SetEtaMax(Double_t val)                     {fMaxEta = val;}
  void          SetPtMin(Double_t val)                      {fMinPt = val;}
  void          SetPtMax(Double_t val)                      {fMaxPt = val;}
  void          SetDistributionV2(TH2D* dist)               {fDistributionV2 = dist;}
  void          ActivateJetRemoval(const char* arrName, Double_t threshold, const char* rhoObj) {fJetRemovalArrayName = arrName; fJetRemovalPtThreshold = threshold; fJetRemovalRhoObj = rhoObj;}

  void          SetInputArrayName(const char* name)         {fInputArrayName = name;}
  void          SetOutputArrayName(const char* name)        {fOutputArrayName = name;}

protected:
  Bool_t        Run();
  void          ExecOnce();
  Double_t      AddFlow(Double_t phi, Double_t pt);

  Bool_t              fRandomizeInPhi;            /// randomize the particle's position in azimuth
  Bool_t              fRandomizeInEta;            /// randomize the particle's position in pseudorap
  Bool_t              fRandomizeInTheta;          /// randomize the particle's position in theta
  Bool_t              fRandomizeInPt;             /// randomize the particle's position in Pt

  Double_t            fMinPhi;                    /// range for phi for randomization
  Double_t            fMaxPhi;                    /// range for phi for randomization
  Double_t            fMinEta;                    /// range for eta for randomization
  Double_t            fMaxEta;                    /// range for eta for randomization
  Double_t            fMinPt;                     /// range for Pt for randomization
  Double_t            fMaxPt;                     /// range for Pt for randomization

  TH2D*               fDistributionV2;            /// Distribution for v2 in bins of pt and centrality

  TString             fInputArrayName;            /// Name of the TClonesArray that will be loaded
  TString             fOutputArrayName;           /// Name of the destination TClonesArray

  TClonesArray*       fInputArray;                //!<! TClonesArray that will be loaded
  TClonesArray*       fOutputArray;               //!<! Destination TClonesArray

  TString             fJetRemovalRhoObj;          /// Name of array to rho object
  TString             fJetRemovalArrayName;       /// Name of the TClonesArray containing jets for removal that will be loaded
  TClonesArray*       fJetRemovalArray;           //!<! TClonesArray containing jets
  Double_t            fJetRemovalPtThreshold;     /// threshold at which jets given in fInputJetArray will be removed

  TRandom3*           fRandom;                    //!<! random number generator

  Bool_t              IsParticleInJet(Int_t part);
  Double_t            GetExternalRho();

  ClassDef(AliAnalysisTaskParticleRandomizer, 4);
};

#endif
