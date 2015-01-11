/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Authors: Henrik Qvigstad, Dmitri Peressounko
// Date   : 20.06.2013
// Adopted for AOD analysis by Boris Polishchuk (10.03.2014)
/* $Id$ */


#ifndef ALIANALYSISTASKPI0FLOWMCAOD_H
#define ALIANALYSISTASKPI0FLOWMCAOD_H

class AliAODMCParticle;
class AliCaloPhoton;

#include "AliAnalysisTaskPi0Flow.h"

class AliAnalysisTaskPi0FlowMCAOD : public AliAnalysisTaskPi0Flow
{
public:
  AliAnalysisTaskPi0FlowMCAOD(const char* name = "AliAnalysisTaskPi0Flow", Period period = kUndefinedPeriod);
  virtual ~AliAnalysisTaskPi0FlowMCAOD();
    
  void SetOffVertexPhotonCut(Bool_t setCut=kTRUE) { kOffVertexCutSet=setCut; }
    
protected: // Override:
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  // Pi0FlowTask
  virtual void SelectPhotonClusters();
  virtual void FillSelectedClusterHistograms();
  virtual void ConsiderPi0s();
  virtual void ConsiderPi0sMix();
  virtual void ProcessMC();

protected: // member functions:
  AliAnalysisTaskPi0FlowMCAOD(const AliAnalysisTaskPi0FlowMCAOD&); // not implemented
  AliAnalysisTaskPi0FlowMCAOD& operator=(const AliAnalysisTaskPi0FlowMCAOD&); // not implemented
  
  TClonesArray* GetMCArray();
  AliAODMCParticle* GetParticle(Int_t); //Returns particle at given position for AOD

protected: // member variables:
  TClonesArray* fMcArray; //mcArray for AOD MC particles
  Bool_t kOffVertexCutSet;
    
  void FillMCHist();
  Double32_t R(AliAODMCParticle* p);
  
  virtual Double_t PrimaryWeight(Int_t primary);
  virtual Double_t PrimaryParticleWeight(AliAODMCParticle * particle);
  void FillSecondaries() ;
  Int_t FindPrimary(AliVCluster* clu,  Bool_t& sure);
  Int_t FindCommonParent(Int_t iPart, Int_t jPart) ;
  Bool_t HaveParent(Int_t iPart, Int_t pdgParent);
  Bool_t InPi0mass(Double_t m, Double_t pt);

  void FillAllHistograms(const char* particleName, AliCaloPhoton* ph1);

  static const Double_t kRCut;
  enum ParticleID {kEta=221};


  ClassDef(AliAnalysisTaskPi0FlowMCAOD, 1); // PHOS analysis task
};

#endif // ALIANALYSISTASKPI0FLOWMCAOD_H
