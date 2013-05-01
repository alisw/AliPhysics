/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Inclusion of AliPHOSHijingEfficiency,
// by Dmitri Peressounko, 05.02.2013
// Authors: Henrik Qvigstad, Dmitri Peressounko
// Date   : 05.04.2013
/* $Id$ */


#ifndef ALIANALYSISTASKPI0FLOWMC_H
#define ALIANALYSISTASKPI0FLOWMC_H

class TParticle;

#include <AliAnalysisTaskPi0Flow.h>


class AliAnalysisTaskPi0FlowMC : public AliAnalysisTaskPi0Flow
{
public:
  AliAnalysisTaskPi0FlowMC(const char* name = "AliAnalysisTaskPi0Flow", Period period = kUndefinedPeriod);
  virtual ~AliAnalysisTaskPi0FlowMC();

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
  AliAnalysisTaskPi0FlowMC(const AliAnalysisTaskPi0FlowMC&); // not implemented
  AliAnalysisTaskPi0FlowMC& operator=(const AliAnalysisTaskPi0FlowMC&); // not implemented
  
  AliStack* GetMCStack();

protected: // member variables:
  AliStack* fStack;
  
  void FillMCHist();
  
  Double_t PrimaryWeight(Int_t primary);
  Double_t PrimaryParticleWeight(TParticle * particle);
  void FillSecondaries() ;
  Int_t FindPrimary(AliVCluster* clu,  Bool_t& sure);
  Int_t FindCommonParent(Int_t iPart, Int_t jPart) ;
  Bool_t HaveParent(Int_t iPart, Int_t pdgParent);
  Bool_t InPi0mass(Double_t m, Double_t pt);

  void FillAllHistograms(const char* particleName, AliCaloPhoton* ph1);

  static const Double_t kRCut = 1.;
  enum ParticleID {kEta=221};


  ClassDef(AliAnalysisTaskPi0FlowMC, 1); // PHOS analysis task
};

#endif // ALIANALYSISTASKPI0FLOWMC_H
