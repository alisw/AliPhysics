#ifndef ALIANALYSISTASKCOMPAREHF_H
#define ALIANALYSISTASKCOMPAREHF_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSECompareHF
// AliAnalysisTaskSE for the comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// Author: A.Dainese, andrea.dainese@lnl.infn.it
//*************************************************************************

#include <TNtuple.h>
#include <TH1F.h>

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSECompareHF : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSECompareHF();
  AliAnalysisTaskSECompareHF(const char *name);
  virtual ~AliAnalysisTaskSECompareHF();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetD0toKpiCuts(Double_t cut0=1000.,Double_t cut1=100000.,
		      Double_t cut2=1.1,Double_t cut3=0.,Double_t cut4=0.,
		      Double_t cut5=100000.,Double_t cut6=100000.,
		      Double_t cut7=100000000.,Double_t cut8=-1.1); 
  // cuts[0] = inv. mass half width [GeV]   
  // cuts[1] = dca [cm]
  // cuts[2] = cosThetaStar 
  // cuts[3] = pTK [GeV/c]
  // cuts[4] = pTPi [GeV/c]
  // cuts[5] = d0K [cm]   upper limit!
  // cuts[6] = d0Pi [cm]  upper limit!
  // cuts[7] = d0d0 [cm^2]
  // cuts[8] = cosThetaPoint
  void SetD0toKpiCuts(const Double_t cuts[9]); 
  
 private:

  AliAnalysisTaskSECompareHF(const AliAnalysisTaskSECompareHF &source);
  AliAnalysisTaskSECompareHF& operator=(const AliAnalysisTaskSECompareHF& source); 
  TList   *fOutput; //! list send on output slot 0
  TNtuple *fNtupleD0Cmp; // output ntuple
  TH1F    *fHistMass;    // output histogram
  Double_t fD0toKpiCuts[9]; // cuts for D0->Kpi selection
  
  ClassDef(AliAnalysisTaskSECompareHF,1); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif

