#ifndef ALIANALYSISTASKSESELECTHF_H
#define ALIANALYSISTASKSESELECTHF_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSESelectHF
// AliAnalysisTaskSE for the selection of heavy-flavour decay candidates
// and creation of a stand-alone AOD
// Author: A.Dainese, andrea.dainese@lnl.infn.it
//*************************************************************************


#include <TClonesArray.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"


class AliAnalysisTaskSESelectHF : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSESelectHF();
  AliAnalysisTaskSESelectHF(const char *name);
  virtual ~AliAnalysisTaskSESelectHF();

  AliAnalysisTaskSESelectHF(const AliAnalysisTaskSESelectHF &source);
  AliAnalysisTaskSESelectHF& operator=(const AliAnalysisTaskSESelectHF& source); 

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

  TClonesArray *fVerticesHFTClArr;     // Array of heavy-flavour vertices
  TClonesArray *fD0toKpiTClArr;        // Array of D0->Kpi
  Double_t fD0toKpiCuts[9];            // cuts for D0->Kpi selection
  
  ClassDef(AliAnalysisTaskSESelectHF,1); // AliAnalysisTaskSE for the reconstruction of heavy-flavour decay candidates
};

#endif

