#ifndef ALIANALYSISTASKSESELECTHF_H
#define ALIANALYSISTASKSESELECTHF_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSESelectHF
// AliAnalysisTaskSE for the selection of heavy-flavour decay candidates
// and creation of a stand-alone AOD
// Author: A.Dainese, andrea.dainese@lnl.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TClonesArray.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"


class AliAnalysisTaskSESelectHF : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSESelectHF();
  AliAnalysisTaskSESelectHF(const char *name);
  virtual ~AliAnalysisTaskSESelectHF();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
 private:

  AliAnalysisTaskSESelectHF(const AliAnalysisTaskSESelectHF &source);
  AliAnalysisTaskSESelectHF& operator=(const AliAnalysisTaskSESelectHF& source); 
  TClonesArray *fVerticesHFTClArr;     //! Array of heavy-flavour vertices
  TClonesArray *fD0toKpiTClArr;        //! Array of D0->Kpi
  AliAnalysisVertexingHF *fVHF; // analysis (used to pass the cuts)
  
  ClassDef(AliAnalysisTaskSESelectHF,2); // AliAnalysisTaskSE for the reconstruction of heavy-flavour decay candidates
};

#endif

