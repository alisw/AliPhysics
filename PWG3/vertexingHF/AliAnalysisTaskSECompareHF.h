#ifndef ALIANALYSISTASKCOMPAREHF_H
#define ALIANALYSISTASKCOMPAREHF_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSECompareHF
// AliAnalysisTaskSE for the comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// Author: A.Dainese, andrea.dainese@ln.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

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

  
 private:

  AliAnalysisTaskSECompareHF(const AliAnalysisTaskSECompareHF &source);
  AliAnalysisTaskSECompareHF& operator=(const AliAnalysisTaskSECompareHF& source); 
  TList   *fOutput; //! list send on output slot 0
  TNtuple *fNtupleCmp; //! output ntuple
  TH1F    *fHistMass;    //! output histogram
  TH1F    *fHistNEvents;    //! output histogram
  AliAnalysisVertexingHF *fVHF;  // Vertexer heavy flavour (used to pass the cuts)
  
  ClassDef(AliAnalysisTaskSECompareHF,4); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif

