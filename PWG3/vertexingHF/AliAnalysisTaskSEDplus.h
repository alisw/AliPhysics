#ifndef ALIANALYSISTASKDPLUS_H
#define ALIANALYSISTASKDPLUS_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEDplus
// AliAnalysisTaskSE for the comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)

//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

class AliAnalysisTaskSEDplus : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEDplus();
  AliAnalysisTaskSEDplus(const char *name, Bool_t fillNtuple=kFALSE);
  virtual ~AliAnalysisTaskSEDplus();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  
 private:

  AliAnalysisTaskSEDplus(const AliAnalysisTaskSEDplus &source);
  AliAnalysisTaskSEDplus& operator=(const AliAnalysisTaskSEDplus& source); 
  TList   *fOutput; //! list send on output slot 0
  TNtuple *fNtupleDplus; //! output ntuple
  TNtuple *fNtupleDplusbackg; //! output ntuple
  Bool_t fFillNtuple;   // flag for filling ntuple
  AliAnalysisVertexingHF *fVHF;  // Vertexer heavy flavour (used to pass the cuts)
  
  ClassDef(AliAnalysisTaskSEDplus,4); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif

