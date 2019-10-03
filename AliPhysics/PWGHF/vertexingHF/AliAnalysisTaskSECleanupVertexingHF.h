#ifndef ALIANALYSISTASKSECLEANUPVERTEXINGHF_H
#define ALIANALYSISTASKSECLEANUPVERTEXINGHF_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSECleanupVertexingHF
// AliAnalysisTaskSE for HF candidates to delete the on-the-fly secondary 
// vertex filled in FillRecoCand
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include "AliAnalysisTaskSE.h"

class AliAODEvent;


class AliAnalysisTaskSECleanupVertexingHF : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSECleanupVertexingHF();
  AliAnalysisTaskSECleanupVertexingHF(const char *name);
  virtual ~AliAnalysisTaskSECleanupVertexingHF();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  //private:
  /// TH1F     *fNentries;            //! histogram with number of events on output slot 3

  ClassDef(AliAnalysisTaskSECleanupVertexingHF,1); // AliAnalysisTaskSE to delete the on-the-fly reconstructed secondary vertex
};

#endif

