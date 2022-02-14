#ifndef ALIANALYSISTASKFilterUPCNanoAOD_H
#define ALIANALYSISTASKFilterUPCNanoAOD_H

/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id$ */

#include "AliAnalysisTaskSE.h"
#include "TString.h"

/// \ingroup rec
/// \class AliAnalysisTaskFilterUPCNanoAOD
/// \brief AliAnalysisTask to convert (filter) full AODs to muon AODs
/// \author Michal Broz

class AliAODBranchReplicator;
class AliAODEvent;

class AliAnalysisTaskFilterUPCNanoAOD : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskFilterUPCNanoAOD(Bool_t withSPDTracklets=kTRUE,Bool_t withMuonTracks=kFALSE,TString extraTriggers="");
  virtual ~AliAnalysisTaskFilterUPCNanoAOD();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

private:
  
  AliAnalysisTaskFilterUPCNanoAOD(const AliAnalysisTaskFilterUPCNanoAOD& rhs); // not implemented on purpose
  AliAnalysisTaskFilterUPCNanoAOD& operator=(const AliAnalysisTaskFilterUPCNanoAOD& rhs); // not implemented on purpose
  
  AliAODBranchReplicator* fBranchReplicator; ///< the class doing the real work
  Bool_t fWithMuonTracks;
  TString fExtraTriggers;
 
  ClassDef(AliAnalysisTaskFilterUPCNanoAOD,3) // class to convert std AOD to muon one
};

#endif
