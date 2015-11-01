#ifndef ALIANALYSISTASKAOD2MUONAOD_H
#define ALIANALYSISTASKAOD2MUONAOD_H

/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id$ */

#include "AliAnalysisTaskSE.h"
#include "TString.h"

/// \ingroup rec
/// \class AliAnalysisTaskAOD2MuonAOD
/// \brief AliAnalysisTask to convert (filter) full AODs to muon AODs
/// \author L. Aphecetche (Subatech)

class AliAODBranchReplicator;
class AliMuonEventCuts;
class AliAODEvent;

class AliAnalysisTaskAOD2MuonAOD : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskAOD2MuonAOD(Int_t mcMode=3, Bool_t withSPDTracklets=kTRUE, AliMuonEventCuts* eventCuts=0x0);
  virtual ~AliAnalysisTaskAOD2MuonAOD();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

private:
  
  AliAnalysisTaskAOD2MuonAOD(const AliAnalysisTaskAOD2MuonAOD& rhs); // not implemented on purpose
  AliAnalysisTaskAOD2MuonAOD& operator=(const AliAnalysisTaskAOD2MuonAOD& rhs); // not implemented on purpose
  
  Bool_t HasMuonInformation(const AliAODEvent& event) const;

  AliAODBranchReplicator* fBranchReplicator; ///< the class doing the real work
  
  AliMuonEventCuts* fMuonEventCuts;
  
  ClassDef(AliAnalysisTaskAOD2MuonAOD,2) // class to convert std AOD to muon one
};

#endif
