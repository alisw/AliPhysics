#ifndef ALIANALYSISTASKMUONREFITVTX_H
#define ALIANALYSISTASKMUONREFITVTX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup muondep
/// \class AliAnalysisTaskMuonRefitVtx
/// \brief intermediate task to extrapolate tracks to another vertex
//Author: Philippe Pillot - SUBATECH Nantes

#include <TString.h>
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMuonRefitVtx : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskMuonRefitVtx();
  AliAnalysisTaskMuonRefitVtx(const char *name);
  virtual ~AliAnalysisTaskMuonRefitVtx();
  
  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fDefaultStorage = ocdbPath; }
  
  /// Use mean SPD vertex (x,y,0) instead of (0,0,0)
  void UseMeanVtxSPD(Bool_t flag = kTRUE) { fUseMeanVtxSPD = flag; }
  
  /// Use MC vertex (x,y,z) instead of (0,0,0)
  void UseMCVtx(Bool_t flag = kTRUE) { fUseMCVtx = flag; }
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   NotifyRun();
  virtual void   Terminate(Option_t *);
  
private:
  
  /// Not implemented
  AliAnalysisTaskMuonRefitVtx(const AliAnalysisTaskMuonRefitVtx& rhs);
  /// Not implemented
  AliAnalysisTaskMuonRefitVtx& operator = (const AliAnalysisTaskMuonRefitVtx& rhs);
  
private:
  
  TString  fDefaultStorage; ///< location of the default OCDB storage
  Bool_t   fOCDBLoaded;     ///< OCDB info properly loaded
  Bool_t   fUseMeanVtxSPD;  ///< use mean SPD vertex (x,y,0) instead of (0,0,0)
  Bool_t   fUseMCVtx;       ///< use MC vertex (x,y,z) instead of (0,0,0)
  Double_t fVtxPos[3];      ///< vertex position
  Double_t fVtxSig[3];      ///< vertex resolution
  
  ClassDef(AliAnalysisTaskMuonRefitVtx, 1);
};

#endif

