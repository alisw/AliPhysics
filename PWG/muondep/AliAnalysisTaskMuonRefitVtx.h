#ifndef ALIANALYSISTASKMUONREFITVTX_H
#define ALIANALYSISTASKMUONREFITVTX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup pwg_muondep_misc
/// \class AliAnalysisTaskMuonRefitVtx
/// \brief intermediate task to extrapolate tracks to another vertex
//Author: Philippe Pillot - SUBATECH Nantes

#include <TString.h>
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMuonRefitVtx : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskMuonRefitVtx();
  AliAnalysisTaskMuonRefitVtx(const char *name);
  virtual ~AliAnalysisTaskMuonRefitVtx();
  
  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fDefaultStorage = ocdbPath; }
  
  /// Use mean SPD vertex instead of (0,0,0)
  void UseMeanVtxSPD(Bool_t flag = kTRUE) {
    fUseMeanVtxSPD = flag;
    AliInfo("!!! Do not mix runs: NotifyRun() which loads the vertex is called only once !!!"); }
  
  /// Use MC vertex instead of (0,0,0)
  void UseMCVtx(Bool_t flag = kTRUE) { fUseMCVtx = flag; }
  
  /// Keep the same vertex (from track position) instead of (0,0,0)
  void UseTrackVtx(Bool_t flag = kTRUE) { fUseTrackVtx = flag; }
  
  // Shift the vertex by (dx,dy,dz)
  void ShiftVtx(Double_t dx, Double_t dy, Double_t dz);
  
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
  Bool_t   fUseMeanVtxSPD;  ///< use mean SPD vertex instead of (0,0,0)
  Bool_t   fUseMCVtx;       ///< use MC vertex instead of (0,0,0)
  Bool_t   fUseTrackVtx;    ///< keep the same vertex (from track position) instead of (0,0,0)
  Double_t fVtxPos[3];      ///< vertex position
  Double_t fVtxSig[3];      ///< vertex resolution
  Double_t fVtxShift[3];    ///< vertex shift
  
  ClassDef(AliAnalysisTaskMuonRefitVtx, 2);
};

//________________________________________________________________________
inline void AliAnalysisTaskMuonRefitVtx::ShiftVtx(Double_t dx, Double_t dy, Double_t dz)
{
  /// Shift the vertex by (dx,dy,dz)
  
  fVtxShift[0] = dx; fVtxShift[1] = dy; fVtxShift[2] = dz;
  
  // shift the default vertex
  for (Int_t i = 0; i < 3; i++) fVtxPos[i] += fVtxShift[i];
  
}

#endif

