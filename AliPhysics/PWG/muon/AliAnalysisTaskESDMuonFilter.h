#ifndef ALIANALYSISTASKESDMUONFILTER_H
#define ALIANALYSISTASKESDMUONFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

///
/// \brief Add the muon tracks to the generic AOD track branch during the
/// filtering of the ESD.
///
/// \author R. Arnaldi 5/5/08 and L. Aphecetche January 2011

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include "TMatrixD.h"

class AliAnalysisFilter;

class AliAnalysisTaskESDMuonFilter : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskESDMuonFilter(Bool_t onlyMuon=kTRUE, Bool_t keepAllEvents=kTRUE, Int_t mcMode=0, Bool_t withSPDtracklets=kFALSE);
  AliAnalysisTaskESDMuonFilter(const char* name, Bool_t onlyMuon=kTRUE, Bool_t keepAllEvents=kTRUE, Int_t mcMode=0, Bool_t withSPDtracklets=kFALSE);
  virtual ~AliAnalysisTaskESDMuonFilter() {;}
  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  
  virtual void ConvertESDtoAOD();
  
  virtual void SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
  void SetWriteMuonAOD(Bool_t enableMuonAOD){fEnableMuonAOD = enableMuonAOD;}
  void SetWithSPDtracklets(Bool_t withSPDtracklets=kTRUE) { fWithSPDTracklets = withSPDtracklets; }
  void SetMCMode(Int_t mcMode) { fMCMode=mcMode; }
  
  void PrintTask(Option_t *option="", Int_t indent=0) const;

  void ConvertCovMatrixMUON2AOD(const TMatrixD& covMUON, Double_t covAOD[21]);
  
private:
  AliAnalysisTaskESDMuonFilter(const AliAnalysisTaskESDMuonFilter&);
  AliAnalysisTaskESDMuonFilter& operator=(const AliAnalysisTaskESDMuonFilter&);
  void AddFilteredAOD(const char* aodfilename, const char* title);
  
  AliAnalysisFilter* fTrackFilter; ///<  Track Filter
  Bool_t fEnableMuonAOD; ///< flag for enabling Muon AOD production
  Bool_t fOnlyMuon; ///< flag for disabling branches irrelevant for (most) muon analyses
  Bool_t fKeepAllEvents; ///< keep even events where there's no muons (to get e.g. unbiased vertex distribution)
  Int_t  fMCMode; ///< whether and how we're filtering MC data
  Bool_t fWithSPDTracklets; ///< whether or not we keep SPD tracklets
  
  ClassDef(AliAnalysisTaskESDMuonFilter, 7); // Analysis task for standard ESD filtering
};

#endif
