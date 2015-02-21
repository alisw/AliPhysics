/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ALIANALYSISTASKTRIGGERRATES_H
#define ALIANALYSISTASKTRIGGERRATES_H

/// \class AliAnalysisTaskTriggerRates
/// \brief task to study online/offline trigger combinations
//Author: Philippe Pillot - SUBATECH Nantes

#include "AliAnalysisTaskSE.h"

class AliCounterCollection;

class AliAnalysisTaskTriggerRates : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskTriggerRates();
  AliAnalysisTaskTriggerRates(const char *name);
  virtual ~AliAnalysisTaskTriggerRates();
  
  virtual void NotifyRun() {}
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  /// separate events according to their centrality
  void SelectCentrality(Bool_t flag = kTRUE) { fnCent = flag ? 10 : 0; }
  
  /// set the range of trigger deviation (15 ± deltaDev) for which the trigger sign is considered as unknown
  void SetZeroDevRange(UInt_t deltaDev) { fDeltaDev = deltaDev; }
  
  /// set the trigger class patterns (e.g. "CMUU7-B-NOPF-")
  void SetTrigClassPatterns(TString mul, TString mll, TString msl, TString msh) {
    fMULPattern = (mul.IsNull()) ? "none" : mul;
    fMLLPattern = (mll.IsNull()) ? "none" : mll;
    fMSLPattern = (msl.IsNull()) ? "none" : msl;
    fMSHPattern = (msh.IsNull()) ? "none" : msh;
  }
  
  /// print raw counts instead of rates
  void PrintCounts(Bool_t flag = kTRUE) { fPrinfCounts = flag; }
  
private:
  AliAnalysisTaskTriggerRates(const AliAnalysisTaskTriggerRates &); // not implemented
  AliAnalysisTaskTriggerRates &operator=(const AliAnalysisTaskTriggerRates &); // not implemented
  
  void InitCentralityBins();
  Int_t TriggerDevSign(AliVParticle *track, UInt_t deltaDev) const;
  void PrintRates(TString ps, TString cent) const;
  
  AliCounterCollection *fTriggerCounters; //! trigger counters
  
  Int_t   fnCent;               /// number of centrality bin used (< 10)
  Float_t fCentBinRange[10][2]; //! centrality bin intervals
  TString fCentBinName[10];     //! centrality bin names
  UInt_t  fDeltaDev;            /// set trgSign = 0 for trigger track with dev = 15 ± fDeltaDev
  TString fMULPattern;          /// MUL trigger class pattern
  TString fMLLPattern;          /// MLL trigger class pattern
  TString fMSLPattern;          /// MSL trigger class pattern
  TString fMSHPattern;          /// MSH trigger class pattern
  Bool_t fPrinfCounts;          /// print raw counts instead of rates
  
  ClassDef(AliAnalysisTaskTriggerRates, 1);
};

#endif

