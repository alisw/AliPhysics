#ifndef ALIANALYSISTASKPILEUP_H
#define ALIANALYSISTASKPILEUP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup muondep
/// \class AliAnalysisTaskPileup
/// \brief Trigger scaler analysis for pileup corrections
/// Based on the work by L. Aphecetche - SUBATECH Nantes
//Author: Diego Stocco - SUBATECH Nantes

#define READOCDB

#include "AliAnalysisTaskSE.h"

class TObjArray;
class TString;
class TArrayI;
class AliCounterCollection;

#ifdef READOCDB
class AliTriggerRunScalers;
#endif

class AliAnalysisTaskPileup : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskPileup();
  AliAnalysisTaskPileup(const char *name);
  virtual ~AliAnalysisTaskPileup();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  virtual void   NotifyRun();

  void SetDefaultStorage(TString dbString);
  void SetSpecificStorage(TString calibType, TString dbString);
  
private:
  
  /// Not implemented
  AliAnalysisTaskPileup(const AliAnalysisTaskPileup& rhs);
  /// Not implemented
  AliAnalysisTaskPileup& operator = (const AliAnalysisTaskPileup& rhs);

  Double_t GetL0Correction(Double_t nCINT1B, Double_t nCBEAMB);

  enum {
    kHevents,  /// Number of events histogram
    kHeventsCorrectL0, /// Number of L0 corrected events histogram
    kNeventHistos  /// Number of trigger histograms
  };

  AliCounterCollection* fEventCounters; //!< Event statistics
  TObjArray* fHistoEventsList;   //!< List of event histograms

  TObjArray* fTriggerClasses; //!< full trigger class name
  TArrayI* fTriggerClassIndex;  //!< Trigger classes mask

  Bool_t fIsInitCDB; //!< Flag telling if CDB is used
  TAxis* fCentralityClasses; //!< Centrality classes

#ifdef READOCDB
  AliTriggerRunScalers* fTriggerRunScalers; //!< Trigger scalers from OCDB
  TObjArray* fStorageList; /// List of storages
#endif
  
  ClassDef(AliAnalysisTaskPileup, 1);
};

#endif

