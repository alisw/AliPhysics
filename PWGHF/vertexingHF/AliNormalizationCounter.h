#ifndef ALINORMALIZATIONCOUNTER_H
#define ALINORMALIZATIONCOUNTER_H

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */ 

//*************************************************************************
// Class AliNormalizationCounter
// Class to store the informations relevant for the normalization in the 
// barrel for each run
// Authors: G. Ortona, ortona@to.infn.it
// D. Caffarri, davide.caffarri@pd.to.infn.it
// with many thanks to P. Pillot
/////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliAODEvent.h>
#include <AliVParticle.h>
#include "AliAnalysisTaskSE.h"
#include "AliCounterCollection.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliRDHFCuts.h"
//#include "AliAnalysisVertexingHF.h"

class AliNormalizationCounter : public TNamed
{
 public:

  AliNormalizationCounter();
  AliNormalizationCounter(const char *name);
  virtual ~AliNormalizationCounter();
  Long64_t Merge(TCollection* list);

  AliCounterCollection* GetCounter(){return &fCounters;}
  void Init();
  void Add(const AliNormalizationCounter*);
  void SetESD(Bool_t flag){fESD=flag;}
  void SetStudyMultiplicity(Bool_t flag){ fMultiplicity=flag; }
  void StoreEvent(AliVEvent*,AliRDHFCuts *,Bool_t mc=kFALSE);
  void StoreCandidates(AliVEvent*, Int_t nCand=0,Bool_t flagFilter=kTRUE);
  TH1D* DrawAgainstRuns(TString candle="candid(filter)",Bool_t drawHist=kTRUE);
  TH1D* DrawRatio(TString candle1="candid(filter)",TString candle2="triggered");
  void PrintRubrics();
  Double_t GetSum(TString candle="triggered");
  TH2F* GetHist(Bool_t filtercuts=kTRUE,Bool_t spdtracklets=kTRUE,Bool_t drawHist=kFALSE);
  Double_t GetNEventsForNorm();
  Double_t GetNEventsForNorm(Int_t runnumber);
  Bool_t GetStudyMultiplicity(){ return fMultiplicity; }
  Double_t GetNEventsForNorm(Int_t minmultiplicity, Int_t maxmultiplicity);
  TH1D* DrawNEventsForNorm(Bool_t drawRatio=kFALSE);

 private:
  AliNormalizationCounter(const AliNormalizationCounter &source);
  AliNormalizationCounter& operator=(const AliNormalizationCounter& source);
  Int_t Multiplicity(AliVEvent* event);
  AliCounterCollection fCounters; //internal counter
  Bool_t fESD; //flag for ESD vs AOD
  Bool_t fMultiplicity; //flag for multiplicity
  TH2F *fHistTrackFilterEvMult; //hist to store no of filter candidates vs no of tracks in the event 
  TH2F *fHistTrackAnaEvMult;//hist to store no of analysis candidates vs no of tracks in the event 
  TH2F *fHistTrackFilterSpdMult; //hist to store no of filter candidates vs  SPD multiplicity 
  TH2F *fHistTrackAnaSpdMult;//hist to store no of analysis candidates vs SPD multiplicity 

  ClassDef(AliNormalizationCounter,5);

};
#endif
