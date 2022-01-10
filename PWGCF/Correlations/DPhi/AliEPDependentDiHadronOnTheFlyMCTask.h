#ifndef ALIVZEROEPONTHEFLYMCTASK_H
#define ALIVZEROEPONTHEFLYMCTASK_H
/* Copyright(c) 1998-1999, AlICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//*****************************************************************
//    Class AliEPDependentDiHadronOnTheFlyMCTask
//    author: Darius Keijdener
//    09/11/2021
//    This analysis task constructs an event plane based on the 
//    VZERO detector range and consecutively fills several dihadron 
//    histograms with the dihadron-pair triggers having a certain 
//    angle to that event plane.
//*****************************************************************

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliUEHistograms.h"
#include "AliUEHist.h"
#include "AliMultSelection.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TList.h>
#include <TObjArray.h>
#include <TArrayF.h>
#include <TString.h>


class AliEPDependentDiHadronOnTheFlyMCTask : public AliAnalysisTaskSE {

 public:
  AliEPDependentDiHadronOnTheFlyMCTask();
  AliEPDependentDiHadronOnTheFlyMCTask(const char* name);
  virtual ~AliEPDependentDiHadronOnTheFlyMCTask();

  void UserCreateOutputObjects();
  void UserExec(Option_t* option);
  void Terminate(Option_t* option);

  void GetEventDetails(AliMCEvent* mcEvent, Double_t& reactionplane, Double_t& zvtx, Double_t& centrality);
  TObjArray* SelectTracks(AliMCEvent* mcEvent);
  void FillEventLevelQA(AliMCEvent* mcEvent);
  void FillEventPlaneHistograms(TObjArray* allTracks, Double_t reactionplane, Double_t& eventplaneV0);
  TH1D* FillVirtualV0(TObjArray* allTracks, Double_t etaMin, Double_t etaMax);
  Double_t ComputeTrueEventPlane(TObjArray* allTracks, Int_t harmonic);
  Double_t ComputeEventPlane(TH1D* sectors, Int_t harmonic);
  Double_t angleBetween(Double_t phi1, Double_t phi2);
  Double_t twoPiToPi(Double_t phi);
  Double_t piToHalfPi(Double_t phi);
  
  TList* fListOutputHistograms;
  TProfile* fProfileFractionPrimaryTracks; // The fraction of tracks that are logged as primary tracks.
  TH1D* fHistNoEvents; // The number of events that triggered the UserExec funtion.
  TH2D* fHistImpactVsMultiplicity; // The impact parameter versus the event multiplicity
  TH1D* fHistReactionPlane; // The distribution of the 2-reactionplane
  TH1D* fHistEventReactionPlane; // The 2-eventplane w.r.t. the 2-reactionplane.
  TH1D* fHistTrueEventReactionPlane; // The distribution of the the 2-eventplane as determined with all particles w.r.t. the 2-reactionplane.
  AliUEHistograms* fHistos; // Handles the dihadron histograms for resp allplane, in-plane, mid-plane and out-of-plane triggers.
  AliUEHistograms* fHistosIn;
  AliUEHistograms* fHistosMid;
  AliUEHistograms* fHistosOut;
  TString fCustomBinning;

  void SetTrackEtaCut(Double_t val) { fTrackEtaCut = val; }
  void SetCustomBinning(const char* binningStr) { fCustomBinning = binningStr; }

 private:
  Double_t fEtaMinV0A; // The V0A resp. C boundaries as used for calculating the event plane
  Double_t fEtaMaxV0A;
  Double_t fEtaMinV0C;
  Double_t fEtaMaxV0C;
  Double_t fTrackEtaCut; // The track cut for particles included in the dihadron histograms.
  Int_t fNoSectorsV0; // The amount of sectors in phi that the simulated V0 detector has
  TString fCentralityEstimator; // The estimator used to determine the centrality. Choices: V0M, TRK, TKL, ZDC or FMD.

  ClassDef(AliEPDependentDiHadronOnTheFlyMCTask, 1);
};

#endif
