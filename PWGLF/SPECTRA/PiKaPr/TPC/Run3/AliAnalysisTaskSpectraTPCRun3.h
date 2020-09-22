/**************************************************************************
*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                         *
*  Author: The ALICE Off-line Project.                                    *
*  Contributors are mentioned in the code where appropriate.              *
*                                                                         *
*  Permission to use, copy, modify and distribute this software and its   *
*  documentation strictly for non-commercial purposes is hereby granted   *
*  without fee, provided that the above copyright notice appears in all   *
*  copies and that both the copyright notice and this permission notice   *
*  appear in the supporting documentation. The authors make no claims     *
*  about the suitability of this software for any purpose. It is          *
*  provided "as is" without express or implied warranty.                  *
**************************************************************************/

///
/// \author Nicolo' Jacazio
/// \file   AliAnalysisTaskSpectraTPCRun3.cxx
/// \brief  Task to produce a TPC spectra with O2 standards
/// \since  14/09/2020
///

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliExternalTrackParam.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLatex.h"
#include "TList.h"
#include "TMath.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#ifndef AliAnalysisTaskSpectraTPCRun3_H
#define AliAnalysisTaskSpectraTPCRun3_H

#define OLD_ALIPHYSICS // flag to compile with older aliphysics

#define AOD_ESD(AOD, ESD, Function) (fAODMode ? AOD->Function : ESD->Function)

class AliAnalysisTaskSpectraTPCRun3 : public AliAnalysisTaskSE {
  public:
  AliAnalysisTaskSpectraTPCRun3();
  AliAnalysisTaskSpectraTPCRun3(const char* name, const Bool_t readAODs, const Bool_t useO2PID);
  virtual ~AliAnalysisTaskSpectraTPCRun3();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  // Parameters
  Bool_t fMCmode = kFALSE;           // MC mode
  Bool_t fUseEventSelection = kTRUE; // Flag for event selection
  Double_t fVtxZMax = 10.f;          // Maximum Z of vertex window
  Double_t fEtaMax = 0.8f;           // Maximum eta window
  // Configuration flags
  const Bool_t fMakeTOFPlots = kFALSE; // Flag to produce TOF plots
  const Bool_t fAODMode;               // Flag to read AODs instead of ESDs
  const Bool_t fUseO2PID;              // Flag to use the PID like in O2
  // PID Parameters
  Double_t bbparam[7];     // Parametrization for the BetheBloch
  Double_t bbresoparam[2]; // Parametrization for the BetheBloch resolution

  private:
  AliAODEvent* fEventAOD = nullptr;       //! input event
  AliAODTrack* fTrackAOD = nullptr;       //! input track
  AliESDEvent* fEventESD = nullptr;       //! input event
  AliESDtrack* fTrackESD = nullptr;       //! input track
  AliPIDResponse* fPIDResponse = nullptr; //! PID response
  TList* fOutputList = nullptr;           //! output list
  AliEventCuts fEventCut;                 //! event selection
  AliESDtrackCuts* fTrackCuts = nullptr;  //! track cuts

  TH1F* hvertexz = nullptr;
  TH2F* htpcsignal = nullptr;
  TH2F* hexpEl = nullptr;
  TH2F* hexpMu = nullptr;
  TH2F* hexpPi = nullptr;
  TH2F* hexpKa = nullptr;
  TH2F* hexpPr = nullptr;
  TH2F* hexpDe = nullptr;
  TH2F* hexpTr = nullptr;
  TH2F* hexpHe = nullptr;
  TH2F* hexpAl = nullptr;
  TH2F* hnsigmaEl = nullptr;
  TH2F* hnsigmaMu = nullptr;
  TH2F* hnsigmaPi = nullptr;
  TH2F* hnsigmaKa = nullptr;
  TH2F* hnsigmaPr = nullptr;
  TH2F* hnsigmaDe = nullptr;
  TH2F* hnsigmaTr = nullptr;
  TH2F* hnsigmaHe = nullptr;
  TH2F* hnsigmaAl = nullptr;
  TH2F* htpcsignalEl = nullptr;
  TH2F* htpcsignalMu = nullptr;
  TH2F* htpcsignalPi = nullptr;
  TH2F* htpcsignalKa = nullptr;
  TH2F* htpcsignalPr = nullptr;
  TH2F* htpcsignalDe = nullptr;
  TH2F* htpcsignalTr = nullptr;
  TH2F* htpcsignalHe = nullptr;
  TH2F* htpcsignalAl = nullptr;

  TH1F* hpt = nullptr;
  TH1F* hpt_El = nullptr;
  TH1F* hpt_Pi = nullptr;
  TH1F* hpt_Ka = nullptr;
  TH1F* hpt_Pr = nullptr;
  TH1F* hp = nullptr;
  TH1F* hp_El = nullptr;
  TH1F* hp_Pi = nullptr;
  TH1F* hp_Ka = nullptr;
  TH1F* hp_Pr = nullptr;

  //
  Double_t tpcExpSignalEl = 0;
  Double_t tpcExpSignalMu = 0;
  Double_t tpcExpSignalPi = 0;
  Double_t tpcExpSignalKa = 0;
  Double_t tpcExpSignalPr = 0;
  Double_t tpcExpSignalDe = 0;
  Double_t tpcExpSignalTr = 0;
  Double_t tpcExpSignalHe = 0;
  Double_t tpcExpSignalAl = 0;
  Double_t tpcNSigmaEl = 0;
  Double_t tpcNSigmaMu = 0;
  Double_t tpcNSigmaPi = 0;
  Double_t tpcNSigmaKa = 0;
  Double_t tpcNSigmaPr = 0;
  Double_t tpcNSigmaDe = 0;
  Double_t tpcNSigmaTr = 0;
  Double_t tpcNSigmaHe = 0;
  Double_t tpcNSigmaAl = 0;
  //
  AliAnalysisTaskSpectraTPCRun3(const AliAnalysisTaskSpectraTPCRun3& t)
      : AliAnalysisTaskSE(t.GetName())
      , fMCmode(t.fMCmode)
      , fUseEventSelection(t.fUseEventSelection)
      , fVtxZMax(t.fVtxZMax)
      , fEtaMax(t.fEtaMax)
      , fMakeTOFPlots(t.fMakeTOFPlots)
      , fAODMode(t.fAODMode)
      , fUseO2PID(t.fUseO2PID)
      , fEventCut()
  {
    for (Int_t i = 0; i < 7; i++)
      bbparam[i] = t.bbparam[i];
    for (Int_t i = 0; i < 2; i++)
      bbresoparam[i] = t.bbresoparam[i];
  }
  AliAnalysisTaskSpectraTPCRun3& operator=(const AliAnalysisTaskSpectraTPCRun3&);
  // PID
  Double_t ExpectedSignal(AliPID::EParticleType id) const
  {
    if (!fUseO2PID) {
#ifdef OLD_ALIPHYSICS
      if (fAODMode)
        return fPIDResponse->GetTPCResponse().GetExpectedSignal(fTrackAOD, id);
      else
        return fPIDResponse->GetTPCResponse().GetExpectedSignal(fTrackESD, id);
#else
      if (fAODMode)
        return fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC, fTrackAOD, id);
      else
        return fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC, fTrackESD, id);
#endif
    }

    const Double_t betaGamma = AOD_ESD(fTrackAOD, fTrackESD, GetTPCmomentum()) / AliPID::ParticleMass(id);
    const Double_t charge = AliPID::ParticleCharge(id);
    const Double_t bb = AliExternalTrackParam::BetheBlochAleph(betaGamma, bbparam[0], bbparam[1], bbparam[2], bbparam[3], bbparam[4]);
    return bbparam[5] * bb * TMath::Power(charge, bbparam[6]);
  }
  Double_t ExpectedReso() const
  {
    if (!fUseO2PID) {
#ifdef OLD_ALIPHYSICS
      if (fAODMode)
        return fPIDResponse->GetTPCResponse().GetExpectedSigma(fTrackAOD, AliPID::kPion);
      else
        return fPIDResponse->GetTPCResponse().GetExpectedSigma(fTrackESD, AliPID::kPion);
#else
      if (fAODMode)
        return fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC, fTrackAOD, AliPID::kPion);
      else
        return fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC, fTrackESD, AliPID::kPion);
#endif
    }

    const Double_t tpcsignal = AOD_ESD(fTrackAOD, fTrackESD, GetTPCsignal());
    const Double_t tpcpoints = AOD_ESD(fTrackAOD, fTrackESD, GetTPCncls());
    return tpcsignal * bbresoparam[0] * (tpcpoints > 0 ? sqrt(1. + bbresoparam[1] / tpcpoints) : 1.f);
  }

  ClassDef(AliAnalysisTaskSpectraTPCRun3, 1);
};

#endif
