#ifndef AliAnalysisTaskEmcalJetPerformance_H
#define AliAnalysisTaskEmcalJetPerformance_H

/**********************************************************************************
 * Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                *
 * All rights reserved.                                                            *
 *                                                                                 *
 * Redistribution and use in source and binary forms, with or without              *
 * modification, are permitted provided that the following conditions are met:     *
 *   * Redistributions of source code must retain the above copyright              *
 *     notice, this list of conditions and the following disclaimer.               *
 *   * Redistributions in binary form must reproduce the above copyright           *
 *     notice, this list of conditions and the following disclaimer in the         *
 *     documentation and/or other materials provided with the distribution.        *
 *   * Neither the name of the <organization> nor the                              *
 *     names of its contributors may be used to endorse or promote products        *
 *     derived from this software without specific prior written permission.       *
 *                                                                                 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY             *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    *
 * *********************************************************************************/

/**
 * \file AliAnalysisTaskEmcalJetPerformance.h
 * \brief Studies of full jets.
 *
 * This header file declares the class AliAnalysisTaskEmcalJetPerformance.
 *
 * \author James Mulligan <james.mulligan@yale.edu>, Yale University
 * \date Oct 20, 2017
 */

#include "THistManager.h"

#include "AliEventCuts.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliEmcalEmbeddingQA.h"

class AliAnalysisTaskEmcalJetPerformance : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetPerformance()                                          ;
  AliAnalysisTaskEmcalJetPerformance(const char *name)                          ;
  virtual ~AliAnalysisTaskEmcalJetPerformance()                                 ;
  
  static AliAnalysisTaskEmcalJetPerformance* AddTaskEmcalJetPerformance(
    const char *ntracks            = "usedefault",
    const char *nclusters          = "usedefault",
    const char *nGenLev            = "mcparticles",
    const Double_t minTrPt         = 0.15,              // Minimum track pT in standard track container
    const Double_t minClPt         = 0.30,              // Minimum cluster E in standard cluster container
    const char *suffix             = "");
  
  // Truth-level particle types
  enum ContributorType {
    kUndefined      = -1,  //!< Undefined
    kPhoton         = 0,   //!< Photon (direct or decay)
    kChargedPion    = 1,   //!< Charged pion
    kProton         = 2,   //!< Proton
    kAntiProton     = 3,   //!< Antiproton
    kChargedKaon    = 4,   //!< Charged Kaon
    kK0L            = 5,   //!< K0L
    kNeutron        = 6,   //!< Neutron
    kAntiNeutron    = 7,   //!< Antineutron
    kElectron       = 8,   //!< Electron
    kMuon           = 9,   //!< Muon
    kOther          = 10   //!< Other
  };
  
  // Detector-level particle types, based on PhysicalPrimary contributors
  enum ParticleType {
    kNotDefined               = -1, //!< Undefined
    kSinglePhoton             = 0,  //!< Photon (direct or decay) is the only contributor
    kSingleElectron           = 1,  //!< Electron is the only contributor
    kSingleChargedPion        = 2,  //!< Charged pion is the only contributor
    kSingleProton             = 3,  //!< Proton is the only contributor
    kSingleAntiProton         = 4,  //!< Antiproton is the only contributor
    kSingleChargedKaon        = 5,  //!< Kaon is the only contributor
    kSingleK0L                = 6,  //!< K0L is the only contributor
    kSingleNeutron            = 7,  //!< Neutron is the only contributor
    kSingleAntiNeutron        = 8,  //!< Antineutron is the only contributor
    kSingleOther              = 9,  //!< One contributor (excluding the above cases)
    kPhotonHadron             = 10, //!< Photon+Hadron are the only contributors, with Photon leading
    kHadronPhoton             = 11, //!< Hadron+Photon are the only contributors, with Hadron leading
    kMergedPi0                = 12, //!< Two particles from merged pi0 are the only contributors
    kPhotonPhotonOther        = 13, //!< Photon+Photon (not from same pi0) are the only contributors
    kHadronHadron             = 14, //!< Hadron+Hadron are the only contributors
    kTwoContributorsOther     = 15, //!< Two contributors (excluding the above cases)
    kMoreThanTwoContributors  = 16  //!< More than two contributors
  };
  
  enum ClusterType {
    kNA       = -1,//!< Undefined
    kEMCal    = 0, //!< EMCal
    kDCal     = 1, //!< DCal
  };

  void UserCreateOutputObjects()                                                ;
  
  // Setters
  void SetMinPt(Double_t d)                                 { fMinPt = d; }
  void SetMaxPt(Double_t d)                                 { fMaxPt = d; }
  void SetUseAliEventCuts(Bool_t b)                         { fUseAliEventCuts = b; }
  void SetUseManualEvtCuts(Bool_t input)                    { fUseManualEventCuts = input;}
  void SetPlotJetHistograms(Bool_t b)                       { fPlotJetHistograms = b; }
  void SetPlotClusterHistograms(Bool_t b)                   { fPlotClusterHistograms = b; }
  void SetPlotParticleCompositionHistograms(Bool_t b)       { fPlotParticleCompositionHistograms = b; }
  void SetComputeBackground(Bool_t b)                       { fComputeBackground = b; }
  void SetDoTriggerSimulation(Bool_t b)                     { fDoTriggerSimulation = b; }
  void SetPlotMatchedJetHistograms(Bool_t b)                { fPlotMatchedJetHistograms = b; }
  void SetComputeMBDownscaling(Bool_t b)                    { fComputeMBDownscaling = b; }
  void SetTrackMatchingDeltaEtaMax(Double_t deta)           { fTrackMatchingDeltaEtaMax = deta; }
  void SetTrackMatchingDeltaPhiMax(Double_t dphi)           { fTrackMatchingDeltaPhiMax = dphi; }
  void SetUseResponseMaker(Bool_t b)                        { fUseResponseMaker = b; }
  void SetPlotDCal(Bool_t b)                                { fPlotDCal = b; }

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;
  Bool_t                      IsEventSelected()                                 ;
  void                        RunChanged(Int_t run)                             ;

  // Analysis and plotting functions
  void                        GenerateHistoBins()                               ;
  void                        AllocateJetHistograms()                           ;
  void                        AllocateClusterHistograms()                       ;
  void                        AllocateParticleCompositionHistograms()           ;
  void                        AllocateBackgroundHistograms()                    ;
  void                        AllocateTriggerSimHistograms()                    ;
  void                        AllocateMatchedJetHistograms()                    ;
  void                        FillJetHistograms()                               ;
  void                        FillClusterHistograms()                           ;
  void                        FillParticleCompositionHistograms()               ;
  void                        FillParticleCompositionClusterHistograms(const AliMCEvent* mcevent);
  void                        FillParticleCompositionJetHistograms(const AliMCEvent* mcevent);
  void                        SetParticleTypeLabels(TAxis* axis)                ;
  void                        ComputeBackground()                               ;
  void                        DoTriggerSimulation()                             ;
  void                        FillTriggerSimHistograms()                        ;
  void                        FillMatchedJetHistograms()                        ;
  
  // Utility functions
  Double_t                    GetJetPt(const AliEmcalJet* jet, Double_t rho);
  Double_t                    GetDeltaR(const AliTLorentzVector* part, Double_t etaRef, Double_t phiRef);
  Double_t                    GetJetType(const AliEmcalJet* jet);
  ContributorType             GetContributorType(const AliVCluster* clus, const AliMCEvent* mcevent, Int_t label);
  Bool_t                      IsHadron(const ContributorType contributor);
  
  // Analysis parameters
  Bool_t                      fPlotJetHistograms;                   ///< Set whether to enable inclusive jet histograms
  Bool_t                      fPlotClusterHistograms;               ///< Set whether to plot cluster histograms
  Bool_t                      fPlotParticleCompositionHistograms;   ///< Set whether to plot jet composition histograms
  Bool_t                      fComputeBackground;                   ///< Set whether to enable study of background
  Bool_t                      fDoTriggerSimulation;                 ///< Set whether to perform a simple trigger simulation
  Bool_t                      fPlotMatchedJetHistograms;            ///< Set whether to plot matched jet histograms (must run ResponseMaker first)
  Bool_t                      fComputeMBDownscaling;                ///< Set whether to compute and plot MB downscaling factors
  Bool_t                      fPlotDCal;                            ///< Set whether to enable several DCal-specific histograms
  
  // Plotting parameters
  Double_t                    fMinPt;                               ///< Histogram min pT limit
  Double_t                    fMaxPt;                               ///< Histogram max pT limit
  Int_t                       fNEtaBins;                            ///< Number of eta bins
  Int_t                       fNPhiBins;                            ///< Number of phi bins
  Int_t                       fNCentHistBins;                       //!<! number of cent bins
  Double_t*                   fCentHistBins;                        //!<! cent bins
  Int_t                       fNPtHistBins;                         //!<! number of variable pt bins
  Double_t*                   fPtHistBins;                          //!<! variable pt bins
  Int_t                       fNM02HistBins;                        //!<! number of variable M02 bins
  Double_t*                   fM02HistBins;                         //!<! variable M02 bins
  Int_t                       fNEoverPBins;                         //!<! number of variable E/p bins
  Double_t*                   fEoverPBins;                          //!<! variable E/p bins
  
  // Track matching parameters (for cluster histogram plots)
  Double_t                    fTrackMatchingDeltaEtaMax;            ///< Maximum delta-eta to consider a track to be matched to a cluster
  Double_t                    fTrackMatchingDeltaPhiMax;            ///< Maximum delta-phi to consider a track to be matched to a cluster
  
  // Trigger parameters
  Double_t                    fMBUpscaleFactor;                     //!<! inverse of downscale factor, for MB trigger
  Double_t                    fMedianEMCal;                         //!<! median patch energy in EMCal, per event
  Double_t                    fMedianDCal;                          //!<! median patch energy in DCal, per event
  Bool_t                      fkEMCEJE;                             //!<! flag telling whether the event is "triggered" or not in "simulation"
  
  // Embedding parameters
  AliEmcalEmbeddingQA         fEmbeddingQA;                         //!<! QA hists for embedding (will only be added if embedding)
  Bool_t                      fUseResponseMaker;                    ///< Flag to use Response Maker rather than JetTagger
  AliJetContainer*            fMCJetContainer;                      //!<!Pointer to jet container of truth-level jets
  
  // Event selection
  Bool_t                      fUseAliEventCuts;                     ///< Flag to use AliEventCuts (otherwise AliAnalysisTaskEmcal will be used)
  AliEventCuts                fEventCuts;                           ///< event selection utility
  TList                      *fEventCutList;                        //!<! Output list for event cut histograms
  Bool_t                      fUseManualEventCuts;                  ///< Flag to use manual event cuts
  
  // MC options
  AliMCParticleContainer*     fGeneratorLevel;                      //!<! generator level container
  
  // Hist manager
  THistManager                fHistManager;                         ///< Histogram manager

 private:
  AliAnalysisTaskEmcalJetPerformance(const AliAnalysisTaskEmcalJetPerformance&)           ; // not implemented
  AliAnalysisTaskEmcalJetPerformance &operator=(const AliAnalysisTaskEmcalJetPerformance&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetPerformance, 14);
  /// \endcond
};
#endif
