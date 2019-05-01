#ifndef AliAnalysisTaskEmcalVsPhos_H
#define AliAnalysisTaskEmcalVsPhos_H

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
 * \file AliAnalysisTaskEmcalVsPhos.h
 * \brief Study of EMCal vs. PHOS clusters.
 *
 * This header file declares the class AliAnalysisTaskEmcalVsPhos.
 *
 * \author James Mulligan <james.mulligan@yale.edu>, Yale University
 * \date Aug 15, 2017
 */

class AliPHOSGeometry;

#include "THistManager.h"

#include "AliEventCuts.h"
#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalVsPhos : public AliAnalysisTaskEmcalJet {
 public:
  
  enum ClusterType {
    kNA       = -1,//!< Undefined
    kEMCal    = 0, //!< EMCal
    kDCal     = 1, //!< DCal
    kPHOS     = 2  //!< PHOS
  };

  AliAnalysisTaskEmcalVsPhos()                                          ;
  AliAnalysisTaskEmcalVsPhos(const char *name)                          ;
  virtual ~AliAnalysisTaskEmcalVsPhos()                                 ;

  void UserCreateOutputObjects()                                        ;
  
  // Setters
  void SetMaxPt(Double_t d)                                 { fMaxPt = d; }
  void SetUseAliEventCuts(Bool_t b)                         { fUseAliEventCuts = b; }
  void SetUseManualEvtCuts(Bool_t input)                    { fUseManualEventCuts = input;}
  void SetPlotNeutralJets(Bool_t b)                         { fPlotNeutralJets = b; }
  void SetPlotClustersInJets(Bool_t b)                      { fPlotClustersInJets = b; }
  void SetPlotClusterHistograms(Bool_t b)                   { fPlotClusterHistograms = b; }
  void SetPlotCellHistograms(Bool_t b)                      { fPlotCellHistograms = b; }
  void SetPlotClusWithoutNonLinCorr(Bool_t b)               { fPlotClusWithoutNonLinCorr = b; }
  void SetPlotExotics(Bool_t b)                             { fPlotExotics = b; }
  void SetPlotStandardClusterTHnSparse(Bool_t b)            { fPlotStandardClusterTHnSparse = b; }
  void SetPlotNearestNeighborDistribution(Bool_t b)         { fPlotNearestNeighborDistribution = b; }
  void SetPlotClusterCone(Bool_t b)                         { fPlotClusterCone = b; }
  void SetPlotCaloCentrality(Bool_t b)                      { fPlotCaloCentrality = b; }
  void SetPlotFineGrainedEtaPhi(Bool_t b)                   { fPlotFineGrainedEtaPhi = b; }
  void SetPlotEvenOddEta(Bool_t b)                          { fPlotEvenOddEta = b; }
  void SetPlotCellSMDensity(Bool_t b)                       { fPlotCellSMDensity = b; }
  void SetExcludeRejectedCells(Bool_t b)                    { fExcludeRejectedCells = b; }
  void SetPlotFineGrainedCentrality(Bool_t b)               { fPlotFineGrainedCentrality = b; }
  void SetPlotEventHistograms(Bool_t b)                     { fPlotEventHistograms = b; }

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;
  Bool_t                      IsEventSelected()                                 ;

  // Analysis and plotting functions
  void                        GenerateHistoBins()                               ;
  void                        AllocateCaloHistograms()                          ;
  void                        AllocateClusterHistograms()                       ;
  void                        AllocateCellHistograms()                          ;
  void                        AllocateNeutralJetHistograms()                    ;
  void                        AllocateClustersInJetsHistograms()                ;
  void                        AllocateEventHistograms()                         ;
  void                        FillCaloHistograms()                              ;
  void                        FillClusterHistograms()                           ;
  void                        FillCellHistograms()                              ;
  void                        FillNeutralJetHistograms()                        ;
  void                        FillClustersInJetsHistograms()                    ;
  void                        FillEventHistograms()                             ;
  void                        FillClusterTHnSparse(TString clustersName, Double_t eta, Double_t phi, Double_t Enonlin, Double_t Ehadcorr, Int_t hasMatchedTrack, Double_t M02, Int_t nCells, Int_t passDispersionCut, Double_t distNN, Int_t isOddEta, Int_t coneType = 0, Double_t R = 0., Double_t Econe = 0.);
  void                        FillClusterTHnSparse(TString clustersName, Double_t eta, Double_t phi, Double_t Enonlin, Double_t eCellCone, Double_t eCellSM, Int_t nCellsCone, Int_t nCellsSM);
  
  // Utility functions
  Double_t                    GetJetPt(const AliEmcalJet* jet, Double_t rho);
  Double_t                    GetDeltaR(AliTLorentzVector part, Double_t etaRef, Double_t phiRef);
  Double_t                    GetDeltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
  Double_t                    GetJetType(const AliEmcalJet* jet);
  Double_t                    GetFcross(const AliVCluster *cluster, AliVCaloCells *cells);
  Double_t                    FindNearestNeighborDistance(AliTLorentzVector cluster);
  Double_t                    GetConeClusterEnergy(Double_t etaRef, Double_t phiRef, Double_t R);
  Double_t                    GetConeCellEnergy(Double_t etaRef, Double_t phiRef, Double_t R, Bool_t returnNcells = kFALSE);
  Double_t                    GetSMCellEnergy(Int_t sm, Int_t clusType, Bool_t returnNcells = kFALSE);
  Bool_t                      IsCellRejected(Int_t absId, Int_t cellType);

  // Analysis configuration and plotting options
  Bool_t                      fPlotClusterHistograms;               ///< Set whether to plot cluster histograms
  Bool_t                      fPlotNeutralJets;                     ///< Set whether to plot neutral jet histo
  Bool_t                      fPlotClustersInJets;                  ///< Set whether to plot histogram of clusters within jets
  Bool_t                      fPlotCellHistograms;                  ///< Set whether to plot cell histograms
  Bool_t                      fPlotClusWithoutNonLinCorr;           ///< If true, use pre-nonlincorr energy in cluster thnsparse
  Bool_t                      fPlotExotics;                         ///< Set whether to plot exotic cluster study
  Bool_t                      fPlotStandardClusterTHnSparse;        ///< Set whether to plot "standard" axes in cluster THnSparse
  Bool_t                      fPlotNearestNeighborDistribution;     ///< Set whether to plot nearest neighbor axis in cluster THnSparse
  Bool_t                      fPlotClusterCone;                     ///< Set whether to plot sum of energy surrounding cluster in THnSparse
  Bool_t                      fPlotCaloCentrality;                  ///< Set whether to bin cluster THnSparse in calorimeter local density
  Bool_t                      fPlotFineGrainedEtaPhi;               ///< Set whether to plot fine-grained eta-phi bins in cluster THnSparse
  Bool_t                      fPlotEvenOddEta;                      ///< Set whether to add axis to THnSparse separating even/odd eta columns
  Bool_t                      fPlotCellSMDensity;                   ///< Set whether to plot SM cell density when computing local density
  Bool_t                      fExcludeRejectedCells;                ///< Set whether to exclude cells from rejected clusters in cone/SM studies
  Bool_t                      fPlotFineGrainedCentrality;           ///< Set whether to plot a more fine grained centrality binning
  Bool_t                      fPlotEventHistograms;                 ///< Set whether to plot some calo event histograms

  // Plotting parameters
  Float_t                     fMaxPt;                               ///< Histogram pt limit
  Int_t                       fNCentHistBins;                       //!<! number of cent bins
  Double_t*                   fCentHistBins;                        //!<! cent bins
  Int_t                       fNPtHistBins;                         //!<! number of variable pt bins
  Double_t*                   fPtHistBins;                          //!<! variable pt bins
  Int_t                       fNM02HistBins;                        //!<! number of variable M02 bins
  Double_t*                   fM02HistBins;                         //!<! variable M02 bins
  
  // Event selection
  Bool_t                      fUseAliEventCuts;                     ///< Flag to use AliEventCuts (otherwise AliAnalysisTaskEmcal will be used)
  AliEventCuts                fEventCuts;                           ///< event selection utility
  TList                      *fEventCutList;                        //!<! Output list for event cut histograms
  Bool_t                      fUseManualEventCuts;                  ///< Flag to use manual event cuts
  
  // Phos geometry (only needed for cluster studies)
  AliPHOSGeometry*            fPHOSGeo;                             //!<! phos geometry
  
  // Hist manager
  THistManager                fHistManager;                         ///< Histogram manager

 private:
  AliAnalysisTaskEmcalVsPhos(const AliAnalysisTaskEmcalVsPhos&)           ; // not implemented
  AliAnalysisTaskEmcalVsPhos &operator=(const AliAnalysisTaskEmcalVsPhos&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalVsPhos, 13);
  /// \endcond
};
#endif
