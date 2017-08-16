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
  void                        FillCaloHistograms()                              ;
  void                        FillClusterHistograms()                           ;
  void                        FillCellHistograms()                              ;
  void                        FillNeutralJetHistograms()                        ;
  void                        FillClustersInJetsHistograms()                    ;
  
  // Utility functions
  Double_t                    GetJetPt(AliJetContainer* jetCont, AliEmcalJet* jet);
  Double_t                    GetDeltaR(AliTLorentzVector* part, Double_t etaRef, Double_t phiRef);
  Double_t                    GetJetType(AliEmcalJet* jet);
  Double_t                    GetFcross(AliVCluster *cluster, AliVCaloCells *cells);
  Double_t                    FindNearestNeighborDistance(AliTLorentzVector cluster, AliClusterContainer* clusters);

  // Analysis configuration and plotting options
  Bool_t                      fPlotClusterHistograms;               ///< Set whether to plot cluster histograms
  Bool_t                      fPlotNeutralJets;                     ///< Set whether to plot neutral jet histo
  Bool_t                      fPlotClustersInJets;                  ///< Set whether to plot histogram of clusters within jets
  Bool_t                      fPlotCellHistograms;                  ///< Set whether to plot cell histograms
  Bool_t                      fPlotClusWithoutNonLinCorr;           ///< If true, use pre-nonlincorr energy in cluster thnsparse
  Bool_t                      fPlotExotics;                         ///< Set whether to plot exotic cluster study
  Bool_t                      fPlotStandardClusterTHnSparse;        ///< Set whether to plot "standard" axes in cluster THnSparse
  Bool_t                      fPlotNearestNeighborDistribution;     ///< Set whether to plot nearest neighbor axis in cluster THnSparse


  // Plotting parameters
  Float_t                     fMaxPt;                               ///< Histogram pt limit
  Int_t                       fNCentHistBins;                       //!<! number of cent bins
  Double_t*                   fCentHistBins;                        //!<! cent bins
  Int_t                       fNPtHistBins;                         //!<! number of variable pt bins
  Double_t*                   fPtHistBins;                          //!<! variable pt bins
  
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
  ClassDef(AliAnalysisTaskEmcalVsPhos, 3);
  /// \endcond
};
#endif
