/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIANALYSISTASKRHODEV_H
#define ALIANALYSISTASKRHODEV_H

#include <utility>

#include "AliAnalysisTaskRhoBaseDev.h"

/** 
 * @class AliAnalysisTaskRhoDev
 * @brief Class for a task that calculates the UE
 * @ingroup PWGJEBASE
 * @author Rosi Reed, Yale University
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * @date June 16, 2017
 *
 * Class for a task that calculates the average background
 * coming from the underlying event (UE) in jet analysis.
 * This task calculates the average background as the median
 * of the pt density of kt clusters. More details at: https://arxiv.org/pdf/0707.1378.pdf.
 * If scale function is given the scaled rho will be exported
 * with the name as "fOutRhoName".Apppend("_Scaled").
 * This is a development version. The stable version of this class
 * is AliAnalysisTaskRho.
 */
class AliAnalysisTaskRhoDev : public AliAnalysisTaskRhoBaseDev {

public:
  /**
   * @brief Default constructor. Needed by ROOT I/O
   */
  AliAnalysisTaskRhoDev();

  /**
   * Standard constructor. Should be used by the user.
   *
   * @param[in] name  Name of the task
   * @param[in] histo If kTRUE, the task will also produce QA histograms
   */
  AliAnalysisTaskRhoDev(const char *name, Bool_t histo=kFALSE);
  
  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskRhoDev() {}

  /**
   * Performing run-independent initialization.
   * Here the histograms should be instantiated.
   */
  void             UserCreateOutputObjects();

  void             SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }
  void             SetRhoSparse(Bool_t b)          { fRhoSparse     = b    ; }
  void             SetExclJetOverlap(TString n)    { fExclJetOverlap= n    ; }

  /**
   * @brief Create an instance of this class and add it to the analysis manager
   * @param trackName name of the track collection
   * @param trackPtCut minimum pt of the tracks
   * @param clusName name of the calorimeter cluster collection
   * @param clusECut minimum energy of the calorimeter clustuers
   * @param nRho name of the output rho object
   * @param jetradius Radius of the kt jets used to calculate the background
   * @param acceptance Fiducial acceptance of the kt jets
   * @param jetType Jet type (full/charged)
   * @param rscheme Recombination scheme
   * @param histo If kTRUE the task will also produce QA histograms
   * @param suffix additional suffix that can be added at the end of the task name
   * @return pointer to the new AliAnalysisTaskRhoDev task
   */
  static AliAnalysisTaskRhoDev* AddTaskRhoDev(
     TString        nTracks                        = "usedefault",
     Double_t       trackPtCut                     = 0.15,
     TString        nClusters                      = "usedefault",
     Double_t       clusECut                       = 0.30,
     TString        nRho                           = "Rho",
     Double_t       jetradius                      = 0.2,
     UInt_t         acceptance                     = AliEmcalJet::kTPCfid,
     AliJetContainer::EJetType_t jetType           = AliJetContainer::kChargedJet,
     AliJetContainer::ERecoScheme_t rscheme        = AliJetContainer::pt_scheme,
     Bool_t         histo                          = kTRUE,
     TString        suffix                         = ""
  );

 protected:
  /**
   * Calculates the average background using the median approach
   * as proposed in https://arxiv.org/pdf/0707.1378.pdf.
   * Rho is stored in fOutRho.
   */
  void          CalculateRho();

  /**
   * Fill histograms.
   */
  Bool_t        FillHistograms();

  /**
   * @brief Verify that the required particle, cluster and jet containers were provided.
   * @return kTRUE if all requirements are satisfied, kFALSE otherwise
   */
  Bool_t        VerifyContainers();

  /**
   * Finds the first two leading jets
   * @return A pair with the leading and sub-leading jets respectively as first and second element
   */
  std::pair<AliEmcalJet*, AliEmcalJet*>
                GetLeadingJets();

  UInt_t           fNExclLeadJets;                 ///< number of leading jets to be excluded from the median calculation
  Bool_t           fRhoSparse;                     ///< flag to run CMS method as described in https://arxiv.org/abs/1207.2392
  TString          fExclJetOverlap;                ///< name of the jet collection that should be used to reject jets that are considered "signal"

  Double_t         fOccupancyFactor;               //!<!occupancy correction factor for sparse events
  TH2F            *fHistOccCorrvsCent;             //!<!occupancy correction vs. centrality

  AliAnalysisTaskRhoDev(const AliAnalysisTaskRhoDev&);             // not implemented
  AliAnalysisTaskRhoDev& operator=(const AliAnalysisTaskRhoDev&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoDev, 2);
};
#endif
