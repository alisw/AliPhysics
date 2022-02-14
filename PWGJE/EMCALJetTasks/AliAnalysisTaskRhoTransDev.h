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
#ifndef ALIANALYSISTASKRHOTRANSDEV_H
#define ALIANALYSISTASKRHOTRANSDEV_H

#include <map>

#include "AliAnalysisTaskRhoBaseDev.h"

/** 
 * @class AliAnalysisTaskRhoTransDev
 * @brief Class for a task that calculates the UE
 * @ingroup PWGJEBASE
 * @author Rosi Reed, Yale University
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * @date June 25, 2017
 *
 * Class for a task that calculates the average background
 * coming from the underlying event (UE) in jet analysis.
 * In each event, tha task finds the leading jet and sum
 * the particle pT in a cone of R = 0.4 perpendicular
 * to the leading jet axis. In addition another histogram is filled
 * with the UE estimated in events that fulfill the back-to-back
 * condition: a back-to-back jet, satisfying |∆φ| > 5/6π and
 * carrying at least a given fraction (default: 60%)
 * of the transverse momentum of the leading jet, is required
 * to be present in the event. Additionally, all the other jets in the event
 * should have pT smaller than a certain threshold (default: 12 GeV/c).
 * This implements the techniques described in the analysis note
 * of the ALICE inclusive jet spectrum in pp collisions at 2.76 TeV,
 * section 7.3: https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/rma/2013-Mar-29-analysis_note-ppJet_note.pdf
 * If scale function is given the scaled rho will be exported
 * with the name as "fOutRhoName".Apppend("_Scaled").
 */
class AliAnalysisTaskRhoTransDev : public AliAnalysisTaskRhoBaseDev {

 public:

  /**
   * Default constructor. Needed by ROOT I/O
   */
  AliAnalysisTaskRhoTransDev();

  /**
   * Standard constructor. Should be used by the user.
   *
   * @param[in] name  Name of the task
   * @param[in] histo If kTRUE, the task will also produce QA histograms
   */
  AliAnalysisTaskRhoTransDev(const char *name, Bool_t histo=kFALSE);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskRhoTransDev() {}

  /**
   * Performing run-independent initialization.
   * Here the histograms should be instantiated.
   */
  void             UserCreateOutputObjects();

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
  static AliAnalysisTaskRhoTransDev* AddTaskRhoTransDev(
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
   * Calculates the average background by summing the transverse momenta
   * of the particles perpendicular to the leading jet axis.
   * This implements the techniques described in the analysis note
   * of the ALICE inclusive jet spectrum in pp collisions at 2.76 TeV,
   * section 7.3: https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/rma/2013-Mar-29-analysis_note-ppJet_note.pdf
   * The background densisty is stored in fOutRho.
   */
  void          CalculateRho();

  /**
   * @brief Fill histograms.
   */
  Bool_t        FillHistograms();

  /**
   * @brief Verify that the required particle, cluster and jet containers were provided.
   * @return kTRUE if all requirements are satisfied, kFALSE otherwise
   */
  Bool_t        VerifyContainers();

  /**
   * Calculates the transervse momentum density in a region that
   * is perpedincular to the leading jet. The perpedincular region
   * spans over an angle pi/4 around pi/4 and -pi/4.
   * This implements the techniques described in the analysis note
   * of the ALICE inclusive jet spectrum in pp collisions at 2.76 TeV,
   * section 7.3: https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/rma/2013-Mar-29-analysis_note-ppJet_note.pdf
   * @param cont Collection of objects (tracks, clusters or particles) used to calculate the momentum density
   * @param leadingJet Leading jet used to define the perpendicular region
   * @return The perpendicular momentum density
   */
  Double_t      GetPerpPtDensity(AliEmcalContainer* cont, AliVParticle* leadingJet);

  TH2                                *fHistB2BRhoVsCent;                 //!<!rho vs. centrality

  std::map<std::string, TH2*>         fHistB2BRhoVsLeadJetPt;            //!<!rho vs. leading jet pt
  TH2                                *fHistB2BRhoVsLeadTrackPt;          //!<!rho vs. leading track pt
  TH2                                *fHistB2BRhoVsNtrack;               //!<!rho vs. no. of tracks
  TH2                                *fHistB2BRhoVsLeadClusterE;         //!<!rho vs. leading cluster energy
  TH2                                *fHistB2BRhoVsNcluster;             //!<!rho vs. no. of clusters
  TH2                                *fHistB2BRhoScaledVsCent;           //!<!rhoscaled vs. centrality
  TH2                                *fHistB2BRhoScaledVsNtrack;         //!<!rhoscaled vs. no. of tracks
  TH2                                *fHistB2BRhoScaledVsNcluster;       //!<!rhoscaled vs. no. of clusters

  AliAnalysisTaskRhoTransDev(const AliAnalysisTaskRhoTransDev&);             // not implemented
  AliAnalysisTaskRhoTransDev& operator=(const AliAnalysisTaskRhoTransDev&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoTransDev, 1);
};
#endif
