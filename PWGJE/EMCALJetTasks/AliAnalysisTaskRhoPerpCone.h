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
#ifndef ALIANALYSISTASKRHOPERPCONE_H
#define ALIANALYSISTASKRHOPERPCONE_H

#include "AliAnalysisTaskRhoBase.h"

/** 
 * @class AliAnalysisTaskRhoPerpCone
 * @brief Class for a task that calculates the UE
 * @ingroup PWGJEBASE
 *
 * Class for a task that calculates the average background
 * coming from the underlying event (UE) in jet analysis.
 * In each event, tha task finds the leading jet and sum
 * the particle pT in a cone of R = 0.4 perpendicular
 * to the leading jet axis. 
 */
class AliAnalysisTaskRhoPerpCone : public AliAnalysisTaskRhoBase
{

public:
    /**
   * Default constructor. Needed by ROOT I/O
   */
    AliAnalysisTaskRhoPerpCone();

    /**
   * Standard constructor. Should be used by the user.
   *
   * @param[in] name  Name of the task
   * @param[in] histo If kTRUE, the task will also produce QA histograms
   */
    AliAnalysisTaskRhoPerpCone(const char *name, Bool_t histo = kFALSE);

    /**
   * @brief Destructor
   */
    virtual ~AliAnalysisTaskRhoPerpCone() {}

    /**
   * Performing run-independent initialization.
   * Here the histograms should be instantiated.
   */
    void UserCreateOutputObjects();

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
  static AliAnalysisTaskRhoPerpCone* AddTaskRhoPerpCone(
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
   * @brief Run the analysis.
   * @return always true
   */
    Bool_t Run();

    TH2 *fHistRhoVsCent;        //!<!rho vs. centrality
    TH2 *fHistRhoVsLeadJetPt;   //!<!rho vs. leading jet pt
    TH2 *fHistRhoVsLeadTrackPt; //!<!rho vs. leading track pt
    TH2 *fHistRhoVsNtrack;      //!<!rho vs. no. of tracks

    AliAnalysisTaskRhoPerpCone(const AliAnalysisTaskRhoPerpCone &);            // not implemented
    AliAnalysisTaskRhoPerpCone &operator=(const AliAnalysisTaskRhoPerpCone &); // not implemented

    ClassDef(AliAnalysisTaskRhoPerpCone, 1);
};
#endif
