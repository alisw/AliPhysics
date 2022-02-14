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
#ifndef ALIEMCALJETTAGGERTASKFAST_H
#define ALIEMCALJETTAGGERTASKFAST_H

//#define JETTAGGERFAST_TEST

class TH1;
class TH2;
class TH3;
class AliJetContainer;

#include "AliAnalysisTaskEmcalJet.h"

namespace PWGJE {
namespace EMCALJetTasks {

/**
 * @class AliEmcalJetTaggerTaskFast
 * @brief Fast jet tagger for geometric matching of jets
 * @ingroup PWGJEBASE
 * @author Martha Verweij
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Nov 8, 2017
 *
 * Class based on AliAnalysisTaskEmcalJetTagger. Navigation finding closest neighbor
 * however is based on a kd-tree.
 *
 */
class AliEmcalJetTaggerTaskFast : public AliAnalysisTaskEmcalJet {
 public:
  enum JetTaggingMethod {
    kGeo      = 0,
    kFraction = 1
  };

  enum JetTaggingType {
    kTag      = 0,
    kClosest  = 1
  };
  /**
   * @enum AcceptanceType
   * @brief Accpetance type used for the two jet containers
   */
  enum AcceptanceType {
    kNoLimit = 0,             ///< No Additional limit compared to jet containers
    kLimitTagEta = 1,         ///< Adding 0.1 in eta limit for tag jets
    kLimitTagEtaPhi = 2,      ///< Adding 0.1 in eta and phi limits for tag jets
    kLimitBaseTagEtaPhi = 3   ///< Adding 0.1 in eta and phi limits for both base and tag jets
  };

  /**
   * @brief Default constructor
   */
  AliEmcalJetTaggerTaskFast();

  /**
   * @brief Standard constructor.
   */
  AliEmcalJetTaggerTaskFast(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalJetTaggerTaskFast() {}

  void                                UserCreateOutputObjects();

  //Setters
  void SetJetContainerBase(Int_t c)                             { fContainerBase = c;}
  void SetJetContainerTag(Int_t c)                              { fContainerTag  = c;}

  void SetJetTaggingType(JetTaggingType t)                      { fJetTaggingType = t;}
  void SetJetTaggingMethod(JetTaggingMethod m)                  { fJetTaggingMethod = m;}

  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f; }

  void SetUseSumw2(Bool_t b)                                    { fUseSumw2 = b;}
  
  void SetTypeAcceptance(AcceptanceType type)                   { fTypeAcc = type; }
  void SetMaxDistance(Double_t dist)                            { fMaxDist = dist; }
  void SetSpecialParticleContainer(Int_t contnumb)              { fSpecPartContTag = contnumb; }


  /**
   * @brief Factory creating new jet matching task
   *
   * The new task is added to the analysis manager, and the containers for
   * input and output are configured. Requires an analysis manager to be
   * present.
   *
   * @param[in] njetsBase Name of the jet container with base jets
   * @param[in] njetsTag
   * @param[in] R resolution parameter
   * @param[in] nrhoBase Name of the \f$\rho\f$-parameter object for base jets
   * @param[in] nrhoTag Name of the \f$\rho\f$-parameter object for tag jets
   * @param[in] ntracks Name of the track container
   * @param[in] nclusters Name of the cluster container
   * @param[in] type Name of the acceptance type (TPC/EMCAL) for base jets
   * @param[in] CentEst Name of the centrality estimator
   * @param[in] pSel Event selection bitmap
   * @param[in] trigClass Name of the trigger class (for output container)
   * @return New jet tagger task with basic configuration
   */
  static AliEmcalJetTaggerTaskFast *AddTaskJetTaggerFast(const char * njetsBase,
      const char * njetsTag,
      const Double_t R,
      const char * nrhoBase,
      const char * nrhoTag,
      const char * ntracks,
      const char * nclusters,
      const char * type,
      const char * CentEst,
      Int_t        pSel,
      const char * trigClass);

 protected:

  /**
   * @brief Run matching
   *
   * In case the task is not yet initialized, initializing first.
   * Before matching the current match status of the jets in both
   * the base and the tag container are reset. For the matching
   * see \ref MatchJetsGeo
   *
   * @return Always true
   */
  Bool_t                              Run();

  /**
   * @brief Filling QA histograms monitoring the quality of the matching
   * @return Always true
   */
  Bool_t                              FillHistograms();

  /**
   * @brief Initializing the matching task by defining the accepted range on both jet containers.
   *
   * Function will not be executed twice.
   */
  void                                Init();

  /**
   * @brief Calculate azimuthal angle between the axes of the jets
   * @param[in] jet1 base jet
   * @param[in] jet2 jet to test
   */
  Double_t GetDeltaPhi(const AliEmcalJet* jet1, const AliEmcalJet* jet2);

  /**
   * @brief Calculate azimuthal angle between the axises of the jets
   * @param[in] phi1 \f$\phi\f$-angle of the first jet
   * @param[in] phi2 \f$\phi\f$-angle of the second jet
   */
  Double_t GetDeltaPhi(Double_t phi1,Double_t phi2);

  /**
   * @brief Match the full jets to the corresponding charged jets
   *
   * For all jets, at both base and tag level, finding the nearest neighbor
   * in distance in the \$\eta\f$-\f$\phi\$ space, accepting only pairs
   * with a distance smaller maxDistance. True jet pairs are accepted only
   * if the base jet is the closest neighbor to the tag jet and vice versa
   * at the same time.
   *
   * @param[in] contBase Container with base jets
   * @param[in] contTag Container with jets to be tagged
   * @param[in] maxDistance Maximum distance allowed in order to accept a pair tag
   */
  Bool_t     MatchJetsGeo(AliJetContainer &contBase, AliJetContainer &contTag, Float_t maxDist = 0.3) const;

  /**
   * @brief Reset tagging for all jets in jet container
   * @param[in] cont Jet container for which to reset the tagging status
   */

  void     ResetTagging(const AliJetContainer &cont) const;
  
 private:
  JetTaggingType                      fJetTaggingType;             ///< jet matching type
  JetTaggingMethod                    fJetTaggingMethod;           ///< jet matching method
  Int_t                               fContainerBase;              ///< jets to be tagged
  Int_t                               fContainerTag;               ///< jets used for tagging
  Int_t                               fSpecPartContTag;            ///< particle container optionally used in AliJetContainer::GetFractionSharedPt(). Set only if needed.
  Double_t                            fMinFractionShared;          ///< only fill histos for jets if shared fraction larger than X
  Bool_t                              fUseSumw2;                   ///< activate sumw2 for output histograms
  Bool_t                              fMatchingDone;               ///< flag to indicate if matching is done or not
  AcceptanceType                      fTypeAcc;                    ///< acceptance cut for the jet containers, see method MatchJetsGeo in .cxx for possibilities
  Double_t                            fMaxDist;                    ///< distance allowed for two jets to match
  Bool_t                              fInit;                       ///< true when the containers are initialized
  TH3            **fh3PtJet1VsDeltaEtaDeltaPhi;  //!<! \f$ p_{t}\f$ jet 1 vs deta vs dphi
  TH2            **fh2PtJet1VsDeltaR;            //!<! \f$ p_{t}\f$ jet 1 vs dR
  TH2            **fh2PtJet2VsFraction;          //!<! \f$ p_{t}\f$ jet 1 vs shared fraction

  TH2            **fh2PtJet1VsLeadPtAllSel;      //!<! all jets after std selection
  TH2            **fh2PtJet1VsLeadPtTagged;      //!<! tagged jets
  TH2            **fh2PtJet1VsPtJet2;            //!<! pT of base jet vs tagged jet
  TH2            **fh2PtJet2VsRelPt;             //!<! pT of tagged jet vs pt base jet / pt tagged jet
  
  TH3             *fh3PtJetDEtaDPhiConst;        //!<! \f$ p_{t}\f$ jet vs delta eta vs delta phi of constituents
  TH3             *fh3PtJetAreaDRConst;          //!<! \f$ p_{t}\f$ jet vs Area vs delta R of constituents
  TH1             *fNAccJets;                    //!<! number of jets per event
#ifdef JETTAGGERFAST_TEST
  TH1             *fIndexErrorRateBase;          //!<! Monitoring number of errors between index in kd-tree and data block for base jets
  TH1             *fIndexErrorRateTag;           //!<! Monitoring number of errors between index in kd-tree and data block for tag jets
  TH1             *fContainerErrorRateBase;      //!<! Monitoring number of errors between index in kd-tree and jet vector for base jets
  TH1             *fContainerErrorRateTag;       //!<! Monitoring number of errors between index in kd-tree and jet vector for tag jets
#endif
  AliEmcalJetTaggerTaskFast(const AliEmcalJetTaggerTaskFast&);            // not implemented
  AliEmcalJetTaggerTaskFast &operator=(const AliEmcalJetTaggerTaskFast&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEmcalJetTaggerTaskFast, 2);
  /// \endcond
};
}
}

#endif

