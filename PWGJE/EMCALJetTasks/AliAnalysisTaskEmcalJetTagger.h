/************************************************************************************
 * Copyright (C) 2013, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIANALYSISTASKEMCALJETTAGGER_H
#define ALIANALYSISTASKEMCALJETTAGGER_H

class TH1;
class TH2;
class TH3;
class AliJetContainer;

#include "AliAnalysisTaskEmcalJet.h"

/**
 * @class AliAnalysisTaskEmcalJetTagger
 * @brief Tagging jets  with jet from another source
 * @ingroup PWGJEBASE
 * @author Martha Verweij
 * @since Nov 29th, 2013
 * 
 * Task tagging jets from one source (the base jet container) with
 * jets from a second source (the tag jet container). Tag jets can be
 * i.e. jets at particle level for base jets at detector level. Two methods
 * can be applied:
 *   * Pure geometric matching, requiring that matched jets have a distance
 *     smaller than a maximum distance
 *   * Fractional particle matching: In addition to geometric matching the
 *     matched jets must also share a minimum amount of particles
 * As matching result either the closest jet or all jets satisfying the matching
 * criteria get matched. 
 * 
 * The matching information is stored in the AliEmcalJet object. The matched jets
 * can be obtained via
 *   - AliEmcalJet::ClosestJet() for the closest jet tagging type
 *   - AliEmcalJet::MatchedJet() for the tag jet tagging type
 */
class AliAnalysisTaskEmcalJetTagger : public AliAnalysisTaskEmcalJet {
 public:
  /**
   * @enum JetTaggingMethod
   * @brief Method used to tag jets from different sources
   */
  enum JetTaggingMethod {
    kGeo      = 0,      ///< Pure geometric tagging (max. distance + acceptance cut)
    kFraction = 1       ///< Geometric tagging + fractional of overlap in particles
  };

  /**
   * @enum JetTaggingType
   * @brief Type of the tagging from different sources
   */
  enum JetTaggingType {
    kTag      = 0,      ///< All jets within a maximum distance are tagged
    kClosest  = 1       ///< Only the closest jet in distance is tagged
  };

  /**
   * @brief Dummy constructor
   */
  AliAnalysisTaskEmcalJetTagger();

  /**
   * @brief Default constructor
   */
  AliAnalysisTaskEmcalJetTagger(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskEmcalJetTagger();

  /**
   * @brief Creating interal QA histograms of the jet tagger task
   */
  void                                UserCreateOutputObjects();

  /**
   * @brief Terminate method - not implemented for the jet tagger task
   */
  void                                Terminate(Option_t *option);

  //Setters

  /**
   * @brief Setting index of the base jet container
   * 
   * The base jet container contains the jets to be tagged (i.e. det. level jets)
   * @param[in] c Index of the jet container
   */
  void SetJetContainerBase(Int_t c)                             { fContainerBase = c;}

  /**
   * @brief Setting index of the tag jet container
   * 
   * The tag jet container contains the jets used for tagging (i.e. part. level jets)
   * @param[in] c Index of the jet container
   */
  void SetJetContainerTag(Int_t c)                              { fContainerTag  = c;}

  /**
   * @brief Set the tagging type
   * 
   * Refer to AliAnalysisTaskEmcalJetTagger::JetTaggingType for the definition
   * of the tagging types
   * @param[in] t Tagging type
   */
  void SetJetTaggingType(JetTaggingType t)                      { fJetTaggingType = t;}

  /**
   * @brief Set the tagging method
   * 
   * Refer to AliAnalysisTaskEmcalJetTagger::JetTaggingMethod for the definition
   * of the tagging method
   * @param[in] m Tagging method
   */
  void SetJetTaggingMethod(JetTaggingMethod m)                  { fJetTaggingMethod = m;}

  /**
   * Set the min. required fraction of shared particles
   * 
   * Only relevant for the method kFraction: In addition 
   * to pure geometric matching jets must have a minimum
   * overlap of particles (particles found in both jets)
   * @param[in] f Min. fraction of overlapping particles
   */
  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f; }

  /**
   * @brief Setting sumw2 option for all QA histograms
   */
  void SetUseSumw2(Bool_t b)                                    { fUseSumw2 = b;}

  /**
   * @brief Setting the acceptance type
   * 
   * Possible acceptance types:
   *   * 0 = use acceptance cuts of container  
   *   * 1 = allow extra margin (default: 0.1) one more for c2 in eta 
   *   * 2 = allow extra margin (default: 0.1) more in eta and phi for c2
   *   * 3 = allow extra margin (default: 0.1) in eta and phi for both containers
   * Use AliAnalysisTaskEmcalJetTagger::SetExtraMargin to define the margin size 
   * (absolute values).
   * 
   * @param[in] type Acceptance type
   */
  void SetTypeAcceptance(Int_t type)                            { fTypeAcc = type; /*see Init()*/}

  /**
   * @brief Set the maximum distance
   * 
   * Jets are only tagged if the closest jet or the matched jets 
   * have a distance smaller than the maximum distance
   * @param[in] dist Maximum distance allowed for the tag jet and the matched jets
   */
  void SetMaxDistance(Double_t dist)                            { fMaxDist = dist; }

  /**
   * @brief Set extra margin for jet acceptance in eta and phi for both the base and the tag jet container 
   * @param marginbase Size of the extra margin in eta and phi for the base jet container (absolute value)
   * @param margintag Size of the extra margin in eta and phi for the tag jet container (absolute value)
   */
  void SetExtraMargins(Double_t marginbase, Double_t margintag) { fExtraMarginAccBase = marginbase; fExtraMarginAccTag = margintag; }

  /**
   * @brief Set extra margin for jet acceptance in eta and phi for both the base and the tag jet container 
   * @param extramargin Size of the extra margin in eta and phi for the tag jet container (absolute value)
   */
  void SetExtraMarginTag(Double_t extramargin)                  { fExtraMarginAccTag = extramargin; }

  /**
   * @brief Set extra margin for jet acceptance in eta and phi for both the base and the tag jet container 
   * @param extramargin Size of the extra margin in eta and phi for the tag jet container (absolute value)
   */
  void SetExtraMarginBase(Double_t extramargin)                 { fExtraMarginAccBase = extramargin; }
  
  void SetSpecialParticleContainer(Int_t contnumb)              { fSpecPartContTag = contnumb; }
 protected:

  /**
   * @brief Retrieving event objects
   */
  Bool_t                              RetrieveEventObjects();

  /**
   * @brief Run jet tagging
   * 
   * Refer to AliAnalysisTaskEmcalJetTagger::MatchJetsGeo for more information about 
   * the jet tagging
   * @return true
   */
  Bool_t                              Run();

  /**
   * @brief Filling internal QA histograms of jet tagger task
   */
  Bool_t                              FillHistograms();

  /**
   * @brief Initialize jet tagger
   * 
   * Setting the acceptance cuts for base and tag jets. Refer to 
   * AliAnalysisTaskEmcalJetTagger::MatchJetsGeo for the definition of the acceptance types
   */
  void                                Init();

  /**
   * @brief Calculate azimuthal angle between the axises of the jets
   * @param[in] jet1 first jet
   * @param[in] jet2 second jet
   */
  Double_t GetDeltaPhi(const AliEmcalJet* jet1, const AliEmcalJet* jet2);

  /**
   * @brief Calculate azimuthal angle between the axises of the jets
   * @param[in] phi1 phi of the first jet
   * @param[in] phi2 phi of the second jet
   */
  Double_t GetDeltaPhi(Double_t phi1,Double_t phi2);


  /**
   * @brief Match the full jets to the corresponding charged jets
   * 
   * Translation of AliAnalysisHelperJetTasks::GetClosestJets to AliEmcalJet objects
   * 
   * Type can be: 
   *   * 0 = use acceptance cuts of container  
   *   * 1 = allow 0.1 one more for c2 in eta 
   *   * 2 = allow 0.1 more in eta and phi for c2
   *   * 3 = allow 0.1 in eta and phi for both containers
   * Jets are only tagged if the closest jet has a distance smaller than the maximum distance
   * 
   * @param[in] c1 Index of the first jet container (base jets)
   * @param[in] c2 Index of the second jet container (tag jets)
   * @param[in] iDebug Debug level of the jet matching method
   * @param[in] maxDist Maximum allowed distance
   * @param[in] type Acceptance type (see definition above)
   * @param[in] bReset If true tags are removed from all jets in the base and tag jet containers
   */
  void     MatchJetsGeo(Int_t c1 = -1, Int_t c2 = -1, Int_t iDebug = 0, Float_t maxDist = 0.3, Int_t type = 2, Bool_t bReset = kTRUE);
  
  /**
   * @brief Remove tag status of all jets in a given jet container
   * @param[in] c Index of the jet container
   */
  void     ResetTagging(const Int_t c);
  
 private:
  JetTaggingType                      fJetTaggingType;             ///< jet matching type
  JetTaggingMethod                    fJetTaggingMethod;           ///< jet matching method
  Int_t                               fContainerBase;              ///< jets to be tagged
  Int_t                               fContainerTag;               ///< jets used for tagging
  Int_t                               fSpecPartContTag;            ///< particle container optionally used in AliJetContainer::GetFractionSharedPt(). Set only if needed.
  Double_t                            fMinFractionShared;          ///< only fill histos for jets if shared fraction larger than X
  Bool_t                              fUseSumw2;                   ///< activate sumw2 for output histograms
  Bool_t                              fMatchingDone;               ///< flag to indicate if matching is done or not
  Int_t                               fTypeAcc;                    ///< acceptance cut for the jet containers, see method MatchJetsGeo in .cxx for possibilities
  Double_t                            fMaxDist;                    ///< distance allowed for two jets to match
  Double_t                            fExtraMarginAccBase;         ///< Extra margin to be added to the acceptance for the different acceptance types (base jet container)
  Double_t                            fExtraMarginAccTag;          ///< Extra margin to be added to the acceptance for the different acceptance types (tag jet container)
  Bool_t                              fInit;                       ///< true when the containers are initialized
  TH3                                 **fh3PtJet1VsDeltaEtaDeltaPhi; //!<! pt jet 1 vs deta vs dphi
  TH2                                 **fh2PtJet1VsDeltaR;         //!<! pt jet 1 vs dR
  TH2                                 **fh2PtJet2VsFraction;       //!<! pt jet 1 vs shared fraction
 
  TH2                                 **fh2PtJet1VsLeadPtAllSel;   //!<! all jets after std selection
  TH2                                 **fh2PtJet1VsLeadPtTagged;   //!<! tagged jets
  TH2                                 **fh2PtJet1VsPtJet2;         //!<! pT of base jet vs tagged jet
  TH2                                 **fh2PtJet2VsRelPt;          //!<! pT of tagged jet vs pt base jet / pt tagged jet
  
  TH3                                 *fh3PtJetDEtaDPhiConst;      //!<! pt jet vs delta eta vs delta phi of constituents
  TH3                                 *fh3PtJetAreaDRConst;        //!<! pt jet vs Area vs delta R of constituents
  TH1                                 *fNAccJets;                  //!<! number of jets per event
  AliAnalysisTaskEmcalJetTagger(const AliAnalysisTaskEmcalJetTagger&);            // not implemented
  AliAnalysisTaskEmcalJetTagger &operator=(const AliAnalysisTaskEmcalJetTagger&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetTagger, 9)
};
#endif

