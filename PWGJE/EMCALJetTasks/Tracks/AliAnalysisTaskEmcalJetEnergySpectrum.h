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
#ifndef __ALINANLYSISTASKJETENERGYSPECTRUM_H__
#define __ALINANLYSISTASKJETENERGYSPECTRUM_H__

#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisEmcalTriggerSelectionHelper.h"
#include "AliVCluster.h"
#include <vector>
#include <TArrayD.h>

class THistManager;
class AliEventCuts;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalJetEnergySpectrum : public AliAnalysisTaskEmcalJet, public AliAnalysisEmcalTriggerSelectionHelperImpl {
public:
  enum EJetTypeOutliers_t {
    kOutlierPartJet,
    kOutlierDetJet
  };
  AliAnalysisTaskEmcalJetEnergySpectrum();
  AliAnalysisTaskEmcalJetEnergySpectrum(EMCAL_STRINGVIEW name);
  virtual ~AliAnalysisTaskEmcalJetEnergySpectrum();

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetNameJetContainer(EMCAL_STRINGVIEW name) { fNameJetContainer = name; }
  void SetNameTriggerDecisionContainer(EMCAL_STRINGVIEW name) { fNameTriggerDecisionContainer = name; }
  void SetTriggerSelection(UInt_t triggerbits, EMCAL_STRINGVIEW triggerstring) { 
    fTriggerSelectionBits = triggerbits;
    fTriggerSelectionString = triggerstring;
  }
  void SetUseDownscaleWeight(bool doUse) { fUseDownscaleWeight = doUse; }
  void SetUseSumw2(Bool_t doUse) { fUseSumw2 = doUse; }
  void SetRangeRun1(Bool_t doUse) { fUseRun1Range = doUse; }
  void SetUseTriggerSelectionForData(bool doUse) { fUseTriggerSelectionForData = doUse; }
  void SetRequireSubsetMB(bool doRequire, ULong_t minbiastrigger = AliVEvent::kAny) { fRequireSubsetMB = doRequire; fMinBiasTrigger = minbiastrigger; }
  void SetUserPtBinning(int nbins, double *binning) { fUserPtBinning.Set(nbins+1, binning); }
  void SetRequestCentrality(bool doRequest) { fRequestCentrality = doRequest; }
  void SetRequestTriggerClusters(bool doRequest) { fRequestTriggerClusters = doRequest; }
  void SetCentralityEstimator(EMCAL_STRINGVIEW centest) { fCentralityEstimator = centest; }
  void SetFillHSparse(Bool_t doFill)               { fFillHSparse = doFill; }
  void SetUseMuonCalo(Bool_t doUse)                { fUseMuonCalo = doUse; }
  void SetEnergyScaleShfit(Double_t scaleshift)    { fScaleShift = scaleshift; } 
  void SetUseStandardOutlierRejection(bool doUse)  { fUseStandardOutlierRejection = doUse; }
  void SetJetTypeOutlierCut(EJetTypeOutliers_t jtype) { fJetTypeOutliers = jtype; }


  static AliAnalysisTaskEmcalJetEnergySpectrum *AddTaskJetEnergySpectrum(
    Bool_t isMC, 
    AliJetContainer::EJetType_t jettype, 
    AliJetContainer::ERecoScheme_t recoscheme, 
    AliVCluster::VCluUserDefEnergy_t energydef, 
    double radius, 
    EMCAL_STRINGVIEW namepartcont, 
    EMCAL_STRINGVIEW trigger, 
    EMCAL_STRINGVIEW suffix = ""
  );

protected:
  virtual void UserCreateOutputObjects();
  virtual bool Run();
  virtual bool IsTriggerSelected();
  virtual Bool_t CheckMCOutliers();
  virtual void RunChanged(Int_t newrun);

private:
  AliAnalysisTaskEmcalJetEnergySpectrum(const AliAnalysisTaskEmcalJetEnergySpectrum &);
  AliAnalysisTaskEmcalJetEnergySpectrum &operator=(const AliAnalysisTaskEmcalJetEnergySpectrum &);

  THistManager                  *fHistos;                       ///< Histogram manager
  Bool_t                        fIsMC;                          ///< Running on simulated events
  Bool_t                        fFillHSparse;                   ///< Fill THnSparses with more information
	UInt_t                        fTriggerSelectionBits;          ///< Trigger selection bits
  TString                       fTriggerSelectionString;        ///< Trigger selection string
  Bool_t                        fRequireSubsetMB;               ///< Require for triggers to be a subset of Min. Bias (for efficiency studies)
  ULong_t                       fMinBiasTrigger;                ///< Min bias trigger for trigger subset (for efficiency studies)
  TString                       fNameTriggerDecisionContainer;  ///< Global trigger decision container
  Bool_t                        fUseTriggerSelectionForData;    ///< Use trigger selection on data (require trigger patch in addition to trigger selection string)
  Bool_t                        fUseDownscaleWeight;            ///< Use 1/downscale as weight
  TString                       fNameJetContainer;              ///< Name of the jet container 
  Bool_t                        fRequestTriggerClusters;        ///< Request distinction of trigger clusters
  Bool_t                        fRequestCentrality;             ///< Request centrality
  Bool_t                        fUseRun1Range;                  ///< Use run1 run range for trending plots     
  Bool_t                        fUseSumw2;                      ///< Switch for sumw2 option in THnSparse (should not be used when a downscale weight is applied)
  Bool_t                        fUseMuonCalo;                   ///< Use events from the (muon)-calo-(fast) cluster
  Bool_t                        fUseStandardOutlierRejection;   ///< Use standard outlier rejection
  EJetTypeOutliers_t            fJetTypeOutliers;               ///< Jet type used for outlier detection
  Double_t                      fScaleShift;                    ///< Artificial jet energy scale shift
  TString                       fCentralityEstimator;           ///< Centrality estimator
  TArrayD                       fUserPtBinning;                 ///< User-defined pt-binning

  ClassDef(AliAnalysisTaskEmcalJetEnergySpectrum, 1);
};

}

}
#endif