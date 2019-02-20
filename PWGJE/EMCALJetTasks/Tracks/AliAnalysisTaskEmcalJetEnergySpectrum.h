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
#include <vector>
#include <TArrayD.h>

class THistManager;
class AliEventCuts;

namespace EmcalTriggerJets {

class AliAnalysisTaskEmcalJetEnergySpectrum : public AliAnalysisTaskEmcalJet {
public:
  enum TriggerCluster_t {
    kTrgClusterANY,
    kTrgClusterCENT,
    kTrgClusterCENTNOTRD,
    kTrgClusterCALO,
    kTrgClusterCALOFAST,
    kTrgClusterCENTBOTH,
    kTrgClusterOnlyCENT,
    kTrgClusterOnlyCENTNOTRD,
    kTrgClusterCALOBOTH,
    kTrgClusterOnlyCALO,
    kTrgClusterOnlyCALOFAST,
    kTrgClusterN
  };
  AliAnalysisTaskEmcalJetEnergySpectrum();
  AliAnalysisTaskEmcalJetEnergySpectrum(const char *name);
  virtual ~AliAnalysisTaskEmcalJetEnergySpectrum();

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetNameJetContainer(const char *name) { fNameJetContainer = name; }
  void SetNameTriggerDecisionContainer(const char *name) { fNameTriggerDecisionContainer = name; }
  void SetTriggerSelection(UInt_t triggerbits, const char *triggerstring) { 
    fTriggerSelectionBits = triggerbits;
    fTriggerSelectionString = triggerstring;
  }
  void SetUseDownscaleWeight(bool doUse) { fUseDownscaleWeight = doUse; }
  void SetUseAliEventCuts(bool doUse) { fUseAliEventCuts = doUse; }
  void SetUseTriggerSelectionForData(bool doUse) { fUseTriggerSelectionForData = doUse; }
  void SetRequireSubsetMB(bool doRequire, ULong_t minbiastrigger = AliVEvent::kAny) { fRequireSubsetMB = doRequire; fMinBiasTrigger = minbiastrigger; }
  void SetUserPtBinning(int nbins, double *binning) { fUserPtBinning.Set(nbins+1, binning); }
  void SetRequestCentrality(bool doRequest) { fRequestCentrality = doRequest; }
  void SetRequestTriggerClusters(bool doRequest) { fRequestTriggerClusters = doRequest; }
  void SetCentralityEstimator(const char *centest) { fCentralityEstimator = centest; }

  static AliAnalysisTaskEmcalJetEnergySpectrum *AddTaskJetEnergySpectrum(Bool_t isMC, AliJetContainer::EJetType_t jettype, double radius, const char *trigger, const char *suffix = "");

protected:
  virtual void UserCreateOutputObjects();
  virtual bool Run();
  virtual bool IsEventSelected();
  virtual bool IsTriggerSelected();
  virtual void RunChanged(Int_t newrun);
  std::vector<TriggerCluster_t> GetTriggerClusterIndices(const TString &triggerstring) const;
  bool IsSelectEmcalTriggers(const std::string &triggerstring) const;
  std::string MatchTrigger(const std::string &striggerstring);

private:
  AliAnalysisTaskEmcalJetEnergySpectrum(const AliAnalysisTaskEmcalJetEnergySpectrum &);
  AliAnalysisTaskEmcalJetEnergySpectrum &operator=(const AliAnalysisTaskEmcalJetEnergySpectrum &);

  THistManager                  *fHistos;                       ///< Histogram manager
  AliEventCuts                  *fEventCuts;                    ///< Event cuts as implemented in AliEventCuts
  Bool_t                        fIsMC;                          ///< Running on simulated events
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
  Bool_t                        fUseAliEventCuts;               ///< Flag switching on AliEventCuts;
  TString                       fCentralityEstimator;           ///< Centrality estimator
  TArrayD                       fUserPtBinning;                 ///< User-defined pt-binning

  ClassDef(AliAnalysisTaskEmcalJetEnergySpectrum, 1);
};

}
#endif