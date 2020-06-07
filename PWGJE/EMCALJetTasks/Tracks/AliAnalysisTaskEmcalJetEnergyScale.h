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
#ifndef ALIANALYSISTASKEMCALJETENERGYSCALE_H
#define ALIANALYSISTASKEMCALJETENERGYSCALE_H

#include <TString.h>
#include "AliAnalysisTaskEmcalJet.h"
#include "AliJetContainer.h"
#include "AliVCluster.h"

class THistManager;
class TRandom;

namespace PWGJE {

namespace EMCALJetTasks{

class AliAnalysisTaskEmcalJetEnergyScale : public AliAnalysisTaskEmcalJet {
public:
  enum EJetTypeOutliers_t {
    kOutlierPartJet,
    kOutlierDetJet
  };
  AliAnalysisTaskEmcalJetEnergyScale();
  AliAnalysisTaskEmcalJetEnergyScale(const char *name);
  virtual ~AliAnalysisTaskEmcalJetEnergyScale();

  void SetNameDetJetContainer(const char *name)  { fNameDetectorJets = name; }
  void SetNamePartJetContainer(const char *name) { fNameParticleJets = name; }
  void SetTriggerName(const char *name)          { fTriggerSelectionString = name; }
  void SetFractionResponseClosure(double fraction) { fFractionResponseClosure = fraction; }
  void SetFillHSparse(Bool_t doFill)             { fFillHSparse = doFill; }
  void SetEnergyScaleShift(Double_t scaleshift)  { fScaleShift = scaleshift; }
  void SetUseStandardOutlierRejection(bool doUse) { fUseStandardOutlierRejection = doUse; }
  void SetDebugMaxJetOutliers(bool doDebug)      { fDebugMaxJetOutliers = doDebug; }
  void SetJetTypeOutlierCut(EJetTypeOutliers_t jtype) { fJetTypeOutliers = jtype; }

  static AliAnalysisTaskEmcalJetEnergyScale *AddTaskJetEnergyScale(
    AliJetContainer::EJetType_t       jetType,
    AliJetContainer::ERecoScheme_t    recoscheme,
    AliVCluster::VCluUserDefEnergy_t  energydef,
    Double_t                          radius,
    Bool_t                            useDCAL,
    const char *                      namepartcont,
    const char *                      trigger,
    const char *                      suffix
  );

protected:
  virtual void UserCreateOutputObjects();
  virtual Bool_t Run(); 
  virtual Bool_t CheckMCOutliers();
  bool IsSelectEmcalTriggers(const TString &triggerstring) const;

private:
  THistManager                *fHistos;                       //!<! Histogram collection
  TString                     fNameDetectorJets;              ///< Name of the data jet container
  TString                     fNameParticleJets;              ///< Name of the MC jet container
  TString                     fTriggerSelectionString;        ///< Trigger selection string
  TString                     fNameTriggerDecisionContainer;  ///< Global trigger decision container
  Double_t                    fFractionResponseClosure;       ///< Fraction of jets used for response in closure test
  Bool_t                      fFillHSparse;                   ///< Fill THnSparses
  Double_t                    fScaleShift;                    ///< Shift of the jet energy scale (fixed)
  Bool_t                      fUseStandardOutlierRejection;   ///< Use standard outlier rejection
  Bool_t                      fDebugMaxJetOutliers;           ///< Debug max jet determination for outlier rejection
  EJetTypeOutliers_t          fJetTypeOutliers;               ///< Jet type used for outlier detection
  TRandom                     *fSampleSplitter;               //!<! Sample splitter

  AliAnalysisTaskEmcalJetEnergyScale(const AliAnalysisTaskEmcalJetEnergyScale &);
  AliAnalysisTaskEmcalJetEnergyScale &operator=(const AliAnalysisTaskEmcalJetEnergyScale &);

  ClassDef(AliAnalysisTaskEmcalJetEnergyScale, 1);
};

}

}
#endif // ALIANALYSISTASKEMCALJETENERGYSCALE_H