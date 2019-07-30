/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIANALYSISTASKEMCALJETITERATIVEDECLUSTERING_H
#define ALIANALYSISTASKEMCALJETITERATIVEDECLUSTERING_H 

#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisEmcalTriggerSelectionHelper.h"
#include <string>

class THistManager;

namespace PWGJE {

namespace EMCALJetTasks {

class AliLundPlaneHelper;

class AliAnalysisTaskEmcalJetIterativeDeclustering : public AliAnalysisTaskEmcalJet, public AliAnalysisEmcalTriggerSelectionHelperImpl {
public:
  AliAnalysisTaskEmcalJetIterativeDeclustering();
  AliAnalysisTaskEmcalJetIterativeDeclustering(EMCAL_STRINGVIEW name);
  virtual ~AliAnalysisTaskEmcalJetIterativeDeclustering();

  void SetUseDownscaleWeight(Bool_t doUse) { fUseDownscaleWeight = doUse; }
  void SetJetPtRange(double ptmin, double ptmax) { fJetPtMin = ptmin; fJetPtMax = ptmax; }
  void SetJetPtMin(double ptmin) { fJetPtMin = ptmin; }
  void SetJetPtMax(double ptmax) { fJetPtMax = ptmax; }
  void SetHardCutoff(double pthard) { fHardCutoff = pthard; }
  void SetSelectTrigger(UInt_t triggerbits, const char * triggerstring) { fTriggerBits = triggerbits; fTriggerString = triggerstring; }

  static AliAnalysisTaskEmcalJetIterativeDeclustering *AddTaskEmcalJetIterativeDeclustering(Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, EMCAL_STRINGVIEW trigger);
protected:
  virtual void UserCreateOutputObjects();
  virtual Bool_t Run();

  virtual bool IsTriggerSelected();
  virtual void RunChanged(Int_t newrun);

  Double_t GetDownscaleWeight() const;

private:
    AliLundPlaneHelper          *fDecluster;                //!<!
    THistManager                *fHistos;                   //!<!
    Bool_t                      fUseDownscaleWeight;
    Bool_t                      fUseChargedConstituents;
    Bool_t                      fUseNeutralConstituents;
    ULong_t                     fTriggerBits;
    std::string                 fTriggerString;
    Double_t                    fJetPtMin;
    Double_t                    fJetPtMax;
    Double_t                    fHardCutoff;

    AliAnalysisTaskEmcalJetIterativeDeclustering(const AliAnalysisTaskEmcalJetIterativeDeclustering &);
    AliAnalysisTaskEmcalJetIterativeDeclustering &operator=(const AliAnalysisTaskEmcalJetIterativeDeclustering &);

    ClassDef(AliAnalysisTaskEmcalJetIterativeDeclustering, 1);
};

}

}
#endif