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
#ifndef __ALINANLYSISTASKJETENERGYSPECTRUMTRIGGERPBPBQA_H__
#define __ALINANLYSISTASKJETENERGYSPECTRUMTRIGGERPBPBQA_H__

#include "AliAnalysisTaskEmcalJet.h"
#include "AliEventCuts.h"
#include <vector>
#include <TArrayD.h>

class THistManager;

namespace EmcalTriggerJets {

class AliAnalysisTaskEmcalJetPbPbTriggerTests : public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskEmcalJetPbPbTriggerTests();
  AliAnalysisTaskEmcalJetPbPbTriggerTests(const char *name);
  virtual ~AliAnalysisTaskEmcalJetPbPbTriggerTests();

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetNameJetContainer(const char *name) { fNameJetContainer = name; }
  void SetUserPtBinning(int nbins, double *binning) { fUserPtBinning.Set(nbins+1, binning); }
  void SetTriggerSelection(UInt_t triggerbits) { fTriggerSelectionBits = triggerbits; }
  void SetUseAliEventCuts(Bool_t b)   { fUseAliEventCuts = b; }
  
  static AliAnalysisTaskEmcalJetPbPbTriggerTests *AddTaskJetPbPbTriggerTests(Bool_t isMC, AliJetContainer::EJetType_t jettype, double radius, const char *suffix = "");

protected:
  virtual void UserCreateOutputObjects();
  virtual bool Run();
  virtual bool IsTriggerSelected();

private:
  AliAnalysisTaskEmcalJetPbPbTriggerTests(const AliAnalysisTaskEmcalJetPbPbTriggerTests &);
  AliAnalysisTaskEmcalJetPbPbTriggerTests &operator=(const AliAnalysisTaskEmcalJetPbPbTriggerTests &);

  THistManager                  *fHistos;                       ///< Histogram manager
  Bool_t                        fIsMC;                          ///< Running on simulated events
  TString                       fNameJetContainer;              ///< Name of the jet container 
  TArrayD                       fUserPtBinning;                 ///< User-defined pt-binning
  // Event selection
  UInt_t                        fTriggerSelectionBits;          ///< Trigger selection bit
  Bool_t                      fUseAliEventCuts;                 ///< Flag to use AliEventCuts (otherwise AliAnalysisTaskEmcal will be used)
  AliEventCuts                fEventCuts;                       ///< event selection utility


  ClassDef(AliAnalysisTaskEmcalJetPbPbTriggerTests, 1);
};

}
#endif
