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
#ifndef ALIANALYSISTASKEMCALJETCONSITUTENTQA_H
#define ALIANALYSISTASKEMCALJETCONSITUTENTQA_H

#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>

#include "AliAnalysisTaskEmcalJet.h"

class THistManager;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalJetConstituentQA : public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskEmcalJetConstituentQA();
  AliAnalysisTaskEmcalJetConstituentQA(const char *name);
  virtual ~AliAnalysisTaskEmcalJetConstituentQA();

  static AliAnalysisTaskEmcalJetConstituentQA *AddTaskEmcalJetConstituentQA(const char *trigger, AliJetContainer::EJetType_t jettype, bool parmode = kFALSE);

  void SetTriggerSelection(const char *trigger) { fTriggerSelectionString = trigger; } 
  void SetNameTrackContainer(const char *name) { fNameTrackContainer = name; }
  void SetNameClusterContainer(const char *name) { fNameClusterContainer = name; }
  void AddNameJetContainer(const char *name) { fNamesJetContainers.Add(new TObjString(name)); }
  void SetUseTriggerSelection(Bool_t doUse) { fUseTriggerSelection = doUse; }
  void SetJetType(AliJetContainer::EJetType_t jettype) { fJetType = jettype; }

  void SetDoHighZClusters(Bool_t doMon)           { fDoHighZClusters = doMon; }
  void SetDoChargedConstituents(Bool_t doMon)     { fDoCharged = doMon; }
  void SetDoNeutralConstituents(Bool_t doMon)     { fDoNeutral = doMon; }
  void SetDoLeadingCharged(Bool_t doMon)          { fDoLeadingCharged = doMon; }
  void SetDoLeadingNeutral(Bool_t doMon)          { fDoLeadingNeutral = doMon; }
  void SetDoLeadingCell(Bool_t doMon)             { fDoLeadingCell = doMon; }

protected:
  virtual void UserCreateOutputObjects();
  virtual bool Run();
  virtual bool IsTriggerSelected();

private:
  THistManager                  *fHistos;                   ///< Histogram manager

  AliJetContainer::EJetType_t   fJetType;                   ///< Jet type
  TString                       fNameTrackContainer;        ///< Name of the track container
  TString                       fNameClusterContainer;      ///< Name of the cluster container
  TObjArray                     fNamesJetContainers;        ///< Names of the connected jet container
  TString                       fTriggerSelectionString;    ///< Trigger selection string
  Bool_t                        fUseTriggerSelection;       ///< Use trigger selection in addition to trigger string
  TString                       fNameTriggerDecisionContainer;  ///< Name of the trigger decision container

  // Switchable options
  Bool_t                        fDoHighZClusters;           ///< Monoitor clusters with a high z
  Bool_t                        fDoCharged;                 ///< Monitor charged constituents
  Bool_t                        fDoNeutral;                 ///< Monitor neutral constituents
  Bool_t                        fDoLeadingCharged;          ///< Monitor leading charged constituents
  Bool_t                        fDoLeadingNeutral;          ///< Monitor leading netural constituents
  Bool_t                        fDoLeadingCell;             ///< Monitor leading cell


  AliAnalysisTaskEmcalJetConstituentQA(const AliAnalysisTaskEmcalJetConstituentQA &);
  AliAnalysisTaskEmcalJetConstituentQA &operator=(const AliAnalysisTaskEmcalJetConstituentQA &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetConstituentQA, 1);
  /// \endcond
};

}

}

#endif