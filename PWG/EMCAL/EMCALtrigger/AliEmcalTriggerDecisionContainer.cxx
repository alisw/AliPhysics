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
#include "AliEmcalTriggerAlias.h"
#include "AliEmcalTriggerDecision.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include <string>
#include <sstream>

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::AliEmcalTriggerDecisionContainer)
/// \endcond

namespace PWG{
namespace EMCAL{

AliEmcalTriggerDecisionContainer::AliEmcalTriggerDecisionContainer():
  TNamed(),
  fContainer()
{
  fContainer.SetOwner(kTRUE);
}

AliEmcalTriggerDecisionContainer::AliEmcalTriggerDecisionContainer(const char* name):
  TNamed(name, ""),
  fContainer()
{
  fContainer.SetOwner(kTRUE);
}

void AliEmcalTriggerDecisionContainer::Reset() {
  if(fContainer.GetEntries()) fContainer.Clear();
}

void AliEmcalTriggerDecisionContainer::AddTriggerDecision(AliEmcalTriggerDecision* const decision) {
  fContainer.Add(decision);
}

const AliEmcalTriggerDecision* AliEmcalTriggerDecisionContainer::FindTriggerDecision(EMCAL_STRINGVIEW decname) const {
  const AliEmcalTriggerDecision* result = nullptr, *tmp = nullptr;
  TIter listiter(&fContainer);
  while((tmp = static_cast<AliEmcalTriggerDecision *>(listiter()))){
    auto alias = tmp->GetTriggerAlias();
    if(alias) {
      // trigger alias present, check for trigger class
      if(alias->HasTriggerClass(decname)) {
        result = tmp;
        break;
      }
    } else {
      // No trigger alias present, check for name of the trigger decision
      if(TString(tmp->GetName()) == TString(decname)) {
        result = tmp;
        break;
      }
    }
  }
  return result;
}

bool AliEmcalTriggerDecisionContainer::IsEventSelected(EMCAL_STRINGVIEW name)  const {
  if(name.find("&") != std::string::npos || name.find("|") != std::string::npos) {
    // combination or intersection of different trigger classes
    char separator;
    enum SeparatorMode {
      kSetCombine,
      kSetIntersect,
      kUnknown
    };
    SeparatorMode mode = SeparatorMode::kUnknown;
    if(name.find("&") != std::string::npos) {
      // separator & has priority
      separator = '&';
      mode = SeparatorMode::kSetIntersect;
    } else {
      separator = '|';
      mode = SeparatorMode::kSetCombine;
    }
    auto selectclasses = parseSelectionString(name, separator);
    bool selectionStatus = (mode == SeparatorMode::kSetIntersect) ? true : false;
    for(const auto trgclass : selectclasses) {
      auto triggerClassResult = IsEventSelected(trgclass.data());
      if((mode == SeparatorMode::kSetIntersect) && !triggerClassResult) {
        // All classes required - event rejection if one class is not present
        selectionStatus = false;
      } else if((mode == SeparatorMode::kSetCombine) && triggerClassResult) {
        // Any class required - event selected if at least one class is present
        selectionStatus = true;
      }
    }
    return selectionStatus;
  } else {
    // single trigger class requested
    const AliEmcalTriggerDecision *trg = FindTriggerDecision(name);
    if(trg) return trg->IsSelected();
    return false;
  }
  return false;
}

std::vector<std::string> AliEmcalTriggerDecisionContainer::parseSelectionString(EMCAL_STRINGVIEW selectionstring, char separator) const {
  std::vector<std::string> result;
  std::stringstream tokenizer(selectionstring.data());
  std::string buffer;
  while(std::getline(tokenizer, buffer, separator)) {
    result.push_back(buffer);
  }
  return result;
}

}
}
