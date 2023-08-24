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
#include <vector>
#include <iostream>
#include <TClonesArray.h>

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerDecision.h"
#include "AliEmcalTriggerSelection.h"
#include "AliEmcalTriggerSelectionCuts.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::AliEmcalTriggerSelection)
/// \endcond

namespace PWG{
namespace EMCAL{

AliEmcalTriggerSelection::AliEmcalTriggerSelection() :
  TNamed(),
  fSelectionCuts(nullptr),
  fTriggerAlias(nullptr)
{
}

AliEmcalTriggerSelection::AliEmcalTriggerSelection(const char *name, const AliEmcalTriggerSelectionCuts * const cuts, const AliEmcalTriggerAlias *alias):
  TNamed(name, ""),
  fSelectionCuts(cuts),
  fTriggerAlias(alias)
{
}

AliEmcalTriggerDecision* AliEmcalTriggerSelection::MakeDecison(const TClonesArray * const inputPatches, const AliEmcalTriggerSelectionCuts::RhoForTrigger &rhocontainer) const {
  AliDebugStream(1) << "Trigger selection " << GetName() << ": Make decision" << std::endl;
  AliEmcalTriggerDecision *result = new AliEmcalTriggerDecision(GetName());
  TIter patchIter(inputPatches);
  AliEMCALTriggerPatchInfo *patch(NULL);
  std::vector<AliEMCALTriggerPatchInfo *> selectedPatches;
  AliDebugStream(1) << "Number of input patches: " << inputPatches->GetEntries() << std::endl;
  while((patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(patchIter()))){
    if(fSelectionCuts->IsSelected(patch, rhocontainer)){
      selectedPatches.push_back(patch);
    }
  }
  AliDebugStream(1) << "Number of patches fulfilling the trigger condition: " << selectedPatches.size() << std::endl;
  // Find the main patch
  AliEMCALTriggerPatchInfo *mainPatch(NULL), *testpatch(NULL);
  for(std::vector<AliEMCALTriggerPatchInfo *>::iterator it = selectedPatches.begin(); it != selectedPatches.end(); ++it){
    testpatch = *it;
    if(!mainPatch) mainPatch = testpatch;
    else if(fSelectionCuts->CompareTriggerPatches(testpatch, mainPatch, rhocontainer) > 0) mainPatch = testpatch;
    result->AddAcceptedPatch(testpatch);
  }
  if(mainPatch) result->SetMainPatch(mainPatch);
  result->SetSelectionCuts(fSelectionCuts);
  if(fTriggerAlias) result->SetTriggerAlias(fTriggerAlias);
  return result;
}

void AliEmcalTriggerSelection::PrintStream(std::ostream &stream) const {
  stream << "  Name of the trigger class: " << GetName() << std::endl;
  if(fTriggerAlias) stream << fTriggerAlias;
  stream << *fSelectionCuts << std::endl;
}

}
}

std::ostream &operator<<(std::ostream &stream, const PWG::EMCAL::AliEmcalTriggerSelection &sel) {
  sel.PrintStream(stream);
  return stream;
}
