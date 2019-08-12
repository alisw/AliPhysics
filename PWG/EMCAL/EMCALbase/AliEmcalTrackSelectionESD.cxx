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
#include <TBits.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObjArray.h>
#include <memory>
#include <iostream>

#include "AliEmcalESDHybridTrackCuts.h"
#include "AliEmcalESDTrackCutsGenerator.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliEmcalTrackSelResultCombined.h"
#include "AliEmcalCutBase.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVCuts.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTrackSelectionESD)
/// \endcond

using namespace PWG::EMCAL;

AliEmcalTrackSelectionESD::AliEmcalTrackSelectionESD():
		AliEmcalTrackSelection()
{
}

AliEmcalTrackSelectionESD::AliEmcalTrackSelectionESD(AliVCuts* cuts):
		AliEmcalTrackSelection()
{
  if(cuts) this->AddTrackCuts(cuts);
}

AliEmcalTrackSelectionESD::AliEmcalTrackSelectionESD(ETrackFilterType_t type, const char* period):
  AliEmcalTrackSelection()
{
  GenerateTrackCuts(type, period);
}

void AliEmcalTrackSelectionESD::GenerateTrackCuts(ETrackFilterType_t type, const char* period)
{
  using PWG::EMCAL::AliEmcalESDHybridTrackCuts;
  if (fListOfCuts) fListOfCuts->Clear();
  fSelectionModeAny = kTRUE;

  switch (type) {
  case kHybridTracks:
    AliDebugStream(1) << "Generate std hybrid track cuts" << std::endl;
    AliEmcalESDTrackCutsGenerator::AddHybridTrackCuts(this, period);
    break;

  case kTPCOnlyTracks:
    AliEmcalESDTrackCutsGenerator::AddTPCOnlyTrackCuts(this, period);
    break;

  case kITSPureTracks:
    AddTrackCuts(AliESDtrackCuts::GetStandardITSPureSATrackCuts2010(kTRUE, kFALSE));
    break;

  case kHybridTracks2010wNoRefit:
  {
    AliDebugStream(1) << "Generate 2010 hybrid track cuts wNoRefit" << std::endl;
    auto hybridcuts = new AliEmcalESDHybridTrackCuts("hybrid_2010_wNoRefit", AliEmcalESDHybridTrackCuts::kDef2010);
    hybridcuts->SetUseNoITSrefitTracks(kTRUE);
    AddTrackCuts(hybridcuts);
    break;
  }

  case kHybridTracks2010woNoRefit:
  {
    AliDebugStream(1) << "Generate 2010 hybrid track cuts woNoRefit" << std::endl;
    auto hybridcuts = new AliEmcalESDHybridTrackCuts("hybrid_2010_woNoRefit", AliEmcalESDHybridTrackCuts::kDef2010);
    hybridcuts->SetUseNoITSrefitTracks(kFALSE);
    AddTrackCuts(hybridcuts);
    break;
  }   

  case kHybridTracks2011wNoRefit:
  {
    AliDebugStream(1) << "Generate 2011 hybrid track cuts wNoRefit" << std::endl;
    auto hybridcuts = new AliEmcalESDHybridTrackCuts("hybrid_2011_wNoRefit", AliEmcalESDHybridTrackCuts::kDef2011);
    hybridcuts->SetUseNoITSrefitTracks(kTRUE);
    AddTrackCuts(hybridcuts);
    break;
  }

  case kHybridTracks2011woNoRefit:
  {
    AliDebugStream(1) << "Generate 2011 hybrid track cuts woNoRefit" << std::endl;
    auto hybridcuts = new AliEmcalESDHybridTrackCuts("hybrid_2011_woNoRefit", AliEmcalESDHybridTrackCuts::kDef2011);
    hybridcuts->SetUseNoITSrefitTracks(kFALSE);
    AddTrackCuts(hybridcuts);
    break;
  }
  case kHybridTracks2018TRD:
  {
    AliDebugStream(1) << "Generate 2018  hybrid track cuts with TRD track length" << std::endl;
    auto hybridcuts = new AliEmcalESDHybridTrackCuts("hybrid_2018_TRD", AliEmcalESDHybridTrackCuts::kDef2018TRD);
    AddTrackCuts(hybridcuts);
    break;
  }

  default:
    break;
  }
}

PWG::EMCAL::AliEmcalTrackSelResultPtr AliEmcalTrackSelectionESD::IsTrackAccepted(AliVTrack* const trk) {
  if (!fListOfCuts){
    AliDebugStream(2) << "No cut array " << std::endl;
    return PWG::EMCAL::AliEmcalTrackSelResultPtr(nullptr, kFALSE);
  } 
  AliESDtrack *esdt = dynamic_cast<AliESDtrack *>(trk);
  if (!esdt) {
    AliPicoTrack *picoTrack = dynamic_cast<AliPicoTrack *>(trk);
    if (picoTrack) {
      esdt = dynamic_cast<AliESDtrack*>(picoTrack->GetTrack());
    }
    else {
      AliError("Neither Pico nor ESD track");
      return PWG::EMCAL::AliEmcalTrackSelResultPtr(nullptr, kFALSE);
    }
  }

  TBits trackbitmap(64);
  trackbitmap.ResetAllBits();
  UInt_t cutcounter = 0;
  AliDebugStream(2) << "Found cut array with " << fListOfCuts->GetEntries() << " cuts\n" << std::endl;
  TClonesArray selectionStatus("PWG::EMCAL::AliEmcalTrackSelResultPtr", fListOfCuts->GetEntries());
  selectionStatus.SetOwner(kTRUE);
  for(auto cutIter : *fListOfCuts){
    AliDebugStream(3) << "executing nect cut: " << static_cast<AliVCuts *>(static_cast<AliEmcalManagedObject *>(cutIter)->GetObject())->GetName() << std::endl;
    PWG::EMCAL::AliEmcalCutBase *mycuts = static_cast<PWG::EMCAL::AliEmcalCutBase *>(static_cast<AliEmcalManagedObject *>(cutIter)->GetObject());
    PWG::EMCAL::AliEmcalTrackSelResultPtr selresult = mycuts->IsSelected(esdt);
    if(selresult) trackbitmap.SetBitNumber(cutcounter);
    new(selectionStatus[selectionStatus.GetEntries()]) PWG::EMCAL::AliEmcalTrackSelResultPtr(selresult);
    cutcounter++;
  }
  // In case of ANY at least one bit has to be set, while in case of ALL all bits have to be set
  PWG::EMCAL::AliEmcalTrackSelResultPtr result(esdt, kFALSE, new PWG::EMCAL::AliEmcalTrackSelResultCombined(&selectionStatus));
  if (fSelectionModeAny){
    result.SetSelectionResult(trackbitmap.CountBits() > 0 || cutcounter == 0);
  } else {
    result.SetSelectionResult(trackbitmap.CountBits() == cutcounter);
  }
  return result;
}

void AliEmcalTrackSelectionESD::SaveQAObjects(TList* outputList) {
  for(auto cutIter : *fListOfCuts){
    AliEmcalManagedObject *ptr = static_cast<AliEmcalManagedObject *>(cutIter);
    if(ptr->GetObject()->IsA() == AliESDtrackCuts::Class()){
      outputList->Add(ptr->GetObject());
      // Remove ownership - taken over by output list
      ptr->SetOwner(false);
    }
  }
}
