/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <TBits.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObjArray.h>
#include <memory>
#include <iostream>

#include "AliEmcalESDHybridTrackCuts.h"
#include "AliEmcalESDTrackCutsGenerator.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVCuts.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTrackSelectionESD)
/// \endcond

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
    auto hybridcuts = new AliEmcalESDHybridTrackCuts("hybrid_2010_wNoRefit", AliEmcalESDHybridTrackCuts::kDef2010);
    hybridcuts->SetUseNoITSrefitTracks(kTRUE);
    AddTrackCuts(hybridcuts);
    break;
  }

  case kHybridTracks2010woNoRefit:
  {
    auto hybridcuts = new AliEmcalESDHybridTrackCuts("hybrid_2010_wNoRefit", AliEmcalESDHybridTrackCuts::kDef2010);
    hybridcuts->SetUseNoITSrefitTracks(kFALSE);
    AddTrackCuts(hybridcuts);
    break;
  }   

  case kHybridTracks2011wNoRefit:
  {
    auto hybridcuts = new AliEmcalESDHybridTrackCuts("hybrid_2010_wNoRefit", AliEmcalESDHybridTrackCuts::kDef2011);
    hybridcuts->SetUseNoITSrefitTracks(kTRUE);
    AddTrackCuts(hybridcuts);
    break;
  }

  case kHybridTracks2011woNoRefit:
  {
    auto hybridcuts = new AliEmcalESDHybridTrackCuts("hybrid_2010_wNoRefit", AliEmcalESDHybridTrackCuts::kDef2011);
    hybridcuts->SetUseNoITSrefitTracks(kFALSE);
    AddTrackCuts(hybridcuts);
    break;
  }

  default:
    break;
  }
}

bool AliEmcalTrackSelectionESD::IsTrackAccepted(AliVTrack* const trk) {
  if (!fListOfCuts){
    AliDebugStream(2) << "No cut array " << std::endl;
    return kTRUE;
  } 
  AliESDtrack *esdt = dynamic_cast<AliESDtrack *>(trk);
  if (!esdt) {
    AliPicoTrack *picoTrack = dynamic_cast<AliPicoTrack *>(trk);
    if (picoTrack) {
      esdt = dynamic_cast<AliESDtrack*>(picoTrack->GetTrack());
    }
    else {
      AliError("Neither Pico nor ESD track");
      return kFALSE;
    }
  }

  fTrackBitmap.ResetAllBits();
  UInt_t cutcounter = 0;
  AliDebugStream(2) << "Found cut array with " << fListOfCuts->GetEntries() << " cuts\n" << std::endl;
  for(auto cutIter : *fListOfCuts){
    AliDebugStream(3) << "executing nect cut: " << static_cast<AliVCuts *>(static_cast<AliEmcalManagedObject *>(cutIter)->GetObject())->GetName() << std::endl;
    if((static_cast<AliVCuts *>(static_cast<AliEmcalManagedObject *>(cutIter)->GetObject()))->IsSelected(esdt)) fTrackBitmap.SetBitNumber(cutcounter);
    cutcounter++;
  }
  // In case of ANY at least one bit has to be set, while in case of ALL all bits have to be set
  if (fSelectionModeAny){
    return fTrackBitmap.CountBits() > 0 || cutcounter == 0;
  } else {
    return fTrackBitmap.CountBits() == cutcounter;
  }
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
