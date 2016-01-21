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
#include <TClonesArray.h>
#include <TBits.h>
#include <TObjArray.h>

#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliEmcalTrackSelectionAOD.h>
#include <AliESDtrack.h>
#include <AliPicoTrack.h>

/// \cond CLASSIMP
ClassImp(AliEmcalTrackSelectionAOD)
/// \endcond

/**
 * Main constructor, initialises fields with 0 (or NULL). For ROOT I/O, not intended
 * to be used by the users.
 */
AliEmcalTrackSelectionAOD::AliEmcalTrackSelectionAOD() :
	AliEmcalTrackSelection(),
	fFilterBits(0),
	fFilterHybridTracks(kFALSE),
	fFilterTPCTracks(kFALSE)
{
}

/**
 * Main Constructor, initalising also track cuts and filter bits. In case the initial cuts
 * is a nullpointer, only filter bits are used for the track selection. This constructor is
 * intended to be used by the users.
 *
 * \param cuts Inital track cut object (of type AliESDtrackCuts, can be a nullpointer)
 * \param filterbits Filter bits required
 */
AliEmcalTrackSelectionAOD::AliEmcalTrackSelectionAOD(AliVCuts* cuts, UInt_t filterbits):
	AliEmcalTrackSelection(),
	fFilterBits(filterbits),
  fFilterHybridTracks(kFALSE),
  fFilterTPCTracks(kFALSE)
{
  AddTrackCuts(cuts);
}

/**
 * Constructor, initalising track cuts depending on the requested type of filtering
 *
 * \param type Track filtering type
 * \param period  Period string (e.g. LHC11h)
 */
AliEmcalTrackSelectionAOD::AliEmcalTrackSelectionAOD(ETrackFilterType_t type, const char* period):
  AliEmcalTrackSelection(),
  fFilterBits(0),
  fFilterHybridTracks(kFALSE),
  fFilterTPCTracks(kFALSE)
{
  GenerateTrackCuts(type, period);
}

/**
 * Automatically generates track cuts depending on the requested type of filtering
 *
 * \param type Track filtering type
 */
void AliEmcalTrackSelectionAOD::GenerateTrackCuts(ETrackFilterType_t type, const char* /*period*/)
{
  switch (type) {
  case kHybridTracks:
    fFilterHybridTracks = kTRUE;
    fFilterTPCTracks = kFALSE;
    break;

  case kTPCOnlyTracks:
    fFilterHybridTracks = kFALSE;
    fFilterTPCTracks = kTRUE;
    break;

  default:
    fFilterHybridTracks = kFALSE;
    fFilterTPCTracks = kFALSE;
    break;
  }
}

/**
 * Function checks whether track is accepted under the given track selection cuts.
 * The function can handle AliAODTrack and AliPicoTrack, while for AliPico track an
 * AliAODTrack is expected to be the underlying structure. If it is not possible to
 * access an AOD track from the input track, the object will not be selected. Otherwise
 * first the status bits are checked (if requested), and if further track cuts (of type
 * AliESDtrackCuts) are provided, the track is converted to an ESD track for further checks.
 *
 * \param trk: Track to check
 * \return true if selected, false otherwise
 */
bool AliEmcalTrackSelectionAOD::IsTrackAccepted(AliVTrack * const trk)
{
  AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(trk);
  if (!aodt){
    AliPicoTrack *picotrack = dynamic_cast<AliPicoTrack*>(trk);
    if(picotrack) {
      aodt = dynamic_cast<AliAODTrack *>(picotrack->GetTrack());
    }
    else {
      AliError("Track neither AOD track nor pico track");
      return kFALSE;
    }
  }
  if(!aodt){
    AliError("Failed getting AOD track");
    return kFALSE;
  }

  fTrackBitmap.ResetAllBits();
  Int_t cutcounter(0);
  if (fFilterBits) {
    if(aodt->TestFilterBit(fFilterBits)) fTrackBitmap.SetBitNumber(cutcounter);
    cutcounter++;
  }
  if(fFilterHybridTracks) {
    if(aodt->IsHybridGlobalConstrainedGlobal()) fTrackBitmap.SetBitNumber(cutcounter++);
  }
  if(fFilterTPCTracks) {
    if(aodt->IsHybridTPCConstrainedGlobal()) fTrackBitmap.SetBitNumber(cutcounter++);
  }
  if (fListOfCuts) {
    for(TIter cutIter = TIter(fListOfCuts).Begin(); cutIter != TIter::End(); ++cutIter){
      AliVCuts *trackCuts = static_cast<AliVCuts *>(*cutIter);
      if(trackCuts->IsA() == AliESDtrackCuts::Class()){
        // If track cuts are AliESDtrackCuts, the track needs to be converted to an AliESDtrack before
        AliESDtrack copyTrack(aodt);
        if(trackCuts->IsSelected(&copyTrack)) fTrackBitmap.SetBitNumber(cutcounter);
      } else{
        if(trackCuts->IsSelected(aodt)) fTrackBitmap.SetBitNumber(cutcounter);
      }
      cutcounter++;
    }
  }

  if (fSelectionModeAny){
    // In case of ANY one of the cuts need to be fulfilled (equivalent to one but set)
    return fTrackBitmap.CountBits() > 0 || cutcounter == 0;
  } else {
    // In case of ALL all of the cuts need to be fulfilled (equivalent to all bits set)
    return fTrackBitmap.CountBits() == cutcounter;
  }
}
