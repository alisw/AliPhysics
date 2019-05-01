/**************************************************************************
* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

#include <TBits.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include "AliVTrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "TParticle.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliCSTrackMaps.h"
#include "AliCSTrackCuts.h"
#include "AliPIDResponse.h"
#include "AliCSPIDCuts.h"
#include "AliLog.h"

/// \file AliCSTrackMaps.cxx
/// \brief Implementation of AOD track mapping class within the correlation studies analysis

TObjArray AliCSTrackMaps::fgAODTracksIdMap = TObjArray(30*1024);

/// Default constructor for serialization
AliCSTrackMaps::AliCSTrackMaps() :
  TNamed()
{
  fgAODTracksIdMap.SetOwner(kFALSE);
}

/// Constructor
/// \param name name of the AOD track maps
/// \param title title of the AOD track maps
AliCSTrackMaps::AliCSTrackMaps(const char *name, const char *title) :
  TNamed(name,title)
{
  fgAODTracksIdMap.SetOwner(kFALSE);
}

/// Destructor
/// We don't own anything, everything we allocate is owned
/// by the output list
AliCSTrackMaps::~AliCSTrackMaps()
{

}


/// Processes the start of a new event
///
/// In case of AOD event builds the AOD tracks id map
/// needed data for track selection
void AliCSTrackMaps::NotifyEvent() {

  /* let's build the potential AOD tracks id map */
  AliInputEventHandler *handler = AliCSAnalysisCutsBase::GetInputEventHandler();
  if (handler != NULL && !handler->InheritsFrom("AliESDInputHandler")) {
    AliInfo("Building the AOD tracks id map");

    fgAODTracksIdMap.Clear();

    AliAODEvent *event = (AliAODEvent *)handler->GetEvent();
    AliAODHeader *header = (AliAODHeader *) event->GetHeader();
    Int_t nesdtracks = header->GetNumberOfESDTracks();
    if (fgAODTracksIdMap.GetSize() < nesdtracks) {
      Int_t newsize = (int(nesdtracks / 1024) * 2) * 1024;
      AliInfo(Form("Expanding the AOD tracks id map from %d to %d", fgAODTracksIdMap.GetSize(), newsize));
      fgAODTracksIdMap.Expand(newsize);
    }

    for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {
      AliVTrack* vtrack = dynamic_cast<AliVTrack*>(event->GetTrack(itrk)); // pointer to reconstructed track
      if(vtrack != NULL) {
        if (vtrack->GetID() < 0) {
          /* check that the original track is already mapped */
          if (!(fgAODTracksIdMap.At(-1 - vtrack->GetID()) != NULL)) {
            AliError(Form("Original track with id: %d is not present for AOD track: %d with id: %d",
                -1 - vtrack->GetID(), itrk, vtrack->GetID()));
          }
        }
        else {
          if (fgAODTracksIdMap.At(vtrack->GetID()) != NULL) {
            AliError(Form("Original track with id: %d already present with id: %d",
                vtrack->GetID(), ((AliVTrack *) fgAODTracksIdMap.At(vtrack->GetID()))->GetID()));
          }
          else
            fgAODTracksIdMap.AddAt(vtrack, vtrack->GetID());
        }
      }
    }
  }
}

/// \cond CLASSIMP
ClassImp(AliCSTrackMaps);
/// \endcond
