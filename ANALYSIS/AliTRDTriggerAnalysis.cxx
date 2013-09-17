/**************************************************************************
 * Copyright(c) 2013, ALICE Experiment at CERN, All rights reserved.      *
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

// evaluate TRD trigger conditions,
// potentially with hardened conditions to remove
// triggers caused by conversions of low-pt photons
// at large radii
//
// Author: Jochen Klein <jochen.klein@cern.ch>

#include "AliLog.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include "AliVTrdTrack.h"

#include "AliTRDTriggerAnalysis.h"

AliTRDTriggerAnalysis::AliTRDTriggerAnalysis() :
  TObject(),
  fTriggerFlags(0),
  fTriggerInputs(0),
  fTriggerClasses(0),
  fRequireMatch(kFALSE),
  fRequireMatchElectron(kTRUE),
  fRequireInTime(kFALSE),
  fTRDptHSE(3.),
  fTRDpidHSE(144),
  fTRDptHQU(2.),
  fTRDpidHQU(164),
  fTRDptHEE(3.),
  fTRDpidHEE(144),
  fTRDminSectorHEE(6),
  fTRDmaxSectorHEE(8),
  fTRDptHJT(3.),
  fTRDnHJT(3)
{
  // ctor
}

AliTRDTriggerAnalysis::~AliTRDTriggerAnalysis()
{
  // dtor
}

Bool_t AliTRDTriggerAnalysis::CalcTriggers(const AliVEvent *event)
{
  // evaluate the TRD trigger conditions,
  // so far HCO, HSE, HQU, HJT, HEE

  ResetTriggers();

  if (!event) {
    AliErrorClass("event pointer is null");
    return kFALSE;
  }

  TString trgClasses = event->GetFiredTriggerClasses();
  // UInt_t  trgInputs  = event->GetHeader()->GetL1TriggerInputs();

  Int_t nTracks[90]      = { 0 }; // stack-wise counted number of tracks above pt threshold

  Int_t nTrdTracks = event->GetNumberOfTrdTracks();

  if (nTrdTracks > 0)
    Fire(kHCO);

  for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
    AliVTrdTrack *trdTrack = event->GetTrdTrack(iTrack);
    if (!trdTrack)
      continue;

    // ignore the track if it was not in time
    // (if required)
    if (fRequireInTime && !trdTrack->GetTrackInTime())
      continue;

    AliVTrack *match = trdTrack->GetTrackMatch();
    AliDebug(2, Form("GTU track %2i with pt = %5.2f has match: %p (pt = %5.2f)",
		     iTrack, trdTrack->Pt(), match, match ? match->Pt() : 0));

    // ignore the track if it does not have a matched global track
    // (if required)
    if (fRequireMatch && !match)
      continue;

    Int_t globalStack = 5*trdTrack->GetSector() + trdTrack->GetStack();

    // stack-wise counting of tracks above pt threshold for jet trigger
    if (TMath::Abs(trdTrack->GetPt()) >= fTRDptHJT) {
      ++nTracks[globalStack];
    }

    // ignore the track for the electron triggers
    // if it does not have a matched global track
    // (if required)
    if (fRequireMatchElectron && !match)
      continue;

    if ((TMath::Abs(trdTrack->Pt()) > fTRDptHQU) && (trdTrack->GetPID() > fTRDpidHQU))
      Fire(kHQU);

    if ((TMath::Abs(trdTrack->Pt()) > fTRDptHSE) && (trdTrack->GetPID() > fTRDpidHSE))
      Fire(kHSE);

    if ((trdTrack->GetSector() >= fTRDminSectorHEE) && (trdTrack->GetSector() <= fTRDmaxSectorHEE) &&
	(TMath::Abs(trdTrack->Pt()) > fTRDptHSE) && (trdTrack->GetPID() > fTRDpidHSE))
      Fire(kHEE);
  }

  // check if HJT condition is fulfilled in any stack
  for (Int_t iStack = 0; iStack < 90; ++iStack) {
    if (nTracks[iStack] >= fTRDnHJT) {
      Fire(kHJT);
      break;
    }
  }

  return kTRUE;
}
