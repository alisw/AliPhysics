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
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVTrdTrack.h"
#include "AliESDTrdTrigger.h"

#include "AliTRDTriggerAnalysis.h"

AliTRDTriggerAnalysis::AliTRDTriggerAnalysis() :
  TObject(),
  fTriggerFlags(),
  fTriggerInputs(0),
  fTriggerClasses(0),
  fVerbosity(0),
  fRequireMatch(kFALSE),
  fRequireMatchElectron(kFALSE),
  fRequireInTime(kTRUE),
  fTRDlayerMaskEl(0x1),
  fTRDnTrackletsEl(5),
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

  memset(fTriggerFlags, 0, sizeof(fTriggerFlags));
  memset(fTriggerContribs, 0, sizeof(fTriggerContribs));
}

AliTRDTriggerAnalysis::~AliTRDTriggerAnalysis()
{
  // dtor
}

void AliTRDTriggerAnalysis::ResetTriggers()
{
  // reset internal cache of trigger status

  memset(fTriggerFlags, 0, sizeof(fTriggerFlags));
  fTriggerInputs = fTriggerClasses = 0;
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

  // GTU information
  UInt_t header = 0x0;

  if (fVerbosity > 0)
    printf("******************************************************************\n");
  const AliESDEvent *esdEvent = dynamic_cast<const AliESDEvent*> (event);
  if (esdEvent) {
    AliESDTrdTrigger* trdTriggerInfo = esdEvent->GetTrdTrigger();

    for (Int_t iSector = 0; iSector < 18; ++iSector) {
      UInt_t  trgFlags = trdTriggerInfo->GetFlags(iSector);

      UInt_t trgContribs = fTriggerContribs[iSector] = trgFlags & 0xfff;
      header |= trgContribs;
      if (fVerbosity > 0)
	printf("sector %2i: %5s %5s %5s %5s %5s %5s %5s %5s (0x%03x)\n", iSector,
	       trgContribs & (1 << 7) ? "TO" : "--",
	       trgContribs & (1 << 6) ? "E2" : "--",
	       trgContribs & (1 << 5) ? "E1" : "--",
	       trgContribs & (1 << 4) ? "J1" : "--",
	       trgContribs & (1 << 3) ? "H2" : "--",
	       trgContribs & (1 << 2) ? "H1" : "--",
	       trgContribs & (1 << 1) ? "E3" : "--",
	       trgContribs & (1 << 0) ? "M1" : "--",
	       trgContribs);

      // trackingDoneTimeSMU = ((trgFlags >> 12) & 0x3ff) * 1./120.;

      // for (Int_t iStack = 0; iStack < 5; ++iStack) {
      //   trackingDoneTMU = ((trgFlags >> 27) & (1 << iStack)) ? kTRUE : kFALSE;  // TMU-level tracking done flag
      //   trackingDoneSMUStack = ((trgFlags >> 22) & (1 << iStack)) ? kTRUE : kFALSE;  // SMU-level stack-related tracking done flag

      //   lmeFlags = trdTriggerInfo->GetLME(iStack) & 0xffffff;
      //   crcErrorFlags = (~(trdTriggerInfo->GetLME(iStack) >> 24)) & 0x3;
      // }
    }
    if (fVerbosity > 0) {
      printf("------------------------------------------------------------------\n");
      printf("total    : %5s %5s %5s %5s %5s %5s %5s %5s (0x%03x)\n",
	     header & (1 << 7) ? "TO" : "--",
	     header & (1 << 6) ? "E2" : "--",
	     header & (1 << 5) ? "E1" : "--",
	     header & (1 << 4) ? "J1" : "--",
	     header & (1 << 3) ? "H2" : "--",
	     header & (1 << 2) ? "H1" : "--",
	     header & (1 << 1) ? "E3" : "--",
	     header & (1 << 0) ? "M1" : "--",
	     header);
    }
  }

  // evaluate trigger classes
  TString trgClasses = event->GetFiredTriggerClasses();
  if (trgClasses.Contains("TRDCO2"))
    MarkClass(kHCO);
  if (trgClasses.Contains("WUHJT"))
    MarkClass(kHJT);
  if (trgClasses.Contains("WUHSE"))
    MarkClass(kHSE);
  if (trgClasses.Contains("WUHQU"))
    MarkClass(kHQU);
  if (trgClasses.Contains("WUHEE"))
    MarkClass(kHEE);

  // evaluate trigger inputs
  UInt_t  trgInputs = 0;
  if (esdEvent)
    trgInputs  = esdEvent->GetHeader()->GetL1TriggerInputs();
  else if (const AliAODEvent *aodEvent = dynamic_cast<const AliAODEvent*> (event))
    trgInputs  = aodEvent->GetHeader()->GetL1TriggerInputs();
  else
    AliError("failed to retrieve L1 trigger inputs");

  if (trgInputs & (1 <<  8))
    MarkInput(kHCO);
  if (trgInputs & (1 <<  9))
    MarkInput(kHJT);
  if (trgInputs & (1 << 10))
    MarkInput(kHSE);
  if (trgInputs & (1 << 12))
    MarkInput(kHQU);
  if (trgInputs & (1 << 13))
    MarkInput(kHEE);

  // evaluate TRD GTU tracks
  Int_t nTracks[90]      = { 0 }; // stack-wise counted number of tracks above pt threshold

  Int_t nTrdTracks = event->GetNumberOfTrdTracks();

  for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
    AliVTrdTrack *trdTrack = event->GetTrdTrack(iTrack);
    if (!trdTrack) {
      AliError(Form("Failed to get track %i", iTrack));
      continue;
    }

    Int_t globalStack = 5*trdTrack->GetSector() + trdTrack->GetStack();

    MarkCondition(kHCO, globalStack);

    for (Int_t iLayer = 0; iLayer < 6; ++iLayer) {
      if (trdTrack->GetLayerMask() & (1 << iLayer)) {
	AliVTrdTracklet *trkl = trdTrack->GetTracklet(iLayer);
	if (!trkl) {
	  AliError(Form("no tracklet in layer %i where one should be for track %i",
		   iLayer, iTrack));
	}
      }
    }

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

    // stack-wise counting of tracks above pt threshold for jet trigger
    if (TMath::Abs(trdTrack->Pt()) >= fTRDptHJT) {
      ++nTracks[globalStack];
    }

    // ignore the track for the electron triggers
    // if it does not have a matched global track
    // (if required)
    if (fRequireMatchElectron && !match)
      continue;

    // ignore the track for the electron triggers
    // if it does not fulfill the tracklet requirement
    if (trdTrack->GetNTracklets() < fTRDnTrackletsEl)
      continue;
    if ((trdTrack->GetLayerMask() & fTRDlayerMaskEl) != fTRDlayerMaskEl)
      continue;

    if ((TMath::Abs(trdTrack->Pt()) >= fTRDptHQU) && (trdTrack->GetPID() >= fTRDpidHQU))
      MarkCondition(kHQU, globalStack);

    if ((TMath::Abs(trdTrack->Pt()) >= fTRDptHSE) && (trdTrack->GetPID() >= fTRDpidHSE))
      MarkCondition(kHSE, globalStack);

    if ((trdTrack->GetSector() >= fTRDminSectorHEE) && (trdTrack->GetSector() <= fTRDmaxSectorHEE) &&
	(TMath::Abs(trdTrack->Pt()) >= fTRDptHEE) && (trdTrack->GetPID() >= fTRDpidHEE))
      MarkCondition(kHEE, globalStack);
  }

  // check if HJT condition is fulfilled in any stack
  for (Int_t iStack = 0; iStack < 90; ++iStack) {
    if (nTracks[iStack] >= fTRDnHJT) {
      MarkCondition(kHJT, iStack);
      break;
    }
  }

  return kTRUE;
}
