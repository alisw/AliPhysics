/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Decoding data from the TRD raw stream                                 //
//  and translation into ADC values, on-line tracklets and tracks         //
//                                                                        //
//  Author: J. Klein (jochen.klein@cern.ch)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdarg>

#include "TClonesArray.h"
#include "TTree.h"

#include "AliLog.h"
#include "AliRawReader.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDdigitsParam.h"
#include "AliTRDtrapConfig.h"
#include "AliTRDarrayADC.h"
#include "AliTRDarrayDictionary.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDtrackletWord.h"
#include "AliESDTrdTrack.h"
#include "AliTreeLoader.h"
#include "AliLoader.h"

#include "AliTRDrawStream.h"

// temporary
#include "AliRunLoader.h"

ClassImp(AliTRDrawStream)

// some static information
Int_t AliTRDrawStream::fgMcmOrder[] = {12, 13, 14, 15,
				       8, 9, 10, 11,
				       4, 5, 6, 7,
				       0, 1, 2, 3};
Int_t  AliTRDrawStream::fgRobOrder [] = {0, 1, 2, 3};
const Int_t  AliTRDrawStream::fgkNlinks = 12;
const Int_t  AliTRDrawStream::fgkNstacks = 5;
const Int_t  AliTRDrawStream::fgkNsectors = 18;
const Int_t  AliTRDrawStream::fgkNtriggers = 12;
const UInt_t AliTRDrawStream::fgkDataEndmarker     = 0x00000000;
const UInt_t AliTRDrawStream::fgkTrackletEndmarker = 0x10001000;

const char* AliTRDrawStream::fgkErrorMessages[] = {
  "Unknown error",
  "Link monitor active",
  "Pretrigger counter mismatch",
  "not a TRD equipment (1024-1041)",
  "Invalid Stack header",
  "Invalid detector number",
  "No digits could be retrieved from the digitsmanager",
  "HC header mismatch",
  "HC check bits wrong",
  "Unexpected position in readout stream",
  "Invalid testpattern mode",
  "Testpattern mismatch",
  "Number of timebins changed",
  "ADC mask inconsistent",
  "ADC check bits invalid",
  "Missing ADC data",
  "Missing expected ADC channels",
  "Missing MCM headers",
  "Missing TP data"
};

Int_t AliTRDrawStream::fgErrorDebugLevel[] = {
  0,
  0,
  2,
  1,
  0,
  1,
  1,
  1,
  1,
  2,
  1,
  0,
  1,
  1,
  2,
  1,
  1,
  1,
  0
};

AliTRDrawStream::ErrorBehav_t AliTRDrawStream::fgErrorBehav[] = {
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kDiscardHC,
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kAbort,
  AliTRDrawStream::kAbort,
  AliTRDrawStream::kAbort,
  AliTRDrawStream::kAbort,
  AliTRDrawStream::kDiscardHC,
  AliTRDrawStream::kDiscardHC,
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kTolerate,
  AliTRDrawStream::kTolerate
};

AliTRDrawStream::AliTRDrawStream(AliRawReader *rawReader) :
  fStoreError(&AliTRDrawStream::ForgetError),
  fRawReader(rawReader),
  fDigitsManager(0x0),
  fDigitsParam(0x0),
  fErrors(0x0),
  fLastError(),
  fErrorFlags(0),
  fStats(),
  fPayloadStart(0x0),
  fPayloadCurr(0x0),
  fPayloadSize(0),
  fNtimebins(-1),
  fLastEvId(-1),
  fCurrSlot(-1),
  fCurrLink(-1),
  fCurrRobPos(-1),
  fCurrMcmPos(-1),
  fCurrEquipmentId(0),
  fCurrSmHeaderSize(0),
  fCurrSmHeaderVersion(0),
  fCurrTrailerReadout(0),
  fCurrTrgHeaderAvail(0),
  fCurrTrgHeaderReadout(0),
  fCurrTrkHeaderAvail(0),
  fCurrEvType(0),
  fCurrTriggerEnable(0),
  fCurrTriggerFired(0),
  fCurrTrackEnable(0),
  fCurrTrackletEnable(0),
  fCurrStackMask(0),
  fCurrTrkHeaderIndexWord(0x0),
  fCurrTrkHeaderSize(0x0),
  fCurrTrkFlags(0x0),
  fCurrTrgHeaderIndexWord(0x0),
  fCurrTrgHeaderSize(0x0),
  fCurrTrgFlags(0x0),
  fCurrStackIndexWord(0x0),
  fCurrStackHeaderSize(0x0),
  fCurrStackHeaderVersion(0x0),
  fCurrLinkMask(0x0),
  fCurrCleanCheckout(0x0),
  fCurrBoardId(0x0),
  fCurrHwRev(-1),
  fCurrHwRevTMU(0x0),
  fCurrLinkMonitorFlags(0x0),
  fCurrLinkDataTypeFlags(0x0),
  fCurrLinkDebugFlags(0x0),
  fCurrSpecial(-1),
  fCurrMajor(-1),
  fCurrMinor(-1),
  fCurrAddHcWords(-1),
  fCurrSm(-1),
  fCurrStack(-1),
  fCurrLayer(-1),
  fCurrSide(-1),
  fCurrHC(-1),
  fCurrCheck(-1),
  fCurrNtimebins(-1),
  fCurrBC(-1),
  fCurrPtrgCnt(-1),
  fCurrPtrgPhase(-1),
  fNDumpMCMs(0),
  fTrackletArray(0x0),
  fAdcArray(0x0),
  fSignalIndex(0x0),
  fTrackletTree(0x0),
  fTracklets(0x0),
  fTracks(0x0),
  fMarkers(0x0)
{
  // default constructor

  fCurrTrkHeaderIndexWord = new UInt_t[fgkNstacks];
  fCurrTrkHeaderSize      = new UInt_t[fgkNstacks];
  fCurrTrkFlags           = new ULong64_t[fgkNsectors*fgkNstacks];
  fCurrTrgHeaderIndexWord = new UInt_t[fgkNtriggers];
  fCurrTrgHeaderSize      = new UInt_t[fgkNtriggers];
  fCurrTrgFlags           = new UInt_t[fgkNsectors];
  fCurrStackIndexWord     = new UInt_t[fgkNstacks];
  fCurrStackHeaderSize    = new UInt_t[fgkNstacks];
  fCurrStackHeaderVersion = new UInt_t[fgkNstacks];
  fCurrLinkMask           = new UInt_t[fgkNstacks];
  fCurrCleanCheckout      = new UInt_t[fgkNstacks];
  fCurrBoardId            = new UInt_t[fgkNstacks];
  fCurrHwRevTMU           = new UInt_t[fgkNstacks];
  fCurrLinkMonitorFlags   = new UInt_t[fgkNstacks * fgkNlinks];
  fCurrLinkDataTypeFlags  = new UInt_t[fgkNstacks * fgkNlinks];
  fCurrLinkDebugFlags     = new UInt_t[fgkNstacks * fgkNlinks];
  for (Int_t iSector = 0; iSector < fgkNsectors; iSector++)
    fCurrTrgFlags[iSector] = 0;
  for (Int_t i = 0; i < 100; i++)
    fDumpMCM[i] = 0;

  // preparing TClonesArray
  fTrackletArray = new TClonesArray("AliTRDtrackletWord", 256);

  // setting up the error tree
  fErrors = new TTree("errorStats", "Error statistics");
  fErrors->SetDirectory(0x0);
  fErrors->Branch("error", &fLastError);
  fErrors->SetCircular(1000);
  for (Int_t i = 0; i < 100; i++) {
    fErrorBuffer[i] = 0;
  }

}

AliTRDrawStream::~AliTRDrawStream()
{
  // destructor

  delete fErrors;

  delete [] fCurrTrkHeaderIndexWord;
  delete [] fCurrTrkHeaderSize;
  delete [] fCurrTrkFlags;
  delete [] fCurrTrgHeaderIndexWord;
  delete [] fCurrTrgHeaderSize;
  delete [] fCurrTrgFlags;
  delete [] fCurrStackIndexWord;
  delete [] fCurrStackHeaderSize;
  delete [] fCurrStackHeaderVersion;
  delete [] fCurrLinkMask;
  delete [] fCurrCleanCheckout;
  delete [] fCurrBoardId;
  delete [] fCurrHwRevTMU;
  delete [] fCurrLinkMonitorFlags;
  delete [] fCurrLinkDataTypeFlags;
  delete [] fCurrLinkDebugFlags;
}

Bool_t AliTRDrawStream::ReadEvent(TTree *trackletTree)
{
  // read the current event from the raw reader and fill it to the digits manager

  if (!fRawReader) {
    AliError("No raw reader available");
    return kFALSE;
  }

  // tracklet output
  ConnectTracklets(trackletTree);

  // some preparations
  fDigitsParam = 0x0;

  // loop over all DDLs
  // data starts with GTU payload, i.e. SM index word
  UChar_t *buffer = 0x0;

  while (fRawReader->ReadNextData(buffer)) {

    fCurrEquipmentId = fRawReader->GetEquipmentId();
    AliDebug(2, Form("equipment: %i", fCurrEquipmentId));

    if (fCurrEquipmentId < kDDLOffset || fCurrEquipmentId > kDDLMax) {
      EquipmentError(kNonTrdEq, "Skipping");
      continue;
    }

    if (fMarkers)
      new ((*fMarkers)[fMarkers->GetEntriesFast()])
	AliTRDrawStreamError(-kSecactive, fCurrEquipmentId - kDDLOffset);

    ReadGTUHeaders((UInt_t*) buffer);

    if (fCurrTrailerReadout)
      ReadGTUTrailer();

    // loop over all active links
    AliDebug(2, Form("Stack mask 0x%02x", fCurrStackMask));
    for (Int_t iStack = 0; iStack < fgkNstacks; iStack++) {
      fCurrSlot = iStack;
      if ((fCurrStackMask & (1 << fCurrSlot)) == 0)
	continue;

      AliDebug(2, Form("Stack %i, Link mask: 0x%02x", fCurrSlot, fCurrLinkMask[fCurrSlot]));
      for (Int_t iLink = 0; iLink < fgkNlinks; iLink++) {
	fCurrLink = iLink;
	fCurrHC   = (fCurrEquipmentId - kDDLOffset) * fgkNstacks * fgkNlinks +
	  fCurrSlot * fgkNlinks + iLink;
	if ((fCurrLinkMask[fCurrSlot] & (1 << fCurrLink)) == 0)
	  continue;

	fErrorFlags = 0;
	// check for link monitor error flag
	if (fCurrLinkMonitorFlags[fCurrSlot*fgkNlinks + fCurrLink] != 0) {
	  LinkError(kLinkMonitor);
	  if (fgErrorBehav[kLinkMonitor] == kTolerate)
	    ReadLinkData();
	}
	else
	  // read the data from one HC
	  ReadLinkData();

	// read all data endmarkers
	SeekNextLink();
      }
    }
  }

  return kTRUE;
}


Bool_t AliTRDrawStream::NextDDL()
{
  // continue reading with the next equipment

  if (!fRawReader)
    return kFALSE;

  fCurrEquipmentId = 0;
  fCurrSlot = 0;
  fCurrLink = 0;

  UChar_t *buffer = 0x0;

  while (fRawReader->ReadNextData(buffer)) {

    fCurrEquipmentId = fRawReader->GetEquipmentId();
    AliDebug(2, Form("equipment: %i", fCurrEquipmentId));

    if (fCurrEquipmentId < kDDLOffset || fCurrEquipmentId > kDDLMax) {
      EquipmentError(kNonTrdEq, "Skipping");
      continue;
    }

    if (fMarkers)
      new ((*fMarkers)[fMarkers->GetEntriesFast()])
	AliTRDrawStreamError(-kSecactive, fCurrEquipmentId - kDDLOffset);

    ReadGTUHeaders((UInt_t*) buffer);

    if (fCurrTrailerReadout)
      ReadGTUTrailer();

    return kTRUE;
  }

  return kFALSE;
}


Int_t AliTRDrawStream::NextChamber(AliTRDdigitsManager *digMgr)
{
  // read the data for the next chamber
  // in case you only want to read the data of a single chamber
  // to read all data ReadEvent(...) is recommended

  fDigitsManager = digMgr;
  fDigitsParam   = 0x0;

  fErrorFlags = 0;

  // tracklet output preparation
  TTree *trklTree = 0x0;
  AliRunLoader *rl = AliRunLoader::Instance();
  AliLoader* trdLoader = rl ? rl->GetLoader("TRDLoader") : NULL;
  AliDataLoader *trklLoader = trdLoader ? trdLoader->GetDataLoader("tracklets") : NULL;
  if (trklLoader) {
    AliTreeLoader *trklTreeLoader = (AliTreeLoader*) trklLoader->GetBaseLoader("tracklets-raw");
    if (trklTreeLoader)
      trklTree = trklTreeLoader->Tree();
    else
      trklTree = trklLoader->Tree();
  }

  if (fTrackletTree != trklTree)
    ConnectTracklets(trklTree);

  if (!fRawReader) {
    AliError("No raw reader available");
    return -1;
  }

  while (fCurrSlot < 0 || fCurrSlot >= fgkNstacks) {
    if (!NextDDL()) {
      fCurrSlot = -1;
      return -1;
    }
    while ((fCurrSlot < fgkNstacks) &&
	   (((fCurrStackMask & (1 << fCurrSlot)) == 0) ||
	    ((fCurrLinkMask[fCurrSlot] & (1 << fCurrLink))) == 0)) {
      fCurrLink++;
      if (fCurrLink >= fgkNlinks) {
	fCurrLink = 0;
	fCurrSlot++;
      }
    }
  }

  AliDebug(2, Form("Stack %i, Link %i, mask: 0x%02x", fCurrSlot, fCurrLink, fCurrLinkMask[fCurrSlot]));
  fCurrHC   = (fCurrEquipmentId - kDDLOffset) * fgkNlinks * fgkNstacks +
    fCurrSlot * fgkNlinks + fCurrLink;

  if (fCurrLinkMonitorFlags[fCurrSlot*fgkNlinks + fCurrLink] != 0) {
    LinkError(kLinkMonitor);
    if (fgErrorBehav[kLinkMonitor] == kTolerate)
      ReadLinkData();
  }
  else
    // read the data from one HC
    ReadLinkData();

  // read all data endmarkers
  SeekNextLink();

  if (fCurrLink % 2 == 0) {
    // if we just read the A-side HC then also check the B-side
    fCurrLink++;
    fCurrHC++;
    if (fCurrLinkMask[fCurrSlot] & (1 << fCurrLink)) {
      if (fCurrLinkMonitorFlags[fCurrSlot*fgkNlinks + fCurrLink] != 0) {
	LinkError(kLinkMonitor);
	if (fgErrorBehav[kLinkMonitor] == kTolerate)
	  ReadLinkData();
      }
      else {
	ReadLinkData();
      }
      SeekNextLink();
    }
  }

  //??? to check
  do {
    fCurrLink++;
    if (fCurrLink >= fgkNlinks) {
      fCurrLink = 0;
      fCurrSlot++;
    }
  } while ((fCurrSlot < fgkNstacks) &&
	   (((fCurrStackMask & (1 << fCurrSlot)) == 0) ||
	    ((fCurrLinkMask[fCurrSlot] & (1 << fCurrLink))) == 0));

  // return chamber information from HC if it is valid
  // otherwise return information from link position
  if (fCurrSm < 0 || fCurrSm >= fgkNsectors || fCurrStack < 0 || fCurrStack >= fgkNstacks || fCurrLayer < 0 || fCurrLayer >= fgkNlinks/2)
    return ((fCurrEquipmentId-kDDLOffset) + fCurrSlot * fgkNlinks/2 + fCurrLink/2);
  else
    return (fCurrSm * fgkNstacks*fgkNlinks/2 + fCurrStack * fgkNlinks/2 + fCurrLayer);
}


Int_t AliTRDrawStream::ReadGTUHeaders(UInt_t *buffer)
{
  // check the data source and read the headers

  if (fCurrEquipmentId >= kDDLOffset && fCurrEquipmentId <= kDDLMax) {
    // this is ROC data

    // setting the pointer to data and current reading position
    fPayloadCurr = fPayloadStart = buffer;
    fPayloadSize = fRawReader->GetDataSize() / sizeof(UInt_t);
    fStats.fStatsSector[fCurrEquipmentId - kDDLOffset].fBytes = fRawReader->GetDataSize();
    AliDebug(2, Form("Read buffer of size: %i", fRawReader->GetDataSize()));

    AliDebug(1, DumpRaw("raw data", fPayloadCurr, TMath::Min(fPayloadSize, 1000)));

    // read SM header
    if (ReadSmHeader() < 0) {
      AliError(Form("Reading SM header failed, skipping this DDL %i", fCurrEquipmentId));
      return -1;
    }

    // read tracking headers (if available)
    if (fCurrTrkHeaderAvail) {
      for (Int_t iStack = 0; iStack < fgkNstacks; iStack++) {
	if ((fCurrStackMask & (1 << iStack)) != 0)
	  ReadTrackingHeader(iStack);
      }
    }

    // read trigger header(s) (if available)
    if (fCurrTrgHeaderAvail)
      ReadTriggerHeaders();

    // read stack header
    for (Int_t iStack = 0; iStack < fgkNstacks; iStack++) {
      if ((fCurrStackMask & (1 << iStack)) != 0)
	ReadStackHeader(iStack);
    }

    return 0;
  }
  else
    return -1;
}

Int_t AliTRDrawStream::ReadSmHeader()
{
  // read the SMU index header at the current reading position
  // and store the information in the corresponding variables

  if (fPayloadCurr - fPayloadStart >= fPayloadSize - 1) {
    EquipmentError(kUnknown, "SM Header incomplete");
    return -1;
  }

  fCurrTrgFlags[fCurrEquipmentId-kDDLOffset] = 0;

  fCurrSmHeaderSize           = ((*fPayloadCurr) >> 16) & 0xffff;
  fCurrSmHeaderVersion        = ((*fPayloadCurr) >> 12) &    0xf;
  fCurrTrackEnable            = ((*fPayloadCurr) >>  6) &    0x1;
  fCurrTrackletEnable         = ((*fPayloadCurr) >>  5) &    0x1;
  fCurrStackMask              = ((*fPayloadCurr)      ) &   0x1f;
  fCurrHwRev                  = (fPayloadCurr[1] >> 12) & 0xffff;

  switch (fCurrSmHeaderVersion) {
  case 0xb:
    fCurrTrailerReadout = 0;
    fCurrTrgHeaderAvail = 0;
    fCurrEvType = 0;
    fCurrTrkHeaderAvail = 0;

    DecodeGTUtracks();
    break;

  case 0xc:
    fCurrTrailerReadout = ((*fPayloadCurr) >> 10) &    0x1;
    fCurrTrgHeaderAvail = 1;
    fCurrTrgHeaderReadout = ((*fPayloadCurr) >>  9) &    0x1;
    fCurrEvType         = ((*fPayloadCurr) >>  7) &    0x3;
    fCurrTrkHeaderAvail = fCurrTrackEnable;
    fCurrTriggerEnable  = (fPayloadCurr[2] >>  8) &  0xfff;
    fCurrTriggerFired   = (fPayloadCurr[2] >>  20) &  0xfff;
    fCurrTrgFlags[fCurrEquipmentId-kDDLOffset] = fCurrTriggerFired;
    break;

  default:
    AliError(Form("unknown SM header version: 0x%x", fCurrSmHeaderVersion));
  }

  AliDebug(5, Form("SM header: size: %i, version: %i, track enable: %i, tracklet enable: %i, stack mask: %2x, trailer: %i, trgheader: %i, trkheader: %i",
		   fCurrSmHeaderSize,
		   fCurrSmHeaderVersion,
		   fCurrTrackEnable,
		   fCurrTrackletEnable,
		   fCurrStackMask,
		   fCurrTrailerReadout,
		   fCurrTrgHeaderAvail,
		   fCurrTrkHeaderAvail ));

  // jump to the first word after the SM header
  fPayloadCurr += fCurrSmHeaderSize + 1;

  return fCurrSmHeaderSize + 1;
}

Int_t AliTRDrawStream::DecodeGTUtracks()
{
  // decode GTU track words
  // this depends on the hardware revision of the SMU

  Int_t sector = fCurrEquipmentId-kDDLOffset;

  if ((sector < 0) || (sector > 17)) {
    AliError(Form("Invalid sector %i for GTU tracks", sector));
    return -1;
  }

  AliDebug(1, DumpRaw(Form("GTU tracks in sector %2i (hw rev %i)", sector, fCurrHwRev),
		      fPayloadCurr + 4, 10, 0xffe0ffff));

  fCurrTrgFlags[sector] = 0;

  if (fCurrHwRev < 1772) {
    UInt_t    fastWord;		// fast trigger word
    ULong64_t trackWord = 0;	// extended track word
    Int_t stack = 0;
    Int_t idx = 0;
    for (UInt_t iWord = 4; iWord < fCurrSmHeaderSize; iWord++) {
      if (fPayloadCurr[iWord] == 0x10000000) { // stack boundary marker
        stack++;
        idx = 0;
      }
      else {
        if ((idx == 0) &&
	    ((fPayloadCurr[iWord] & 0xfffff0f0) == 0x13370000)) {
	  fastWord = fPayloadCurr[iWord];
	  fCurrTrgFlags[sector] |= 1 << (stack+11); // assume tracking done if fast word sent
	  AliDebug(1, Form("stack %i: fast trigger word: 0x%08x", stack, fastWord));
	  continue;
        }
        else if ((idx & 0x1) == 0x1) {
	  trackWord |= ((ULong64_t) fPayloadCurr[iWord]) << 32;
	  AliDebug(1,Form("track debug word: 0x%016llx", trackWord));
	  if (fTracks) {
	    AliESDTrdTrack *trk = new ((*fTracks)[fTracks->GetEntriesFast()])
	      AliESDTrdTrack();

	    trk->SetSector(sector);
	    trk->SetStack((trackWord >> 60) & 0x7);
	    trk->SetA(0);
	    trk->SetB(0);
	    trk->SetPID(0);
	    trk->SetLayerMask((trackWord >> 16) & 0x3f);
	    trk->SetTrackletIndex((trackWord >> 22) & 0x3f, 0);
	    trk->SetTrackletIndex((trackWord >> 28) & 0x3f, 1);
	    trk->SetTrackletIndex((trackWord >> 34) & 0x3f, 2);
	    trk->SetTrackletIndex((trackWord >> 40) & 0x3f, 3);
	    trk->SetTrackletIndex((trackWord >> 46) & 0x3f, 4);
	    trk->SetTrackletIndex((trackWord >> 52) & 0x3f, 5);

	    trk->SetFlags(0);
	    trk->SetReserved(0);

	    Float_t pt = (((Int_t) (trackWord & 0xffff) ^ 0x8000) - 0x8000)/128.;
	    if (TMath::Abs(pt) > 0.1) {
	      trk->SetA((Int_t) (0.15*51625./100./pt / 160e-4 * 2));
	    }
	  }
	}
        else {
	  trackWord = fPayloadCurr[iWord];
        }
        idx++;
      }
    }
  }
  else if (fCurrHwRev < 1804) {
    UInt_t    fastWord;		// fast trigger word
    ULong64_t trackWord = 0;	// extended track word
    Int_t stack = 0;
    Int_t idx = 0;
    for (UInt_t iWord = 4; iWord < fCurrSmHeaderSize; iWord++) {
      if (fPayloadCurr[iWord] == 0xffe0ffff) { // stack boundary marker
        stack++;
        idx = 0;
      }
      else {
        if ((idx == 0) &&
	    ((fPayloadCurr[iWord] & 0xfffff0f0) == 0x13370000)) {
	  fastWord = fPayloadCurr[iWord];
	  fCurrTrgFlags[sector] |= 1 << (stack+11); // assume tracking done if fast word sent
	  AliDebug(1, Form("stack %i: fast trigger word: 0x%08x", stack, fastWord));
	  continue;
        }
        else if ((idx & 0x1) == 0x1) {
	  trackWord |= ((ULong64_t) fPayloadCurr[iWord]) << 32;
	  AliDebug(1, Form("track debug word: 0x%016llx", trackWord));
	  if (fTracks) {
	    AliESDTrdTrack *trk = new ((*fTracks)[fTracks->GetEntriesFast()])
	      AliESDTrdTrack();

	    trk->SetSector(fCurrEquipmentId-kDDLOffset);
	    trk->SetStack((trackWord >> 60) & 0x7);
	    trk->SetA(0);
	    trk->SetB(0);
	    trk->SetPID(0);
	    trk->SetLayerMask((trackWord >> 16) & 0x3f);
	    trk->SetTrackletIndex((trackWord >> 22) & 0x3f, 0);
	    trk->SetTrackletIndex((trackWord >> 28) & 0x3f, 1);
	    trk->SetTrackletIndex((trackWord >> 34) & 0x3f, 2);
	    trk->SetTrackletIndex((trackWord >> 40) & 0x3f, 3);
	    trk->SetTrackletIndex((trackWord >> 46) & 0x3f, 4);
	    trk->SetTrackletIndex((trackWord >> 52) & 0x3f, 5);

	    trk->SetFlags(0);
	    trk->SetReserved(0);

	    Float_t pt = (((Int_t) (trackWord & 0xffff) ^ 0x8000) - 0x8000)/128.;
	    if (TMath::Abs(pt) > 0.1) {
	      trk->SetA((Int_t) (-0.15*51625./100./pt / 160e-4 * 2));
	    }
	  }
        }
        else {
	  trackWord = fPayloadCurr[iWord];
        }
        idx++;
      }
    }
  }
  else if (fCurrHwRev < 1819) {
    UInt_t    fastWord;		// fast trigger word
    ULong64_t trackWord = 0;	// extended track word
    Int_t stack = 0;
    Int_t idx = 0;
    for (UInt_t iWord = 4; iWord < fCurrSmHeaderSize; iWord++) {
      if (fPayloadCurr[iWord] == 0xffe0ffff) { // stack boundary marker
	stack++;
	idx = 0;
      }
      else {
	if ((idx == 0) &&
	    ((fPayloadCurr[iWord] & 0xfffff0f0) == 0x13370000)) {
	  fastWord = fPayloadCurr[iWord];
	  if (fastWord & (1 << 13))
	    fCurrTrgFlags[sector] |= 1 << (stack+11);
	  AliDebug(1, Form("stack %i: fast trigger word: 0x%08x", stack, fastWord));
	  continue;
	}
	else if ((idx & 0x1) == 0x1) {
	  trackWord |= ((ULong64_t) fPayloadCurr[iWord]) << 32;
	  AliDebug(1, Form("track debug word: 0x%016llx", trackWord));

	  if (fTracks) {
	    AliESDTrdTrack *trk = new ((*fTracks)[fTracks->GetEntriesFast()])
	      AliESDTrdTrack();

	    trk->SetSector(fCurrEquipmentId-kDDLOffset);
	    trk->SetStack((trackWord >> 60) & 0x7);
	    trk->SetA(0);
	    trk->SetB(0);
	    // trk->SetPt(((trackWord & 0xffff) ^ 0x8000) - 0x8000);
	    trk->SetPID(0);
	    trk->SetLayerMask((trackWord >> 16) & 0x3f);
	    trk->SetTrackletIndex((trackWord >> 22) & 0x3f, 0);
	    trk->SetTrackletIndex((trackWord >> 28) & 0x3f, 1);
	    trk->SetTrackletIndex((trackWord >> 34) & 0x3f, 2);
	    trk->SetTrackletIndex((trackWord >> 40) & 0x3f, 3);
	    trk->SetTrackletIndex((trackWord >> 46) & 0x3f, 4);
	    trk->SetTrackletIndex((trackWord >> 52) & 0x3f, 5);

	    trk->SetFlags(0);
	    trk->SetReserved(0);

	    Float_t pt = (((Int_t) (trackWord & 0xffff) ^ 0x8000) - 0x8000)/128.;
	    if (TMath::Abs(pt) > 0.1) {
	      trk->SetA((Int_t) (0.15*51625./100./trk->Pt() / 160e-4 * 2));
	    }
	  }
	}
	else {
	  trackWord = fPayloadCurr[iWord];
	}
	idx++;
      }
    }
  }
  else if (fCurrHwRev < 1860) {
    UInt_t    fastWord;		// fast trigger word
    ULong64_t trackWord = 0;	// extended track word
    Int_t stack = 0;
    Int_t idx = 0;
    Bool_t upperWord = kFALSE;
    Int_t word = 0;
    for (UInt_t iWord = 4; iWord < fCurrSmHeaderSize; iWord++) {
      if (fPayloadCurr[iWord] == 0xffe0ffff) { // stack boundary marker
        stack++;
        idx = 0;
	upperWord = kFALSE;
      }
      else {
	// assemble the 32-bit words out of 16-bit blocks
	if (upperWord) {
	  word |= (fPayloadCurr[iWord] & 0xffff0000);
	  upperWord = kFALSE;
	}
	else {
	  // lower word is read first
	  word = (fPayloadCurr[iWord] & 0xffff0000) >> 16;
	  upperWord = kTRUE;
	  continue;
	}

        if ((word & 0xffff0008) == 0x13370008) {
	  fastWord = word;
	  AliDebug(1, Form("stack %i: fast track word: 0x%08x", stack, fastWord));
	  if (fastWord & (1 << 13))
	    fCurrTrgFlags[sector] |= 1 << (stack+11);
	  continue;
        }
        else if ((idx & 0x1) == 0x1) {
	  trackWord |= ((ULong64_t) word) << 32;
	  AliDebug(1, Form("track debug word: 0x%016llx", trackWord));
	  if (fTracks) {
	    AliESDTrdTrack *trk = new ((*fTracks)[fTracks->GetEntriesFast()])
	      AliESDTrdTrack();

	    trk->SetSector(fCurrEquipmentId-kDDLOffset);
	    trk->SetStack((trackWord >> 60) & 0x7);
	    trk->SetA(0);
	    trk->SetB(0);
	    trk->SetPID(0);
	    trk->SetLayerMask((trackWord >> 16) & 0x3f);
	    trk->SetTrackletIndex((trackWord >> 22) & 0x3f, 0);
	    trk->SetTrackletIndex((trackWord >> 28) & 0x3f, 1);
	    trk->SetTrackletIndex((trackWord >> 34) & 0x3f, 2);
	    trk->SetTrackletIndex((trackWord >> 40) & 0x3f, 3);
	    trk->SetTrackletIndex((trackWord >> 46) & 0x3f, 4);
	    trk->SetTrackletIndex((trackWord >> 52) & 0x3f, 5);

	    trk->SetFlags(0);
	    trk->SetReserved(0);

	    Float_t pt = (((Int_t) (trackWord & 0xffff) ^ 0x8000) - 0x8000)/128.;
	    if (TMath::Abs(pt) > 0.1) {
	      trk->SetA((Int_t) (0.15*51625./100./pt / 160e-4 * 2));
	    }
	  }
        }
        else {
	  trackWord = word;
        }
        idx++;
      }
    }

  }
  else {
    ULong64_t trackWord = 0; // this is the debug word
    Int_t stack = 0;
    Int_t idx = 0;
    Bool_t upperWord = kFALSE;
    Int_t word = 0;
    for (UInt_t iWord = 4; iWord < fCurrSmHeaderSize; iWord++) {
      if (fPayloadCurr[iWord] == 0xffe0ffff) {
        stack++;
        idx = 0;
	upperWord = kFALSE;
      }
      else {
	// assemble the 32-bit words out of 16-bit blocks
	if (upperWord) {
	  word |= (fPayloadCurr[iWord] & 0xffff0000);
	  upperWord = kFALSE;
	}
	else {
	  // lower word is read first
	  word = (fPayloadCurr[iWord] & 0xffff0000) >> 16;
	  upperWord = kTRUE;
	  continue;
	}

        if ((word & 0xffff0008) == 0x13370008) {
	  AliDebug(1, Form("stack %i: fast track word: 0x%08x", stack, word));
	  continue;
        }
        else if ((word & 0xffff0010) == 0x13370010) {
	  AliDebug(1, Form("stack %i: tracking done word: 0x%08x", stack, word));
	  fCurrTrgFlags[sector] |= 1 << (stack+11);
	  continue;
	}
        else if ((idx & 0x1) == 0x1) {
	  trackWord |= ((ULong64_t) word) << 32;
	  AliDebug(1, Form("track debug word: 0x%16llx", trackWord));
	  if (fTracks) {
	    AliESDTrdTrack *trk = new ((*fTracks)[fTracks->GetEntriesFast()])
	      AliESDTrdTrack();
	    trk->SetSector(fCurrEquipmentId-kDDLOffset);
	    trk->SetStack((trackWord >> 60) & 0x7);
	    trk->SetA(0);
	    trk->SetB(0);
	    trk->SetPID(0);
	    trk->SetLayerMask((trackWord >> 16) & 0x3f);
	    trk->SetTrackletIndex((trackWord >> 22) & 0x3f, 0);
	    trk->SetTrackletIndex((trackWord >> 28) & 0x3f, 1);
	    trk->SetTrackletIndex((trackWord >> 34) & 0x3f, 2);
	    trk->SetTrackletIndex((trackWord >> 40) & 0x3f, 3);
	    trk->SetTrackletIndex((trackWord >> 46) & 0x3f, 4);
	    trk->SetTrackletIndex((trackWord >> 52) & 0x3f, 5);

	    trk->SetFlags(0);
	    trk->SetReserved(0);

	    Float_t pt = (((Int_t) (trackWord & 0xffff) ^ 0x8000) - 0x8000)/128.;
	    if (TMath::Abs(pt) > 0.1) {
	      trk->SetA(-(Int_t) (0.15*51625./100./pt / 160e-4 * 2));
	    }
	  }
        }
        else {
	  trackWord = word;
        }
        idx++;
      }
    }
  }
  return 0;
}

Int_t AliTRDrawStream::ReadTrackingHeader(Int_t stack)
{
  // read the tracking information and store it for the given stack

  // index word

  fCurrTrkHeaderIndexWord[stack] = *fPayloadCurr;
  fCurrTrkHeaderSize[stack]      = ((*fPayloadCurr) >> 16) & 0x3ff;

  AliDebug(1, Form("tracking header index word: 0x%08x, size: %i (hw rev: %i)",
		   fCurrTrkHeaderIndexWord[stack], fCurrTrkHeaderSize[stack], fCurrHwRev));
  Int_t trackingTime = *fPayloadCurr & 0x3ff;

  fCurrTrgFlags[fCurrEquipmentId-kDDLOffset] |= ((fCurrTrkHeaderIndexWord[stack] >> 10) & 0x1) << (22 + stack);
  fPayloadCurr++;

  // data words
  ULong64_t trackWord = 0;
  Int_t idx = 0;
  Int_t trackIndex = fTracks ? fTracks->GetEntriesFast() : -1;

  for (UInt_t iWord = 0; iWord < fCurrTrkHeaderSize[stack]; iWord++) {

    if (!(idx & 0x1)) {
      // first part of 64-bit word
      trackWord = fPayloadCurr[iWord];
    }
    else {
      trackWord |= ((ULong64_t) fPayloadCurr[iWord]) << 32;

      if (trackWord & (1ul << 63)) {
	if ((trackWord & (0x3ful << 56)) != 0) {
	  // track word
	  AliDebug(2, Form("track word: 0x%016llx", trackWord));

	  if (fTracks) {
	    AliESDTrdTrack *trk = new ((*fTracks)[fTracks->GetEntriesFast()])
	      AliESDTrdTrack();

	    trk->SetSector(fCurrEquipmentId-kDDLOffset);
	    trk->SetLayerMask((trackWord >> 56) & 0x3f);
	    trk->SetA( (((trackWord >> 38) & 0x3ffff) ^ 0x20000) - 0x20000);
	    trk->SetB( (((trackWord >> 20) & 0x3ffff) ^ 0x20000) - 0x20000);
	    trk->SetC( (((trackWord >> 8)  &  0xffff) ^  0x8000) -  0x8000);
	    trk->SetPID((trackWord >>  0) & 0xff);
	    trk->SetStack(stack);

	    // now compare the track word with the one generated from the ESD information
	    if (trackWord != trk->GetTrackWord(0)) {
	      AliError(Form("track word 0x%016llx does not match the read one 0x%016llx",
			    trk->GetTrackWord(0), trackWord));
	    }
	  }
	}
	else {
	  // done marker (so far only used to set trigger flag)
	  fCurrTrgFlags[fCurrEquipmentId-kDDLOffset] |= 1 << (27 + stack);
	  fCurrTrkFlags[(fCurrEquipmentId-kDDLOffset)*fgkNstacks + stack] = trackWord;

	  AliDebug(2, Form("tracking done marker: 0x%016llx, trigger flags: 0x%08x",
			   trackWord, fCurrTrgFlags[fCurrEquipmentId-kDDLOffset]));
	  AliDebug(2, Form("seg / stack / first / last / done / index : %i %i %lli %lli %lli %i",
			   fCurrEquipmentId - kDDLOffset, stack,
			   (trackWord >> 20) & 0x3ff,
			   (trackWord >> 10) & 0x3ff,
			   (trackWord >>  0) & 0x3ff,
			   trackingTime));
	}
      }
      else {
	// extended track word
	AliDebug(2, Form("extended track word: 0x%016llx", trackWord));

	if (fTracks) {
	  AliESDTrdTrack *trk = (AliESDTrdTrack*) (*fTracks)[trackIndex];

	  trk->SetFlags((trackWord >> 52) & 0x7ff);
	  trk->SetReserved((trackWord >> 49) & 0x7);
	  trk->SetY((trackWord >> 36) & 0x1fff);
	  trk->SetTrackletIndex((trackWord >>  0) & 0x3f, 0);
	  trk->SetTrackletIndex((trackWord >>  6) & 0x3f, 1);
	  trk->SetTrackletIndex((trackWord >> 12) & 0x3f, 2);
	  trk->SetTrackletIndex((trackWord >> 18) & 0x3f, 3);
	  trk->SetTrackletIndex((trackWord >> 24) & 0x3f, 4);
	  trk->SetTrackletIndex((trackWord >> 30) & 0x3f, 5);

	  if (trackWord != trk->GetExtendedTrackWord(0)) {
	    AliError(Form("extended track word 0x%016llx does not match the read one 0x%016llx",
			  trk->GetExtendedTrackWord(0), trackWord));
	    }

	  trackIndex++;
	}
      }
    }
    idx++;
  }

  fPayloadCurr += fCurrTrkHeaderSize[stack];

  return fCurrTrkHeaderSize[stack];
}

Int_t AliTRDrawStream::ReadTriggerHeaders()
{
  // read all trigger headers present

  AliDebug(1, Form("trigger mask: 0x%03x, fired: 0x%03x\n",
		   fCurrTriggerEnable, fCurrTriggerFired));
  // loop over potential trigger blocks
  for (Int_t iTrigger = 0; iTrigger < fgkNtriggers; iTrigger++) {
    // check for trigger enable
    if (fCurrTriggerEnable & (1 << iTrigger)) {
      // check for readout mode and trigger fired
      if ((fCurrTrgHeaderReadout == 0) || (fCurrTriggerFired & (1 << iTrigger))) {
	// index word
	AliDebug(1, Form("trigger index word %i: 0x%08x\n", iTrigger, *fPayloadCurr));
	fCurrTrgHeaderIndexWord[iTrigger] = *fPayloadCurr;
	fCurrTrgHeaderSize[iTrigger]      = ((*fPayloadCurr) >> 16) & 0x3ff;
	if (iTrigger == 7) {
	  // timeout trigger, use to extract tracking time
	  fCurrTrgFlags[fCurrEquipmentId-kDDLOffset] |= (*fPayloadCurr & 0x3ff) << 12;
	}

	fPayloadCurr++;
	// data words
	fPayloadCurr += fCurrTrgHeaderSize[iTrigger];
      }
    }
  }

  return 0;
}

Int_t AliTRDrawStream::ReadStackHeader(Int_t stack)
{
  // read the stack header
  // and store the information in the corresponding variables

  fCurrStackIndexWord[stack]     = *fPayloadCurr;
  fCurrStackHeaderSize[stack]    = (((*fPayloadCurr) >> 16) & 0xffff) + 1;
  fCurrStackHeaderVersion[stack] = ((*fPayloadCurr) >> 12) & 0xf;
  fCurrLinkMask[stack]           = (*fPayloadCurr) & 0xfff;

  // dumping stack header
  AliDebug(1, DumpRaw(Form("stack %i header", stack), fPayloadCurr, fCurrStackHeaderSize[stack]));

  if (fPayloadCurr - fPayloadStart >= fPayloadSize - (Int_t) fCurrStackHeaderSize[stack]) {
    LinkError(kStackHeaderInvalid, "Stack index header aborted");
    return -1;
  }

  switch (fCurrStackHeaderVersion[stack]) {
  case 0xa:
    if (fCurrStackHeaderSize[stack] < 8) {
      LinkError(kStackHeaderInvalid, "Stack header smaller than expected!");
      return -1;
    }

    fCurrCleanCheckout[stack] = fPayloadCurr[1] & 0x1;
    fCurrBoardId[stack]       = (fPayloadCurr[1] >> 8) & 0xff;
    fCurrHwRevTMU[stack]      = (fPayloadCurr[1] >> 16) & 0xffff;

    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
      // A side
      fCurrLinkMonitorFlags  [stack*fgkNlinks + iLayer*2]      = fPayloadCurr[iLayer+2] & 0xf;
      fCurrLinkDataTypeFlags [stack*fgkNlinks + iLayer*2]      = (fPayloadCurr[iLayer+2] >> 4) & 0x3;
      fCurrLinkDebugFlags    [stack*fgkNlinks + iLayer*2]      = (fPayloadCurr[iLayer+2] >> 12) & 0xf;
      // B side
      fCurrLinkMonitorFlags  [stack*fgkNlinks + iLayer*2 + 1]  = (fPayloadCurr[iLayer+2] >> 16) & 0xf;
      fCurrLinkDataTypeFlags [stack*fgkNlinks + iLayer*2 + 1]  = (fPayloadCurr[iLayer+2] >> 20) & 0x3;
      fCurrLinkDebugFlags    [stack*fgkNlinks + iLayer*2 + 1]  = (fPayloadCurr[iLayer+2] >> 28) & 0xf;
    }
    break;

  default:
    LinkError(kStackHeaderInvalid, "Invalid Stack Header version %x", fCurrStackHeaderVersion[stack]);
  }

  fPayloadCurr += fCurrStackHeaderSize[stack];

  return fCurrStackHeaderSize[stack];
}

Int_t AliTRDrawStream::ReadGTUTrailer()
{
  // read the SM trailer containing CRCs from various stages

  UInt_t* trailer = fPayloadStart + fPayloadSize -1;

  // look for the trailer index word from the end
  for (Int_t iWord = 0; iWord < fPayloadSize; iWord++) {
    if ((fPayloadStart[fPayloadSize-1-iWord] & 0xffff) == 0x1f51) {
      trailer = fPayloadStart + fPayloadSize - 1 - iWord;
      break;
    }
  }

  if (((*trailer) & 0xffff) == 0x1f51) {
    UInt_t trailerIndexWord = (*trailer);
    Int_t trailerSize = (trailerIndexWord >> 16) & 0xffff;
    AliDebug(2, DumpRaw("GTU trailer", trailer, trailerSize+1));
    // parse the trailer
  }
  else
    EquipmentError(kUnknown, "trailer index marker mismatch");

  return 0;
}

Int_t AliTRDrawStream::ReadLinkData()
{
  // read the data in one link (one HC) until the data endmarker is reached
  // returns the number of words read!

  Int_t count = 0;
  UInt_t* startPosLink = fPayloadCurr;

  AliDebug(1, DumpRaw(Form("link data from seg %2i slot %i link %2i", fCurrEquipmentId-kDDLOffset, fCurrSlot, fCurrLink),
		      fPayloadCurr, TMath::Min((Int_t) (fPayloadSize - (fPayloadCurr-fPayloadStart)), 100), 0x00000000));

  if (fMarkers)
    new ((*fMarkers)[fMarkers->GetEntriesFast()])
      AliTRDrawStreamError(-kHCactive, fCurrEquipmentId-kDDLOffset, fCurrSlot, fCurrLink);

  if (fErrorFlags & kDiscardHC)
    return count;

  //??? add check whether tracklets are enabled
  count += ReadTracklets();
  if (fErrorFlags & kDiscardHC)
    return count;

  AliDebug(1, DumpRaw("HC header", fPayloadCurr, 4, 0x00000000));
  count += ReadHcHeader();
  if (fErrorFlags & kDiscardHC)
    return count;

  Int_t det = fCurrSm * 30 + fCurrStack * 6 + fCurrLayer;

  if (det > -1 && det < 540) {

    // ----- check which kind of data -----
    if (fCurrMajor & 0x40) {
      if ((fCurrMajor & 0x7) == 0x7) {
	AliDebug(1, "This is a config event");
	UInt_t *startPos = fPayloadCurr;
	while (fPayloadCurr - fPayloadStart < fPayloadSize &&
	       *fPayloadCurr != fgkDataEndmarker)
	  fPayloadCurr++;
	count += fPayloadCurr - startPos;

	// feeding TRAP config
	AliTRDtrapConfig *trapcfg = AliTRDtrapConfig::Instance();
	trapcfg->ReadPackedConfig(fCurrHC, startPos, fPayloadCurr - startPos);
      }
      else {
	Int_t tpmode = fCurrMajor & 0x7;
	AliDebug(1, Form("Checking testpattern (mode %i) data", tpmode));
	ReadTPData(tpmode);
      }
    }
    else {
      // reading real data
      if (fDigitsManager) {
	if ((fAdcArray = fDigitsManager->GetDigits(det))) {
	  //fAdcArray->Expand();
	  if (fAdcArray->GetNtime() != fCurrNtimebins)
	    fAdcArray->Allocate(16, 144, fCurrNtimebins);
	}
	else {
	  LinkError(kNoDigits);
	}

	if (!fDigitsParam) {
	  fDigitsParam = fDigitsManager->GetDigitsParam();
	}
	if (fDigitsParam) {
	  fDigitsParam->SetPretriggerPhase(det, fCurrPtrgPhase);
	  fDigitsParam->SetNTimeBins(det, fCurrNtimebins);
	  fDigitsParam->SetADCbaseline(det, 10);
	}

	if (fDigitsManager->UsesDictionaries()) {
	  fDigitsManager->GetDictionary(det, 0)->Reset();
	  fDigitsManager->GetDictionary(det, 1)->Reset();
	  fDigitsManager->GetDictionary(det, 2)->Reset();
	}

	if ((fSignalIndex = fDigitsManager->GetIndexes(det))) {
	  fSignalIndex->SetSM(fCurrSm);
	  fSignalIndex->SetStack(fCurrStack);
	  fSignalIndex->SetLayer(fCurrLayer);
	  fSignalIndex->SetDetNumber(det);
	  if (!fSignalIndex->IsAllocated())
	    fSignalIndex->Allocate(16, 144, fCurrNtimebins);
	}

	if (fCurrMajor & 0x20) {
	  AliDebug(1, "This is a zs event");
	  count += ReadZSData();
	}
	else {
	  AliDebug(1, "This is a nozs event");
	  count += ReadNonZSData();
	}
      }
      else {
	// just read until data endmarkers
	while (fPayloadCurr - fPayloadStart < fPayloadSize &&
	       *fPayloadCurr != fgkDataEndmarker)
	  fPayloadCurr++;
      }
    }
  }
  else {
    LinkError(kInvalidDetector, "%i", det);
    while (fPayloadCurr - fPayloadStart < fPayloadSize &&
	   *fPayloadCurr != fgkDataEndmarker)
      fPayloadCurr++;
  }

  if (fCurrSm > -1 && fCurrSm < 18) {
    fStats.fStatsSector[fCurrSm].fStatsHC[fCurrHC%60].fBytes     += (fPayloadCurr - startPosLink) * sizeof(UInt_t);
    fStats.fStatsSector[fCurrSm].fStatsHC[fCurrHC%60].fBytesRead += count * sizeof(UInt_t);
    fStats.fStatsSector[fCurrSm].fBytesRead                      += count * sizeof(UInt_t);
    fStats.fBytesRead                                            += count * sizeof(UInt_t);
  }

  return count;
}

Int_t AliTRDrawStream::ReadTracklets()
{
  // read the tracklets from one HC

  fTrackletArray->Clear();

  UInt_t *start = fPayloadCurr;
  while (*(fPayloadCurr) != fgkTrackletEndmarker &&
	 fPayloadCurr - fPayloadStart < fPayloadSize) {
    new ((*fTrackletArray)[fTrackletArray->GetEntriesFast()]) AliTRDtrackletWord(*(fPayloadCurr), fCurrHC);

    fPayloadCurr++;
  }

  if (fTrackletArray->GetEntriesFast() > 0) {
    AliDebug(1, Form("Found %i tracklets in %i %i %i (ev. %i)", fTrackletArray->GetEntriesFast(),
		     (fCurrEquipmentId-kDDLOffset), fCurrSlot, fCurrLink, fRawReader->GetEventIndex()));
    if (fCurrSm > -1 && fCurrSm < 18) {
      fStats.fStatsSector[fCurrSm].fStatsHC[fCurrHC%60].fNTracklets += fTrackletArray->GetEntriesFast();
      fStats.fStatsSector[fCurrSm].fNTracklets                      += fTrackletArray->GetEntriesFast();
    }
    if (fTrackletTree)
      fTrackletTree->Fill();
    if (fTracklets)
      for (Int_t iTracklet = 0; iTracklet < fTrackletArray->GetEntriesFast(); iTracklet++) {
	new ((*fTracklets)[fTracklets->GetEntriesFast()]) AliTRDtrackletWord(*((AliTRDtrackletWord*)(*fTrackletArray)[iTracklet]));
      }
  }

  // loop over remaining tracklet endmarkers
  while ((*(fPayloadCurr) == fgkTrackletEndmarker &&
	  fPayloadCurr - fPayloadStart < fPayloadSize))
    fPayloadCurr++;

  return fPayloadCurr - start;
}

Int_t AliTRDrawStream::ReadHcHeader()
{
  // read and parse the HC header of one HC
  // and store the information in the corresponding variables

  AliDebug(1, Form("HC header: 0x%08x", *fPayloadCurr));
  UInt_t *start = fPayloadCurr;
  // check not to be at the data endmarker
  if (*fPayloadCurr == fgkDataEndmarker)
    return 0;

  fCurrSpecial    = (*fPayloadCurr >> 31) & 0x1;
  fCurrMajor      = (*fPayloadCurr >> 24) & 0x7f;
  fCurrMinor      = (*fPayloadCurr >> 17) & 0x7f;
  fCurrAddHcWords = (*fPayloadCurr >> 14) & 0x7;
  fCurrSm         = (*fPayloadCurr >> 9) & 0x1f;
  fCurrLayer      = (*fPayloadCurr >> 6) & 0x7;
  fCurrStack      = (*fPayloadCurr >> 3) & 0x7;
  fCurrSide       = (*fPayloadCurr >> 2) & 0x1;
  fCurrCheck      = (*fPayloadCurr) & 0x3;

  if ((fCurrSm != (((Int_t) fCurrEquipmentId) - kDDLOffset)) ||
      (fCurrStack != fCurrSlot) ||
      (fCurrLayer != fCurrLink / 2) ||
      (fCurrSide != fCurrLink % 2)) {
    LinkError(kHCmismatch,
	      "HC: %i, %i, %i, %i\n 0x%08x 0x%08x 0x%08x 0x%08x",
	      fCurrSm, fCurrStack, fCurrLayer, fCurrSide,
	      fPayloadCurr[0], fPayloadCurr[1], fPayloadCurr[2], fPayloadCurr[3]);
  }
  if (fCurrCheck != 0x1) {
    LinkError(kHCcheckFailed);
  }

  if (fCurrAddHcWords > 0) {
    fCurrNtimebins = (fPayloadCurr[1] >> 26) & 0x3f;
    fCurrBC = (fPayloadCurr[1] >> 10) & 0xffff;
    fCurrPtrgCnt = (fPayloadCurr[1] >> 6) & 0xf;
    fCurrPtrgPhase = (fPayloadCurr[1] >> 2) & 0xf;
  }

  fPayloadCurr += 1 + fCurrAddHcWords;

  return (fPayloadCurr - start);
}

Int_t AliTRDrawStream::ReadTPData(Int_t mode)
{
  // testing of testpattern 1 to 3 (hardcoded), 0 missing
  // evcnt checking missing
  Int_t cpu = 0;
  Int_t cpufromchannel[] = {0, 0, 0, 0, 0,  1, 1, 1, 1, 1,  2, 2, 2, 2, 2,  3, 3, 3, 3, 3, 3};
  Int_t evcnt = 0;
  Int_t count = 0;
  Int_t mcmcount = -1;
  Int_t wordcount = 0;
  Int_t channelcount = 0;
  UInt_t expword = 0;
  UInt_t expadcval = 0;
  UInt_t diff = 0;
  Int_t lastmcmpos = -1;
  Int_t lastrobpos = -1;

  UInt_t* start = fPayloadCurr;

  while (*(fPayloadCurr) != fgkDataEndmarker &&
	 fPayloadCurr - fPayloadStart < fPayloadSize - 1) {

    // ----- Checking MCM Header -----
    AliDebug(2, DumpMcmHeader("MCM header: ", *fPayloadCurr));
    UInt_t *startPosMCM = fPayloadCurr;
    mcmcount++;

    // ----- checking for proper readout order - ROB -----
    if (GetROBReadoutPos(ROB(*fPayloadCurr) / 2) >= lastrobpos) {
      lastrobpos = GetROBReadoutPos(ROB(*fPayloadCurr) / 2);
    }
    else {
      ROBError(kPosUnexp, Form("#%i after #%i in readout order", GetROBReadoutPos(ROB(*fPayloadCurr) / 2), lastrobpos));
    }
    fCurrRobPos = ROB(*fPayloadCurr);

    // ----- checking for proper readout order - MCM -----
    if (GetMCMReadoutPos(MCM(*fPayloadCurr)) >= (lastmcmpos + 1) % 16) {
      lastmcmpos = GetMCMReadoutPos(MCM(*fPayloadCurr));
    }
    else {
      MCMError(kPosUnexp, Form("#%i after #%i in readout order", GetMCMReadoutPos(MCM(*fPayloadCurr)), lastmcmpos));
    }
    fCurrMcmPos = MCM(*fPayloadCurr);


    fPayloadCurr++;

    evcnt = 0x3f & *fPayloadCurr >> 26;
    cpu = -1;
    channelcount = 0;
    while (channelcount < 21) {
      count = 0;
      if (cpu != cpufromchannel[channelcount]) {
	cpu = cpufromchannel[channelcount];
	expadcval = (1 << 9) | (fCurrRobPos << 6) | (fCurrMcmPos << 2) | cpu;
	wordcount = 0;
      }

      while (count < 10) {
	if (*fPayloadCurr == fgkDataEndmarker) {
	  MCMError(kMissTpData);
	  return (fPayloadCurr - start);
	}

	if (channelcount % 2 == 0)
	  expword = 0x3;
	else
	  expword = 0x2;

	if (mode == 1) {
	  // ----- TP 1 -----
	  expword |= expadcval << 2;
	  expadcval = ( (expadcval << 1) | ( ( (expadcval >> 9) ^ (expadcval >> 6) ) & 1) ) & 0x3FF;
	  expword |= expadcval << 12;
	  expadcval = ( (expadcval << 1) | ( ( (expadcval >> 9) ^ (expadcval >> 6) ) & 1) ) & 0x3FF;
	  expword |= expadcval << 22;
	  expadcval = ( (expadcval << 1) | ( ( (expadcval >> 9) ^ (expadcval >> 6) ) & 1) ) & 0x3FF;
	}
	else if (mode == 2) {
	  // ----- TP 2 ------
	  expword = ((0x3f & evcnt) << 26) | ((fCurrSm + 1) << 21) | ((fCurrLayer + 1) << 18) |
	    ((fCurrStack + 1) << 15) |
	    (fCurrRobPos << 12) | (fCurrMcmPos << 8) | (cpu << 6) | (wordcount + 1);
	}
	else if (mode == 3) {
	  // ----- TP 3 -----
	  expword = ((0xfff & evcnt) << 20) | (fCurrSm << 15) | (fCurrLink/2 << 12) | (fCurrStack << 9) |
	    (fCurrRobPos << 6) | (fCurrMcmPos << 2) | (cpu << 0);
	}
	else {
	  expword = 0;
	  LinkError(kTPmodeInvalid, "Just reading");
	}

	diff = *fPayloadCurr ^ expword;
	AliDebug(11, Form("Comparing ch %2i, word %2i (cpu %i): 0x%08x <-> 0x%08x",
			  channelcount, wordcount, cpu, *fPayloadCurr, expword));

	if (diff != 0) {
	  MCMError(kTPmismatch,
		   "Seen 0x%08x, expected 0x%08x, diff: 0x%08x (0x%02x, 0x%04x) - word %2i (cpu %i, ch %i)",
		   *fPayloadCurr, expword, diff,
		   0xff & (diff | diff >> 8 | diff >> 16 | diff >> 24),
		   0xffff & (diff | diff >> 16),
		   wordcount, cpu, channelcount);;
	}
	fPayloadCurr++;
	count++;
	wordcount++;
      }
      channelcount++;
    }
    // continue with next MCM

    if (IsDumping() && DumpingMCM(fCurrHC/2, fCurrRobPos, fCurrMcmPos)) {
      AliInfo(DumpRaw(Form("Event %i: Det %3i ROB %i MCM %2i", fRawReader->GetEventIndex(), fCurrHC/2, fCurrRobPos, fCurrMcmPos),
		      startPosMCM, fPayloadCurr - startPosMCM));
    }

  }
  return fPayloadCurr - start;
}


Int_t AliTRDrawStream::ReadZSData()
{
  // read the zs data from one link from the current reading position

  UInt_t *start = fPayloadCurr;

  Int_t mcmcount = 0;
  Int_t mcmcountExp = fCurrStack == 2 ? 48 : 64;
  Int_t channelcount = 0;
  Int_t channelcountExp = 0;
  Int_t channelcountMax = 0;
  Int_t timebins;
  Int_t currentTimebin = 0;
  Int_t adcwc = 0;
  Int_t evno = -1;
  Int_t lastmcmpos = -1;
  Int_t lastrobpos = -1;

  if (fCurrNtimebins != fNtimebins) {
    if (fNtimebins > 0)
      LinkError(kNtimebinsChanged,
		"No. of timebins changed from %i to %i", fNtimebins, fCurrNtimebins);
    fNtimebins = fCurrNtimebins;
  }

  timebins = fNtimebins;

  while (*(fPayloadCurr) != fgkDataEndmarker &&
	 fPayloadCurr - fPayloadStart < fPayloadSize) {

    // ----- Checking MCM Header -----
    AliDebug(2, DumpMcmHeader("MCM header: ", *fPayloadCurr));
    UInt_t *startPosMCM = fPayloadCurr;

    // ----- checking for proper readout order - ROB -----
    if (GetROBReadoutPos(ROB(*fPayloadCurr) / 2) >= lastrobpos) {
      if (GetROBReadoutPos(ROB(*fPayloadCurr) / 2) > lastrobpos)
	lastmcmpos = -1;
      lastrobpos = GetROBReadoutPos(ROB(*fPayloadCurr) / 2);
    }
    else {
      ROBError(kPosUnexp, Form("#%i after #%i and #%i in readout order",
			       GetROBReadoutPos(ROB(*fPayloadCurr) / 2), lastrobpos, GetROBReadoutPos(fCurrRobPos)));
    }
    fCurrRobPos = ROB(*fPayloadCurr);

    // ----- checking for proper readout order - MCM -----
    if (GetMCMReadoutPos(MCM(*fPayloadCurr)) > lastmcmpos) {
      lastmcmpos = GetMCMReadoutPos(MCM(*fPayloadCurr));
    }
    else {
      MCMError(kPosUnexp, Form("#%i after #%i and #%i in readout order",
			       GetMCMReadoutPos(MCM(*fPayloadCurr)), lastmcmpos, GetMCMReadoutPos(fCurrMcmPos)));
    }
    fCurrMcmPos = MCM(*fPayloadCurr);

    if (EvNo(*fPayloadCurr) != evno) {
      if (evno == -1)
	evno = EvNo(*fPayloadCurr);
      else {
	MCMError(kPtrgCntMismatch, "%i <-> %i", evno, EvNo(*fPayloadCurr));
      }
    }
    Int_t adccoloff = AdcColOffset(*fPayloadCurr);
    Int_t padcoloff = PadColOffset(*fPayloadCurr);
    Int_t row = Row(*fPayloadCurr);
    fPayloadCurr++;

    // ----- Reading ADC channels -----
    AliDebug(2, DumpAdcMask("ADC mask: ", *fPayloadCurr));

    // ----- analysing the ADC mask -----
    channelcount = 0;
    channelcountExp = GetNActiveChannelsFromMask(*fPayloadCurr);
    channelcountMax = GetNActiveChannels(*fPayloadCurr);
    Int_t channelmask = GetActiveChannels(*fPayloadCurr);
    Int_t channelno = -1;
    fPayloadCurr++;

    if (channelcountExp != channelcountMax) {
      if (channelcountExp > channelcountMax) {
	Int_t temp = channelcountExp;
	channelcountExp = channelcountMax;
	channelcountMax = temp;
      }
      while (channelcountExp < channelcountMax && channelcountExp < 21 &&
	     fPayloadCurr - fPayloadStart < fPayloadSize - 10 * channelcountExp - 1) {
	MCMError(kAdcMaskInconsistent,
		 "Possible MCM-H: 0x%08x, possible ADC-mask: 0x%08x",
		 *(fPayloadCurr + 10 * channelcountExp),
		 *(fPayloadCurr + 10 * channelcountExp + 1) );
	if (!CouldBeMCMhdr( *(fPayloadCurr + 10 * channelcountExp)) && !CouldBeADCmask( *(fPayloadCurr + 10 * channelcountExp + 1)))
	  channelcountExp++;
	else {
	  break;
	}
      }
      MCMError(kAdcMaskInconsistent,
	       "Inconsistency in no. of active channels: Counter: %i, Mask: %i, chosen: %i!",
	       GetNActiveChannels(fPayloadCurr[-1]), GetNActiveChannelsFromMask(fPayloadCurr[-1]), channelcountExp);
    }
    AliDebug(2, Form("expecting %i active channels, %i timebins", channelcountExp, fCurrNtimebins));

    // ----- reading marked ADC channels -----
    while (channelcount < channelcountExp && *(fPayloadCurr) != fgkDataEndmarker) {
      if (channelno < 20)
	channelno++;
      while (channelno < 20 && (channelmask & 1 << channelno) == 0)
	channelno++;

      if (fCurrNtimebins > 30) {
	currentTimebin = ((*fPayloadCurr >> 2) & 0x3f);
	timebins = ((*fPayloadCurr >> 8) & 0xf) * 3;
      }
      else {
	currentTimebin = 0;
      }

      adcwc = 0;
      Int_t nADCwords = (timebins + 2) / 3;
      AliDebug(3, Form("Now reading %i words for channel %2i", nADCwords, channelno));
      Int_t adccol = adccoloff - channelno;
      Int_t padcol = padcoloff - channelno;
//      if (adccol < 3 || adccol > 165)
//	AliInfo(Form("writing channel %i of det %3i %i:%2i to adcrow/-col: %i/%i padcol: %i",
//		     channelno, fCurrHC/2, fCurrRobPos, fCurrMcmPos, row, adccol, padcol));

      while ((adcwc < nADCwords) &&
	     (*(fPayloadCurr) != fgkDataEndmarker) &&
	     (fPayloadCurr - fPayloadStart < fPayloadSize)) {
	int check = 0x3 & *fPayloadCurr;
	if (channelno % 2 != 0)	{ // odd channel
	  if (check != 0x2 && channelno < 21) {
	    MCMError(kAdcCheckInvalid,
		     "%i for %2i. ADC word in odd channel %i",
		     check, adcwc+1, channelno);
	  }
	}
	else {			// even channel
	  if (check != 0x3 && channelno < 21) {
	    MCMError(kAdcCheckInvalid,
		     "%i for %2i. ADC word in even channel %i",
		     check, adcwc+1, channelno);
	  }
	}

	// filling the actual timebin data
	int tb2 = 0x3ff & *fPayloadCurr >> 22;
	int tb1 = 0x3ff & *fPayloadCurr >> 12;
	int tb0 = 0x3ff & *fPayloadCurr >> 2;
	if (adcwc != 0 || fCurrNtimebins <= 30)
	  fAdcArray->SetDataByAdcCol(row, adccol, currentTimebin++, tb0);
	else
	  tb0 = -1;
	fAdcArray->SetDataByAdcCol(row, adccol, currentTimebin++, tb1);
	fAdcArray->SetDataByAdcCol(row, adccol, currentTimebin++, tb2);

	adcwc++;
	fPayloadCurr++;
      }

      if (adcwc != nADCwords)
	MCMError(kAdcDataAbort);

      // adding index
      if (padcol > 0 && padcol < 144) {
	fSignalIndex->AddIndexRC(row, padcol);
      }

      channelcount++;
    }

    if (fCurrSm > -1 && fCurrSm < 18) {
      fStats.fStatsSector[fCurrSm].fStatsHC[fCurrHC%60].fNChannels += channelcount;
      fStats.fStatsSector[fCurrSm].fNChannels                      += channelcount;
    }
    if (channelcount != channelcountExp)
      MCMError(kAdcChannelsMiss);

    mcmcount++;
    if (fCurrSm > -1 && fCurrSm < 18) {
      fStats.fStatsSector[fCurrSm].fStatsHC[fCurrHC%60].fNMCMs++;
      fStats.fStatsSector[fCurrSm].fNMCMs++;
    }

    if (IsDumping() && DumpingMCM(fCurrHC/2, fCurrRobPos, fCurrMcmPos)) {
      AliInfo(DumpRaw(Form("Event %i: Det %3i ROB %i MCM %2i", fRawReader->GetEventIndex(), fCurrHC/2, fCurrRobPos, fCurrMcmPos),
		      startPosMCM, fPayloadCurr - startPosMCM));
    }

    // continue with next MCM
  }

  // check for missing MCMs (if header suppression is inactive)
  if (((fCurrMajor & 0x1) == 0) && (mcmcount != mcmcountExp)) {
    LinkError(kMissMcmHeaders,
	      "No. of MCM headers %i not as expected: %i",
	      mcmcount, mcmcountExp);
  }

  return (fPayloadCurr - start);
}

Int_t AliTRDrawStream::ReadNonZSData()
{
  // read the non-zs data from one link from the current reading position

  UInt_t *start = fPayloadCurr;

  Int_t mcmcount = 0;
  Int_t mcmcountExp = fCurrStack == 2 ? 48 : 64;
  Int_t channelcount = 0;
  Int_t channelcountExp = 0;
  Int_t timebins;
  Int_t currentTimebin = 0;
  Int_t adcwc = 0;
  Int_t evno = -1;
  Int_t lastmcmpos = -1;
  Int_t lastrobpos = -1;

  if (fCurrNtimebins != fNtimebins) {
    if (fNtimebins > 0)
      LinkError(kNtimebinsChanged,
		"No. of timebins changed from %i to %i", fNtimebins, fCurrNtimebins);
    fNtimebins = fCurrNtimebins;
  }

  timebins = fNtimebins;

  while (*(fPayloadCurr) != fgkDataEndmarker &&
	 fPayloadCurr - fPayloadStart < fPayloadSize - 2) {

    // ----- Checking MCM Header -----
    AliDebug(2, Form("MCM header: 0x%08x", *fPayloadCurr));

    // ----- checking for proper readout order - ROB -----
    if (GetROBReadoutPos(ROB(*fPayloadCurr) / 2) >= lastrobpos) {
      lastrobpos = GetROBReadoutPos(ROB(*fPayloadCurr) / 2);
    }
    else {
      ROBError(kPosUnexp, Form("#%i after #%i in readout order", GetROBReadoutPos(ROB(*fPayloadCurr) / 2), lastrobpos));
    }
    fCurrRobPos = ROB(*fPayloadCurr);

    // ----- checking for proper readout order - MCM -----
    if (GetMCMReadoutPos(MCM(*fPayloadCurr)) >= (lastmcmpos + 1) % 16) {
      lastmcmpos = GetMCMReadoutPos(MCM(*fPayloadCurr));
    }
    else {
      MCMError(kPosUnexp, Form("#%i after #%i in readout order", GetMCMReadoutPos(MCM(*fPayloadCurr)), lastmcmpos));
    }
    fCurrMcmPos = MCM(*fPayloadCurr);

    if (EvNo(*fPayloadCurr) != evno) {
      if (evno == -1)
	evno = EvNo(*fPayloadCurr);
      else {
	MCMError(kPtrgCntMismatch, "%i <-> %i", evno, EvNo(*fPayloadCurr));
      }
    }

    channelcount = 0;
    channelcountExp = 21;
    int channelno = -1;

    Int_t adccoloff = AdcColOffset(*fPayloadCurr);
    Int_t padcoloff = PadColOffset(*fPayloadCurr);
    Int_t row = Row(*fPayloadCurr);

    fPayloadCurr++;

    // ----- reading marked ADC channels -----
    while (channelcount < channelcountExp &&
	   *(fPayloadCurr) != fgkDataEndmarker) {
      if (channelno < 20)
	channelno++;

      currentTimebin = 0;

      adcwc = 0;
      Int_t nADCwords = (timebins + 2) / 3;
      AliDebug(2, Form("Now looking %i words", nADCwords));
      Int_t adccol = adccoloff - channelno;
      Int_t padcol = padcoloff - channelno;
      while ((adcwc < nADCwords) &&
	     (*(fPayloadCurr) != fgkDataEndmarker) &&
	     (fPayloadCurr - fPayloadStart < fPayloadSize)) {
	int check = 0x3 & *fPayloadCurr;
	if (channelno % 2 != 0)	{ // odd channel
	  if (check != 0x2 && channelno < 21) {
	    MCMError(kAdcCheckInvalid,
		     "%i for %2i. ADC word in odd channel %i",
		     check, adcwc+1, channelno);
	  }
	}
	else {			// even channel
	  if (check != 0x3 && channelno < 21) {
	    MCMError(kAdcCheckInvalid,
		     "%i for %2i. ADC word in even channel %i",
		     check, adcwc+1, channelno);
	  }
	}

	// filling the actual timebin data
	int tb2 = 0x3ff & *fPayloadCurr >> 22;
	int tb1 = 0x3ff & *fPayloadCurr >> 12;
	int tb0 = 0x3ff & *fPayloadCurr >> 2;
	if (adcwc != 0 || fCurrNtimebins <= 30)
	  fAdcArray->SetDataByAdcCol(row, adccol, currentTimebin++, tb0);
	else
	  tb0 = -1;
	fAdcArray->SetDataByAdcCol(row, adccol, currentTimebin++, tb1);
	fAdcArray->SetDataByAdcCol(row, adccol, currentTimebin++, tb2);

	adcwc++;
	fPayloadCurr++;
      }

      if (adcwc != nADCwords)
	MCMError(kAdcDataAbort);

      // adding index
      if (padcol > 0 && padcol < 144) {
	fSignalIndex->AddIndexRC(row, padcol);
      }

      channelcount++;
    }

    if (channelcount != channelcountExp)
      MCMError(kAdcChannelsMiss);
    mcmcount++;
    // continue with next MCM
  }

  // check for missing MCMs (if header suppression is inactive)
  if (mcmcount != mcmcountExp) {
    LinkError(kMissMcmHeaders,
	      "%i not as expected: %i", mcmcount, mcmcountExp);
  }

  return (fPayloadCurr - start);
}

Int_t AliTRDrawStream::SeekNextLink()
{
  // proceed in raw data stream till the next link

  UInt_t *start = fPayloadCurr;

  // read until data endmarkers
  while (fPayloadCurr - fPayloadStart < fPayloadSize &&
	 *fPayloadCurr != fgkDataEndmarker)
    fPayloadCurr++;

  // read all data endmarkers
  while (fPayloadCurr - fPayloadStart < fPayloadSize &&
	 *fPayloadCurr == fgkDataEndmarker)
    fPayloadCurr++;

  return (fPayloadCurr - start);
}

Bool_t AliTRDrawStream::ConnectTracklets(TTree *trklTree)
{
  // connect the tracklet tree used to store the tracklet output

  fTrackletTree = trklTree;
  if (!fTrackletTree)
    return kTRUE;

  if (!fTrackletTree->GetBranch("hc"))
    fTrackletTree->Branch("hc", &fCurrHC, "hc/I");
  else
    fTrackletTree->SetBranchAddress("hc", &fCurrHC);

  if (!fTrackletTree->GetBranch("trkl"))
    fTrackletTree->Branch("trkl", &fTrackletArray);
  else
    fTrackletTree->SetBranchAddress("trkl", &fTrackletArray);

  return kTRUE;
}


void AliTRDrawStream::EquipmentError(ErrorCode_t err, const char *const msg, ...)
{
  // register error according to error code on equipment level
  // and return the corresponding error message

  fLastError.fSector = fCurrEquipmentId - kDDLOffset;
  fLastError.fStack  = -1;
  fLastError.fLink   = -1;
  fLastError.fRob    = -1;
  fLastError.fMcm    = -1;
  fLastError.fError  = err;
  (this->*fStoreError)();

  va_list ap;
  if (fgErrorDebugLevel[err] > 10)
    AliDebug(fgErrorDebugLevel[err],
	     Form("Event %6i: Eq. %2d - %s : %s",
		  fRawReader->GetEventIndex(), fCurrEquipmentId, fgkErrorMessages[err],
		  (va_start(ap, msg), vsprintf(fErrorBuffer, msg, ap), va_end(ap), fErrorBuffer) ));
  else
    AliError(Form("Event %6i: Eq. %2d - %s : %s",
		  fRawReader->GetEventIndex(), fCurrEquipmentId, fgkErrorMessages[err],
		  (va_start(ap, msg), vsprintf(fErrorBuffer, msg, ap), va_end(ap), fErrorBuffer) ));
  fErrorFlags |= fgErrorBehav[err];
}


void AliTRDrawStream::StackError(ErrorCode_t err, const char *const msg, ...)
{
  // register error according to error code on stack level
  // and return the corresponding error message

  fLastError.fSector = fCurrEquipmentId - kDDLOffset;
  fLastError.fStack  = fCurrSlot;
  fLastError.fLink   = -1;
  fLastError.fRob    = -1;
  fLastError.fMcm    = -1;
  fLastError.fError  = err;
  (this->*fStoreError)();

  va_list ap;
  if (fgErrorDebugLevel[err] > 0)
    AliDebug(fgErrorDebugLevel[err],
	     Form("Event %6i: Eq. %2d S %i - %s : %s",
		  fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fgkErrorMessages[err],
		  (va_start(ap, msg), vsprintf(fErrorBuffer, msg, ap), va_end(ap), fErrorBuffer) ));
  else
    AliError(Form("Event %6i: Eq. %2d S %i - %s : %s",
		  fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fgkErrorMessages[err],
		  (va_start(ap, msg), vsprintf(fErrorBuffer, msg, ap), va_end(ap), fErrorBuffer) ));
  fErrorFlags |= fgErrorBehav[err];
}


void AliTRDrawStream::LinkError(ErrorCode_t err, const char *const msg, ...)
{
  // register error according to error code on link level
  // and return the corresponding error message

  fLastError.fSector = fCurrEquipmentId - kDDLOffset;
  fLastError.fStack  = fCurrSlot;
  fLastError.fLink   = fCurrLink;
  fLastError.fRob    = -1;
  fLastError.fMcm    = -1;
  fLastError.fError  = err;
  (this->*fStoreError)();

  va_list ap;
  if (fgErrorDebugLevel[err] > 0)
    AliDebug(fgErrorDebugLevel[err],
	     Form("Event %6i: Eq. %2d S %i l %2i - %s : %s",
		  fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fCurrLink, fgkErrorMessages[err],
		  (va_start(ap, msg), vsprintf(fErrorBuffer, msg, ap), va_end(ap), fErrorBuffer) ));
  else
    AliError(Form("Event %6i: Eq. %2d S %i l %2i - %s : %s",
		  fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fCurrLink, fgkErrorMessages[err],
		  (va_start(ap, msg), vsprintf(fErrorBuffer, msg, ap), va_end(ap), fErrorBuffer) ));
  fErrorFlags |= fgErrorBehav[err];
}


void AliTRDrawStream::ROBError(ErrorCode_t err, const char *const msg, ...)
{
  // register error according to error code on ROB level
  // and return the corresponding error message

  fLastError.fSector = fCurrEquipmentId - kDDLOffset;
  fLastError.fStack  = fCurrSlot;
  fLastError.fLink   = fCurrLink;
  fLastError.fRob    = fCurrRobPos;
  fLastError.fMcm    = -1;
  fLastError.fError  = err;
  (this->*fStoreError)();

  va_list ap;
  if (fgErrorDebugLevel[err] > 0)
    AliDebug(fgErrorDebugLevel[err],
	     Form("Event %6i: Eq. %2d S %i l %2i ROB %i - %s : %s",
		  fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fCurrLink, fCurrRobPos, fgkErrorMessages[err],
		  (va_start(ap, msg), vsprintf(fErrorBuffer, msg, ap), va_end(ap), fErrorBuffer) ));
  else
    AliError(Form("Event %6i: Eq. %2d S %i l %2i ROB %i - %s : %s",
		  fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fCurrLink, fCurrRobPos, fgkErrorMessages[err],
		  (va_start(ap, msg), vsprintf(fErrorBuffer, msg, ap), va_end(ap), fErrorBuffer) ));
  fErrorFlags |= fgErrorBehav[err];
}


void AliTRDrawStream::MCMError(ErrorCode_t err, const char *const msg, ...)
{
  // register error according to error code on MCM level
  // and return the corresponding error message

  fLastError.fSector = fCurrEquipmentId - kDDLOffset;
  fLastError.fStack  = fCurrSlot;
  fLastError.fLink   = fCurrLink;
  fLastError.fRob    = fCurrRobPos;
  fLastError.fMcm    = fCurrMcmPos;
  fLastError.fError  = err;
  (this->*fStoreError)();

  va_list ap;
  if (fgErrorDebugLevel[err] > 0)
    AliDebug(fgErrorDebugLevel[err],
	     Form("Event %6i: Eq. %2d S %i l %2i ROB %i MCM %2i - %s : %s",
		  fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fCurrLink, fCurrRobPos, fCurrMcmPos, fgkErrorMessages[err],
		  (va_start(ap, msg), vsprintf(fErrorBuffer, msg, ap), va_end(ap), fErrorBuffer) ));
  else
    AliError(Form("Event %6i: Eq. %2d S %i l %2i ROB %i MCM %2i - %s : %s",
		  fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fCurrLink, fCurrRobPos, fCurrMcmPos, fgkErrorMessages[err],
		  (va_start(ap, msg), vsprintf(fErrorBuffer, msg, ap), va_end(ap), fErrorBuffer) ));
  fErrorFlags |= fgErrorBehav[err];
}

const char* AliTRDrawStream::GetErrorMessage(ErrorCode_t errCode)
{
  // return the error message for the given error code

  if (errCode > 0 && errCode < kLastErrorCode)
    return fgkErrorMessages[errCode];
  else
    return "";
}

void AliTRDrawStream::AliTRDrawStats::ClearStats()
{
  // clear statistics (includes clearing sector-wise statistics)

  fBytesRead = 0;
  for (Int_t iSector = 0; iSector < 18; iSector++) {
    fStatsSector[iSector].ClearStats();
  }

}

void AliTRDrawStream::AliTRDrawStats::AliTRDrawStatsSector::ClearStats()
{
  // clear statistics (includes clearing HC-wise statistics)

  fBytes = 0;
  fBytesRead = 0;
  fNTracklets = 0;
  fNMCMs = 0;
  fNChannels = 0;

  for (Int_t iHC = 0; iHC < 60; iHC++) {
    fStatsHC[iHC].ClearStats();
  }
}

void AliTRDrawStream::AliTRDrawStats::AliTRDrawStatsSector::AliTRDrawStatsHC::ClearStats()
{
  // clear statistics

  fBytes = 0;
  fBytesRead = 0;
  fNTracklets = 0;
  fNMCMs = 0;
  fNChannels = 0;
}

void AliTRDrawStream::SetDumpMCM(Int_t det, Int_t rob, Int_t mcm, Bool_t dump)
{
  // mark MCM for dumping of raw data

  if (dump) {
    fDumpMCM[fNDumpMCMs++] = (det << 7) | (rob << 4) | mcm;
  }
  else {
    Int_t iMCM;
    for (iMCM = 0; iMCM < fNDumpMCMs; iMCM++) {
      if (fDumpMCM[iMCM] == ((det << 7) | (rob << 4) | mcm)) {
	fNDumpMCMs--;
	break;
      }
    }
    for ( ; iMCM < fNDumpMCMs; iMCM++) {
      fDumpMCM[iMCM] = fDumpMCM[iMCM+1];
    }
  }
}

Bool_t AliTRDrawStream::DumpingMCM(Int_t det, Int_t rob, Int_t mcm)  const
{
  // check if MCM data should be dumped

  for (Int_t iMCM = 0; iMCM < fNDumpMCMs; iMCM++) {
    if (fDumpMCM[iMCM] == ((det << 7) | (rob << 4) | mcm)) {
      return kTRUE;
    }
  }
  return kFALSE;
}

TString AliTRDrawStream::DumpRaw(TString title, UInt_t *start, Int_t length, UInt_t endmarker)
{
  // dump raw data

  title += "\n";
  for (Int_t pos = 0; pos < length; pos += 4) {
    if ((start[pos+0] != endmarker) && pos+0 < length)
      if ((start[pos+1] != endmarker && pos+1 < length))
	if ((start[pos+2] != endmarker && pos+2 < length))
	  if ((start[pos+3] != endmarker && pos+3 < length))
	    title += Form("   0x%08x 0x%08x 0x%08x 0x%08x\n",
			  start[pos+0], start[pos+1], start[pos+2], start[pos+3]);
	  else {
	    title += Form("   0x%08x 0x%08x 0x%08x 0x%08x\n",
			  start[pos+0], start[pos+1], start[pos+2], start[pos+3]);
	    return title;
	  }
	else {
	  title += Form("   0x%08x 0x%08x 0x%08x\n",
			start[pos+0], start[pos+1], start[pos+2]);
	  return title;
	}
      else {
	title += Form("   0x%08x 0x%08x\n",
		      start[pos+0], start[pos+1]);
	return title;
      }
    else {
      title += Form("   0x%08x\n",
		    start[pos+0]);
      return title;
    }
  }
  return title;
}

TString AliTRDrawStream::DumpMcmHeader(TString title, UInt_t word)
{
  title += Form("0x%08x -> ROB: %i, MCM: %2i",
		word, ROB(word), MCM(word));
  return title;
}

TString AliTRDrawStream::DumpAdcMask(TString title, UInt_t word)
{
  title += Form("0x%08x -> #ch : %2i, 0x%06x (%2i ch)",
		word, GetNActiveChannels(word), GetActiveChannels(word), GetNActiveChannelsFromMask(word));
  return title;
}

AliTRDrawStream::AliTRDrawStreamError::AliTRDrawStreamError(Int_t error, Int_t sector, Int_t stack, Int_t link, Int_t rob, Int_t mcm) :
  fError(error),
  fSector(sector),
  fStack(stack),
  fLink(link),
  fRob(rob),
  fMcm(mcm)
{
  // ctor

}
