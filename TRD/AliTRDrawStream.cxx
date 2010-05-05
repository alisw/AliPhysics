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
//  and translation into ADC values                                       //
//                                                                        //
//  Author: J. Klein (jochen.klein@cern.ch)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

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
#include "AliTreeLoader.h"

#include "AliTRDrawStream.h"

// temporary
#include "AliRunLoader.h"

ClassImp(AliTRDrawStream)

// some static information 
const Int_t AliTRDrawStream::fgkMcmOrder[] = {12, 13, 14, 15, 
					      8, 9, 10, 11, 
					      4, 5, 6, 7, 
					      0, 1, 2, 3};
const Int_t  AliTRDrawStream::fgkRobOrder [] = {0, 1, 2, 3};
const Int_t  AliTRDrawStream::fgkNlinks = 12;
const Int_t  AliTRDrawStream::fgkNstacks = 5;
const UInt_t AliTRDrawStream::fgkDataEndmarker     = 0x00000000;
const UInt_t AliTRDrawStream::fgkTrackletEndmarker = 0x10001000;

const char* AliTRDrawStream::fgErrorMessages[] = {
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
  "Missing MCM headers"
};

AliTRDrawStream::AliTRDrawStream(AliRawReader *rawReader) :
  fStats(), 
  fRawReader(rawReader),
  fDigitsManager(0x0),
  fDigitsParam(0x0),
  fErrors(0x0),
  fLastError(),
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
  fCurrSmuIndexHeaderSize(0),
  fCurrSmuIndexHeaderVersion(0),
  fCurrTrackEnable(0),
  fCurrTrackletEnable(0),
  fCurrStackMask(0),
  fCurrStackIndexWord(0x0),
  fCurrStackHeaderSize(0x0),
  fCurrStackHeaderVersion(0x0),
  fCurrLinkMask(0x0),
  fCurrCleanCheckout(0x0),
  fCurrBoardId(0x0),
  fCurrHwRev(0x0),
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
  fTrackletTree(0x0)
{
  // default constructor

  fCurrStackIndexWord     = new UInt_t[fgkNstacks];	 
  fCurrStackHeaderSize    = new UInt_t[fgkNstacks];	 
  fCurrStackHeaderVersion = new UInt_t[fgkNstacks];
  fCurrLinkMask           = new UInt_t[fgkNstacks];		 
  fCurrCleanCheckout      = new UInt_t[fgkNstacks];	 
  fCurrBoardId            = new UInt_t[fgkNstacks];		 
  fCurrHwRev              = new UInt_t[fgkNstacks];             
  fCurrLinkMonitorFlags   = new UInt_t[fgkNstacks * fgkNlinks];
  fCurrLinkDataTypeFlags  = new UInt_t[fgkNstacks * fgkNlinks];
  fCurrLinkDebugFlags     = new UInt_t[fgkNstacks * fgkNlinks];

  // preparing TClonesArray
  fTrackletArray = new TClonesArray("AliTRDtrackletWord", 256);

  // setting up the error tree
  fErrors = new TTree("errorStats", "Error statistics");
  fErrors->SetDirectory(0x0);
  fErrors->Branch("error", &fLastError, "sector/I:stack:link:error:rob:mcm");
  fErrors->SetCircular(1000);
}

AliTRDrawStream::~AliTRDrawStream()
{
  // destructor

  delete fErrors;

  delete [] fCurrStackIndexWord;
  delete [] fCurrStackHeaderSize;
  delete [] fCurrStackHeaderVersion;
  delete [] fCurrLinkMask;
  delete [] fCurrCleanCheckout;
  delete [] fCurrBoardId;
  delete [] fCurrHwRev;
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
  // data starts with GTU payload, i.e. SMU index word
  UChar_t *buffer = 0x0;

  while (fRawReader->ReadNextData(buffer)) {

    fCurrEquipmentId = fRawReader->GetEquipmentId();
    AliDebug(2, Form("equipment: %i", fCurrEquipmentId));

    if (fCurrEquipmentId < 1024 || fCurrEquipmentId > 1041) {
      AliError(EquipmentError(kNonTrdEq, "Skipping"));
      continue;
    }

    // setting the pointer to data and current reading position
    fPayloadCurr = fPayloadStart = (UInt_t*) (buffer);
    fPayloadSize = fRawReader->GetDataSize() / sizeof(UInt_t);
    fStats.fStatsSector[fCurrEquipmentId - 1024].fBytes = fRawReader->GetDataSize();
    AliDebug(2, Form("Read buffer of size: %i", fRawReader->GetDataSize()));

    // read SMU index header
    if (ReadSmHeader() < 0) {
      AliError(Form("Reading SMU header failed, skipping this DDL %i", fCurrEquipmentId));
      continue;
    }

    // read stack index header
    for (Int_t iStack = 0; iStack < 5; iStack++) {
      if ((fCurrStackMask & (1 << iStack)) != 0) 
	ReadStackIndexHeader(iStack);
    }

    for (Int_t iStack = 0; iStack < 5; iStack++) {
      fCurrSlot = iStack;
      if ((fCurrStackMask & (1 << fCurrSlot)) == 0) 
	continue;

      AliDebug(2, Form("Stack %i, Link mask: 0x%02x", fCurrSlot, fCurrLinkMask[fCurrSlot]));
      for (Int_t iLink = 0; iLink < 12; iLink++) {
	fCurrLink = iLink;
	fCurrHC   = fCurrSm * 60 + fCurrSlot * 12 + iLink;
	if ((fCurrLinkMask[fCurrSlot] & (1 << fCurrLink)) == 0)
	  continue;
	
	if (fCurrLinkMonitorFlags[fCurrSlot*fgkNlinks + fCurrLink] != 0)
	  LinkError(kLinkMonitor);
	AliDebug(2, Form("Payload for S%il%i", fCurrSlot, fCurrLink));
	for (Int_t i = 0; i < 10; i++) //???
	  AliDebug(5, Form("%3i: 0x%08x 0x%08x 0x%08x 0x%08x", i*4,
		 fPayloadCurr[4*i], fPayloadCurr[4*i+1], 
		 fPayloadCurr[4*i+2], fPayloadCurr[4*i+3]));

	// read the data from one HC
	ReadLinkData();

	// read all data endmarkers
	while (fPayloadCurr - fPayloadStart < fPayloadSize &&
	       *fPayloadCurr == fgkDataEndmarker)
	  fPayloadCurr++;

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
    
    if (fCurrEquipmentId < 1024 || fCurrEquipmentId > 1041) {
      AliError(EquipmentError(kNonTrdEq, "Skipping"));
      continue;
    }

    // setting the pointer to data and current reading position
    fPayloadCurr = fPayloadStart = (UInt_t*) (buffer);
    fPayloadSize = fRawReader->GetDataSize() / sizeof(UInt_t);
    AliDebug(2, Form("Read buffer of size: %i", fRawReader->GetDataSize()));

    // read SMU index header
    if (ReadSmHeader() < 0) {
      AliError(Form("Reading SMU header failed, skipping this DDL %i", fCurrEquipmentId));
      continue;
    }
    
    // read stack index header
    for (Int_t iStack = 0; iStack < 5; iStack++) {
      if ((fCurrStackMask & (1 << iStack)) != 0) {
	ReadStackIndexHeader(iStack);
      }
    }
    return kTRUE;
  }

  return kFALSE;
}


Int_t AliTRDrawStream::NextChamber(AliTRDdigitsManager *digMgr, UInt_t ** /* trackletContainer */, UShort_t ** /* errorContainer */)
{
  // read the data for the next chamber
  // in case you only want to read the data of a single chamber
  // to read all data ReadEvent(...) is recommended

  fDigitsManager = digMgr; 
  fDigitsParam   = 0x0;

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

  if (fCurrSlot < 0 || fCurrSlot >= 5) {
    if (!NextDDL()) {
      fCurrSlot = -1;
      return -1;
    }
  }

  AliDebug(2, Form("Stack %i, Link %i, mask: 0x%02x", fCurrSlot, fCurrLink, fCurrLinkMask[fCurrSlot]));
  fCurrHC   = (fCurrEquipmentId - 1024) * 60 + fCurrSlot * 12 + fCurrLink;

  if (fCurrLinkMonitorFlags[fCurrSlot*fgkNlinks + fCurrLink] != 0)
    LinkError(kLinkMonitor);

  // read the data from one HC
  ReadLinkData();
  
  // read all data endmarkers
  while (fPayloadCurr - fPayloadStart < fPayloadSize && 
	 *fPayloadCurr == fgkDataEndmarker) 
    fPayloadCurr++;

  if (fCurrLink % 2 == 0) {
    // if we just read the A-side HC then also check the B-side
    fCurrLink++;
    if (fCurrLinkMask[fCurrSlot] & (1 << fCurrLink)) {
      ReadLinkData();
      while (fPayloadCurr - fPayloadStart < fPayloadSize && 
	     *fPayloadCurr == fgkDataEndmarker) 
	fPayloadCurr++;
    }
  }

  //??? to check 
  do {
    fCurrLink++; 
    if (fCurrLink > 11) {
      fCurrLink = 0;
      fCurrSlot++;
    }
  } while ((fCurrSlot < 5) && 
	   (((fCurrStackMask & (1 << fCurrSlot)) == 0) || 
	    ((fCurrLinkMask[fCurrSlot] & (1 << fCurrLink))) == 0));

  return (fCurrSm * 30 + fCurrStack * 6 + fCurrLayer);
}


Int_t AliTRDrawStream::ReadSmHeader()
{
  // read the SMU index header at the current reading position 
  // and store the information in the corresponding variables

  if (fPayloadCurr - fPayloadStart >= fPayloadSize - 1) {
    AliError(EquipmentError(kUnknown, "SM Header incomplete"));
    return -1;
  }

  fCurrSmuIndexHeaderSize     = ((*fPayloadCurr) >> 16) & 0xffff;
  fCurrSmuIndexHeaderVersion  = ((*fPayloadCurr) >> 12) &    0xf;
  fCurrTrackEnable            = ((*fPayloadCurr) >>  6) &    0x1;
  fCurrTrackletEnable         = ((*fPayloadCurr) >>  5) &    0x1;
  fCurrStackMask              = ((*fPayloadCurr)      ) &   0x1f;

  AliDebug(5, Form("SMU header: size: %i, version: %i, track enable: %i, tracklet enable: %i, stack mask: %2x",
		   fCurrSmuIndexHeaderSize, 
		   fCurrSmuIndexHeaderVersion, 
		   fCurrTrackEnable, 
		   fCurrTrackletEnable,
		   fCurrStackMask));
  
  fPayloadCurr += fCurrSmuIndexHeaderSize + 1;

  return fCurrSmuIndexHeaderSize + 1;
}

Int_t AliTRDrawStream::ReadStackIndexHeader(Int_t stack)
{
  // read the stack index header 
  // and store the information in the corresponding variables

  fCurrStackIndexWord[stack]     = *fPayloadCurr;
  fCurrStackHeaderSize[stack]    = (((*fPayloadCurr) >> 16) & 0xffff) + 1;
  fCurrStackHeaderVersion[stack] = ((*fPayloadCurr) >> 12) & 0xf;
  fCurrLinkMask[stack]           = (*fPayloadCurr) & 0xfff;

  if (fPayloadCurr - fPayloadStart >= fPayloadSize - (Int_t) fCurrStackHeaderSize[stack]) {
    AliError(StackError(kStackHeaderInvalid, "Stack index header aborted"));
    return -1;
  }

  switch (fCurrStackHeaderVersion[stack]) {
  case 0xa: 
    if (fCurrStackHeaderSize[stack] < 8) {
      AliError(StackError(kStackHeaderInvalid, "Stack header smaller than expected!"));
      return -1;
    }
    
    fCurrCleanCheckout[stack] = fPayloadCurr[1] & 0x1;
    fCurrBoardId[stack]       = (fPayloadCurr[1] >> 8) & 0xff;
    fCurrHwRev[stack]         = (fPayloadCurr[1] >> 16) & 0xffff;
    
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
    AliError(StackError(kStackHeaderInvalid, Form("Invalid Stack Index Header version %x", 
						  fCurrStackHeaderVersion[stack])));
  }
  
  fPayloadCurr += fCurrStackHeaderSize[stack];

  return fCurrStackHeaderSize[stack];
}

Int_t AliTRDrawStream::ReadLinkData()
{
  // read the data in one link (one HC) until the data endmarker is reached
  // returns the number of words read!

  Int_t count = 0;
  UInt_t* startPosLink = fPayloadCurr;

//  printf("----- HC: %i -----\n", fCurrHC);
//  for (Int_t i = 0; i < 3; i++) {
//    printf("0x%08x 0x%08x 0x%08x 0x%08x\n", 
//	   fPayloadCurr[i*4+0], fPayloadCurr[i*4+1], fPayloadCurr[i*4+2], fPayloadCurr[i*4+3]);
//  }

  count += ReadTracklets();

  count += ReadHcHeader();

  Int_t det = fCurrSm * 30 + fCurrStack * 6 + fCurrLayer;

  if (det > -1 && det < 540) {
    
    if ((fAdcArray = fDigitsManager->GetDigits(det))) {
      //fAdcArray->Expand();
      if (fAdcArray->GetNtime() != fCurrNtimebins)
	fAdcArray->Allocate(16, 144, fCurrNtimebins);
    }
    else {
      AliError(LinkError(kNoDigits));
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
    else if (fCurrMajor & 0x20) {
      AliDebug(1, "This is a zs event");
      count += ReadZSData();
    }
    else {
      AliDebug(1, "This is a nozs event");
      count += ReadNonZSData();
    }
  }
  else {
    AliError(LinkError(kInvalidDetector, Form("%i", det)));
    while (fPayloadCurr - fPayloadStart < fPayloadSize &&
	   *fPayloadCurr != fgkDataEndmarker)
      fPayloadCurr++;
  }
  
  fStats.fStatsSector[fCurrSm].fStatsHC[fCurrHC%60].fBytes     += (fPayloadCurr - startPosLink) * sizeof(UInt_t);
  fStats.fStatsSector[fCurrSm].fStatsHC[fCurrHC%60].fBytesRead += count * sizeof(UInt_t);
  fStats.fStatsSector[fCurrSm].fBytesRead                      += count * sizeof(UInt_t);
  fStats.fBytesRead                                            += count * sizeof(UInt_t);

  return count;
}

Int_t AliTRDrawStream::ReadTracklets()
{
  // read the tracklets from one HC

  fTrackletArray->Clear();

  UInt_t *start = fPayloadCurr;
  while (*(fPayloadCurr) != fgkTrackletEndmarker && 
	 fPayloadCurr - fPayloadStart < fPayloadSize) {

    new ((*fTrackletArray)[fTrackletArray->GetEntriesFast()]) AliTRDtrackletWord(*(fPayloadCurr));

    fPayloadCurr++;
  }

  if (fTrackletArray->GetEntriesFast() > 0) {
    AliDebug(1, Form("Found %i tracklets in %i %i %i (ev. %i)", fTrackletArray->GetEntriesFast(), 
		     fCurrSm, fCurrSlot, fCurrLink, fRawReader->GetEventIndex()));
    fStats.fStatsSector[fCurrSm].fStatsHC[fCurrHC%60].fNTracklets += fTrackletArray->GetEntriesFast();
    fStats.fStatsSector[fCurrSm].fNTracklets                      += fTrackletArray->GetEntriesFast();
    if (fTrackletTree)
      fTrackletTree->Fill();
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

  if (fCurrSm != (((Int_t) fCurrEquipmentId) - 1024) || 
      fCurrStack != fCurrSlot || 
      fCurrLayer != fCurrLink / 2 || 
      fCurrSide != fCurrLink % 2) {
    AliError(LinkError(kHCmismatch,
		       Form("HC: %i, %i, %i, %i\n 0x%08x 0x%08x 0x%08x 0x%08x", 
			    fCurrSm, fCurrStack, fCurrLayer, fCurrSide,
			    fPayloadCurr[0], fPayloadCurr[1], fPayloadCurr[2], fPayloadCurr[3]) ));
  }
  if (fCurrCheck != 0x1) {
    AliError(LinkError(kHCcheckFailed));
  }
  
  if (fCurrAddHcWords > 0) {
    fCurrNtimebins = (fPayloadCurr[1] >> 26) & 0x3f;
    fCurrBC = (fPayloadCurr[1] >> 10) & 0xffff;
    fCurrPtrgCnt = (fPayloadCurr[1] >> 6) & 0xf;
    fCurrPtrgPhase = (fPayloadCurr[1] >> 2) & 0xf;
  }
  
  fPayloadCurr += 1 + fCurrAddHcWords;
  
  return (fPayloadCurr - start) / sizeof(UInt_t);
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
    AliDebug(2, Form("MCM header: 0x%08x", *fPayloadCurr));
    mcmcount++;
    
    // ----- checking for proper readout order - ROB -----
    if (GetROBReadoutPos(ROB(*fPayloadCurr) / 2) >= lastrobpos) {
      lastrobpos = GetROBReadoutPos(ROB(*fPayloadCurr) / 2);
    }
    else {
      AliError(ROBError(kPosUnexp));
    }
    fCurrRobPos = ROB(*fPayloadCurr);
    
    // ----- checking for proper readout order - MCM -----
    if (GetMCMReadoutPos(MCM(*fPayloadCurr)) >= (lastmcmpos + 1) % 16) {
      lastmcmpos = GetMCMReadoutPos(MCM(*fPayloadCurr));
    }
    else {
      AliError(MCMError(kPosUnexp));
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
	  AliError(LinkError(kTPmodeInvalid, "Just reading"));
	}

	diff = *fPayloadCurr ^ expword;
	if (diff != 0) {
	  AliError(MCMError(kTPmismatch,
			    Form("Seen 0x%08x, expected 0x%08x, diff: 0x%08x (0x%02x)", 
				 *fPayloadCurr, expword, diff, 0xff & (diff | diff >> 8 | diff >> 16 | diff >> 24)) ));
	}
	fPayloadCurr++;
	count++;
	wordcount++;
      }
      channelcount++;
    }
    // continue with next MCM
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
      AliError(LinkError(kNtimebinsChanged,
			 Form("No. of timebins changed from %i to %i", fNtimebins, fCurrNtimebins) ));
    fNtimebins = fCurrNtimebins;
  }
  
  timebins = fNtimebins;
  
  while (*(fPayloadCurr) != fgkDataEndmarker && 
	 fPayloadCurr - fPayloadStart < fPayloadSize) {
    
    // ----- Checking MCM Header -----
    AliDebug(2, Form("MCM header: 0x%08x", *fPayloadCurr));
    UInt_t *startPosMCM = fPayloadCurr;
    
    // ----- checking for proper readout order - ROB -----
    if (GetROBReadoutPos(ROB(*fPayloadCurr) / 2) >= lastrobpos) {
      lastrobpos = GetROBReadoutPos(ROB(*fPayloadCurr) / 2);
    }
    else {
      AliError(ROBError(kPosUnexp));
    }
    fCurrRobPos = ROB(*fPayloadCurr);
    
    // ----- checking for proper readout order - MCM -----
    if (GetMCMReadoutPos(MCM(*fPayloadCurr)) > lastmcmpos) {
      lastmcmpos = GetMCMReadoutPos(lastmcmpos);
    }
    else {
      AliError(MCMError(kPosUnexp));
    }
    fCurrMcmPos = MCM(*fPayloadCurr);
    
    if (EvNo(*fPayloadCurr) != evno) {
      if (evno == -1)
	evno = EvNo(*fPayloadCurr);
      else {
	AliError(MCMError(kPtrgCntMismatch,
			  Form("%i <-> %i", evno, EvNo(*fPayloadCurr)) ));
      }
    }
    Int_t adccoloff = AdcColOffset(*fPayloadCurr);
    Int_t padcoloff = PadColOffset(*fPayloadCurr);
    Int_t row = Row(*fPayloadCurr);
    fPayloadCurr++;
    
    // ----- Reading ADC channels -----
    AliDebug(2, Form("ADC mask: 0x%08x", *fPayloadCurr));
    
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
	AliError(MCMError(kAdcMaskInconsistent,
			  Form("Possible MCM-H: 0x%08x, possible ADC-mask: 0x%08x", 
			       *(fPayloadCurr + 10 * channelcountExp), 
			       *(fPayloadCurr + 10 * channelcountExp + 1) ) ));
	if (!CouldBeMCMhdr( *(fPayloadCurr + 10 * channelcountExp)) && !CouldBeADCmask( *(fPayloadCurr + 10 * channelcountExp + 1))) 
	  channelcountExp++;
	else {
	  break;
	}
      }
      AliError(MCMError(kAdcMaskInconsistent,
			Form("Inconsistency in no. of active channels: Counter: %i, Mask: %i, chosen: %i!", 
			     GetNActiveChannels(fPayloadCurr[-1]), GetNActiveChannelsFromMask(fPayloadCurr[-1]), channelcountExp) ));
    }
    AliDebug(2, Form("expecting %i active channels, timebins: %i", channelcountExp, fCurrNtimebins));
    
    // ----- reading marked ADC channels -----
    while (channelcount < channelcountExp && *(fPayloadCurr) != fgkDataEndmarker) {
      if (channelno < 21)
	channelno++;
      while (channelno < 21 && (channelmask & 1 << channelno) == 0)
	channelno++;
      
      if (fCurrNtimebins > 30) {
	currentTimebin = ((*fPayloadCurr >> 2) & 0x3f);
	timebins = ((*fPayloadCurr >> 8) & 0xf) * 3;
      } 
      else {
	currentTimebin = 0;
      }
      
      adcwc = 0;
      AliDebug(2, Form("Now looking %i words", timebins / 3));
      Int_t adccol = adccoloff - channelno;
      Int_t padcol = padcoloff - channelno;
//      if (adccol < 3 || adccol > 165) 
//	AliInfo(Form("writing channel %i of det %3i %i:%2i to adcrow/-col: %i/%i padcol: %i", 
//		     channelno, fCurrHC/2, fCurrRobPos, fCurrMcmPos, row, adccol, padcol));

      while (adcwc < timebins / 3 && 
	     *(fPayloadCurr) != fgkDataEndmarker && 
	     fPayloadCurr - fPayloadStart < fPayloadSize) {
	int check = 0x3 & *fPayloadCurr;
	if (channelno % 2 != 0)	{ // odd channel
	  if (check != 0x2 && channelno < 21) {
	    AliError(MCMError(kAdcCheckInvalid,
			      Form("%i for %2i. ADC word in odd channel %i", 
				   check, adcwc+1, channelno) ));
	  }
	}
	else {			// even channel
	  if (check != 0x3 && channelno < 21) {
	    AliError(MCMError(kAdcCheckInvalid,
			      Form("%i for %2i. ADC word in even channel %i", 
				   check, adcwc+1, channelno) ));
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
      
      if (adcwc != timebins / 3) 
	AliError(MCMError(kAdcDataAbort));
      
      // adding index 
      if (padcol > 0 && padcol < 144) {
	fSignalIndex->AddIndexRC(row, padcol);
      }
      
      channelcount++;
    }

    fStats.fStatsSector[fCurrSm].fStatsHC[fCurrHC%60].fNChannels += channelcount;
    fStats.fStatsSector[fCurrSm].fNChannels                      += channelcount;
    if (channelcount != channelcountExp)
      AliError(MCMError(kAdcChannelsMiss));
    
    mcmcount++;
    fStats.fStatsSector[fCurrSm].fStatsHC[fCurrHC%60].fNMCMs++;
    fStats.fStatsSector[fCurrSm].fNMCMs++;

    if (IsDumping() && DumpingMCM(fCurrHC/2, fCurrRobPos, fCurrMcmPos)) {
      DumpRaw(Form("Det %3i ROB %i MCM %2i", fCurrHC/2, fCurrRobPos, fCurrMcmPos),
	      startPosMCM, fPayloadCurr - startPosMCM);
    }

    // continue with next MCM
  }

  // check for missing MCMs (if header suppression is inactive)
  if (((fCurrMajor & 0x1) == 0) && (mcmcount != mcmcountExp)) {
    AliError(LinkError(kMissMcmHeaders,
		       Form("No. of MCM headers %i not as expected: %i", 
			    mcmcount, mcmcountExp) ));
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
      AliError(LinkError(kNtimebinsChanged,
			 Form("No. of timebins changed from %i to %i", fNtimebins, fCurrNtimebins) ));
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
      AliError(ROBError(kPosUnexp));
    }
    fCurrRobPos = ROB(*fPayloadCurr);
    
    // ----- checking for proper readout order - MCM -----
    if (GetMCMReadoutPos(MCM(*fPayloadCurr)) >= (lastmcmpos + 1) % 16) {
      lastmcmpos = GetMCMReadoutPos(*fPayloadCurr);
    }
    else {
      AliError(MCMError(kPosUnexp));
    }
    fCurrMcmPos = MCM(*fPayloadCurr);
    
    if (EvNo(*fPayloadCurr) != evno) {
      if (evno == -1)
	evno = EvNo(*fPayloadCurr);
      else {
	AliError(MCMError(kPtrgCntMismatch,
			  Form("%i <-> %i", evno, EvNo(*fPayloadCurr)) ));
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
      if (channelno < 21)
	channelno++;
      
      currentTimebin = 0;
      
      adcwc = 0;
      AliDebug(2, Form("Now looking %i words", timebins / 3));
      Int_t adccol = adccoloff - channelno;
      Int_t padcol = padcoloff - channelno;
      while (adcwc < timebins / 3 && 
	     *(fPayloadCurr) != fgkDataEndmarker && 
	     fPayloadCurr - fPayloadStart < fPayloadSize) {
	int check = 0x3 & *fPayloadCurr;
	if (channelno % 2 != 0)	{ // odd channel
	  if (check != 0x2 && channelno < 21) {
	    AliError(MCMError(kAdcCheckInvalid,
			      Form("%i for %2i. ADC word in odd channel %i", 
				   check, adcwc+1, channelno) ));
	  }
	}
	else {			// even channel
	  if (check != 0x3 && channelno < 21) {
	    AliError(MCMError(kAdcCheckInvalid,
			      Form("%i for %2i. ADC word in even channel %i", 
				   check, adcwc+1, channelno) ));
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

      if (adcwc != timebins / 3) 
	AliError(MCMError(kAdcDataAbort));
      
      // adding index 
      if (padcol > 0 && padcol < 144) {
	fSignalIndex->AddIndexRC(row, padcol);
      }

      channelcount++;
    }

    if (channelcount != channelcountExp)
      AliError(MCMError(kAdcChannelsMiss));
    mcmcount++;
    // continue with next MCM
  }

  // check for missing MCMs (if header suppression is inactive)
  if (mcmcount != mcmcountExp) {
    AliError(LinkError(kMissMcmHeaders,
		       Form("%i not as expected: %i", mcmcount, mcmcountExp) ));
  }

  return (fPayloadCurr - start);
}

Bool_t AliTRDrawStream::ConnectTracklets(TTree *trklTree) 
{
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


TString AliTRDrawStream::EquipmentError(ErrorCode_t err, TString msg)
{ 
  // register error according to error code on equipment level 
  // and return the corresponding error message

  fLastError.fSector = fCurrEquipmentId - 1024;
  fLastError.fStack  = -1;
  fLastError.fLink   = -1;
  fLastError.fRob    = -1;
  fLastError.fMcm    = -1;
  fLastError.fError  = err;
  fErrors->Fill();

  return (Form("Event %6i: Eq. %2d - %s : ", 
	       fRawReader->GetEventIndex(), fCurrEquipmentId, fgErrorMessages[err]) + msg); 
}										


TString AliTRDrawStream::StackError(ErrorCode_t err, TString msg)
{ 
  // register error according to error code on stack level 
  // and return the corresponding error message

  fLastError.fSector = fCurrEquipmentId - 1024;
  fLastError.fStack  = fCurrSlot;
  fLastError.fLink   = -1;
  fLastError.fRob    = -1;
  fLastError.fMcm    = -1;
  fLastError.fError  = err;
  fErrors->Fill();

  return (Form("Event %6i: Eq. %2d S %i - %s : ", 
	       fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fgErrorMessages[err]) + msg); 
} 


TString AliTRDrawStream::LinkError(ErrorCode_t err, TString msg)
{ 
  // register error according to error code on link level 
  // and return the corresponding error message

  fLastError.fSector = fCurrEquipmentId - 1024;
  fLastError.fStack  = fCurrSlot;
  fLastError.fLink   = fCurrLink;
  fLastError.fRob    = -1;
  fLastError.fMcm    = -1;
  fLastError.fError  = err;
  fErrors->Fill();

  return (Form("Event %6i: Eq. %2d S %i l %2i - %s : ", 
	       fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fCurrLink, fgErrorMessages[err]) + msg); 
} 


TString AliTRDrawStream::ROBError(ErrorCode_t err, TString msg)
{ 
  // register error according to error code on ROB level 
  // and return the corresponding error message

  fLastError.fSector = fCurrEquipmentId - 1024;
  fLastError.fStack  = fCurrSlot;
  fLastError.fLink   = fCurrLink;
  fLastError.fRob    = fCurrRobPos;
  fLastError.fMcm    = -1;
  fLastError.fError  = err;
  fErrors->Fill();

  return (Form("Event %6i: Eq. %2d S %i l %2i ROB %i - %s : ", 
	       fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fCurrLink, fCurrRobPos, fgErrorMessages[err]) + msg); 
} 


TString AliTRDrawStream::MCMError(ErrorCode_t err, TString msg)
{ 
  // register error according to error code on MCM level 
  // and return the corresponding error message

  fLastError.fSector = fCurrEquipmentId - 1024;
  fLastError.fStack  = fCurrSlot;
  fLastError.fLink   = fCurrLink;
  fLastError.fRob    = fCurrRobPos;
  fLastError.fMcm    = fCurrMcmPos;
  fLastError.fError  = err;
  fErrors->Fill();

  return (Form("Event %6i: Eq. %2d S %i l %2i ROB %i MCM %2i - %s : ", 
	       fRawReader->GetEventIndex(), fCurrEquipmentId, fCurrSlot, fCurrLink, fCurrRobPos, fCurrMcmPos, fgErrorMessages[err]) + msg); 
}

const char* AliTRDrawStream::GetErrorMessage(ErrorCode_t errCode)
{ 
  // return the error message for the given error code

  if (errCode > 0 && errCode < kLastErrorCode) 
    return fgErrorMessages[errCode];
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

Bool_t AliTRDrawStream::DumpingMCM(Int_t det, Int_t rob, Int_t mcm) 
{
  // check if MCM data should be dumped

  for (Int_t iMCM = 0; iMCM < fNDumpMCMs; iMCM++) {
    if (fDumpMCM[iMCM] == ((det << 7) | (rob << 4) | mcm)) {
      return kTRUE;
    }
  }
  return kFALSE;
}

void AliTRDrawStream::DumpRaw(TString title, UInt_t *start, Int_t length)
{
  // dump raw data

  title += "\n";
  Int_t pos = 0;
  for ( ; pos+3 < length; pos += 4) {
    title += Form("0x%08x 0x%08x 0x%08x 0x%08x\n", 
		  start[pos+0], start[pos+1], start[pos+2], start[pos+3]);
  }
  for ( ; pos < length; pos++) {
    title += Form("0x%08x ", start[pos]);
  }
  AliInfo(title);
}
