// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCData.h"

#include <Alieve/AliEveTPCSectorData.h>

#include <AliSimDigits.h>
#include <AliTPCParam.h>
#include <AliTPCRawStream.h>
#include <TTree.h>


//______________________________________________________________________________
// AliEveTPCData
//
// A central manager for TPC data of an event.  Can read digits (from
// a tree: LoadDigits()) and raw-data (via AliRawReader: LoadRaw()).
//
// The sector data is stored in 36 AliEveTPCSectorData objects.
// Sectors 0 - 17: +z side, 18 - 35: -z side.
// No separation of inner/outer segments, use row numbers for addressing.
//
// Threshold application and pedestal subtraction can be performed at
// load time: use SetLoadThreshold(thresh) and SetLoadPedestal(ped).
//
// For raw-data (loaded using LoadRaw) pedestals can be calculated
// automatically per pad. Use SetAutoPedestal(kTRUE) to activate it.
// You might still want to set load threshold (default iz zero).
//

ClassImp(AliEveTPCData)

AliEveTPCData::AliEveTPCData() :
  fSectors(36), fSectorBlockSize(65536),
  fLoadThreshold(0), fLoadPedestal(0), fAutoPedestal(kFALSE)
{
  AliEveTPCSectorData::InitStatics();
}

AliEveTPCData::~AliEveTPCData()
{
  DeleteAllSectors();
}

/******************************************************************************/

void AliEveTPCData::CreateSector(Int_t sector)
{
  if(fSectors[sector] == 0)
    fSectors[sector] = new AliEveTPCSectorData(sector, fSectorBlockSize);
}

void AliEveTPCData::CreateAllSectors()
{
  for(Int_t s=0; s<36; ++s)
    CreateSector(s);
}

void AliEveTPCData::DropAllSectors()
{
  for(Int_t s=0; s<36; ++s) {
    if(fSectors[s] != 0)
      fSectors[s]->DropData();
  }
}

void AliEveTPCData::DeleteAllSectors()
{
  for(Int_t s=0; s<36; ++s) {
    delete fSectors[s];
    fSectors[s] = 0;
  }
}

/******************************************************************************/

AliEveTPCSectorData* AliEveTPCData::GetSectorData(Int_t sector, Bool_t spawnSectors)
{
  if(sector < 0 || sector > 35) return 0;
  if(fSectors[sector] == 0 && spawnSectors)
    CreateSector(sector);
  return fSectors[sector];
}

/******************************************************************************/

void AliEveTPCData::LoadDigits(TTree* tree, Bool_t spawnSectors)
{
  // Load data from TTree of AliSimDigits.
  // If spawnSectors is false only sectors that have been created previously
  // via CreateSector() are loaded.
  // If spawnSectors is true sectors are created if data for them is encountered.

  AliSimDigits digit, *digitPtr = &digit;
  tree->GetBranch("Segment")->SetAddress(&digitPtr);

  Int_t sector, row, pad, curPad;
  Short_t time, signal;
  Bool_t  inFill = kFALSE;
  AliEveTPCSectorData* secData = 0;

  Int_t numEnt = (Int_t) tree->GetEntries();
  for (Int_t ent=0; ent<numEnt; ent++) {
    tree->GetEntry(ent);
    AliEveTPCSectorData::GetParam().AdjustSectorRow(digit.GetID(), sector, row);
    if(sector >= 36) {
      sector -= 36;
      row    += AliEveTPCSectorData::GetInnSeg().GetNRows();
    }
    secData = GetSectorData(sector, spawnSectors);
    if(secData == 0)
      continue;

    if(digit.First() == kFALSE)
      continue;
    curPad = -1;
    do {
      pad    = digit.CurrentColumn();
      time   = digit.CurrentRow();
      signal = digit.CurrentDigit();

      if(pad != curPad) {
	if(inFill)
	  secData->EndPad(fAutoPedestal, fLoadThreshold);
	secData->BeginPad(row, pad, kFALSE);
	curPad = pad;
	inFill = kTRUE;
      }
      if(fAutoPedestal) {
	secData->RegisterData(time, signal);
      } else {
	signal -= fLoadPedestal;
	if(signal >= fLoadThreshold)
	  secData->RegisterData(time, signal);
      }

    } while (digit.Next());
    if(inFill) {
      secData->EndPad(fAutoPedestal, fLoadThreshold);
      inFill = kFALSE;
    }
  }
}

void AliEveTPCData::LoadRaw(AliTPCRawStream& input, Bool_t spawnSectors, Bool_t warn)
{
  // Load data from AliTPCRawStream.
  // If spawnSectors is false only sectors that have been created previously
  // via CreateSector() are loaded.
  // If spawnSectors is true sectors are created if data for them is encountered.

  static const TEveException eH("AliEveTPCData::LoadRaw ");

  Int_t   sector = -1, row = -1, pad = -1, rowOffset = 0;
  Short_t time,  signal;
  Bool_t  inFill   = kFALSE;
  Short_t lastTime = 9999;
  Bool_t  lastTimeWarn = kFALSE;
  AliEveTPCSectorData* secData = 0;

  Short_t threshold = fLoadThreshold;

  while (input.Next()) {
    if (input.IsNewSector()) {
      if(inFill) {
	secData->EndPad(fAutoPedestal, threshold);
	inFill = kFALSE;
      }
      sector = input.GetSector();
      if(sector >= 36) {
	sector -= 36;
	rowOffset = AliEveTPCSectorData::GetInnSeg().GetNRows();
      } else {
	rowOffset = 0;
      }
      secData = GetSectorData(sector, spawnSectors);
    }
    if (secData == 0)
      continue;

    if (input.IsNewPad()) {
      if(inFill) {
	secData->EndPad(fAutoPedestal, threshold);
	inFill = kFALSE;
      }
      row = input.GetRow() + rowOffset;
      pad = input.GetPad();

      if(pad >= AliEveTPCSectorData::GetNPadsInRow(row)) {
	if(warn) {
	  Warning(eH.Data(), "pad out of range (row=%d, pad=%d, maxpad=%d).",
		  row, pad, AliEveTPCSectorData::GetNPadsInRow(row));
	}
	continue;
      }

      AliEveTPCSectorData::PadRowHack* prh = secData->GetPadRowHack(row, pad);
      if(prh != 0) {
	threshold = prh->fThrExt + Short_t(prh->fThrFac*fLoadThreshold);
      } else {
	threshold = fLoadThreshold;
      }

      secData->BeginPad(row, pad, kTRUE);
      inFill   = kTRUE;
      lastTime = 1024;  lastTimeWarn = kFALSE;
    }

    time   = input.GetTime();
    signal = input.GetSignal();
    if(time >= lastTime) {
      if(lastTimeWarn == kFALSE) {
	if(warn)
	  Warning(eH.Data(), "time out of order (row=%d, pad=%d, time=%d, lastTime=%d).",
		  row, pad, time, lastTime);
        lastTimeWarn = kTRUE;
      }
      continue;
    }
    lastTime = time;
    if(fAutoPedestal) {
      secData->RegisterData(time, signal);
    } else {
      signal -= fLoadPedestal;
      if(signal > threshold)
	secData->RegisterData(time, signal);
    }
  }

  if(inFill) {
    secData->EndPad(fAutoPedestal, threshold);
    inFill = kFALSE;
  }
}
