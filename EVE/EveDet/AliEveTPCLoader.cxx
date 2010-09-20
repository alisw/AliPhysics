// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCLoader.h"
#include "AliEveTPCData.h"
#include <EveDet/AliEveTPCSector2D.h>
#include <EveDet/AliEveTPCSector3D.h>
#include <TEveManager.h>
#include <TEveGedEditor.h>

#include <AliRawReaderRoot.h>
#include <AliTPCRawStreamV3.h>

#include <TSystem.h>

//==============================================================================
//==============================================================================
// AliEveTPCLoader
//==============================================================================

//______________________________________________________________________________
//
// Front-end for stand-alone inspection of TPC raw-data.
// GUI is implemented in AliEveTPCLoaderEditor class.

ClassImp(AliEveTPCLoader)

AliEveTPCLoader::AliEveTPCLoader(const Text_t* n, const Text_t* t) :
  TEveElementList(n, t),

  fFile(),
  fEvent(-1),
  fDoubleSR(kFALSE),

  fTPCEquipementMap(),
  fReader(0),
  fData(0),

  fSec2Ds(36),
  fSec3Ds(36),

  fSetInitSectorParams(kFALSE),
  fInitMinTime(0), fInitMaxTime(460), fInitThreshold(5), fInitMaxVal(128),
  fCutOnEta(kFALSE), fEtaMax(1.5), fEtaMin(-1.5)
{
  // Constructor.

  fData = new AliEveTPCData;
}

AliEveTPCLoader::~AliEveTPCLoader()
{
  // Destructor.

  delete fReader;
  delete fData;
}

/******************************************************************************/

void AliEveTPCLoader::RemoveElementLocal(TEveElement* el)
{
  // Local removal an element, virtual from TEveElement.
  // Need to search for it in visual sector representations and remove
  // it there as well.

  for (Int_t i=0; i<36; ++i) {
    if (fSec2Ds[i] == el) fSec2Ds[i] = 0;
    if (fSec3Ds[i] == el) fSec3Ds[i] = 0;
  }
}

void AliEveTPCLoader::RemoveElementsLocal()
{
  // Local removal of all elements, virtual from TEveElement.
  // Must remove all visual sector representations.

  for (Int_t i=0; i<36; ++i) {
    fSec2Ds[i] = 0;
    fSec3Ds[i] = 0;
  }
}

/******************************************************************************/

void AliEveTPCLoader::SetData(AliEveTPCData* d)
{
  // Set external TPC-data container. The old one is deleted.

  delete fData;
  fData = d;
}

/******************************************************************************/

void AliEveTPCLoader::OpenFile()
{
  // Open the raw-data file, as set in fFile data-member.
  // Raw-reader is also created and attached to the file.
  // First event is loaded and all sectors for which the data-exists
  // are created.

  static const TEveException kEH("AliEveTPCLoader::OpenFile ");

  if (gSystem->AccessPathName(fFile, kReadPermission))
      throw(kEH + "can not read '" + fFile + "'.");

  fData->DeleteAllSectors();

  delete fReader;
  fReader =  0;
  fEvent  = -1;

  fReader = new AliRawReaderRoot(fFile);
  if (fTPCEquipementMap != "")
    fReader->LoadEquipmentIdsMap
      (gSystem->ExpandPathName(fTPCEquipementMap.Data()));

  NextEvent();
  LoadEvent();
  UpdateSectors(kTRUE);
}

void AliEveTPCLoader::LoadEvent()
{
  // Load an event.

  static const TEveException kEH("AliEveTPCLoader::LoadEvent ");

  if (fReader == 0)
    throw(kEH + "data file not opened.");

  printf("Now loading event %d\n", fEvent);
  fReader->Reset();
  AliTPCRawStreamV3 input(fReader);
  fReader->Select("TPC");

  fData->DropAllSectors();
  fData->LoadRaw(input, kTRUE, kTRUE);
}

void AliEveTPCLoader::NextEvent(Bool_t rewindOnEnd)
{
  // Go o the next event.
  // When the last event is reached and rewindOnEnd is true, the file
  // is rewound back to the first event. Otherwise an exception is thrown.

  static const TEveException kEH("AliEveTPCLoader::NextEvent ");

  if (fReader == 0)
    throw(kEH + "data file not opened.");

  if (fReader->NextEvent() == kTRUE) {
    ++fEvent;
  } else {
    if (fEvent == -1)
      throw(kEH + "no events available.");
    if (rewindOnEnd) {
      printf("Reached end of stream (event=%d), rewinding to first event.\n", fEvent);
      fReader->RewindEvents();
      fReader->NextEvent();
      fEvent = 0;
    } else {
      throw(kEH + "last event reached.");
    }
  }
}

void AliEveTPCLoader::GotoEvent(Int_t event)
{
  // Go to specified event.

  static const TEveException kEH("AliEveTPCLoader::GotoEvent ");

  if (fReader == 0)
    throw(kEH + "data file not opened.");

  if (event == fEvent)
    return;
  Bool_t checkEnd;
  if (event < fEvent) {
    fReader->RewindEvents();
    fEvent = -1;
    checkEnd = kFALSE;
  } else {
    checkEnd = kTRUE;
  }
  do {
    NextEvent();
  } while (fEvent != event && !(checkEnd == kTRUE && fEvent == 0));
  LoadEvent();
  UpdateSectors();
}

void* AliEveTPCLoader::LoopEvent(AliEveTPCLoader* loader)
{
  // Loop to next event on given loader. Static member to be used for
  // external control during animations.
  // Calls NextEvent(), LoadEvent() and UpdateSectors().

  loader->NextEvent();
  loader->LoadEvent();
  loader->UpdateSectors();
  if (gEve->GetEditor()->GetModel() == loader)
    gEve->EditElement(loader);
  return 0;
}

/******************************************************************************/

void AliEveTPCLoader::UpdateSectors(Bool_t dropNonPresent)
{
  // Update visual representations of sectors.
  // If dropNonPresent is true, the sectors for which there is no data in
  // the current event are removed.

  gEve->DisableRedraw();
  for (Int_t i=0; i<=35; ++i)
  {
    AliEveTPCSectorData* sd = fData->GetSectorData(i);

    // 2D sectors
    if (fSec2Ds[i] != 0)
    {
      if (dropNonPresent && sd == 0) {
	gEve->RemoveElement(fSec2Ds[i], this);
	fSec2Ds[i] = 0;
      } else {
	fSec2Ds[i]->IncRTS();
        fSec2Ds[i]->ElementChanged();
      }
    }
    else
    {
      if (sd != 0) {
	AliEveTPCSector2D* s = new AliEveTPCSector2D(Form("Sector2D %d", i));
	fSec2Ds[i] = s;
	s->SetSectorID(i);
	s->SetDataSource(fData);

	if (fDoubleSR)
	  s->SetMaxTime(1023);

        if (fSetInitSectorParams) {
          s->SetMinTime(fInitMinTime);
          s->SetMaxTime(fInitMaxTime);
          s->SetThreshold(fInitThreshold);
	  s->SetMaxVal(fInitMaxVal);
        }

	s->SetAutoTrans(kTRUE);
	s->SetFrameColor(36);

	gEve->AddElement(s, this);
      }
    }

    // 3D sectors
    if (fSec3Ds[i] != 0)
    {
      if (dropNonPresent && sd == 0) {
	gEve->RemoveElement(fSec3Ds[i], this);
	fSec3Ds[i] = 0;
      } else {
        fSec3Ds[i]->SetCutOnEta(fCutOnEta);
        fSec3Ds[i]->SetEtaMax(fEtaMax);
        fSec3Ds[i]->SetEtaMin(fEtaMin);
	fSec3Ds[i]->IncRTS();
        fSec3Ds[i]->ElementChanged();
      }
    }
  }
  gEve->Redraw3D(kTRUE, kFALSE);
  gEve->EnableRedraw();
}

void AliEveTPCLoader::ReloadSectors()
{
  // Reload current event and update sectors.

  LoadEvent();
  UpdateSectors();
}

void AliEveTPCLoader::CreateSectors3D()
{
  // Create 3D representations of sectors.

  gEve->DisableRedraw();
  for (Int_t i=0; i<=35; ++i) {
    AliEveTPCSectorData* sd = fData->GetSectorData(i);
    if (sd != 0 && fSec3Ds[i] == 0) {
      AliEveTPCSector3D* s = new AliEveTPCSector3D(Form("Sector3D %d", i));
      fSec3Ds[i] = s;
      s->SetSectorID(i);
      s->SetDataSource(fData);
      s->SetCutOnEta(fCutOnEta);
      s->SetEtaMax(fEtaMax);
      s->SetEtaMin(fEtaMin);

      if (fDoubleSR)
	s->SetDriftVel(2.273);
      if (fSec2Ds[i] != 0)
	s->CopyVizParams(fSec2Ds[i]);

      s->SetAutoTrans(kTRUE);
      s->SetFrameColor(36);

      gEve->AddElement(s, this);
    }
  }
  gEve->EnableRedraw();
}

void AliEveTPCLoader::DeleteSectors3D()
{
  // Delete 3D representations of sectors.

  gEve->DisableRedraw();
  for (Int_t i=0; i<=35; ++i) {
    TEveElement* re = fSec3Ds[i];
    if (re != 0) {
      gEve->RemoveElement(re, this);
      // delete re; // Done automatically.
      fSec3Ds[i] = 0;
    }
  }
  gEve->EnableRedraw();
}

/******************************************************************************/

void AliEveTPCLoader::SetInitParams(Int_t mint, Int_t maxt, Int_t thr, Int_t maxval)
{
  // Set initial viualization parameters for 2D and 3D sector representations.

  fSetInitSectorParams = kTRUE;
  fInitMinTime   = mint;
  fInitMaxTime   = maxt;
  fInitThreshold = thr;
  fInitMaxVal    = maxval;
}
