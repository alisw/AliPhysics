// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveEventManager.h"
#include <TEveManager.h>

#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliRawReaderRoot.h>
#include <AliRawReaderFile.h>
#include <AliRawReaderDate.h>
#include <AliMagFMaps.h>
#include <AliCDBManager.h>
#include <AliHeader.h>
#include <AliGeomManager.h>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

//==============================================================================
//==============================================================================
// AliEveEventManager
//==============================================================================

//______________________________________________________________________________
//
// Provide interface for loading and navigating standard AliRoot data
// (AliRunLoader) and ESDs.
//
// Missing support for raw-data. For now this is handled individually
// by each sub-detector.
//
// Also provides interface to magnetic-field and geometry. Mostly
// intended as wrappers over standard AliRoot functionality for
// convenient use from visualizateion macros.

ClassImp(AliEveEventManager)

AliEveEventManager* gAliEveEvent = 0;

Bool_t AliEveEventManager::fgAssertRunLoader = kFALSE;
Bool_t AliEveEventManager::fgAssertESD       = kFALSE;
Bool_t AliEveEventManager::fgAssertRaw       = kFALSE;

TString  AliEveEventManager::fgESDFileName("AliESDs.root");
TString  AliEveEventManager::fgRawFileName("raw.root");
TString  AliEveEventManager::fgCdbUri("local://$ALICE_ROOT");

AliMagF* AliEveEventManager::fgMagField = 0;


AliEveEventManager::AliEveEventManager() :
  TEveEventManager(),

  fPath      ( ), fEventId (-1),
  fRunLoader (0),
  fESDFile   (0), fESDTree (0), fESD (0),
  fESDfriend (0), fESDfriendExists(kFALSE),
  fRawReader (0)
{
  // Default constructor.
}

AliEveEventManager::AliEveEventManager(TString path, Int_t ev) :
  TEveEventManager("AliEVE AliEveEventManager"),

  fPath   (path), fEventId(-1),
  fRunLoader (0),
  fESDFile   (0), fESDTree (0), fESD (0),
  fESDfriend (0), fESDfriendExists(kFALSE),
  fRawReader (0)
{
  // Constructor with event-directory URL and event-id.

  Open();
  if (ev >= 0) GotoEvent(ev);
}

AliEveEventManager::~AliEveEventManager()
{
  // Destructor.

  // Somewhat unclear what to do here.
  // In principle should close all data sources and deregister from
  // TEveManager.
}

/******************************************************************************/

void AliEveEventManager::SetESDFileName(const Text_t* esd)
{
  // Set file-name for opening ESD, default "AliESDs.root".

  if (esd) fgESDFileName = esd;
}

void AliEveEventManager::SetRawFileName(const Text_t* raw)
{
  // Set file-name for opening of raw-data, default "raw.root"
  if (raw) fgRawFileName = raw;
}

void AliEveEventManager::SetCdbUri(const Text_t* cdb)
{
  // Set path to CDB, default "local://$ALICE_ROOT".

  if (cdb) fgCdbUri = cdb;
}

void AliEveEventManager::SetAssertElements(Bool_t assertRunloader,
                                           Bool_t assertEsd,
                                           Bool_t assertRaw)
{
  // Set global flags that detrmine which parts of the event-data must
  // be present when the event is opened.

  fgAssertRunLoader = assertRunloader;
  fgAssertESD = assertEsd;
  fgAssertRaw = assertRaw;
}

/******************************************************************************/

void AliEveEventManager::Open()
{
  // Open event-data from URL specified in fPath.
  // Attempts to create AliRunLoader() and to open ESD with ESDfriends.
  // Warning is reported if run-loader or ESD is not found.
  // Global data-members fgAssertRunLoader and fgAssertESD can be set
  // to throw exceptions instead.

  static const TEveException kEH("AliEveEventManager::Open ");

  gSystem->ExpandPathName(fPath);
  if (fPath[0] != '/')
    fPath = Form("%s/%s", gSystem->WorkingDirectory(), fPath.Data());

  Int_t runNo = -1;

  TString gaPath(Form("%s/galice.root", fPath.Data()));
  // If i use open directly, we get fatal.
  // Is this (AccessPathName check) ok for xrootd / alien?
  if (gSystem->AccessPathName(gaPath, kReadPermission) == kFALSE)
  {
    fRunLoader = AliRunLoader::Open(gaPath);
    if (fRunLoader)
    {
      TString alicePath = fPath + "/";
      fRunLoader->SetDirName(alicePath);

      if (fRunLoader->LoadgAlice() != 0)
        Warning(kEH, "failed loading gAlice via run-loader.");

      if (fRunLoader->LoadHeader() == 0)
      {
        runNo = fRunLoader->GetHeader()->GetRun();
      }
      else
      {
        Warning(kEH, "failed loading run-loader's header.");
        delete fRunLoader;
        fRunLoader = 0;
      }
    }
    else // run-loader open failed
    {
      Warning(kEH, "failed opening ALICE run-loader from '%s'.", gaPath.Data());
    }
  }
  else // galice not readable
  {
    Warning(kEH, "can not read '%s'.", gaPath.Data());
  }
  if (fRunLoader == 0)
  {
    if (fgAssertRunLoader)
      throw (kEH + "Bootstraping of run-loader failed. Its precence was requested.");
    else
      Warning(kEH, "Bootstraping of run-loader failed.");
  }


  TString esdPath(Form("%s/%s", fPath.Data(), fgESDFileName.Data()));
  if ((fESDFile = TFile::Open(esdPath)))
  {
    fESD = new AliESDEvent();
    fESDTree = (TTree*) fESDFile->Get("esdTree");
    if (fESDTree != 0)
    {
      fESD->ReadFromTree(fESDTree);
      runNo = fESD->GetESDRun()->GetRunNumber();

      // Check if ESDfriends exists and attach the branch
      TString p = Form("%s/AliESDfriends.root", fPath.Data());
      if (gSystem->AccessPathName(p, kReadPermission) == kFALSE)
      {
        fESDfriendExists = kTRUE;
        fESDTree->SetBranchStatus ("ESDfriend*", 1);
        fESDTree->SetBranchAddress("ESDfriend.", &fESDfriend);
      }
    }
    else // esdtree == 0
    {
      delete fESDFile; fESDFile = 0;
      Warning(kEH, "failed getting the esdTree.");
    }
  }
  else // esd not readable
  {
    Warning(kEH, "can not read ESD file '%s'.", esdPath.Data());
  }
  if (fESDTree == 0)
  {
    if (fgAssertESD)
    {
      throw (kEH + "ESD not initialized. Its precence was requested.");
    } else {
      Warning(kEH, "ESD not initialized.");
    }
  }

  TString rawPath(Form("%s/%s", fPath.Data(), fgRawFileName.Data()));
  // If i use open directly, raw-reader reports an error but i have
  // no way to detect it.
  // Is this (AccessPathName check) ok for xrootd / alien?
  if (gSystem->AccessPathName(rawPath, kReadPermission) == kFALSE)
  {
    if (fgRawFileName.EndsWith("/"))
    {
      fRawReader = new AliRawReaderFile(rawPath);
    }
    else if (fgRawFileName.EndsWith(".root"))
    {
      fRawReader = new AliRawReaderRoot(rawPath);
    }
    else if (!fgRawFileName.IsNull())
    {
      fRawReader = new AliRawReaderDate(rawPath);
    } 
  }

  if (fRawReader == 0)
  {
    if (fgAssertRaw)
    {
      throw (kEH + "raw-data not initialized. Its precence was requested.");
    } else {
      Warning(kEH, "raw-data not initialized.");
    }
  }

  if (runNo < 0)
  {
    if (fRawReader)
    {
      fRawReader->NextEvent();
      runNo = fRawReader->GetRunNumber();
      printf("Determining run-no from raw ... run=%d\n", runNo);
      fRawReader->RewindEvents();
    } else {
      throw (kEH + "unknown run number.");
    }
  }

  {
    AliCDBManager* cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage(fgCdbUri);
    if (cdb->IsDefaultStorageSet() == kFALSE)
      throw (kEH + "CDB initialization failed.");
    cdb->SetRun(runNo);
  }

  SetName(Form("Event %d", fEventId));
  SetTitle(fPath);
}

void AliEveEventManager::GotoEvent(Int_t event)
{
  // Load data for specified event.
  // If event is out of range an exception is thrown and old state
  // is preserved.
  // After successful loading of event, the virtual function
  // AfterNewEventLoaded() is called. This executes commands that
  // were registered via TEveEventManager::AddNewEventCommand().

  static const TEveException kEH("AliEveEventManager::GotoEvent ");

  if (event < 0) {
    Error(kEH, "event must be non-negative.");
    return;
  }

  Int_t maxEvent = 0;
  if (fRunLoader) {
    maxEvent = fRunLoader->GetNumberOfEvents() - 1;
  } else if (fESDTree) {
    maxEvent = fESDTree->GetEntries() - 1;
  } else if (fRawReader) {
    maxEvent = 10000000;
    Info(kEH, "number of events unknown for raw-data, setting max-event id to 10M.");
  } else {
    throw (kEH + "neither RunLoader, ESD nor Raw loaded.");
  }
  if (event < 0 || event > maxEvent)
    throw (kEH + Form("event %d not present, available range [%d, %d].",
                      event, 0, maxEvent));

  TEveManager::TRedrawDisabler rd(gEve);
  gEve->Redraw3D(kFALSE, kTRUE); // Enforce drop of all logicals.

  // !!! MT this is somewhat brutal; at least optionally, one could be
  // a bit gentler, checking for objs owning their external refs and having
  // additinal parents.
  DestroyElements();

  if (fRunLoader) {
    if (fRunLoader->GetEvent(event) != 0)
      throw (kEH + "failed getting required event.");
  }

  if (fESDTree) {
    if (fESDTree->GetEntry(event) <= 0)
      throw (kEH + "failed getting required event from ESD.");

    if (fESDfriendExists)
      fESD->SetESDfriend(fESDfriend);
  }

  if (fRawReader)
  {
    Int_t rawEv = fEventId;
    if (event < rawEv)
    {
      fRawReader->RewindEvents();
      rawEv = -1;
    }

    while (rawEv < event)
    {
      if ( ! fRawReader->NextEvent())
      {
        fRawReader->RewindEvents();
        fEventId = -1;
        throw (kEH + Form("Error going to next raw-event from event %d.", rawEv));
      }
      ++rawEv;
    }

    printf ("Loaded raw-event %d.\n", rawEv);
  }

  fEventId = event;
  SetName(Form("Event %d", fEventId));
  UpdateItems();

  AfterNewEventLoaded();
}

void AliEveEventManager::Close()
{
  // Close the event files.
  // For the moment only ESD is closed. Needs to be investigated for
  // AliRunLoader and Raw.

  if (fESDTree) {
    delete fESD;       fESD       = 0;
    delete fESDfriend; fESDfriend = 0;

    delete fESDTree; fESDTree = 0;
    delete fESDFile; fESDFile = 0;
  }
}


/******************************************************************************/
// Static convenience functions, mainly used from macros.
/******************************************************************************/

AliRunLoader* AliEveEventManager::AssertRunLoader()
{
  // Make sure AliRunLoader is initialized and return it.
  // Throws exception in case run-loader is not available.
  // Static utility for macros.

  static const TEveException kEH("AliEveEventManager::AssertRunLoader ");

  if (gAliEveEvent == 0)
    throw (kEH + "ALICE event not ready.");
  if (gAliEveEvent->fRunLoader == 0)
    throw (kEH + "AliRunLoader not initialised.");
  return gAliEveEvent->fRunLoader;
}

AliESDEvent* AliEveEventManager::AssertESD()
{
  // Make sure AliESDEvent is initialized and return it.
  // Throws exception in case ESD is not available.
  // Static utility for macros.

  static const TEveException kEH("AliEveEventManager::AssertESD ");

  if (gAliEveEvent == 0)
    throw (kEH + "ALICE event not ready.");
  if (gAliEveEvent->fESD == 0)
    throw (kEH + "AliESD not initialised.");
  return gAliEveEvent->fESD;
}

AliESDfriend* AliEveEventManager::AssertESDfriend()
{
  // Make sure AliESDfriend is initialized and return it.
  // Throws exception in case ESDfriend-loader is not available.
  // Static utility for macros.

  static const TEveException kEH("AliEveEventManager::AssertESDfriend ");

  if (gAliEveEvent == 0)
    throw (kEH + "ALICE event not ready.");
  if (gAliEveEvent->fESDfriend == 0)
    throw (kEH + "AliESDfriend not initialised.");
  return gAliEveEvent->fESDfriend;
}

AliRawReader* AliEveEventManager::AssertRawReader()
{
  // Make sure raw-reader is initialized and return it.

  static const TEveException kEH("AliEveEventManager::AssertRawReader ");

  if (gAliEveEvent == 0)
    throw (kEH + "ALICE event not ready.");
  if (gAliEveEvent->fRawReader == 0)
    throw (kEH + "RawReader not ready.");

  return gAliEveEvent->fRawReader;
}

AliMagF* AliEveEventManager::AssertMagField()
{
  // Make sure AliMagF is initialized and return it.
  // Throws exception in case magnetic field is not available.
  // Static utility for macros.

  if (fgMagField == 0)
  {
    if (gAliEveEvent && gAliEveEvent->fRunLoader && gAliEveEvent->fRunLoader->GetAliRun())
      fgMagField = gAliEveEvent->fRunLoader->GetAliRun()->Field();
    else
      fgMagField = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
  }
  return fgMagField;
}

TGeoManager* AliEveEventManager::AssertGeometry()
{
  // Make sure AliGeomManager is initialized and returns the
  // corresponding TGeoManger.
  // gGeoManager is not set, maybe it should be.
  // Throws exception in case run-loader is not available.
  // Static utility for macros.

  // !!!! Should we set gGeoManager here?

  static const TEveException kEH("AliEveEventManager::AssertGeometry ");

  if (AliGeomManager::GetGeometry() == 0)
  {
    gGeoManager = 0;
    AliGeomManager::LoadGeometry();
    if ( ! AliGeomManager::GetGeometry())
    {
      throw (kEH + "can not load geometry.");
    }
    if ( ! AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO ACORDE"))
    {
      ::Warning(kEH, "mismatch of alignable volumes. Proceeding.");
      // throw (kEH + "could not apply align objs.");
    }

    // Temporary fix.
    // In AliEve several simplified geometries can be loaded at a later stage.
    // Should handle this inTEveManager::GetGeometry() by properly
    // unlocking and re-locking the geo-manager.
    TGeoManager::UnlockGeometry();
  }

  return AliGeomManager::GetGeometry();
}
