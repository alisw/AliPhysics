// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveEventManager.h"
#include "AliEveMacroExecutor.h"
#include <TEveManager.h>

#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliDAQ.h>
#include <AliRawEventHeaderBase.h>
#include <AliRawReaderRoot.h>
#include <AliRawReaderFile.h>
#include <AliRawReaderDate.h>
#include <AliMagFMaps.h>
#include <AliCDBManager.h>
#include <AliHeader.h>
#include <AliGeomManager.h>

#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TSystem.h>
#include <TTimeStamp.h>

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
  fRawReader (0),
  fAutoLoad(kFALSE),
  fAutoLoadTime(5.),
  fAutoLoadTimer(0),
  fIsOnline(kFALSE),
  fExecutor(new AliEveMacroExecutor)
{
  // Default constructor.
}

AliEveEventManager::AliEveEventManager(TString path, Int_t ev) :
  TEveEventManager("AliEVE AliEveEventManager"),

  fPath   (path), fEventId(-1),
  fRunLoader (0),
  fESDFile   (0), fESDTree (0), fESD (0),
  fESDfriend (0), fESDfriendExists(kFALSE),
  fRawReader (0),
  fAutoLoad(kFALSE),
  fAutoLoadTime(5.),
  fAutoLoadTimer(0),
  fIsOnline(kFALSE),
  fExecutor(new AliEveMacroExecutor)
{
  // Constructor with event-directory URL and event-id.

  Open();
  if (ev >= 0) GotoEvent(ev);
}

AliEveEventManager::~AliEveEventManager()
{
  // Destructor.

  if (fAutoLoadTimer) delete fAutoLoadTimer;
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
  // The following magick is required for ESDriends to be loaded properly
  // from non-current directory.
  if (fPath.IsNull() || fPath == ".")
  {
    fPath = gSystem->WorkingDirectory();
  }
  else if ( ! fPath.BeginsWith("file:/"))
  {
    TUrl    url(fPath, kTRUE);
    TString protocol(url.GetProtocol());
    if (protocol == "file" && fPath[0] != '/')
      fPath = Form("%s/%s", gSystem->WorkingDirectory(), fPath.Data());
  }

  Int_t runNo = -1;

  // Open ESD and ESDfriends

  TString esdPath(Form("%s/%s", fPath.Data(), fgESDFileName.Data()));
  if ((fESDFile = TFile::Open(esdPath)))
  {
    fESD = new AliESDEvent();
    fESDTree = (TTree*) fESDFile->Get("esdTree");
    if (fESDTree != 0)
    {
      // Check if ESDfriends exists and attach the branch
      TString p(Form("%s/AliESDfriends.root", fPath.Data()));
      TFile *esdFriendFile = TFile::Open(p);
      if (esdFriendFile) {
	if (!esdFriendFile->IsZombie())
	  {
	    esdFriendFile->Close();
	    delete esdFriendFile;
	    fESDfriendExists = kTRUE;
	    fESDTree->SetBranchStatus ("ESDfriend*", 1);
	    fESDTree->SetBranchAddress("ESDfriend.", &fESDfriend);
	  }
	else
	  {
	    esdFriendFile->Close();
	    delete esdFriendFile;
	  }
      }

      fESD->ReadFromTree(fESDTree);
      if (!fESDfriendExists) fESDTree->SetBranchStatus ("ESDfriend*", 0);
      if (fESDTree->GetEntry(0) <= 0)
	{
	  delete fESDFile; fESDFile = 0;
	  delete fESD; fESD = 0;
	  Warning(kEH, "failed getting the first entry from esdTree.");
	}
      else
	{
	  if (runNo < 0)
	    runNo = fESD->GetESDRun()->GetRunNumber();
	}
    }
    else // esdtree == 0
    {
      delete fESDFile; fESDFile = 0;
      delete fESD; fESD = 0;
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

  // Open RunLoader from galice.root

  TString gaPath(Form("%s/galice.root", fPath.Data()));
  // If i use open directly, we get fatal.
  // Is AccessPathName check ok for xrootd / alien? Yes, not for http.
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
        if (runNo < 0)
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

  // Open raw-data file

  TString rawPath(Form("%s/%s", fPath.Data(), fgRawFileName.Data()));
  // If i use open directly, raw-reader reports an error but i have 	 
  // no way to detect it. 	 
  // Is this (AccessPathName check) ok for xrootd / alien? Yes, not for http. 	 
  if (gSystem->AccessPathName(rawPath, kReadPermission) == kFALSE) 	 
  {
    fRawReader = AliRawReader::Create(rawPath);
  }
  else
  {
    fRawReader = AliRawReader::Create(fgRawFileName);
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
      Info(kEH, "Determining run-no from raw ... run=%d.", runNo);
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

void AliEveEventManager::SetEvent(AliRunLoader *runLoader, AliRawReader *rawReader, AliESDEvent *esd)
{
  // Set an event from an external source
  // The method is used in the online visualisation
  fRunLoader = runLoader;
  fRawReader = rawReader;
  fESD = esd;
  fIsOnline = kTRUE;
  SetTitle("Online event in memory");
  SetName("Online Event");

  ElementChanged();
  AfterNewEventLoaded();
}

Int_t AliEveEventManager::GetMaxEventId(Bool_t /*refreshESD*/) const
{
  // Returns maximum available event id.
  // If raw-data is the only data-source this can not be known
  // and 10,000,000 is returned.
  // If neither data-source is initialised an exception is thrown.
  // If refresh_esd is true and ESD is the primary event-data source
  // its header is re-read from disk.

  static const TEveException kEH("AliEveEventManager::GetMaxEventId ");

  if (fESDTree)
  {
    // Refresh crashes with root-5.21.1-alice.
    // Fixed by Philippe 5.8.2008 r25053, can be reactivated
    // when we move to a newer root.
    // if (refreshESD)
    //   fESDTree->Refresh();
    return fESDTree->GetEntries() - 1;
  }
  else if (fRunLoader)
  {
    return fRunLoader->GetNumberOfEvents() - 1;
  }
  else if (fRawReader)
  {
    Int_t n = fRawReader->GetNumberOfEvents() - 1;
    return n > -1 ? n : 10000000;
  }
  else
  {
    throw (kEH + "neither RunLoader, ESD nor Raw loaded.");
  }
}

void AliEveEventManager::GotoEvent(Int_t event)
{
  // Load data for specified event.
  // If event is out of range an exception is thrown and old state
  // is preserved.
  // After successful loading of event, the virtual function
  // AfterNewEventLoaded() is called. This executes commands that
  // were registered via TEveEventManager::AddNewEventCommand().
  //
  // If event is negative, it is subtracted from the number of
  // available events, thus passing -1 will load the last event.
  // This is not supported when raw-data is the only data-source
  // as the number of events is not known.

  static const TEveException kEH("AliEveEventManager::GotoEvent ");

  Int_t maxEvent = 0;
  if (fESDTree)
  {
    // Refresh crashes with root-5.21.1-alice.
    // Fixed by Philippe 5.8.2008 r25053, can be reactivated
    // when we move to a newer root.
    // fESDTree->Refresh();
    maxEvent = fESDTree->GetEntries() - 1;
    if (event < 0)
      event = fESDTree->GetEntries() + event;
  }
  else if (fRunLoader)
  {
    maxEvent = fRunLoader->GetNumberOfEvents() - 1;
    if (event < 0)
      event = fRunLoader->GetNumberOfEvents() + event;
  }
  else if (fRawReader)
  {
    maxEvent = fRawReader->GetNumberOfEvents() - 1;
    if (maxEvent < 0)
    {
      maxEvent = 10000000;
      if (event < 0) {
        Error(kEH, "current raw-data source does not support direct event access.");
        return;
      }
      Info(kEH, "number of events unknown for current raw-data source, setting max-event id to 10M.");
    }
    else
    {
      if (event < 0)
        event = fRawReader->GetNumberOfEvents() + event;
    }
  }
  else
  {
    throw (kEH + "neither RunLoader, ESD nor Raw loaded.");
  }
  if (event < 0 || event > maxEvent)
  {
    throw (kEH + Form("event %d not present, available range [%d, %d].",
                      event, 0, maxEvent));
  }

  TEveManager::TRedrawDisabler rd(gEve);
  gEve->Redraw3D(kFALSE, kTRUE); // Enforce drop of all logicals.

  // !!! MT this is somewhat brutal; at least optionally, one could be
  // a bit gentler, checking for objs owning their external refs and having
  // additinal parents.
  DestroyElements();

  if (fESDTree) {
    if (fESDTree->GetEntry(event) <= 0)
      throw (kEH + "failed getting required event from ESD.");

    if (fESDfriendExists)
      fESD->SetESDfriend(fESDfriend);
  }

  if (fRunLoader) {
    if (fRunLoader->GetEvent(event) != 0)
      throw (kEH + "failed getting required event.");
  }

  if (fRawReader)
  {
    // AliRawReader::GotoEvent(Int_t) works for AliRawReaderRoot/Chain.
    if (fRawReader->GotoEvent(event) == kFALSE)
    {
      // Use fallback method - iteration with NextEvent().
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
      Warning(kEH, "Loaded raw-event %d with fallback method.\n", rawEv);
    }
  }

  fEventId = event;
  SetName(Form("Event %d", fEventId));
  ElementChanged();

  AfterNewEventLoaded();
}

void AliEveEventManager::NextEvent()
{
  // Loads next event
  // either in automatic (online) or
  // manual mode
  
  if (fIsOnline)
  {
    if (fAutoLoadTimer) fAutoLoadTimer->Stop();

    DestroyElements();

    gSystem->ExitLoop();
  }
  else
  {
    if (fEventId < GetMaxEventId(kTRUE))
      GotoEvent(fEventId + 1);
    else
      GotoEvent(0);
    StartStopAutoLoadTimer();
  }
}

void AliEveEventManager::PrevEvent()
{
  // Loads previous event
  // only in case of manual mode
  if (!fIsOnline) {
    GotoEvent(fEventId - 1);
    StartStopAutoLoadTimer();
  }
}

void AliEveEventManager::Close()
{
  // Close the event data-files and delete ESD, ESDfriend, run-loader
  // and raw-reader.

  if (fESDTree) {
    delete fESD;       fESD       = 0;
    delete fESDfriend; fESDfriend = 0;

    delete fESDTree;   fESDTree = 0;
    delete fESDFile;   fESDFile = 0;
  }

  if (fRunLoader) {
    delete fRunLoader; fRunLoader = 0;
  }

  if (fRawReader) {
    delete fRawReader; fRawReader = 0;
  }
}


/******************************************************************************/
// Static convenience functions, mainly used from macros.
/******************************************************************************/

Bool_t AliEveEventManager::HasRunLoader()
{
  // Check if AliRunLoader is initialized.

  return gAliEveEvent && gAliEveEvent->fRunLoader;
}

Bool_t AliEveEventManager::HasESD()
{
  // Check if AliESDEvent is initialized.

  return gAliEveEvent && gAliEveEvent->fESD;
}

Bool_t AliEveEventManager::HasESDfriend()
{
  // Check if AliESDfriend is initialized.

  return gAliEveEvent && gAliEveEvent->fESDfriend;
}

Bool_t AliEveEventManager::HasRawReader()
{
  // Check if raw-reader is initialized.

  return gAliEveEvent && gAliEveEvent->fRawReader;
}

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
  // gGeoManager is set to the return value.
  // Throws exception if geometry can not be loaded or if it is not
  // available and the TGeoManager is locked.
  // Static utility for macros.

  static const TEveException kEH("AliEveEventManager::AssertGeometry ");

  if (AliGeomManager::GetGeometry() == 0)
  {
    if (TGeoManager::IsLocked())
      throw (kEH + "geometry is not loaded but TGeoManager is locked.");

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
    AliGeomManager::GetGeometry()->DefaultColors();
  }

  gGeoManager = AliGeomManager::GetGeometry();
  return gGeoManager;
}

void AliEveEventManager::SetAutoLoad(Bool_t autoLoad)
{
  // Set the automatic event loading mode
  fAutoLoad = autoLoad;
  StartStopAutoLoadTimer();
}

void AliEveEventManager::SetAutoLoadTime(Double_t time)
{
  // Set the auto-load time in seconds
  fAutoLoadTime = time;
  StartStopAutoLoadTimer();
}

void AliEveEventManager::StartStopAutoLoadTimer()
{
  // Create if needed and start
  // the automatic event loading timer
  if (fAutoLoad)
  {
    if (!fAutoLoadTimer)
    {
      fAutoLoadTimer = new TTimer;
      fAutoLoadTimer->Connect("Timeout()","AliEveEventManager",this,"NextEvent()");
    }
    fAutoLoadTimer->Start((Long_t)fAutoLoadTime*1000,kTRUE);
  }
  else
  {
    if (fAutoLoadTimer) fAutoLoadTimer->Stop();
  }
}

void AliEveEventManager::AfterNewEventLoaded()
{
  // Execute registered macros and commands.
  // At the end emit NewEventLoaded signal.
  //
  // Virtual from TEveEventManager.

  if (fExecutor)
    fExecutor->ExecMacros();

  TEveEventManager::AfterNewEventLoaded();

  NewEventLoaded();
}

void AliEveEventManager::NewEventLoaded()
{
  // Emit NewEventLoaded signal.

  Emit("NewEventLoaded()");
}

//==============================================================================

TString AliEveEventManager::GetEventInfoHorizontal() const
{
  // Dumps the event-header contents in vertical formatting.

  TString rawInfo, esdInfo;

  if (!fRawReader)
  {
    rawInfo = "No raw-data event info is available!\n";
  }
  else
  {
    const UInt_t* attr = fRawReader->GetAttributes();
    TTimeStamp ts(fRawReader->GetTimestamp());
    rawInfo.Form("RAW event info: Run#: %d  Event type: %d (%s)  Period: %x  Orbit: %x  BC: %x\n"
		 "Trigger: %llx\nDetectors: %x (%s)\nAttributes:%x-%x-%x  Timestamp: %s\n",
		 fRawReader->GetRunNumber(),fRawReader->GetType(),AliRawEventHeaderBase::GetTypeName(fRawReader->GetType()),
		 fRawReader->GetPeriod(),fRawReader->GetOrbitID(),fRawReader->GetBCID(),
		 fRawReader->GetClassMask(),
		 *fRawReader->GetDetectorPattern(),AliDAQ::ListOfTriggeredDetectors(*fRawReader->GetDetectorPattern()),
		 attr[0],attr[1],attr[2], ts.AsString("s"));
  }

  if (!fESD)
  {
    esdInfo = "No ESD event info is available!";
  }
  else
  {
    TString acttrclasses   = fESD->GetESDRun()->GetActiveTriggerClasses();
    TString firedtrclasses = fESD->GetFiredTriggerClasses();
    TTimeStamp ts(fESD->GetTimeStamp());
    esdInfo.Form("ESD event info: Run#: %d  Event type: %d (%s)  Period: %x  Orbit: %x  BC: %x\n"
		 "Active trigger classes: %s\nTrigger: %llx (%s)\nEvent# in file: %d  Timestamp: %s",
		 fESD->GetRunNumber(),
		 fESD->GetEventType(),AliRawEventHeaderBase::GetTypeName(fESD->GetEventType()),
		 fESD->GetPeriodNumber(),fESD->GetOrbitNumber(),fESD->GetBunchCrossNumber(),
		 acttrclasses.Data(),
		 fESD->GetTriggerMask(),firedtrclasses.Data(),
		 fESD->GetEventNumberInFile(), ts.AsString("s"));
  }

  return rawInfo + esdInfo;
}

TString AliEveEventManager::GetEventInfoVertical() const
{
  // Dumps the event-header contents in vertical formatting.

  TString rawInfo, esdInfo;

  if (!fRawReader)
  {
    rawInfo = "No raw-data event info is available!\n";
  }
  else
  {
    const UInt_t* attr = fRawReader->GetAttributes();
    rawInfo.Form("Raw-data event info:\nRun#: %d\nEvent type: %d (%s)\nPeriod: %x\nOrbit: %x   BC: %x\nTrigger: %llx\nDetectors: %x (%s)\nAttributes:%x-%x-%x\nTimestamp: %x\n",
		 fRawReader->GetRunNumber(),fRawReader->GetType(),AliRawEventHeaderBase::GetTypeName(fRawReader->GetType()),
		 fRawReader->GetPeriod(),fRawReader->GetOrbitID(),fRawReader->GetBCID(),
		 fRawReader->GetClassMask(),
		 *fRawReader->GetDetectorPattern(),AliDAQ::ListOfTriggeredDetectors(*fRawReader->GetDetectorPattern()),
		 attr[0],attr[1],attr[2],
		 fRawReader->GetTimestamp());
  }

  if (!fESD)
  {
    esdInfo = "No ESD event info is available!\n";
  }
  else
  {
    TString acttrclasses   = fESD->GetESDRun()->GetActiveTriggerClasses();
    TString firedtrclasses = fESD->GetFiredTriggerClasses();
    esdInfo.Form("ESD event info:\nRun#: %d\nActive trigger classes: %s\nEvent type: %d (%s)\nPeriod: %x\nOrbit: %x   BC: %x\nTrigger: %llx (%s)\nEvent# in file:%d\nTimestamp: %x\n",
		 fESD->GetRunNumber(),
		 acttrclasses.Data(),
		 fESD->GetEventType(),AliRawEventHeaderBase::GetTypeName(fESD->GetEventType()),
		 fESD->GetPeriodNumber(),fESD->GetOrbitNumber(),fESD->GetBunchCrossNumber(),
		 fESD->GetTriggerMask(),firedtrclasses.Data(),
		 fESD->GetEventNumberInFile(),
		 fESD->GetTimeStamp());
  }

  return rawInfo + "\n" + esdInfo;
}
  
