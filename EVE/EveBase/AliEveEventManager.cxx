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
#include <AliESDRun.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliAODEvent.h>

#include <AliDAQ.h>
#include <AliRawEventHeaderBase.h>
#include <AliRawReaderRoot.h>
#include <AliRawReaderFile.h>
#include <AliRawReaderDate.h>
#include <AliMagF.h>
#include <AliCDBManager.h>
#include <AliCDBStorage.h>
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
// Provides interface for loading and navigating standard AliRoot data
// (AliRunLoader), ESD, AOD and RAW.
//
// ESDfriend is attached automatically, if the file is found.
//
// AODfriends are not attached automatically as there are several
// possible files involved. To have a specific AODfriend attached, call
// static method
//   AliEveEventManager::AddAODfriend("AliAOD.VertexingHF.root");
// before initializing the event-manager.
//
// Also provides interface to magnetic-field and geometry. Mostly
// intended as wrappers over standard AliRoot functionality for
// convenient use from visualizateion macros.
//
// There can be a single main event-manger, it is stored in private
// data member fgMaster and can be accessed via static member function
// GetMaster().
//
// For event overlaying and embedding one can instantiate additional
// event-managers via static method AddDependentManager(const TString& path).
// This interface is under development.

ClassImp(AliEveEventManager)

Bool_t AliEveEventManager::fgAssertRunLoader = kFALSE;
Bool_t AliEveEventManager::fgAssertESD       = kFALSE;
Bool_t AliEveEventManager::fgAssertAOD       = kFALSE;
Bool_t AliEveEventManager::fgAssertRaw       = kFALSE;

TString  AliEveEventManager::fgESDFileName("AliESDs.root");
TString  AliEveEventManager::fgAODFileName("AliAOD.root");
TString  AliEveEventManager::fgRawFileName("raw.root");
TString  AliEveEventManager::fgCdbUri("local://$ALICE_ROOT/OCDB");

TList*   AliEveEventManager::fgAODfriends = 0;

AliEveEventManager* AliEveEventManager::fgMaster  = 0;
AliEveEventManager* AliEveEventManager::fgCurrent = 0;

void AliEveEventManager::InitInternals()
{
  // Initialize internal members.

  static const TEveException kEH("AliEveEventManager::InitInternals ");

  if (fgCurrent != 0)
  {
    throw(kEH + "Dependent event-managers should be created via static method AddDependentManager().");
  }

  if (fgMaster == 0)
  {
    fgMaster = this;
  }

  fgCurrent = this;

  fAutoLoadTimer = new TTimer;
  fAutoLoadTimer->Connect("Timeout()", "AliEveEventManager", this, "AutoLoadNextEvent()");

  fExecutor = new AliEveMacroExecutor;
}

AliEveEventManager::AliEveEventManager(const TString& name) :
  TEveEventManager(name),

  fPath      ( ), fEventId (-1),
  fRunLoader (0),
  fESDFile   (0), fESDTree (0), fESD (0),
  fESDfriend (0), fESDfriendExists(kFALSE),
  fAODFile   (0), fAODTree (0), fAOD (0),
  fRawReader (0),
  fAutoLoad  (kFALSE), fAutoLoadTime (5.),     fAutoLoadTimer(0),
  fIsOpen    (kFALSE), fHasEvent     (kFALSE), fExternalCtrl (kFALSE),
  fSelectOnTriggerType(kFALSE), fTriggerType(""),
  fExecutor    (0),
  fSubManagers (0),
  fAutoLoadTimerRunning(kFALSE)
{
  // Default constructor.

  InitInternals();
}

AliEveEventManager::AliEveEventManager(const TString& name, const TString& path, Int_t ev) :
  TEveEventManager(name, path),

  fPath   (path), fEventId(-1),
  fRunLoader (0),
  fESDFile   (0), fESDTree (0), fESD (0),
  fESDfriend (0), fESDfriendExists(kFALSE),
  fAODFile   (0), fAODTree (0), fAOD (0),
  fRawReader (0),
  fAutoLoad  (kFALSE), fAutoLoadTime (5),      fAutoLoadTimer(0),
  fIsOpen    (kFALSE), fHasEvent     (kFALSE), fExternalCtrl (kFALSE),
  fSelectOnTriggerType(kFALSE), fTriggerType(""),
  fExecutor    (0),
  fSubManagers (0),
  fAutoLoadTimerRunning(kFALSE)
{
  // Constructor with event-directory URL and event-id.

  InitInternals();

  Open();
  if (ev >= 0)
  {
    GotoEvent(ev);
  }
}

AliEveEventManager::~AliEveEventManager()
{
  // Destructor.

  delete fSubManagers;

  if (fIsOpen)
  {
    Close();
  }
}

/******************************************************************************/

void AliEveEventManager::SetESDFileName(const TString& esd)
{
  // Set file-name for opening ESD, default "AliESDs.root".

  if ( ! esd.IsNull()) fgESDFileName = esd;
}

void AliEveEventManager::SetAODFileName(const TString& aod)
{
  // Set file-name for opening AOD, default "AliAOD.root".

  if ( ! aod.IsNull()) fgAODFileName = aod;
}

void AliEveEventManager::AddAODfriend(const TString& friendFileName)
{
  // Add new AOD friend file-name to be attached when opening AOD.
  // This should include '.root', as in 'AliAOD.VertexingHF.root'.

  if (fgAODfriends == 0)
  {
    fgAODfriends = new TList;
    fgAODfriends->SetOwner(kTRUE);
  }
  fgAODfriends->Add(new TObjString(friendFileName));
}

void AliEveEventManager::SetRawFileName(const TString& raw)
{
  // Set file-name for opening of raw-data, default "raw.root"
  if ( ! raw.IsNull()) fgRawFileName = raw;
}

void AliEveEventManager::SetCdbUri(const TString& cdb)
{
  // Set path to CDB, default "local://$ALICE_ROOT/OCDB".

  if ( ! cdb.IsNull()) fgCdbUri = cdb;
}

void AliEveEventManager::SetAssertElements(Bool_t assertRunloader, Bool_t assertEsd,
					   Bool_t assertAod, Bool_t assertRaw)
{
  // Set global flags that detrmine which parts of the event-data must
  // be present when the event is opened.

  fgAssertRunLoader = assertRunloader;
  fgAssertESD = assertEsd;
  fgAssertAOD = assertAod;
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

  if (fExternalCtrl)
  {
    throw (kEH + "Event-loop is under external control.");
  }
  if (fIsOpen)
  {
    throw (kEH + "Event-files already opened.");
  }

  gSystem->ExpandPathName(fPath);
  // The following magick is required for ESDfriends to be loaded properly
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
      // Check if ESDfriends exists and attach the branch.
      // We use TFile::Open() instead of gSystem->AccessPathName
      // as it seems to work better when attachine alieve to a
      // running reconstruction process with auto-save on.
      // There was also a problem with TTree::Refresh() - it didn't
      // save the friend branch on a separate file, fixed in 5.22.2 -
      // so we might want to try the old way again soon.
      TString p(Form("%s/AliESDfriends.root", fPath.Data()));
      TFile *esdFriendFile = TFile::Open(p);
      if (esdFriendFile)
      {
	if (!esdFriendFile->IsZombie())
	{
	  esdFriendFile->Close();
	  fESDfriendExists = kTRUE;
	  fESDTree->SetBranchStatus ("ESDfriend*", 1);
	}
	delete esdFriendFile;
      }

      fESD->ReadFromTree(fESDTree);
      if (fESDfriendExists)
      {
	fESDfriend = (AliESDfriend*) fESD->FindListObject("AliESDfriend");
	Info(kEH, "found and attached ESD friend.");
      }
      else
      {
	Warning(kEH, "ESDfriend not found.");
      }

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

  // Open AOD and registered friends

  TString aodPath(Form("%s/%s", fPath.Data(), fgAODFileName.Data()));
  if ((fAODFile = TFile::Open(aodPath)))
  {
    fAOD = new AliAODEvent();
    fAODTree = (TTree*) fAODFile->Get("aodTree");
    if (fAODTree != 0)
    {
      // Check if AODfriends exist and attach them.
      TIter       friends(fgAODfriends);
      TObjString *name;
      while ((name = (TObjString*) friends()) != 0)
      {
	TString p(Form("%s/%s", fPath.Data(), name->GetName()));
	if (gSystem->AccessPathName(p, kReadPermission) == kFALSE)
	{
	  fAODTree->AddFriend("aodTree", name->GetName());
	}
      }

      fAOD->ReadFromTree(fAODTree);

      if (fAODTree->GetEntry(0) <= 0)
      {
	delete fAODFile; fAODFile = 0;
	delete fAOD;     fAOD     = 0;
	Warning(kEH, "failed getting the first entry from addTree.");
      }
      else
      {
	if (runNo < 0)
	  runNo = fAOD->GetRunNumber();
      }
    }
    else // aodtree == 0
    {
      delete fAODFile; fAODFile = 0;
      delete fAOD;     fAOD     = 0;
      Warning(kEH, "failed getting the aodTree.");
    }
  }
  else // aod not readable
  {
    Warning(kEH, "can not read AOD file '%s'.", aodPath.Data());
  }
  if (fAODTree == 0)
  {
    if (fgAssertAOD)
    {
      throw (kEH + "AOD not initialized. Its precence was requested.");
    } else {
      Warning(kEH, "AOD not initialized.");
    }
  }

  // Open RunLoader from galice.root

  TString gaPath(Form("%s/galice.root", fPath.Data()));
  // If i use open directly, we get fatal.
  // Is AccessPathName check ok for xrootd / alien? Yes, not for http.
  if (gSystem->AccessPathName(gaPath, kReadPermission) == kFALSE)
  {
    fRunLoader = AliRunLoader::Open(gaPath, GetName());
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
  AliLog::EType_t oldLogLevel = (AliLog::EType_t) AliLog::GetGlobalLogLevel();
  if (fgAssertRaw == kFALSE)
  {
    AliLog::SetGlobalLogLevel(AliLog::kFatal);
  }
  if (gSystem->AccessPathName(rawPath, kReadPermission) == kFALSE)
  {
    fRawReader = AliRawReader::Create(rawPath);
  }
  else
  {
    fRawReader = AliRawReader::Create(fgRawFileName);
  }
  if (fgAssertRaw == kFALSE)
  {
    AliLog::SetGlobalLogLevel(oldLogLevel);
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

  // Initialize OCDB ... only in master event-manager

  if (this == fgMaster)
  {
    AliCDBManager* cdb = AliCDBManager::Instance();
    if (cdb->IsDefaultStorageSet() == kTRUE)
    {
      Warning(kEH, "CDB already set - using the old storage:\n  '%s'",
	      cdb->GetDefaultStorage()->GetURI().Data());
    }
    else
    {
      cdb->SetDefaultStorage(fgCdbUri);
      if (cdb->IsDefaultStorageSet() == kFALSE)
	throw (kEH + "CDB initialization failed.");
    }
    cdb->SetRun(runNo);
  }

  fIsOpen = kTRUE;
}

void AliEveEventManager::SetEvent(AliRunLoader *runLoader, AliRawReader *rawReader, AliESDEvent *esd, AliESDfriend *esdf)
{
  // Set an event from an external source.
  // The method is used in the online visualisation.
  // AOD is not supported.

  static const TEveException kEH("AliEveEventManager::SetEvent ");

  if (fIsOpen)
  {
    Warning(kEH, "Event-files were open. Closing and switching to external control.");
    Close();
  }

  fRunLoader = runLoader;
  fRawReader = rawReader;
  fESD       = esd;
  fESDfriend = esdf;
  fAOD       = 0;

  fEventId++;
  fHasEvent     = kTRUE;
  fExternalCtrl = kTRUE;

  SetTitle("Online event in memory");
  SetName ("Online Event");
  ElementChanged();

  AfterNewEventLoaded();

  if (fAutoLoad) StartAutoLoadTimer();
}

Int_t AliEveEventManager::GetMaxEventId(Bool_t /*refreshESD*/) const
{
  // Returns maximum available event id.
  // If under external control or event is not opened -1 is returned.
  // If raw-data is the only data-source this can not be known
  // and 10,000,000 is returned.
  // If neither data-source is initialised an exception is thrown.
  // If refresh_esd is true and ESD is the primary event-data source
  // its header is re-read from disk.

  static const TEveException kEH("AliEveEventManager::GetMaxEventId ");

  if (fExternalCtrl || fIsOpen == kFALSE)
  {
    return -1;
  }

  if (fESDTree)
  {
    // Refresh crashes with root-5.21.1-alice.
    // Fixed by Philippe 5.8.2008 r25053, can be reactivated
    // when we move to a newer root.
    // if (refreshESD)
    //   fESDTree->Refresh();
    return fESDTree->GetEntries() - 1;
  }
  else if (fAODTree)
  {
    return fAODTree->GetEntries() - 1;
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
    throw (kEH + "neither ESD, AOD, RunLoader nor Raw loaded.");
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

  if (fAutoLoadTimerRunning)
  {
    throw (kEH + "Event auto-load timer is running.");
  }
  if (fExternalCtrl)
  {
    throw (kEH + "Event-loop is under external control.");
  }
  else if (!fIsOpen)
  {
    throw (kEH + "Event-files not opened.");
  }

  fHasEvent = kFALSE;

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
  else if (fAODTree)
  {
    maxEvent = fAODTree->GetEntries() - 1;
    if (event < 0)
      event = fAODTree->GetEntries() + event;
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

  if (fAODTree) {
    if (fAODTree->GetEntry(event) <= 0)
      throw (kEH + "failed getting required event from AOD.");
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

  fHasEvent = kTRUE;
  fEventId  = event;
  if (this == fgMaster)
  {
    SetName(Form("Event %d", fEventId));
    ElementChanged();
  }

  AfterNewEventLoaded();
}

void AliEveEventManager::NextEvent()
{
  // Loads next event.
  // Does magick needed for online display when under external event control.

  static const TEveException kEH("AliEveEventManager::NextEvent ");

  if (fAutoLoadTimerRunning)
  {
    throw (kEH + "Event auto-load timer is running.");
  }

  if (fExternalCtrl)
  {
    // !!! This should really go somewhere else. It is done in GotoEvent(),
    // so here we should do it in SetEvent().
    DestroyElements();

    gSystem->ExitLoop();
  }
  else
  {
    Int_t nexteventbytrigger=0;
    if (fSelectOnTriggerType)
    {
      if (FindNextByTrigger(nexteventbytrigger)) //if not found do nothing
        GotoEvent(nexteventbytrigger);
    }
    else if (fEventId < GetMaxEventId(kTRUE))
      GotoEvent(fEventId + 1);
    else
      GotoEvent(0);
  }
}

void AliEveEventManager::PrevEvent()
{
  // Loads previous event.

  static const TEveException kEH("AliEveEventManager::PrevEvent ");

  if (fAutoLoadTimerRunning)
  {
    throw (kEH + "Event auto-load timer is running.");
  }
  if (fExternalCtrl)
  {
    throw (kEH + "Event-loop is under external control.");
  }
  Int_t nexteventbytrigger=0;
  if (fSelectOnTriggerType)
  {
    if (FindPrevByTrigger(nexteventbytrigger))
      GotoEvent(nexteventbytrigger);
  }
  else
    GotoEvent(fEventId - 1);
}

void AliEveEventManager::Close()
{
  // Close the event data-files and delete ESD, ESDfriend, run-loader
  // and raw-reader.

  static const TEveException kEH("AliEveEventManager::Close ");

  if (!fIsOpen)
  {
    throw (kEH + "Event-files not opened.");
  }

  if (fAutoLoadTimerRunning)
    StopAutoLoadTimer();

  if (fESDTree) {
    delete fESD;       fESD       = 0;
    delete fESDfriend; fESDfriend = 0;
    fESDfriendExists = kFALSE;

    delete fESDTree;   fESDTree = 0;
    delete fESDFile;   fESDFile = 0;
  }

  if (fAODTree) {
    delete fAOD;       fAOD       = 0;

    delete fAODTree;   fAODTree = 0;
    delete fAODFile;   fAODFile = 0;
  }

  if (fRunLoader) {
    delete fRunLoader; fRunLoader = 0;
  }

  if (fRawReader) {
    delete fRawReader; fRawReader = 0;
  }

  fEventId  = -1;
  fIsOpen   = kFALSE;
  fHasEvent = kFALSE;
}


//------------------------------------------------------------------------------
// Static convenience functions, mainly used from macros.
//------------------------------------------------------------------------------

Bool_t AliEveEventManager::HasRunLoader()
{
  // Check if AliRunLoader is initialized.

  return fgCurrent && fgCurrent->fHasEvent && fgCurrent->fRunLoader;
}

Bool_t AliEveEventManager::HasESD()
{
  // Check if AliESDEvent is initialized.

  return fgCurrent && fgCurrent->fHasEvent && fgCurrent->fESD;
}

Bool_t AliEveEventManager::HasESDfriend()
{
  // Check if AliESDfriend is initialized.

  return fgCurrent && fgCurrent->fHasEvent && fgCurrent->fESDfriend;
}

Bool_t AliEveEventManager::HasAOD()
{
  // Check if AliESDEvent is initialized.

  return fgCurrent && fgCurrent->fHasEvent && fgCurrent->fAOD;
}

Bool_t AliEveEventManager::HasRawReader()
{
  // Check if raw-reader is initialized.

  return fgCurrent && fgCurrent->fHasEvent && fgCurrent->fRawReader;
}

AliRunLoader* AliEveEventManager::AssertRunLoader()
{
  // Make sure AliRunLoader is initialized and return it.
  // Throws exception in case run-loader is not available.
  // Static utility for macros.

  static const TEveException kEH("AliEveEventManager::AssertRunLoader ");

  if (fgCurrent == 0 || fgCurrent->fHasEvent == kFALSE)
    throw (kEH + "ALICE event not ready.");
  if (fgCurrent->fRunLoader == 0)
    throw (kEH + "AliRunLoader not initialised.");
  return fgCurrent->fRunLoader;
}

AliESDEvent* AliEveEventManager::AssertESD()
{
  // Make sure AliESDEvent is initialized and return it.
  // Throws exception in case ESD is not available.
  // Static utility for macros.

  static const TEveException kEH("AliEveEventManager::AssertESD ");

  if (fgCurrent == 0 || fgCurrent->fHasEvent == kFALSE)
    throw (kEH + "ALICE event not ready.");
  if (fgCurrent->fESD == 0)
    throw (kEH + "AliESD not initialised.");
  return fgCurrent->fESD;
}

AliESDfriend* AliEveEventManager::AssertESDfriend()
{
  // Make sure AliESDfriend is initialized and return it.
  // Throws exception in case ESDfriend-loader is not available.
  // Static utility for macros.

  static const TEveException kEH("AliEveEventManager::AssertESDfriend ");

  if (fgCurrent == 0 || fgCurrent->fHasEvent == kFALSE)
    throw (kEH + "ALICE event not ready.");
  if (fgCurrent->fESDfriend == 0)
    throw (kEH + "AliESDfriend not initialised.");
  return fgCurrent->fESDfriend;
}

AliAODEvent* AliEveEventManager::AssertAOD()
{
  // Make sure AliAODEvent is initialized and return it.
  // Throws exception in case AOD is not available.
  // Static utility for macros.

  static const TEveException kEH("AliEveEventManager::AssertAOD ");

  if (fgCurrent == 0 || fgCurrent->fHasEvent == kFALSE)
    throw (kEH + "ALICE event not ready.");
  if (fgCurrent->fAOD == 0)
    throw (kEH + "AliAOD not initialised.");
  return fgCurrent->fAOD;
}

AliRawReader* AliEveEventManager::AssertRawReader()
{
  // Make sure raw-reader is initialized and return it.

  static const TEveException kEH("AliEveEventManager::AssertRawReader ");

  if (fgCurrent == 0 || fgCurrent->fHasEvent == kFALSE)
    throw (kEH + "ALICE event not ready.");
  if (fgCurrent->fRawReader == 0)
    throw (kEH + "RawReader not ready.");

  return fgCurrent->fRawReader;
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

//------------------------------------------------------------------------------

AliEveEventManager* AliEveEventManager::AddDependentManager(const TString& name, const TString& path)
{
  // Create and attach a dependent event-manager.
  // It is not added into eve list tree.

  static const TEveException kEH("AliEveEventManager::AddDependentManager ");

  if (fgMaster == 0)
    throw(kEH + "Master event-manager must be instantiated first.");

  if (fgMaster->fSubManagers == 0)
  {
    fgMaster->fSubManagers = new TList;
    fgMaster->fSubManagers->SetOwner(kTRUE);
  }

  AliEveEventManager* new_mgr = 0;
  fgCurrent = 0;
  try
  {
    new_mgr = new AliEveEventManager(name, path, fgMaster->fEventId);
    fgMaster->fSubManagers->Add(new_mgr);
  }
  catch (TEveException& exc)
  {
    ::Error(kEH, "Creation of new event-manager failed: '%s'.", exc.Data());
  }
  fgCurrent = fgMaster;

  return new_mgr;
}

AliEveEventManager* AliEveEventManager::GetDependentManager(const TString& name)
{
  // Get a dependant manager by name.
  // This will not change the current manager, use helper class
  // AliEveEventManager::CurrentChanger for that.

  static const TEveException kEH("AliEveEventManager::GetDependentManager ");

  if (fgMaster == 0)
    throw(kEH + "Master event-manager must be instantiated first.");

  if (fgMaster->fSubManagers == 0)
    return 0;

  return dynamic_cast<AliEveEventManager*>(fgMaster->fSubManagers->FindObject(name));
}

AliEveEventManager* AliEveEventManager::GetMaster()
{
  // Get master event-manager.

  return fgMaster;
}

AliEveEventManager* AliEveEventManager::GetCurrent()
{
  // Get current event-manager.

  return fgCurrent;
}

//------------------------------------------------------------------------------
// Autoloading of events
//------------------------------------------------------------------------------

void AliEveEventManager::SetAutoLoadTime(Float_t time)
{
  // Set the auto-load time in seconds

  fAutoLoadTime = time;
}

void AliEveEventManager::SetAutoLoad(Bool_t autoLoad)
{
  // Set the automatic event loading mode

  static const TEveException kEH("AliEveEventManager::SetAutoLoad ");

  if (fAutoLoad == autoLoad)
  {
    Warning(kEH, "Setting autoload to the same value as before - %s. Ignoring.", fAutoLoad ? "true" : "false");
    return;
  }

  fAutoLoad = autoLoad;
  if (fAutoLoad)
  {
    StartAutoLoadTimer();
  }
  else
  {
    StopAutoLoadTimer();
  }
}

void AliEveEventManager::StartAutoLoadTimer()
{
  // Start the auto-load timer.

  fAutoLoadTimer->SetTime((Long_t)(1000*fAutoLoadTime));
  fAutoLoadTimer->Reset();
  fAutoLoadTimer->TurnOn();
  fAutoLoadTimerRunning = kTRUE;
}

void AliEveEventManager::StopAutoLoadTimer()
{
  // Stop the auto-load timer.

  fAutoLoadTimerRunning = kFALSE;
  fAutoLoadTimer->TurnOff();
}

void AliEveEventManager::AutoLoadNextEvent()
{
  // Called from auto-load timer, so it has to be public.
  // Do NOT call it directly.

  static const TEveException kEH("AliEveEventManager::AutoLoadNextEvent ");

  if ( ! fAutoLoadTimerRunning || ! fAutoLoadTimer->HasTimedOut())
  {
    Warning(kEH, "Called unexpectedly - ignoring the call. Should ONLY be called from an internal timer.");
    return;
  }

  StopAutoLoadTimer();
  NextEvent();
  if (fAutoLoad && !fExternalCtrl)
    StartAutoLoadTimer();
}


//------------------------------------------------------------------------------
// Event selection by trigger
//------------------------------------------------------------------------------

Bool_t AliEveEventManager::FindNextByTrigger(Int_t& event)
{
  // Find next event that matches the trigger.
  // If a matching event is not found, we loop around and eventually
  // end up at the same event.

  static const TEveException kEH("AliEveEventManager::FindNextByTrigger ");

  if (!fESDTree) return kFALSE;
  TString firedtrclasses;
  for (Int_t i = fEventId+1; i<GetMaxEventId(kTRUE)+1; i++)
  {
    if (fESDTree->GetEntry(i) <= 0)
      throw (kEH + "failed getting required event from ESD.");
    firedtrclasses = fESD->GetFiredTriggerClasses();
    if (firedtrclasses.Contains(fTriggerType))
    {
      event=i;
      return kTRUE;
    }
  }
  for (Int_t i = 0; i<fEventId+1; i++)
  {
    if (fESDTree->GetEntry(i) <= 0)
      throw (kEH + "failed getting required event from ESD.");
    firedtrclasses = fESD->GetFiredTriggerClasses();
    if (firedtrclasses.Contains(fTriggerType))
    {
      event=i;
      return kTRUE;
    }
  }
  return kFALSE;
}

Bool_t AliEveEventManager::FindPrevByTrigger(Int_t& event)
{
  // Find previous event that matches the trigger.

  static const TEveException kEH("AliEveEventManager::FindPrevByTrigger ");

  if (!fESDTree) return kFALSE;
  TString firedtrclasses;
  for (Int_t i = fEventId-1; i>=0; i--)
  {
    if (fESDTree->GetEntry(i) <= 0)
      throw (kEH + "failed getting required event from ESD.");
    firedtrclasses = fESD->GetFiredTriggerClasses();
    if (firedtrclasses.Contains(fTriggerType))
    {
      event=i;
      return kTRUE;
    }
  }
  for (Int_t i = GetMaxEventId(kTRUE); i>fEventId-1; i--)
  {
    if (fESDTree->GetEntry(i) <= 0)
      throw (kEH + "failed getting required event from ESD.");
    firedtrclasses = fESD->GetFiredTriggerClasses();
    if (firedtrclasses.Contains(fTriggerType))
    {
      event=i;
      return kTRUE;
    }
  }
  return kFALSE;
}


//------------------------------------------------------------------------------
// Post event-loading functions
//------------------------------------------------------------------------------

void AliEveEventManager::AfterNewEventLoaded()
{
  // Execute registered macros and commands.
  // At the end emit NewEventLoaded signal.
  //
  // Virtual from TEveEventManager.

  static const TEveException kEH("AliEveEventManager::AfterNewEventLoaded ");

  if (fExecutor)
    fExecutor->ExecMacros();

  TEveEventManager::AfterNewEventLoaded();

  NewEventLoaded();

  if (this == fgMaster && fSubManagers != 0)
  {
    TIter next(fSubManagers);
    while ((fgCurrent = dynamic_cast<AliEveEventManager*>(next())) != 0)
    {
      gEve->SetCurrentEvent(fgCurrent);
      try
      {
	fgCurrent->GotoEvent(fEventId);
      }
      catch (TEveException& exc)
      {
	// !!! Should somehow tag / disable / remove it?
	Error(kEH, "Getting event %d for sub-event-manager '%s' failed: '%s'.",
	      fEventId, fgCurrent->GetName(), exc.Data());
      }
    }
    fgCurrent = fgMaster;
    gEve->SetCurrentEvent(fgMaster);
  }
}

void AliEveEventManager::NewEventLoaded()
{
  // Emit NewEventLoaded signal.

  Emit("NewEventLoaded()");
}


//------------------------------------------------------------------------------
// Event info dumpers
//------------------------------------------------------------------------------

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
		 "Active trigger classes: %s\nTrigger: %llx (%s)\nEvent# in file: %d  Timestamp: %s, MagField: %.2e",
		 fESD->GetRunNumber(),
		 fESD->GetEventType(),AliRawEventHeaderBase::GetTypeName(fESD->GetEventType()),
		 fESD->GetPeriodNumber(),fESD->GetOrbitNumber(),fESD->GetBunchCrossNumber(),
		 acttrclasses.Data(),
		 fESD->GetTriggerMask(),firedtrclasses.Data(),
		 fESD->GetEventNumberInFile(), ts.AsString("s"), fESD->GetMagneticField());
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
