// $Header$

#include "EventAlieve.h"
#include <Reve/Reve.h>
#include <Reve/ReveManager.h>

#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliMagFMaps.h>
#include <AliCDBManager.h>
#include <AliHeader.h>
#include <AliGeomManager.h>

#include <TFile.h>
#include <TTree.h>
#include <TError.h>

#include <TROOT.h>
#include <TSystem.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// Event
//
// Provide interface for loading and navigating standard AliRoot data
// (AliRunLoader) and ESDs.
//
// Missing support for raw-data. For now this is handled individually
// by each sub-detector.
// 

ClassImp(Event)

Event* Alieve::gEvent = 0;

Bool_t Alieve::Event::fgAssertRunLoader = kFALSE;
Bool_t Alieve::Event::fgAssertESDTree   = kFALSE;

TString  Alieve::Event::fgCdbUri("local://$ALICE_ROOT");

AliMagF* Alieve::Event::fgMagField = 0;


Event::Event() :
  EventBase(),

  fPath (), fEventId   (0),
  fRunLoader (0),
  fESDFile   (0), fESDTree (0), fESD (0),
  fESDfriend (0), fESDfriendExists(kFALSE)
{}

Event::Event(TString path, Int_t ev) :
  EventBase("AliEVE Event"),

  fPath (path), fEventId(-1),
  fRunLoader (0),
  fESDFile   (0), fESDTree (0), fESD (0),
  fESDfriend (0), fESDfriendExists(kFALSE)
{
  Open();
  if (ev >= 0) GotoEvent(ev);
}

/**************************************************************************/

void Event::Open()
{
  static const Exc_t eH("Event::Open ");

  gSystem->ExpandPathName(fPath);
  if(fPath[0] != '/')
    fPath = Form("%s/%s", gSystem->WorkingDirectory(), fPath.Data());

  Int_t runNo = -1;

  TString ga_path(Form("%s/galice.root", fPath.Data()));
  if(gSystem->AccessPathName(ga_path, kReadPermission) == kFALSE)
  {
    fRunLoader = AliRunLoader::Open(ga_path);
    if (fRunLoader)
    {
      TString alice_path = fPath + "/";
      fRunLoader->SetDirName(alice_path);

      if (fRunLoader->LoadgAlice() != 0)
	Warning(eH, "failed loading gAlice via run-loader.");

      if (fRunLoader->LoadHeader() == 0)
      {
	runNo = fRunLoader->GetHeader()->GetRun();
      }
      else
      {
	Warning(eH, "failed loading run-loader's header.");
	delete fRunLoader;
	fRunLoader = 0;
      }
    }
    else // run-loader open failed
    {
      Warning(eH, "failed opening ALICE run-loader from '%s'.", ga_path.Data());
    }
  }
  else // galice not readable
  {
    Warning(eH, "can not read '%s'.", ga_path.Data());
  }
  if (fRunLoader == 0)
  {
    if(fgAssertRunLoader)
      throw(eH + "Bootstraping of run-loader failed. Its precence was requested.");
    else
      Warning(eH, "Bootstraping of run-loader failed.");
  }
  

  TString esd_path(Form("%s/AliESDs.root", fPath.Data()));
  if(gSystem->AccessPathName(esd_path, kReadPermission) == kFALSE)
  {
    fESDFile = new TFile(esd_path);
    if(fESDFile->IsZombie() == kFALSE)
    {
      fESD = new AliESDEvent();
      fESDTree = (TTree*) fESDFile->Get("esdTree");
      if (fESDTree != 0)
      {
	fESD->ReadFromTree(fESDTree);
	runNo = fESD->GetESDRun()->GetRunNumber();

	// Check if ESDfriends exists and attach the branch
	TString p = Form("%s/AliESDfriends.root", fPath.Data());
	if(gSystem->AccessPathName(p, kReadPermission) == kFALSE)
	{
	  fESDfriendExists = kTRUE;
	  fESDTree->SetBranchStatus ("ESDfriend*", 1);
	  fESDTree->SetBranchAddress("ESDfriend.", &fESDfriend);
	}
      }
      else // esdtree == 0
      {
	delete fESDFile; fESDFile = 0;
	Warning(eH, "failed getting the esdTree.");
      }
    }
    else // esd tfile is zombie
    {
      delete fESDFile; fESDFile = 0;
      Warning(eH, "failed opening ESD from '%s'.", esd_path.Data());
    }
  }
  else // esd not readable
  {
    Warning(eH, "can not read ESD file '%s'.", esd_path.Data());
  }
  if (fESDTree == 0)
  {
    if (fgAssertESDTree)
    {
      throw(eH + "ESD not initialized. Its precence was requested.");
    } else { 
      Warning(eH, "ESD not initialized.");
    }
  }

  if (runNo < 0)
    throw(eH + "invalid run number.");

  {
    AliCDBManager* cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage(fgCdbUri);
    if (cdb->IsDefaultStorageSet() == kFALSE)
      throw(eH + "CDB initialization failed.");
    cdb->SetRun(runNo);
  }

  SetName(Form("Event %d", fEventId));
  SetTitle(fPath);
}

void Event::GotoEvent(Int_t event)
{
  static const Exc_t eH("Event::GotoEvent ");

  Int_t maxEvent = 0;
  if(fRunLoader)
    maxEvent = fRunLoader->GetNumberOfEvents() - 1;
  else if(fESDTree)
    maxEvent = fESDTree->GetEntries() - 1;
  else
    throw(eH + "neither RunLoader nor ESD loaded.");

  if(event < 0 || event > maxEvent)
    throw(eH + Form("event %d not present, available range [%d, %d].",
		    event, 0, maxEvent));

  ReveManager::RedrawDisabler rd(gReve);
  gReve->Redraw3D(kFALSE, kTRUE); // Enforce drop of all logicals.

  // !!! MT this is somewhat brutal; at least optionally, one could be
  // a bit gentler, checking for objs owning their external refs and having
  // additinal parents.
  DestroyElements();
  fEventId = event;
  SetName(Form("Event %d", fEventId));
  UpdateItems();

  if(fRunLoader) {
    if(fRunLoader->GetEvent(fEventId) != 0)
      throw(eH + "failed getting required event.");
  }

  if(fESDTree) {
    if(fESDTree->GetEntry(fEventId) <= 0)
      throw(eH + "failed getting required event from ESD.");

    if (fESDfriendExists)
      fESD->SetESDfriend(fESDfriend);
  }

  AfterNewEventLoaded();
}

void Event::Close()
{
  if (fESDTree) {
    delete fESD;       fESD       = 0;
    delete fESDfriend; fESDfriend = 0;

    delete fESDTree; fESDTree = 0;
    delete fESDFile; fESDFile = 0;
  }
}


/**************************************************************************/
/**************************************************************************/

// Static convenience functions.

AliRunLoader* Event::AssertRunLoader()
{
  static const Exc_t eH("Event::AssertRunLoader ");

  if(gEvent == 0)
    throw(eH + "ALICE event not ready.");
  if(gEvent->fRunLoader == 0)
    throw(eH + "AliRunLoader not initialised.");
  return gEvent->fRunLoader;
}

AliESDEvent* Event::AssertESD()
{
  static const Exc_t eH("Event::AssertESD ");

  if(gEvent == 0)
    throw(eH + "ALICE event not ready.");
  if(gEvent->fESD == 0)
    throw(eH + "AliESD not initialised.");
  return gEvent->fESD;
}

AliESDfriend* Event::AssertESDfriend()
{
  static const Exc_t eH("Event::AssertESDfriend ");

  if(gEvent == 0)
    throw(eH + "ALICE event not ready.");
  if(gEvent->fESDfriend == 0)
    throw(eH + "AliESDfriend not initialised.");
  return gEvent->fESDfriend;
}

AliMagF* Event::AssertMagField()
{
  if (fgMagField == 0)
  {
    if (gEvent && gEvent->fRunLoader && gEvent->fRunLoader->GetAliRun())
      fgMagField = gEvent->fRunLoader->GetAliRun()->Field();
    else
      fgMagField = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
  }
  return fgMagField;
}

TGeoManager* Event::AssertGeometry()
{
  static const Exc_t eH("Event::AssertGeometry ");

  if (AliGeomManager::GetGeometry() == 0)
  {
    gGeoManager = 0;
    AliGeomManager::LoadGeometry();
    if ( ! AliGeomManager::GetGeometry())
    {
      throw(eH + "can not load geometry.");
    }
    if ( ! AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO ACORDE"))
    {
      ::Warning(eH, "mismatch of alignable volumes. Proceeding.");
      // throw(eH + "could not apply align objs.");
    }
  }

  return AliGeomManager::GetGeometry();
}
