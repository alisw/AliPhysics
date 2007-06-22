// $Header$

#include "EventAlieve.h"
#include <Reve/Reve.h>
#include <Reve/RGTopFrame.h>

#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliESD.h>
#include <AliESDfriend.h>
#include <AliMagFMaps.h>
#include <AliCDBManager.h>
#include <AliHeader.h>
#include <AliGeomManager.h>

#include <TFile.h>
#include <TTree.h>
#include <TObjString.h>
#include <TError.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TCint.h>

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

Bool_t Alieve::Event::fgUseRunLoader   = kTRUE;
Bool_t Alieve::Event::fgUseESDTree     = kTRUE;
Bool_t Alieve::Event::fgAvoidExcOnOpen = kTRUE;

TString  Alieve::Event::fgCdbUri("local://$ALICE_ROOT");

AliMagF* Alieve::Event::fgMagField = 0;


void Event::Initialize(Bool_t use_runloader, Bool_t use_esd,
		       Bool_t avoid_exc_on_open)
{
  static const Exc_t eH("Event::Initialize ");

  fgUseRunLoader   = use_runloader;
  fgUseESDTree     = use_esd;
  fgAvoidExcOnOpen = avoid_exc_on_open;

  /*
  if(fgUseRunLoader == false && fgUseESDTree == false)
    throw(eH + "should use at least one data source.");

  if(fgUseRunLoader) {
    AssertMacro("loadlibs.C");
  }
  else if(fgUseESDTree) {
    gSystem->Load("libESD.so");
  }
  */
}

/**************************************************************************/

Event::Event() :
  EventBase(),

  fPath (), fEventId   (0),
  fRunLoader (0),
  fESDFile   (0), fESDTree (0), fESD (0),
  fESDfriend (0), fESDfriendExists(kFALSE),
  fNewEventCommands()
{}

Event::Event(TString path, Int_t ev) :
  EventBase("AliEVE Event"),

  fPath (path), fEventId(-1),
  fRunLoader (0),
  fESDFile   (0), fESDTree (0), fESD (0),
  fESDfriend (0), fESDfriendExists(kFALSE),
  fNewEventCommands()
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

  if(fgUseRunLoader)
  {
    TString ga_path(Form("%s/galice.root", fPath.Data()));
    if(gSystem->AccessPathName(ga_path, kReadPermission))
    {
      if (fgAvoidExcOnOpen) {
	Warning(eH, "RunLoader not initialized.");
	goto end_run_loader;
      } else {
	throw(eH + "can not read '" + ga_path + "'.");
      }
    }
    fRunLoader = AliRunLoader::Open(ga_path);
    if(!fRunLoader)
      throw(eH + "failed opening ALICE run loader from '" + ga_path + "'.");
    {
      TString alice_path = fPath + "/";
      fRunLoader->SetDirName(alice_path);
    }
    if(fRunLoader->LoadgAlice() != 0)
      throw(eH + "failed loading gAlice.");
    if(fRunLoader->LoadHeader() != 0)
      throw(eH + "failed loading header.");
    runNo = fRunLoader->GetHeader()->GetRun();
  }
end_run_loader:

  if(fgUseESDTree)
  {
    TString p(Form("%s/AliESDs.root", fPath.Data()));
    if(gSystem->AccessPathName(p, kReadPermission))
    {
      if (fgAvoidExcOnOpen) {
	Warning(eH, "ESD not initialized.");
	goto end_esd_loader;
      } else { 
	throw(eH + "can not read '" + p + "'.");
      }
    }
    fESDFile = new TFile(p);
    if(fESDFile->IsZombie()) {
      delete fESDFile; fESDFile = 0;
      throw(eH + "failed opening ALICE ESD from '" + p + "'.");
    }

    fESD = new AliESD();
    fESDTree = (TTree*) fESDFile->Get("esdTree");
    if(fESDTree == 0)
      throw(eH + "failed getting the esdTree.");
    fESD->ReadFromTree(fESDTree);
    runNo = fESD->GetESDRun()->GetRunNumber();

    // Check if ESDfriends exists and attach the branch
    p = Form("%s/AliESDfriends.root", fPath.Data());
    if(gSystem->AccessPathName(p, kReadPermission) == kFALSE)
    {
      fESDfriendExists = kTRUE;
      fESDTree->SetBranchStatus ("ESDfriend*", 1);
      fESDTree->SetBranchAddress("ESDfriend.", &fESDfriend);
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

end_esd_loader:
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

  RGTopFrame::RedrawDisabler rd(gReve);
  gReve->Redraw3D(kFALSE, kTRUE); // Enforce drop of all logicals.

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

void Event::AfterNewEventLoaded()
{
  TIter next(&fNewEventCommands);
  TObject* o;
  while ((o = next())) {
    TObjString* s = dynamic_cast<TObjString*>(o);
    if (s)
      gInterpreter->ProcessLine(s->String());
  }
}

void Event::AddNewEventCommand(const Text_t* cmd)
{
  fNewEventCommands.Add(new TObjString(cmd));
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

AliESD* Event::AssertESD()
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
