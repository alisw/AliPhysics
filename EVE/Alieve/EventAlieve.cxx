// $Header$

#include "EventAlieve.h"
#include <Reve/Reve.h>

#include <AliRunLoader.h>
#include <AliESD.h>
#include <AliESDfriend.h>

#include <TFile.h>
#include <TTree.h>

#include <TROOT.h>
#include <TSystem.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// Event
//

ClassImp(Event)

Event* Alieve::gEvent = 0;

Bool_t Alieve::Event::fgUseRunLoader   = kTRUE;
Bool_t Alieve::Event::fgUseESDTree     = kTRUE;
Bool_t Alieve::Event::fgAvoidExcOnOpen = kTRUE;

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
  fESDFile       (0), fESDTree       (0), fESD       (0),
  /* fESDfriendFile (0), fESDfriendTree (0), */
  fESDfriend (0), fESDfriendExists(kFALSE)
{}

Event::Event(TString path, Int_t ev) :
  EventBase("AliEVE Event"),

  fPath (path), fEventId(ev),
  fRunLoader (0),
  fESDFile       (0), fESDTree       (0), fESD       (0),
  /* fESDfriendFile (0), fESDfriendTree (0), */
  fESDfriend (0), fESDfriendExists(kFALSE)
{
  Open();
}

/**************************************************************************/

void Event::Open()
{
  static const Exc_t eH("Event::Open ");

  gSystem->ExpandPathName(fPath);
  if(fPath[0] != '/')
    fPath = Form("%s/%s", gSystem->WorkingDirectory(), fPath.Data());

  if(fgUseRunLoader) {
    TString ga_path(Form("%s/galice.root", fPath.Data()));
    if(gSystem->AccessPathName(ga_path, kReadPermission)) {
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

    if(fRunLoader->GetEvent(fEventId) != 0)
      throw(eH + "failed getting required event.");
  }
end_run_loader:

  if(fgUseESDTree) {
    TString p(Form("%s/AliESDs.root", fPath.Data()));
    if(gSystem->AccessPathName(p, kReadPermission)) {
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

    fESDTree = (TTree*) fESDFile->Get("esdTree");
    if(fESDTree == 0)
      throw(eH + "failed getting the esdTree.");
    fESDTree->SetBranchAddress("ESD", &fESD);

    // Check if ESDfriends exists and attach the branch
    p = Form("%s/AliESDfriends.root", fPath.Data());
    if(gSystem->AccessPathName(p, kReadPermission) == kFALSE) {
      //fESDfriendFile = new TFile(p);
      //if(fESDfriendFile->IsZombie()) {
      //delete fESDfriendFile; fESDfriendFile = 0;
      //throw(eH + "failed opening ALICE ESDfriend from '" + p + "'.");
      //}

      //fESDfriendTree = (TTree*) fESDfriendFile->Get("esdFriendTree");
      //if(fESDfriendTree == 0)
      //  throw(eH + "failed getting the esdFriendTree.");
      //fESDfriendTree->SetBranchAddress("ESDfriend", &fESDfriend);
      //if(fESDfriendTree->GetEntry(fEventId) <= 0)
      //throw(eH + "failed getting required event from ESDfriend.");

      fESDfriendExists = kTRUE;
      fESDTree->SetBranchStatus ("ESDfriend*", 1);
      fESDTree->SetBranchAddress("ESDfriend.", &fESDfriend);

    }

    if(fESDTree->GetEntry(fEventId) <= 0)
      throw(eH + "failed getting required event from ESD.");

    if (fESDfriendExists)
      fESD->SetESDfriend(fESDfriend);
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

  DestroyElements();
  fEventId = event;
  SetName(Form("Event %d", fEventId));
  UpdateItems();

  if(fRunLoader) {
    if(fRunLoader->GetEvent(fEventId) != 0)
      throw(eH + "failed getting required event.");
  }

  if(fESDTree) {
    delete fESD;       fESD       = 0;
    delete fESDfriend; fESDfriend = 0;

    if(fESDTree->GetEntry(fEventId) <= 0)
      throw(eH + "failed getting required event from ESD.");

    //if(fESDfriendTree != 0) {
    //  if(fESDfriendTree->GetEntry(fEventId) <= 0)
    //	throw(eH + "failed getting required event from ESDfriend.");

    if (fESDfriendExists)
      fESD->SetESDfriend(fESDfriend);
  }
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
