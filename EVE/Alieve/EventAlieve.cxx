// $Header$

#include "EventAlieve.h"
#include <Reve/Reve.h>

#include <AliRunLoader.h>
#include <AliESD.h>

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

Bool_t Alieve::Event::fgUseRunLoader = true;
Bool_t Alieve::Event::fgUseESDTree   = true;

void Event::Initialize(Bool_t use_runloader, Bool_t use_esd)
{
  static const Exc_t eH("Event::Initialize ");

  fgUseRunLoader = use_runloader;
  fgUseESDTree   = use_esd;

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

void Event::Init()
{
  fRunLoader = 0;
  fESDFile   = 0;
  fESDTree   = 0;
  fESD       = 0;
}

Event::Event() : TNamed(), fEventId(0)
{
  Init();
}

Event::Event(TString path, Int_t ev) : fPath(path), fEventId(ev)
{
  Init();
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
    if(gSystem->AccessPathName(ga_path, kReadPermission))
      throw(eH + "can not read '" + ga_path + "'.");
    fRunLoader = AliRunLoader::Open(ga_path);
    if(!fRunLoader)
      throw(eH + "failed opening ALICE run loader from '" + ga_path + "'.");
    {
      TString alice_path = fPath + "/";
      fRunLoader->SetDirName(alice_path);
    }
    if(fRunLoader->LoadgAlice() != 0) {
      throw(eH + "failed loading gAlice.");
    }

    if(fRunLoader->GetEvent(fEventId) != 0) {
      throw(eH + "failed getting required event.");
    }
  }

  if(fgUseESDTree) {
    TString p(Form("%s/AliESDs.root", fPath.Data()));
    fESDFile = new TFile(p);
    if(fESDFile->IsZombie()) {
      delete fESDFile; fESDFile = 0;
      throw(eH + "failed opening ALICE ESD from '" + p + "'.");
    }

    fESDTree = (TTree*) fESDFile->Get("esdTree");
    if(fESDTree == 0)
      throw(eH + "failed getting the esdTree.");
    fESDTree->SetBranchAddress("ESD", &fESD);
    if(fESDTree->GetEntry(fEventId) <= 0)
      throw(eH + "failed getting required event.");
  }

  SetName(Form("Event%d", fEventId));
  SetTitle(fPath);
}

void Event::Close()
{

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
