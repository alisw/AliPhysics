// $Header$

#include "VSDEvent.h"

#include <TClass.h>

using namespace Reve;

//______________________________________________________________________
// VSDEvent
//

ClassImp(VSDEvent)

/**************************************************************************/

const TString VSDEvent::sVSDHeaderName("VSDheader");

VSDEvent::VSDEvent() :
  EventBase("VSDEvent"),
  fRun(0), fEvent(0), fUrl(),
  fTags(),
  fDirectory(0), fCreator(0), fSelector(0)
{}

VSDEvent::VSDEvent(Int_t run, Int_t evt, const TString& url) :
  EventBase("VSDEvent"),
  fRun(run), fEvent(evt), fUrl(url),
  fTags(),
  fDirectory(0), fCreator(0), fSelector(0)
{}

VSDEvent::VSDEvent(const TString& url) :
  EventBase("VSDEvent"),
  fRun(0), fEvent(0), fUrl(url),
  fTags(),
  fDirectory(0), fCreator(0), fSelector(0)
{}

/**************************************************************************/

VSDEvent* VSDEvent::OpenDirectory(const TString& dir_name)
{
  static const Exc_t eH("VSDEvent::OpenDirectory ");

  TDirectory* dir = dynamic_cast<TDirectory*>(gDirectory->Get(dir_name));
  if (!dir)
    throw(eH + "direcotry '" + dir_name + "' not found.");
  return OpenDirectory(dir);
}

VSDEvent* VSDEvent::OpenDirectory(TDirectory* dir)
{
  static const Exc_t eH("VSDEvent::OpenDirectory ");

  VSDEvent* evt = dynamic_cast<VSDEvent*>(dir->Get(sVSDHeaderName));
  if (!evt)
    throw(eH + "VSD header '" + sVSDHeaderName + "' not found.");
  evt->fDirectory = dir;
  return evt;
}

/**************************************************************************/
/**************************************************************************/

void VSDEvent::SetDirectory(TDirectory* dir)
{
  if(fDirectory)
    fDirectory->RecursiveRemove(this);
  fDirectory = dir;
  if(fDirectory)
    fDirectory->Append(this);
}

TDirectory* VSDEvent::MakeDirectory(const Text_t* name, const Text_t* title)
{
  TDirectory* dir = new TDirectory(name, title);
  SetDirectory(dir);
  return fDirectory;
}

/**************************************************************************/

void VSDEvent::Print(Option_t* ) const
{
  printf("%s: '%s', '%s'\n", IsA()->GetName(), GetName(), GetTitle());
  printf("  run=%d, event=%d, url='%s'\n", fRun, fEvent, fUrl.Data());
  if(fDirectory)
    printf("  directory: '%s', '%s'\n", fDirectory->GetName(), fDirectory->GetTitle());
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

#include <TBrowser.h>

void EvTree::Browse(TBrowser* b)
{
  // fFolder.Browse(b); // This adds all elements to top-level.
  b->Add(&fFolder);
  TTree::Browse(b);
}
