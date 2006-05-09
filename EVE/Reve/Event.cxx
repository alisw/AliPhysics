// $Header$

#include "Event.h"

#include <TClass.h>

using namespace Reve;

//______________________________________________________________________
// Event
//

ClassImp(Event)

/**************************************************************************/

const TString Event::sVSDHeaderName("VSDheader");

void Event::Init()
{
  SetName(sVSDHeaderName);
  fDirectory = 0;
  fCreator   = 0;
  fSelector  = 0;
}

Event::Event() :
  fRun(0), fEvent(0)
{ Init(); }

Event::Event(Int_t run, Int_t evt, const TString& url) :
  fRun(run), fEvent(evt), fUrl(url)
{ Init(); }

Event::Event(const TString& url) :
  fEvent(0), fUrl(url)
{ Init(); }

/**************************************************************************/

Event* Event::OpenDirectory(const TString& dir_name)
{
  static const Exc_t eH("Event::OpenDirectory ");

  TDirectory* dir = dynamic_cast<TDirectory*>(gDirectory->Get(dir_name));
  if (!dir)
    throw(eH + "direcotry '" + dir_name + "' not found.");
  return OpenDirectory(dir);
}

Event* Event::OpenDirectory(TDirectory* dir)
{
  static const Exc_t eH("Event::OpenDirectory ");

  Event* evt = dynamic_cast<Event*>(dir->Get(sVSDHeaderName));
  if (!evt)
    throw(eH + "VSD header '" + sVSDHeaderName + "' not found.");
  evt->fDirectory = dir;
  return evt;
}

/**************************************************************************/
/**************************************************************************/

void Event::SetDirectory(TDirectory* dir)
{
  if(fDirectory)
    fDirectory->RecursiveRemove(this);
  fDirectory = dir;
  if(fDirectory)
    fDirectory->Append(this);
}

TDirectory* Event::MakeDirectory(const Text_t* name, const Text_t* title)
{
  TDirectory* dir = new TDirectory(name, title);
  SetDirectory(dir);
  return fDirectory;
}

/**************************************************************************/

void Event::Print(Option_t* ) const
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
