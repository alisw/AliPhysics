// $Header$

#include "EventBase.h"

#include <TObjString.h>
#include <TCint.h>

using namespace Reve;

//______________________________________________________________________
// Reve::EventBase
//

ClassImp(EventBase)

EventBase::EventBase(const Text_t* n, const Text_t* t) :
  RenderElementList(n, t),
  fNewEventCommands()
{}

/**************************************************************************/

void EventBase::AfterNewEventLoaded()
{
  TIter next(&fNewEventCommands);
  TObject* o;
  while ((o = next())) {
    TObjString* s = dynamic_cast<TObjString*>(o);
    if (s)
      gInterpreter->ProcessLine(s->String());
  }
}

void EventBase::AddNewEventCommand(const Text_t* cmd)
{
  fNewEventCommands.Add(new TObjString(cmd));
}
