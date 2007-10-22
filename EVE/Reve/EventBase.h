// $Header$

#ifndef REVE_EventBase_H
#define REVE_EventBase_H

#include <Reve/RenderElement.h>

#include <TList.h>

namespace Reve {

class EventBase : public RenderElementList
{
protected:
  TList        fNewEventCommands;

public:
  EventBase(const Text_t* n="EventBase", const Text_t* t="");
  virtual ~EventBase() {}

  TList& GetNewEventCommands() { return fNewEventCommands; }

  virtual void Open() {}
  virtual void GotoEvent(Int_t /*event*/) {}
  virtual void NextEvent() {}
  virtual void PrevEvent() {}
  virtual void Close() {}

  virtual void AfterNewEventLoaded();
  virtual void AddNewEventCommand(const Text_t* cmd);


  ClassDef(EventBase, 1);
}; // endclass Event

}

#endif
