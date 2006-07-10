// $Header$

#ifndef REVE_EventBase_H
#define REVE_EventBase_H

#include <Reve/RenderElement.h>

namespace Reve {

class EventBase : public RenderElementList
{
protected:

public:
  EventBase(const Text_t* n="EventBase", const Text_t* t="");
  virtual ~EventBase() {}

  ClassDef(EventBase, 1);
}; // endclass Event

}

#endif
