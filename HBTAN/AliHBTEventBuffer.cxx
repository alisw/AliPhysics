#include "AliHBTEventBuffer.h"

ClassImp(AliHBTEventBuffer)

//______________________________________________________
////////////////////////////////////////////////////////
//
// class AliHBTEventBuffer
//
// FIFO type event buffer
//
// Piotr.Skowronski@cern.ch
//
////////////////////////////////////////////////////////

AliHBTEventBuffer::AliHBTEventBuffer():
 fSize(-1),fEvents(),fIter(&fEvents)
{
  //ctor
}
/***********************************************************/

AliHBTEventBuffer::AliHBTEventBuffer(Int_t size):
 fSize(size),fEvents(),fIter(&fEvents)
{
  //ctor
}
/***********************************************************/

AliHBTEventBuffer::~AliHBTEventBuffer()
{
  //dtor -- TList::IsOwner(1) does not work - Valgrind says that there is mem leak
  //take care owerseves
  if (fEvents.IsOwner())
   { 
     AliHBTEvent* e=0x0;
     while (( e=RemoveLast() )) delete e;
   }
}
/***********************************************************/

AliHBTEvent* AliHBTEventBuffer::Push(AliHBTEvent* event)
{
  //adds a new event, and returns old of do not fit in size
  if (fSize == 0) return event;
  
  AliHBTEvent* ret = 0x0;
  
  if (fSize == fEvents.GetSize()) 
    ret = dynamic_cast<AliHBTEvent*>(fEvents.Remove(fEvents.Last()));
  if (event) fEvents.AddFirst(event);
  return ret;
}

