#include "AliHBTEventBuffer.h"

ClassImp(AliHBTEventBuffer)

//______________________________________________________
////////////////////////////////////////////////////////
//
// class AliHBTEventBuffer
//
// FIFO type event buffer

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

