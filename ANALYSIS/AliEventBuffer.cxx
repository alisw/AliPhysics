#include "AliEventBuffer.h"

ClassImp(AliEventBuffer)

//______________________________________________________
////////////////////////////////////////////////////////
//
// class AliEventBuffer
//
// FIFO type event buffer
//
// Piotr.Skowronski@cern.ch
//
////////////////////////////////////////////////////////

AliEventBuffer::AliEventBuffer():
 fSize(-1),fEvents(),fIter(&fEvents)
{
  //ctor
}
/***********************************************************/

AliEventBuffer::AliEventBuffer(Int_t size):
 fSize(size),fEvents(),fIter(&fEvents)
{
  //ctor
}
/***********************************************************/

AliEventBuffer::~AliEventBuffer()
{
  //dtor -- TList::IsOwner(1) does not work - Valgrind says that there is mem leak
  //take care owerseves
  if (fEvents.IsOwner())
   { 
     AliAOD* e=0x0;
     while (( e=RemoveLast() )) delete e;
   }
}
/***********************************************************/

void AliEventBuffer::Reset()
{
  //Resets the queue
  if (fEvents.IsOwner())
   { 
     AliAOD* e=0x0;
     while (( e=RemoveLast() )) delete e;
   }
  else
   {
     fEvents.RemoveAll();
   } 
}
/***********************************************************/

AliAOD* AliEventBuffer::Push(AliAOD* event)
{
  //adds a new event, and returns old of do not fit in size
  if (fSize == 0) return event;
  
  AliAOD* ret = 0x0;
  
  if (fSize == fEvents.GetSize()) 
    ret = dynamic_cast<AliAOD*>(fEvents.Remove(fEvents.Last()));
  if (event) fEvents.AddFirst(event);
  return ret;
}

