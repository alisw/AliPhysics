#ifndef ALIHBTEVENTBUFFER_H
#define ALIHBTEVENTBUFFER_H

#include <TObject.h>
#include <TList.h>
#include "AliHBTEvent.h"

class AliHBTEventBuffer: public TObject
{
  public:
    AliHBTEventBuffer();
    AliHBTEventBuffer(Int_t size);
    virtual ~AliHBTEventBuffer(){}
    
    AliHBTEvent* Push(AliHBTEvent* event);//adds a new event, and returns old of do not fit in size
    AliHBTEvent* RemoveLast(){return dynamic_cast<AliHBTEvent*>(fEvents.Remove(fEvents.Last()));}
    void         ResetIter(){fIter.Reset();}
    AliHBTEvent* Next(){return dynamic_cast<AliHBTEvent*>( fIter.Next() );}
    void         SetSize(Int_t size){fSize = size;}
    Int_t        GetSize() const {return fSize;}
    void         SetOwner(Bool_t flag) {fEvents.SetOwner(flag);}
  protected:
  private:
    Int_t  fSize;//size of buffer; if 0 infinite size
    TList  fEvents;//list with arrays
    TIter  fIter;//iterator
    ClassDef(AliHBTEventBuffer,1)
};


#endif
