#ifndef AliEventBuffer_H
#define AliEventBuffer_H
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

#include <TObject.h>
#include <TList.h>
#include "AliAOD.h"

class AliEventBuffer: public TObject
{
  public:
    AliEventBuffer();
    AliEventBuffer(Int_t size);
    virtual ~AliEventBuffer();
    
    AliAOD* Push(AliAOD* event);//adds a new event, and returns old of do not fit in size
    AliAOD* RemoveLast(){return dynamic_cast<AliAOD*>(fEvents.Remove(fEvents.Last()));}
    void    ResetIter(){fIter.Reset();}
    AliAOD* Next(){return dynamic_cast<AliAOD*>( fIter.Next() );}
    void    SetSize(Int_t size){fSize = size;}
    Int_t   GetSize() const {return fSize;}
    void    SetOwner(Bool_t flag) {fEvents.SetOwner(flag);}
    void    Reset();
  protected:
  private:
    Int_t   fSize;//size of buffer; if -1 infinite size
    TList   fEvents;//list with arrays
    TIter   fIter;//iterator
    ClassDef(AliEventBuffer,1)
};


#endif
