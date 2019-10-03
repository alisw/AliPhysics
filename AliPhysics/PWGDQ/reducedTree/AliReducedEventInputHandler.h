#ifndef ALIREDUCEDEVENTINPUTHANDLER_H
#define ALIREDUCEDEVENTINPUTHANDLER_H

//
//     AliInputEventHandler implementation for handling AliReducedEvent information
//     Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no
//

#include "AliInputEventHandler.h"
#include "AliReducedBaseEvent.h"
//#include "AliReducedEventInfo.h"
class TTree;

class AliReducedEventInputHandler : public AliInputEventHandler {
  public:
     
   enum EReducedEventInputType {
      kReducedBaseEvent=0,     // minimal event information (AliReducedBaseEvent)
      kReducedEventInfo            // extended event information (AliReducedEventInfo)
   };
    AliReducedEventInputHandler();
    AliReducedEventInputHandler(const char* name, const char* title);
    virtual ~AliReducedEventInputHandler();
    
    virtual Bool_t                             Init(Option_t* /*opt*/) {return kTRUE;}
    virtual Bool_t                             Init(TTree* tree, Option_t* opt);
                 AliReducedBaseEvent* GetReducedEvent() const {return fReducedEvent;}
//                  AliReducedEventInfo*   GetReducedEvent() const {return fReducedEvent;}
    virtual Bool_t                             BeginEvent(Long64_t entry);
    virtual Bool_t                             Notify() { return AliVEventHandler::Notify();};
    virtual Bool_t                             Notify(const char* path);
    virtual Bool_t                             FinishEvent();
             
                 void                                SetInputEventType(Int_t type) {fEventInputOption = type;} ;
                 Int_t                               GetInputEventType() const {return fEventInputOption;};
                 
 private:
    AliReducedEventInputHandler(const AliReducedEventInputHandler& handler);             
    AliReducedEventInputHandler& operator=(const AliReducedEventInputHandler& handler);      
    
    Int_t  fEventInputOption;                          // one of the options listed in EReducedEventInputType
    AliReducedBaseEvent* fReducedEvent;   //! Pointer to the event
    //AliReducedEventInfo* fReducedEvent;   //! Pointer to the event
    
    ClassDef(AliReducedEventInputHandler, 2);
};

#endif
