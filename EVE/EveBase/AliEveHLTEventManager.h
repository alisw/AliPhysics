//
//  AliEveHLTEventManager.h
//
//  blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//
//

#ifndef __AliEveHLTEventManager__
#define __AliEveHLTEventManager__

#include "AliEveHLTEventManager.h"
class AliEveHLTEventManager : public AliEveEventManager
{
public:
    AliEveHLTEventManager(Int_t ev=0, bool storageManager=false);
    ~AliEveHLTEventManager();
    
private:
    static void* DispatchEventListenerHLT(void *arg){static_cast<AliEveEventManager*>(arg)->GetNextEventHLT();return nullptr;}
    void GetNextEvent();
    void InitOCDB(int runNo);
    void GotoEvent(Int_t event);
    void NextEvent();
    
    void* fZMQContext;
    void* fZMQeventQueue; //this is the ONLY queue for threads!
    TString fHLTPublisherAddress;
    ClassDef(AliEveHLTEventManager, 0); // Interface for getting all event components in a uniform way.
};

#endif
