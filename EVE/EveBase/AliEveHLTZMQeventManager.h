//
//  AliEveHLTZMQeventManager.h
//
//  blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//
//

#ifndef __AliEveHLTZMQeventManager__
#define __AliEveHLTZMQeventManager__

#include "AliEveEventManager.h"
#include "AliHLTDataTypes.h"

#include "TThread.h"

class AliEveHLTZMQeventManager : public AliEveEventManager
{
public:
    AliEveHLTZMQeventManager(bool storageManager=false);
    ~AliEveHLTZMQeventManager();
    
private:
    static void* DispatchEventListenerHLT(void *arg){static_cast<AliEveHLTZMQeventManager*>(arg)->PullEventFromHLT();return nullptr;}
    void PullEventFromHLT();
    void InitOCDB(int runNo);
    void GotoEvent(Int_t event);
    void NextEvent();
    
    TThread *fEventListenerThreadHLT;
    
    int fCurrentRun;
    
    void* fZMQContext;
    void* fZMQeventQueue; //this is the ONLY queue for threads!
    TString fHLTPublisherAddress;
    ClassDef(AliEveHLTZMQeventManager, 0); // Interface for getting all event components in a uniform way.
};

const int kAliHLTComponentDataTypeTopicSize = kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize;

inline bool Topicncmp(const char* topic, const char* reference, int topicSize=kAliHLTComponentDataTypeTopicSize, int referenceSize=kAliHLTComponentDataTypeTopicSize)
{
    for (int i=0; i<((topicSize<referenceSize)?topicSize:referenceSize); i++)
    {
        if (!(topic[i]=='*' || reference[i]=='*' || topic[i]==reference[i])) {return false;}
    }
    return true;
}




#endif
