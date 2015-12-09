//
//  AliEveDataSourceHLTZMQ.h
//
//  blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//
//

#ifndef __AliEveDataSourceHLTZMQ__
#define __AliEveDataSourceHLTZMQ__

#include "AliEveEventManager.h"
#include "AliEveDataSource.h"
#include "AliHLTDataTypes.h"

#include "TThread.h"

class AliEveDataSourceHLTZMQ : public AliEveDataSource
{
public:
    AliEveDataSourceHLTZMQ(bool storageManager=false);
    ~AliEveDataSourceHLTZMQ();

private:
    static void* DispatchEventListenerHLT(void *arg){static_cast<AliEveDataSourceHLTZMQ*>(arg)->PullEventFromHLT();return nullptr;}
    void PullEventFromHLT();
    void InitOCDB(int runNo);
    void GotoEvent(Int_t event);
    void Init();
    void NextEvent();
    
    TThread *fEventListenerThreadHLT;
    
    void* fZMQContext;
    void* fZMQeventQueue; //this is the ONLY queue for threads!
    TString fHLTPublisherAddress;

    AliEveDataSourceHLTZMQ(const AliEveDataSourceHLTZMQ&);
    AliEveDataSourceHLTZMQ& operator=(const AliEveDataSourceHLTZMQ&);

    ClassDef(AliEveDataSourceHLTZMQ, 0); // Interface for getting all event components in a uniform way.
};

//const int kAliHLTComponentDataTypeTopicSize = kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize;
//
//inline bool Topicncmp(const char* topic, const char* reference, int topicSize=kAliHLTComponentDataTypeTopicSize, int referenceSize=kAliHLTComponentDataTypeTopicSize)
//{
//    for (int i=0; i<((topicSize<referenceSize)?topicSize:referenceSize); i++)
//    {
//        if (!(topic[i]=='*' || reference[i]=='*' || topic[i]==reference[i])) {return false;}
//    }
//    return true;
//}




#endif
