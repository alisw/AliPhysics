//
//  AliEveDataSourceHLTZMQ.h
//
//  blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//
//

#ifndef __AliEveDataSourceHLTZMQ__
#define __AliEveDataSourceHLTZMQ__

#include "AliEveDataSource.h"

class AliEveDataSourceHLTZMQ : public AliEveDataSource
{
public:
    AliEveDataSourceHLTZMQ(bool storageManager=false);
    ~AliEveDataSourceHLTZMQ();

private:
    void GotoEvent(Int_t event);
    void Init();
    void NextEvent();
    void RequestData();
    Bool_t ReceivePromptRecoParameters(Int_t runNo);
    
    void* fZMQContext;
    void* fZMQin;
    int fZMQtimeout;

    AliEveDataSourceHLTZMQ(const AliEveDataSourceHLTZMQ&);
    AliEveDataSourceHLTZMQ& operator=(const AliEveDataSourceHLTZMQ&);

    ClassDef(AliEveDataSourceHLTZMQ, 0);
};

#endif
