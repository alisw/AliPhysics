
#ifndef ALIHLTPHOSONLINEDISPLAYTAB_H
#define ALIHLTPHOSONLINEDISPLAYTAB_H

#include "TGTab.h"
#include "AliHLTHOMERReader.h"
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSBase.h"

//#define X_RANGE_START 120
//#define X_RANGE_LENGTH 80

#define X_RANGE_START 128
#define X_RANGE_LENGTH 64

#define X_RANGE_END  X_RANGE_START + X_RANGE_LENGTH


using namespace PhosHLTConst;

class AliHLTHOMERReader;


class AliHLTPHOSOnlineDisplayTab : public TGTab, public AliHLTPHOSBase
{
 public:
  virtual ~AliHLTPHOSOnlineDisplayTab();
  AliHLTPHOSOnlineDisplayTab();
  virtual void InitDisplay(TGTab *tabPtr) = 0;
  void PrintBlockInfo(AliHLTHOMERReader *homeReaderPtr, int i);
  int GetEventInfo(AliHLTHOMERReader *homeReaderPtr, int i);
  virtual void ReadBlockData(AliHLTHOMERReader *homeReaderPtr) = 0;
 
  
  void SetRunNumber(const int runnumber) 
  {
    
   fRunNumber = runnumber ;
   cout << __FILE__ <<":"<< __LINE__ << "RunNumber was set to "<< fRunNumber  <<endl; ;
  };
  


 protected:
  Bool_t fgAccumulate;
  Bool_t fgSyncronize;
  AliHLTHOMERReader* fgHomerReaderPtr;
  AliHLTHOMERReader* fgHomerReadersPtr[MAX_HOSTS];
  int DoGetNextEvent();
  int fgEvntCnt;
  int fgNHosts;

  int fRunNumber;
};


#endif
