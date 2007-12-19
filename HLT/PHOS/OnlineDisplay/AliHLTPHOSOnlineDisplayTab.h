#ifndef ALIHLTPHOSONLINEDISPLAYTAB_H
#define ALIHLTPHOSONLINEDISPLAYTAB_H

#include "TGTab.h"
#include "HOMERReader.h"
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSBase.h"

//#define X_RANGE_START 120
//#define X_RANGE_LENGTH 80

#define X_RANGE_START 128
#define X_RANGE_LENGTH 64

#define X_RANGE_END  X_RANGE_START + X_RANGE_LENGTH


using namespace PhosHLTConst;

class HOMERReader;


class AliHLTPHOSOnlineDisplayTab : public TGTab, public AliHLTPHOSBase
{
 public:
  virtual ~AliHLTPHOSOnlineDisplayTab();
  AliHLTPHOSOnlineDisplayTab();
  virtual void InitDisplay(TGTab *tabPtr) = 0;
  void PrintBlockInfo(HOMERReader *homeReaderPtr, int i);
  int GetEventInfo(HOMERReader *homeReaderPtr, int i);
  virtual void ReadBlockData(HOMERReader *homeReaderPtr) = 0;

 protected:
  Bool_t fgAccumulate;
  Bool_t fgSyncronize;
  HOMERReader* fgHomerReaderPtr;
  HOMERReader* fgHomerReadersPtr[MAX_HOSTS];
  int DoGetNextEvent();
  int fgEvntCnt;
  int fgNHosts;
};


#endif
