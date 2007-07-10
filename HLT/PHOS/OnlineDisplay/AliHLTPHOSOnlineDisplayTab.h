#ifndef ALIHLTPHOSONLINEDISPLAYTAB_H
#define ALIHLTPHOSONLINEDISPLAYTAB_H

#include "TGTab.h"
#include "HOMERReader.h"
#include "AliHLTPHOSCommonDefs.h"

#include "AliHLTPHOSConstants.h"
using namespace PhosHLTConst;

class HOMERReader;


class AliHLTPHOSOnlineDisplayTab : public TGTab
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
