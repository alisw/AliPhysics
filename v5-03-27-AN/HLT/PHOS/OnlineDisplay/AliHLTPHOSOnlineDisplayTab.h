//-*- Mode: C++ -*-
// $Id$


#ifndef ALIHLTPHOSONLINEDISPLAYTAB_H
#define ALIHLTPHOSONLINEDISPLAYTAB_H

#include "TGTab.h"
#include "AliHLTHOMERReader.h"
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSConstants.h"

#include <iostream>
using namespace std;

//#include "AliHLTPHOSBase.h"

//#define XRANGESTART 120
//#define XRANGELENGTH 80

#define XRANGESTART 128
#define XRANGELENGTH 64*3

#define XRANGEEND  XRANGESTART + XRANGELENGTH


// using namespace PhosHLTConst;


class AliHLTHOMERReader;


//class AliHLTPHOSOnlineDisplayTab : public TGTab, public AliHLTPHOSBase
class AliHLTPHOSOnlineDisplayTab : public TGTab
{
 public:
  virtual ~AliHLTPHOSOnlineDisplayTab();
  AliHLTPHOSOnlineDisplayTab();

  void PrintBlockInfo(AliHLTHOMERReader *homeReaderPtr, int i);
  int GetEventInfo(AliHLTHOMERReader *homeReaderPtr, int i);
  virtual void ReadBlockData(AliHLTHOMERReader *homeReaderPtr) = 0;
  virtual void InitDisplay(TGTab *tabPtr) = 0; 
  
  void SetRunNumber(const int runnumber) 
  {
    fRunNumber = runnumber ;
   cout << __FILE__ <<":"<< __LINE__ << "RunNumber was set to "<< fRunNumber  <<endl; ;
  };
  



 protected:
  Bool_t fgAccumulate;
  Bool_t fgSyncronize;
  AliHLTHOMERReader* fgHomerReaderPtr;
  AliHLTHOMERReader* fgHomerReadersPtr[MAXHOSTS];
  int DoGetNextEvent();
  int fgEvntCnt;
  int fgNHosts;

  int fRunNumber;




};


#endif
