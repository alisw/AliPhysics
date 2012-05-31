//-*- Mode: C++ -*-
// $Id: AliHLTEMCALOnlineDisplayTab.h 35108 2009-09-30 01:58:37Z phille $


#ifndef ALIHLTEMCALONLINEDISPLAYTAB_H
#define ALIHLTEMCALONLINEDISPLAYTAB_H

#include "TGTab.h"
#include "AliHLTHOMERReader.h"
//#include "AliHLTEMCALCommonDefs.h"
//#include "AliHLTEMCALConstants.h"
#include "AliHLTCaloConstants.h"


#include <iostream>
using namespace std;

//#include "AliHLTEMCALBase.h"

//#define XRANGESTART 120
//#define XRANGELENGTH 80

#define XRANGESTART 128
#define XRANGELENGTH 64*3

//#define XRANGESTART 0
//#define XRANGELENGTH 150

#define XRANGEEND  XRANGESTART + XRANGELENGTH


//using namespace EmcalHLTConst;
//using namespace CaloHLTConst;


class AliHLTHOMERReader;

using CALO::MAXHOSTS;


//class AliHLTEMCALOnlineDisplayTab : public TGTab, public AliHLTEMCALBase
class AliHLTEMCALOnlineDisplayTab : public TGTab
{
 public:
  virtual ~AliHLTEMCALOnlineDisplayTab();
  AliHLTEMCALOnlineDisplayTab();

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
