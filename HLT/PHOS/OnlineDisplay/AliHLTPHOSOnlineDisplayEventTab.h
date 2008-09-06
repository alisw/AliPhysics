#ifndef ALIHLTPHOSONLINEDISPLAYEVENTTAB_H
#define ALIHLTPHOSONLINEDISPLAYEVENTTAB_H

#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include "AliHLTPHOSOnlineDisplayTab.h"
#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>
#include "AliHLTPHOSOnlineDisplayTH2D.h"
#include "AliHLTPHOSConstants.h"

#define N_ZRCU_COORD 2
#define N_XRCU_COORD 2

using namespace PhosHLTConst;

class AliHLTPHOSGetEventButton;
class AliHLTHOMERReader;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSOnlineDisplay;
class AliHLTPHOSSharedMemoryInterface;


class AliHLTPHOSOnlineDisplayEventTab : public AliHLTPHOSOnlineDisplayTab
{
 public:
 
  virtual ~AliHLTPHOSOnlineDisplayEventTab();
 
  AliHLTPHOSOnlineDisplayEventTab(AliHLTPHOSOnlineDisplay *onlineDisplayPtr, TGTab *tabPtr, 
				  AliHLTHOMERReader *fgHomerReaderPtr, 
				  AliHLTHOMERReader *fgHomerReadersPtr[MAX_HOSTS], 
				  int nHosts, const int runnumber = -1);
    //    {

 
  
 

/* 
  void SetRunNumber(const int runnumber) 
  {
    
    fRunNumber = runnumber ;
    cout << __FILE__ <<":"<< __LINE__ << "RunNumber was set to "<< fRunNumber  <<endl; ;
  };
  */

  Int_t GetRawData(TH1D *histPtr, int x, int z, int gain);
  void UpdateDisplay();
  int GetNextEvent();
  virtual void ReadBlockData(AliHLTHOMERReader *homeReaderPtr);
  void FindFourierBlocks(AliHLTHOMERReader *homeReaderPtr);

  void ResetDisplay();
  TGTab               *fTab;
  TGTab               *fSubTab1;
  TRootEmbeddedCanvas *fEc1, *fEc2, *fEc3, *fEc4, *fEc5, *fEc6;
  TGCompositeFrame    *fSubF1, *fSubF2, *fSubF3;
  TCanvas *fgCanvasPtr[N_GAINS];
  AliHLTPHOSOnlineDisplayTH2D *fgLegoPlotPtr[N_GAINS];
  int *fChannelData[N_MODULES][N_XRCU_COORD][N_ZRCU_COORD][N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  Int_t fNChannelSamples[N_MODULES][N_XRCU_COORD][N_ZRCU_COORD][N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  Int_t fChannelEnergy[N_MODULES][N_XRCU_COORD][N_ZRCU_COORD][N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];

 protected:
  Bool_t fgAccumulate;

 private:
  AliHLTPHOSOnlineDisplayEventTab();
  AliHLTPHOSGetEventButton* fgEventButtPtr; 
  void InitDisplay(TGTab *tabPtr) {};
  void InitDisplay(TGTab *tabPtr, const int runnumber);

  AliHLTPHOSOnlineDisplay *fOnlineDisplayPtr;
  AliHLTPHOSSharedMemoryInterface *fShmPtr;   

  ///int fEvent

};


#endif
