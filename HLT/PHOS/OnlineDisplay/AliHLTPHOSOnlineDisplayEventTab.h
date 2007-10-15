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
//#include "AliHLTPHOSOnlineDisplay.h"

//#define N_SAMPLES 70 //BAD, someone is going to pay for this
//#define N_SAMPLES 140 //BAD, someone is going to pay for this
#define N_ZRCU_COORD 2
#define N_XRCU_COORD 2

using namespace PhosHLTConst;

 
class AliHLTPHOSGetEventButton;
class HOMERReader;
//class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSOnlineDisplay;

//      AliHLTPHOSRcuCellEnergyDataStruct.h 
class AliHLTPHOSOnlineDisplayEventTab : public AliHLTPHOSOnlineDisplayTab
{
 public:
  virtual ~AliHLTPHOSOnlineDisplayEventTab();
  AliHLTPHOSOnlineDisplayEventTab(AliHLTPHOSOnlineDisplay *onlineDisplayPtr, TGTab *tabPtr, HOMERReader *fgHomerReaderPtr, HOMERReader *fgHomerReadersPtr[MAX_HOSTS], int nHosts);
  //  void GetRawData(TH1D *histPtr);
  //AliHLTPHOSOnlineDisplayEventTab::GetRawData(TH1D *histPtr, int mod, int rcuX, int rcuZ, int x, int z, int gain)
  void GetRawData(TH1D *histPtr, int mod, int rcuX, int rcuZ, int x, int z, int gain);
  void GetRawData(TH1D *histPtr, int x, int z, int gain);

  void UpdateDisplay();
  int GetNextEvent();
  virtual void ReadBlockData(HOMERReader *homeReaderPtr);
  void ResetDisplay();
  TGTab               *fTab;
  TGTab               *fSubTab1;
  TRootEmbeddedCanvas *fEc1, *fEc2, *fEc3, *fEc4, *fEc5, *fEc6;
  TGCompositeFrame    *fSubF1, *fSubF2, *fSubF3;
  TCanvas *fgCanvasHGPtr;
  TCanvas *fgCanvasLGPtr;

  //  TH2D *fgLegoPlotLGPtr;
  // TH2D *fgLegoPlotHGPtr;

  AliHLTPHOSOnlineDisplayTH2D *fgLegoPlotLGPtr;
  AliHLTPHOSOnlineDisplayTH2D *fgLegoPlotHGPtr;


  //  int *fChannelData[N_MODULES][N_RCUS_PER_MODULE][N_ZROWS_RCU][N_XCOLUMNS_RCU][N_GAINS];
  int *fChannelData[N_MODULES][N_XRCU_COORD][N_ZRCU_COORD][N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];

 protected:
  Bool_t fgAccumulate;

 private:
  AliHLTPHOSOnlineDisplayEventTab();
  AliHLTPHOSGetEventButton* fgEventButtPtr; 
  void InitDisplay(TGTab *tabPtr);

  //  AliHLTPHOSOnlineDisplay. 
  AliHLTPHOSOnlineDisplay *fOnlineDisplayPtr;
    
};


#endif
