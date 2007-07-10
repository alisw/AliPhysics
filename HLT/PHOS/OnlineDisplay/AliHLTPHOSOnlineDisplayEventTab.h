#ifndef ALIHLTPHOSONLINEDISPLAYEVENTTAB_H
#define ALIHLTPHOSONLINEDISPLAYEVENTTAB_H

#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include "AliHLTPHOSOnlineDisplayTab.h"
#include <TCanvas.h>
#include <TH2D.h>

#include "AliHLTPHOSConstants.h"
using namespace PhosHLTConst;


class AliHLTPHOSGetEventButton;
class HOMERReader;
class AliHLTPHOSRcuCellEnergyDataStruct;

class AliHLTPHOSOnlineDisplayEventTab : public AliHLTPHOSOnlineDisplayTab
{
 public:
  virtual ~AliHLTPHOSOnlineDisplayEventTab();
  AliHLTPHOSOnlineDisplayEventTab(TGTab *tabPtr, HOMERReader *fgHomerReaderPtr, HOMERReader *fgHomerReadersPtr[MAX_HOSTS], int nHosts);
  AliHLTPHOSOnlineDisplayEventTab();
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
  TH2D *fgLegoPlotLGPtr;
  TH2D *fgLegoPlotHGPtr;

 protected:
  Bool_t fgAccumulate;

 private:
  AliHLTPHOSGetEventButton* fgEventButtPtr; 
  void InitDisplay(TGTab *tabPtr);
};


#endif
