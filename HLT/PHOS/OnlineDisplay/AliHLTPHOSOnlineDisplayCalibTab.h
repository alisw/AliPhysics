//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSONLINEDISPLAYCALIBTAB_H
#define ALIHLTPHOSONLINEDISPLAYCALIBTAB_H

#include "TGTab.h"
#include "AliHLTPHOSOnlineDisplayTab.h"
#include "TH2.h"
#include "AliHLTPHOSCommonDefs.h"
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include "AliHLTDataTypes.h"


#include "AliHLTPHOSConstants.h"
using namespace PhosHLTConst;


class AliHLTPHOSGetEventButton;
class AliHLTPHOSOnlineDisplayCalibTab : public AliHLTPHOSOnlineDisplayTab
{
 public:
  AliHLTPHOSOnlineDisplayCalibTab();
  AliHLTPHOSOnlineDisplayCalibTab(TGTab  *tabPtr, HOMERReader *homerSyncPtr, HOMERReader *homerPtrs[MAXHOSTS], int nHosts);
  virtual ~AliHLTPHOSOnlineDisplayCalibTab();

  void InitDisplay(TGTab *tabPtr);
  void UpdateDisplay();
  int GetNextEvent();
  virtual void ReadBlockData(HOMERReader *homeReaderPtr);
  void ResetDisplay();

  TH2D *fgCalibHistPtr[NGAINS];
  TH2I *fgHitsHistPtr[NGAINS]; 
  TH2D *fgAveragePtr[NGAINS];
  TH2D *fgDCSViewPtr[NGAINS];
       
  TH2D *fDeadCannelMapPtr[NGAINS];
  TGTab               *fTab;
  TRootEmbeddedCanvas *fEc7, *fEc8, *fEc9, *fEc10, *fEc11, *fEc12, *fEc13, *fEc14, *fEc15, *fEc16, *fEc17, *fEc18;
  TGTab               *fSubTab2;
  TGCompositeFrame    *fSubF4, *fSubF5, *fSubF6, *fSubF7,*fSubF8, *fSubF9;
  TCanvas *fgCanvasHGPtr;
  TCanvas *fgCanvasLGPtr;
  TH2D *fgLegoPlotLGPtr;
  TH2D *fgLegoPlotHGPtr;
  AliHLTPHOSGetEventButton* fgEventButtPtr; 
};


#endif
