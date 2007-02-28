#ifndef ALIHLTPHOSONLINEDISPLAY
#define ALIHLTPHOSONLINEDISPLAY

/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#include "HOMERData.h"
#include "HOMERReader.h"
#include "HOMERWriter.h"
#include "Rtypes.h"
#include <TGFrame.h>
#include "TH2.h"
#include "AliHLTPHOSGetEventButton.h" 
#include "TGTab.h"
#include <TRootEmbeddedCanvas.h>
#include "TGFrame.h"
#include "AliHLTPHOSCommonDefs.h"

#define MAX_HOSTS 10
#define MAX_HOSTNAME_LENGTH 64
#define DEFAULT_EVENT_PORT 42001 
#define DEFAULT_HISTO_PORT 42002 

class TCanvas;
class AliHLTPHOSRcuCellEnergyDataStruct;

class AliHLTPHOSOnlineDisplay : public  TGMainFrame
{
 public:
  ~AliHLTPHOSOnlineDisplay();
  int GetNextEvent();
  int GetHistogram();
  void InitDisplay();
  void UpdateDisplay();
  void UpdateHistograms();
  static int ScanArguments(int argc, char** argv);
  static AliHLTPHOSOnlineDisplay* Instance();  
  
 private:
  AliHLTPHOSOnlineDisplay();
  static TGCompositeFrame    *fFrame1, *fF1, *fF2, *fF3, *fF4, *fF5, *fSubF1, *fSubF2, *fSubF3, *fSubF4, *fSubF5, *fSubF6, *fSubF7;
  static TGTab               *fTab;
  static TGTab               *fSubTab1;
  static TGTab               *fSubTab2;
  static TRootEmbeddedCanvas *fEc1, *fEc2, *fEc3, *fEc4, *fEc5, *fEc6, *fEc7, *fEc8, *fEc9, *fEc10, *fEc11, *fEc12;
  static AliHLTPHOSGetEventButton* fgEventButtPtr; 
  static AliHLTPHOSOnlineDisplay* fgInstancePtr;
  static TH2D *fgLegoPlotLGPtr;
  static TH2D *fgLegoPlotHGPtr;
  static TH2D *fgCalibHistPtr[N_GAINS];
  static TH2I *fgHitsHistPtr[N_GAINS]; 
  static char *fgDefaultDet;        
  static char *fgDefaultDataType;   
  static int fgEvntCnt;
  static TCanvas *fgCanvasHGPtr;
  static TCanvas *fgCanvasLGPtr;
  static unsigned int fgNHosts;
  static unsigned int fgNPorts;
  static HOMERReader* fgHomerReaderPtr;
  static HOMERReader* fgHomerReadersPtr[MAX_HOSTS];
  static HOMERReader* fgCalibReadersPtr[MAX_HOSTS];
  static char  *fgHosts[MAX_HOSTS];
  static short unsigned    *fgPorts;
  static Bool_t fgAccumulate;
  static Bool_t fgSyncronize;
};


#endif
