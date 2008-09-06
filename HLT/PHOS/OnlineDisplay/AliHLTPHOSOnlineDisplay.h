#ifndef ALIHLTPHOSONLINEDISPLAY
#define ALIHLTPHOSONLINEDISPLAY

/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */


#include "AliHLTHOMERData.h"
#include "AliHLTHOMERReader.h"
#include "AliHLTHOMERWriter.h"
#include "Rtypes.h"
#include <TGFrame.h>
#include "TH2.h"
#include "AliHLTPHOSGetEventButton.h" 
#include "TGTab.h"
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include "TGFrame.h"
#include "AliHLTPHOSOnlineDisplayEventTab.h"
#include "AliHLTPHOSOnlineDisplayCalibTab.h"

#include "AliHLTPHOSConstants.h"

#include "AliHLTPHOSOnlineDisplayEventTab.h"
#include "AliHLTPHOSOnlineDisplayFourierTab.h"
#include "AliHLTPHOSBase.h"

using namespace PhosHLTConst;

//#define N_SAMPLES 70
//#define N_SAMPLES 140
//#define MAX_HISTOGRAMS 25
#define MAX_HISTOGRAMS 320

class TCanvas;

class AliHLTPHOSOnlineDisplay : public  TGMainFrame, public AliHLTPHOSBase
{
 public:
  ~AliHLTPHOSOnlineDisplay();

  
  int GetNextEvent();
  int GetHistogram();
  
  void InitDisplay();
  void EvaluateAverage();
  int ScanArguments(int argc, char** argv);
  static AliHLTPHOSOnlineDisplay* Instance(int argc, char** argv);  
  static AliHLTPHOSOnlineDisplayEventTab  *fgEventTabPtr;
  static AliHLTPHOSOnlineDisplayFourierTab  *fgFourierTabPtr;
  void Gain2Text(const int gain,  char *txt) const;


 private:
  AliHLTPHOSOnlineDisplay();
  AliHLTPHOSOnlineDisplay(int argc, char** argv);
  static AliHLTPHOSOnlineDisplayCalibTab  *fgCalibTabPtr;
  static TGTab               *fTab;
  static AliHLTPHOSOnlineDisplay* fgInstancePtr;
  static unsigned int fgNHosts;
  static unsigned int fgNPorts;
  static AliHLTHOMERReader* fgHomerReaderPtr;
  static AliHLTHOMERReader* fgHomerReadersPtr[MAX_HOSTS];
  static char  *fgHosts[MAX_HOSTS];
  static short unsigned    *fgPorts;
  static Bool_t fgAccumulate;
  static Bool_t fgSyncronize;
  TCanvas  *fgRawDataCanvas;
  TH1D     *fgRawDataPlotsPtr[MAX_HISTOGRAMS];
};


#endif
