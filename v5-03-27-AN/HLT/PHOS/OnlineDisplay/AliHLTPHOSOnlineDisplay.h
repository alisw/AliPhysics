//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSONLINEDISPLAY_H
#define ALIHLTPHOSONLINEDISPLAY_H

/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#include "Rtypes.h"
#include <TGFrame.h>
#include "TH2.h"
#include "AliHLTPHOSConstants.h"


class TStyle;
class TSystem;
class AliHLTPHOSOnlineDisplayFourierTab;
class AliHLTPHOSOnlineDisplayEventTab;
class AliHLTPHOSConstants;
class AliHLTPHOSOnlineDisplayCalibTab;
class AliHLTPHOSOnlineDisplayEventTab;
class TGFrame;
class TCanvas;
class TRootEmbeddedCanvas;
class TGTab;
class AliHLTPHOSGetEventButton;
class TH2;
class TGFrame;

class AliHLTHOMERWriter;
class AliHLTHOMERReader;
class AliHLTHOMERData;

// using namespace PhosHLTConst;

//#define N_SAMPLES 70
//#define N_SAMPLES 140
//#define MAX_HISTOGRAMS 25


#define MAXHISTOGRAMS 320
#define MAXHOSTS 20

class TCanvas;

//class AliHLTPHOSOnlineDisplay : public  TGMainFrame, public AliHLTPHOSBase
class AliHLTPHOSOnlineDisplay : public  TGMainFrame
{
 public:
  ~AliHLTPHOSOnlineDisplay();

  /** Copy constructor */  
  AliHLTPHOSOnlineDisplay(const AliHLTPHOSOnlineDisplay &) : 
    TGMainFrame(),
    //   AliHLTPHOSBase(),
    fRunNumber(0),
    fgRawDataCanvas(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSOnlineDisplay & operator = (const AliHLTPHOSOnlineDisplay)
  {
    //Assignment
    return *this; 
  }

  
  int GetNextEvent();
  int GetHistogram();
  
  void InitDisplay();
  void EvaluateAverage();
  int ScanArguments(int argc, char** argv);
  static AliHLTPHOSOnlineDisplay* Instance(int argc, char** argv);  
  static AliHLTPHOSOnlineDisplayEventTab  *fgEventTabPtr; //COMMENT
  static AliHLTPHOSOnlineDisplayFourierTab  *fgFourierTabPtr; //COMMENT
  void Gain2Text(const int gain,  char *txt) const;

 protected:
    
  int fRunNumber; //COMMENT
  //  bool fIsSetRunNumber;


 private:
  AliHLTPHOSOnlineDisplay();
  AliHLTPHOSOnlineDisplay(int argc, char** argv);
  static AliHLTPHOSOnlineDisplayCalibTab  *fgCalibTabPtr; //COMMENT
  static TGTab               *fgTab; //COMMENT
  static AliHLTPHOSOnlineDisplay* fgInstancePtr; //COMMENT
  static unsigned int fgNHosts; //COMMENT
  static unsigned int fgNPorts; //COMMENT 
  static AliHLTHOMERReader* fgHomerReaderPtr; //COMMENT 
  static AliHLTHOMERReader* fgHomerReadersPtr[MAXHOSTS]; //COMMENT
  static char  *fgHosts[MAXHOSTS]; //COMMENT
  static short unsigned    *fgPorts;  //COMMENT 
  static Bool_t fgAccumulate; //COMMENT 
  static Bool_t fgSyncronize; //COMMENT
  TCanvas  *fgRawDataCanvas; //COMMENT
  TH1D     *fgRawDataPlotsPtr[MAXHISTOGRAMS]; //COMMENT
  
  // int fRunNumber;
  //  bool fIsSetRunNumber;


};


#endif
