//-*- Mode: C++ -*-
// $Id: AliHLTEMCALOnlineDisplay.h 35108 2009-09-30 01:58:37Z phille $

#ifndef ALIHLTEMCALONLINEDISPLAY_H
#define ALIHLTEMCALONLINEDISPLAY_H

/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#include "Rtypes.h"
#include <TGFrame.h>
#include "TH2.h"
#include "AliHLTCaloConstants.h"

using CALO::MAXHOSTS;
#define DEFAULTEVENTPORT 42001

class TStyle;
class TSystem;
class AliHLTEMCALOnlineDisplayFourierTab;
class AliHLTEMCALOnlineDisplayEventTab;
class AliHLTEMCALConstants;
class AliHLTEMCALOnlineDisplayCalibTab;
class AliHLTEMCALOnlineDisplayEventTab;
class TGFrame;
class TCanvas;
class TRootEmbeddedCanvas;
class TGTab;
class AliHLTEMCALGetEventButton;
class TH2;
class TGFrame;

class AliHLTHOMERWriter;
class AliHLTHOMERReader;
class AliHLTHOMERData;


//using namespace PhosHLTConst;

 
using namespace  CaloHLTConst;

//#define N_SAMPLES 70
//#define N_SAMPLES 140
//#define MAX_HISTOGRAMS 25
#define MAXHISTOGRAMS 320

class TCanvas;

//class AliHLTEMCALOnlineDisplay : public  TGMainFrame, public AliHLTEMCALBase
class AliHLTEMCALOnlineDisplay : public  TGMainFrame
{
 public:
  ~AliHLTEMCALOnlineDisplay();

  /** Copy constructor */  
  AliHLTEMCALOnlineDisplay(const AliHLTEMCALOnlineDisplay &) : 
    TGMainFrame(),
    //   AliHLTEMCALBase(),
    fRunNumber(0),
    fgRawDataCanvas(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTEMCALOnlineDisplay & operator = (const AliHLTEMCALOnlineDisplay)
  {
    //Assignment
    return *this; 
  }

  
  int GetNextEvent();
  int GetHistogram();
  
  void InitDisplay();
  void EvaluateAverage();
  int ScanArguments(int argc, char** argv);
  static AliHLTEMCALOnlineDisplay* Instance(int argc, char** argv);  
  static AliHLTEMCALOnlineDisplayEventTab  *fgEventTabPtr; //COMMENT
  static AliHLTEMCALOnlineDisplayFourierTab  *fgFourierTabPtr; //COMMENT
  void Gain2Text(const int gain,  char *txt) const;

 protected:
    
  int fRunNumber; //COMMENT
  //  bool fIsSetRunNumber;


 private:
  AliHLTEMCALOnlineDisplay();
  AliHLTEMCALOnlineDisplay(int argc, char** argv);
  static AliHLTEMCALOnlineDisplayCalibTab  *fgCalibTabPtr; //COMMENT
  static TGTab               *fgTab; //COMMENT
  static AliHLTEMCALOnlineDisplay* fgInstancePtr; //COMMENT
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
