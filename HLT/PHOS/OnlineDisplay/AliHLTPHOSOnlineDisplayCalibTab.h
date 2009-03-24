//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Per Thomas Hille for the ALICE                                *
 * offline/HLT Project. Contributors are mentioned in the code where      *
 * appropriate.                                                           *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#ifndef ALIHLTPHOSONLINEDISPLAYCALIBTAB_H
#define ALIHLTPHOSONLINEDISPLAYCALIBTAB_H

#include "AliHLTPHOSCommonDefs.h"
#include <TCanvas.h>
#include "AliHLTDataTypes.h"
#include "AliHLTPHOSOnlineDisplayTab.h"
#include "AliHLTPHOSConstants.h"

using namespace PhosHLTConst;

class TGTab;

class TRootEmbeddedCanvas;
class TH2D;
class TH2I;
class TCanvas;
class AliHLTPHOSGetEventButton;

class AliHLTPHOSOnlineDisplayCalibTab : public AliHLTPHOSOnlineDisplayTab
{
 public:
  AliHLTPHOSOnlineDisplayCalibTab(); // Default Constructor
  AliHLTPHOSOnlineDisplayCalibTab(TGTab  *tabPtr, HOMERReader *homerSyncPtr, HOMERReader *homerPtrs[MAXHOSTS], int nHosts); // Constructor
  virtual ~AliHLTPHOSOnlineDisplayCalibTab(); // Destructor

  void InitDisplay(TGTab *tabPtr); // InitDisplay
  void UpdateDisplay(); // UpdateDisplay
  int GetNextEvent(); // GetNextEvent
  virtual void ReadBlockData(HOMERReader *homeReaderPtr); // ReadBlockData
  void ResetDisplay() const; // ResetDisplay

  TH2D *fgCalibHistPtr[NGAINS]; // fgCalibHistPtr
  TH2I *fgHitsHistPtr[NGAINS];  // fgHitsHistPtr
  TH2D *fgAveragePtr[NGAINS]; // fgAveragePtr
  TH2D *fgDCSViewPtr[NGAINS]; // fgDCSViewPtr
       
  TH2D *fDeadCannelMapPtr[NGAINS]; // fDeadCannelMapPtr
  TGTab               *fTab; // fTab
  TRootEmbeddedCanvas *fEc7, *fEc8, *fEc9, *fEc10, *fEc11, *fEc12, *fEc13, *fEc14, *fEc15, *fEc16, *fEc17, *fEc18; // Canvases
  TGTab               *fSubTab2; // fSubTab2
  TGCompositeFrame    *fSubF4, *fSubF5, *fSubF6, *fSubF7,*fSubF8, *fSubF9; // Composite frames
  TCanvas *fgCanvasHGPtr; // fgCanvasHGPtr
  TCanvas *fgCanvasLGPtr; // fgCanvasLGPtr
  TH2D *fgLegoPlotLGPtr; // fgLegoPlotLGPtr
  TH2D *fgLegoPlotHGPtr; // fgLegoPlotHGPtr
  AliHLTPHOSGetEventButton* fgEventButtPtr; // fgEventButtPtr

private:
  
  AliHLTPHOSOnlineDisplayCalibTab(const AliHLTPHOSOnlineDisplayCalibTab&); //copy constructor
  AliHLTPHOSOnlineDisplayCalibTab & operator=(const AliHLTPHOSOnlineDisplayCalibTab); //assignement operator

};


#endif
