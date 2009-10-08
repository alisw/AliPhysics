//-*- Mode: C++ -*-
// $Id: AliHLTEMCALOnlineDisplayCalibTab.h 31683 2009-03-24 21:17:03Z odjuvsla $

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


#ifndef ALIHLTEMCALONLINEDISPLAYCALIBTAB_H
#define ALIHLTEMCALONLINEDISPLAYCALIBTAB_H

//#include "AliHLTEMCALCommonDefs.h"
#include <TCanvas.h>
#include "AliHLTDataTypes.h"
#include "AliHLTEMCALOnlineDisplayTab.h"
#include "AliHLTEMCALConstants.h"

using namespace EmcalHLTConst;

class TGTab;

class TRootEmbeddedCanvas;
class TH2D;
class TH2I;
class TCanvas;
class AliHLTEMCALGetEventButton;

class AliHLTEMCALOnlineDisplayCalibTab : public AliHLTEMCALOnlineDisplayTab
{
 public:
  AliHLTEMCALOnlineDisplayCalibTab(); // Default Constructor
  AliHLTEMCALOnlineDisplayCalibTab(TGTab  *tabPtr, HOMERReader *homerSyncPtr, HOMERReader *homerPtrs[MAXHOSTS], int nHosts); // Constructor
  virtual ~AliHLTEMCALOnlineDisplayCalibTab(); // Destructor

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
  AliHLTEMCALGetEventButton* fgEventButtPtr; // fgEventButtPtr

private:
  
  AliHLTEMCALOnlineDisplayCalibTab(const AliHLTEMCALOnlineDisplayCalibTab&); //copy constructor
  AliHLTEMCALOnlineDisplayCalibTab & operator=(const AliHLTEMCALOnlineDisplayCalibTab); //assignement operator

};


#endif
