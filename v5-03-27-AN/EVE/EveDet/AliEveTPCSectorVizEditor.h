// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTPCSectorVizEditor_H
#define AliEveTPCSectorVizEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;
class TGHSlider;

class TEveGValuator;
class TEveGDoubleValuator;

class AliEveTPCSectorViz;

//------------------------------------------------------------------------------
// AliEveTPCSectorVizEditor
//
// Editor for AliEveTPCSectorViz.

class AliEveTPCSectorVizEditor : public TGedFrame
{
public:
  AliEveTPCSectorVizEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		     UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTPCSectorVizEditor() {}

  virtual void SetModel(TObject* obj);

  void DoSectorID();
  void DoAutoTrans();

  void DoRnrInn();
  void DoRnrOut1();
  void DoRnrOut2();

  void DoThreshold();
  void DoMaxVal();

  void DoTime();

protected:
  AliEveTPCSectorViz   *fM;          // Model dynamic-casted to AliEveTPCSectorVizEditor

  TEveGValuator        *fSectorID;   // Widget for SectorID.
  TGCheckButton        *fAutoTrans;  // Widget for AutoTrans.

  TGCheckButton        *fRnrInn;     // Widget for RnrInn.
  TGCheckButton        *fRnrOut1;    // Widget for RnrOut1.
  TGCheckButton        *fRnrOut2;    // Widget for RnrOut2.

  TEveGValuator        *fThreshold;  // Widget for Threshold.
  TEveGValuator        *fMaxVal;     // Widget for MaxVal.

  TEveGDoubleValuator  *fTime;       // Widget for time-range.

private:
  AliEveTPCSectorVizEditor(const AliEveTPCSectorVizEditor&);            // Not implemented
  AliEveTPCSectorVizEditor& operator=(const AliEveTPCSectorVizEditor&); // Not implemented

  ClassDef(AliEveTPCSectorVizEditor, 0); // Editor for AliEveTPCSectorViz.
};

#endif
