// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

#ifndef ALIEVE_TPCSectorVizEditor_H
#define ALIEVE_TPCSectorVizEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;
class TGHSlider;

class TEveGValuator;
class TEveGDoubleValuator;
class TEveTransSubEditor;


class AliEveTPCSectorViz;

class AliEveTPCSectorVizEditor : public TGedFrame
{
  AliEveTPCSectorVizEditor(const AliEveTPCSectorVizEditor&);            // Not implemented
  AliEveTPCSectorVizEditor& operator=(const AliEveTPCSectorVizEditor&); // Not implemented

protected:
  AliEveTPCSectorViz* fM; // fModel dynamic-casted to AliEveTPCSectorVizEditor

  TEveTransSubEditor* fHMTrans;

  TEveGValuator* fSectorID;
  TGCheckButton*    fAutoTrans;

  TGCheckButton*    fRnrInn;
  TGCheckButton*    fRnrOut1;
  TGCheckButton*    fRnrOut2;

  TEveGValuator* fThreshold;
  TEveGValuator* fMaxVal;   

  TEveGDoubleValuator* fTime;

public:
  AliEveTPCSectorVizEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		     UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~AliEveTPCSectorVizEditor();

  virtual void SetModel(TObject* obj);

  void DoSectorID();
  void DoAutoTrans();

  void DoRnrInn();
  void DoRnrOut1();
  void DoRnrOut2();

  void DoThreshold();
  void DoMaxVal();

  void DoTime();
 
  ClassDef(AliEveTPCSectorVizEditor, 0); // Editor for AliEveTPCSectorViz
}; // endclass AliEveTPCSectorVizEditor

#endif
