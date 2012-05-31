// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTrackFitterEditor_H
#define AliEveTrackFitterEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class AliEveTrackFitter;

class AliEveTrackFitterEditor : public TGedFrame
{
public:
  AliEveTrackFitterEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTrackFitterEditor() {}

  virtual void SetModel(TObject* obj);

  void DoStart();
  void DoFit();
  void DoReset();
  void DoStop();
  void DoGraph();

protected:
  AliEveTrackFitter* fM; // fModel dynamic-casted to AliEveTrackFitterEditor

  TGTextButton* fFit;    // button to fit selection
  TGTextButton* fReset;  // button to reset selection
  TGTextButton* fStart;  // button to connect to signal
  TGTextButton* fStop;   // button to disconnect from signal
  TGTextButton* fGraph;  // button to draw graph

private:
  AliEveTrackFitterEditor(const AliEveTrackFitterEditor&);            // Not implemented
  AliEveTrackFitterEditor& operator=(const AliEveTrackFitterEditor&); // Not implemented

  ClassDef(AliEveTrackFitterEditor, 0); // Editor for AliEveTrackFitter class.
};

#endif
