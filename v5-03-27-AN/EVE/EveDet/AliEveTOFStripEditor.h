#ifndef ALIEVETOFSTRIPEDITOR_H
#define ALIEVETOFSTRIPEDITOR_H

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

//
// 
//

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class TGHSlider;

class TEveGValuator;
class TEveGDoubleValuator;
 

class AliEveTOFStrip;

class AliEveTOFStripEditor : public TGedFrame
{
public:
  AliEveTOFStripEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTOFStripEditor() {}

  virtual void SetModel(TObject* obj);

  void DoThreshold();
  void DoMaxVal();

protected:
  AliEveTOFStrip *fM; // Model object.

  TEveGValuator  *fThreshold; // Value widget for threshold.
  TEveGValuator  *fMaxVal;    // Value widget for maximal value.

private:
  AliEveTOFStripEditor(const AliEveTOFStripEditor&);            // Not implemented
  AliEveTOFStripEditor& operator=(const AliEveTOFStripEditor&); // Not implemented

  ClassDef(AliEveTOFStripEditor, 0); // Editor for AliEveTOFStrip
}; // endclass AliEveTOFStripEditor


#endif
