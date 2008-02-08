// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_TOFStripEditor_H
#define ALIEVE_TOFStripEditor_H

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
private:
  AliEveTOFStripEditor(const AliEveTOFStripEditor&);            // Not implemented
  AliEveTOFStripEditor& operator=(const AliEveTOFStripEditor&); // Not implemented

protected:
  AliEveTOFStrip* fM; // fModel dynamic-casted to AliEveTOFStripEditor

  // Declare widgets
  // TGSomeWidget*   fXYZZ;

  TEveGValuator*    fThreshold;
  TEveGValuator*    fMaxVal;

public:
  AliEveTOFStripEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTOFStripEditor();

  virtual void SetModel(TObject* obj);
  void DoThreshold();
  void DoMaxVal();

  // Declare callback/slot methods
  // void DoXYZZ();

  ClassDef(AliEveTOFStripEditor, 0); // Editor for AliEveTOFStrip
}; // endclass AliEveTOFStripEditor


#endif
