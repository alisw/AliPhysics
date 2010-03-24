#ifndef ALIEVETOFDIGITSINFOEDITOR_H
#define ALIEVETOFDIGITSINFOEDITOR_H

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class AliEveTOFDigitsInfo;

class AliEveTOFDigitsInfoEditor : public TGedFrame
{
public:
  AliEveTOFDigitsInfoEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTOFDigitsInfoEditor() {}

  virtual void SetModel(TObject* obj);

protected:
  AliEveTOFDigitsInfo* fM; // Model object.

private:
  AliEveTOFDigitsInfoEditor(const AliEveTOFDigitsInfoEditor&);            // Not implemented
  AliEveTOFDigitsInfoEditor& operator=(const AliEveTOFDigitsInfoEditor&); // Not implemented

  ClassDef(AliEveTOFDigitsInfoEditor, 0); // Editor for AliEveTOFDigitsInfo
};

#endif
