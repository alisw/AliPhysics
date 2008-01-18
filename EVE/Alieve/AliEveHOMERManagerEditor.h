// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

#ifndef ALIEVE_AliEVEHOMERManagerEditor_H
#define ALIEVE_AliEVEHOMERManagerEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGTextButton;
class TGNumberEntry;
class TGColorSelect;

class AliEveHOMERManager;

class AliEveHOMERManagerEditor : public TGedFrame
{
private:
  AliEveHOMERManagerEditor(const AliEveHOMERManagerEditor&);            // Not implemented
  AliEveHOMERManagerEditor& operator=(const AliEveHOMERManagerEditor&); // Not implemented

protected:
  AliEveHOMERManager* fM; // fModel dynamic-casted to AliEveHOMERManagerEditor

  // Declare widgets
  // TGSomeWidget*   fXYZZ;
  TGTextButton  *fButt;

public:
  AliEveHOMERManagerEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveHOMERManagerEditor();

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();
  void DoButt();

  ClassDef(AliEveHOMERManagerEditor, 0); // Editor for AliEveHOMERManager
}; // endclass AliEveHOMERManagerEditor

#endif
