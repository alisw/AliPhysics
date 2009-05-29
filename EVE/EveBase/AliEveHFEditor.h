// @(#)root/eve:$Id$
// Main author: Davide Caffarri 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveHFEditor_H
#define AliEveHFEditor_H

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class AliEveHF;

//______________________________________________________________________________
// Short description of AliEveHFEditor
//

class AliEveHFEditor : public TGedFrame
{
public:
  AliEveHFEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                 UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveHFEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();
  void DisplayDetailed();

protected:
  AliEveHF  *fM; // Model object.

  Int_t     fnProng;

  TGLabel   *fInfoLabel0; // label
  TGLabel   *fInfoLabel1; // label
  TGLabel   *fInfoLabel2; // label

  TGButton  *fXButton;

private:
  AliEveHFEditor(const AliEveHFEditor&);            // Not implemented
  AliEveHFEditor& operator=(const AliEveHFEditor&); // Not implemented

  ClassDef(AliEveHFEditor, 0); // GUI editor for AliEveHF.
};

#endif
