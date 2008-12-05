// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVECASCADEEDITOR_H
#define ALIEVECASCADEEDITOR_H

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class AliEveCascade;

//______________________________________________________________________________
// Short description of AliEveCascadeEditor
//

class AliEveCascadeEditor : public TGedFrame
{
public:
  AliEveCascadeEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                 UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveCascadeEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();
  void DisplayDetailed();

protected:
  AliEveCascade  *fM; // Model object.

  TGLabel   *fInfoLabel0; // label
  TGLabel   *fInfoLabel1; // label

  TGButton  *fXButton;

private:
  AliEveCascadeEditor(const AliEveCascadeEditor&);            // Not implemented
  AliEveCascadeEditor& operator=(const AliEveCascadeEditor&); // Not implemented

  ClassDef(AliEveCascadeEditor, 0); // GUI editor for AliEveCascade.
};

#endif
