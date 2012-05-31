// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveV0Editor_H
#define AliEveV0Editor_H

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class AliEveV0;

//______________________________________________________________________________
// Short description of AliEveV0Editor
//

class AliEveV0Editor : public TGedFrame
{
public:
  AliEveV0Editor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                 UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveV0Editor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();
  void DisplayDetailed();

protected:
  AliEveV0  *fM; // Model object.

  TGLabel   *fInfoLabel0; // label
  TGLabel   *fInfoLabel1; // label
  TGLabel   *fInfoLabelNegDaughter; // label
  TGLabel   *fInfoLabelPosDaughter; // label

  TGButton  *fXButton;

private:
  AliEveV0Editor(const AliEveV0Editor&);            // Not implemented
  AliEveV0Editor& operator=(const AliEveV0Editor&); // Not implemented

  ClassDef(AliEveV0Editor, 0); // GUI editor for AliEveV0.
};

#endif
