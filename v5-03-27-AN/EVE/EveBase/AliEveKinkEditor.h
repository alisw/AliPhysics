// $Id$
// Author: Paraskevi Ganoti: 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveKinkEditor_H
#define AliEveKinkEditor_H

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class AliEveKink;

//______________________________________________________________________________
// Short description of AliEveKinkEditor
//

class AliEveKinkEditor : public TGedFrame
{
public:
  AliEveKinkEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                 UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveKinkEditor() {}

  virtual void SetModel(TObject* obj);
  void DisplayDetailed();

protected:
  AliEveKink  *fM; // Model object.

  TGLabel   *fInfoLabel0; // label
  TGLabel   *fInfoLabel1; // label
  TGLabel   *fInfoLabelDaughter; // label  
  
  TGButton  *fXButton;

private:
  AliEveKinkEditor(const AliEveKinkEditor&);            // Not implemented
  AliEveKinkEditor& operator=(const AliEveKinkEditor&); // Not implemented

  ClassDef(AliEveKinkEditor, 0); // GUI editor for AliEveKink.
};

#endif
