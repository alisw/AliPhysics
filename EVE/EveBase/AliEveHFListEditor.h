// @(#)root/eve:$Id$
// Main author: Davide Caffarri 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveHFListEditor_H
#define AliEveHFListEditor_H

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TEveGValuator;
class TEveGDoubleValuator;
class TGComboBox;

class AliEveHFList;

//______________________________________________________________________________
// Short description of AliEveHFListEditor
//

class AliEveHFListEditor : public TGedFrame
{
public:
  AliEveHFListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                     UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveHFListEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  void DoMinMaxPt();
  void DoMinMaxCosPointingAngle();
  void DoMinMaxInvMass();


protected:
  AliEveHFList            *fM; // Model object.

  // Declare widgets

  TEveGDoubleValuator* fMinMaxPt;
  TEveGDoubleValuator* fMinMaxCosPointingAngle;
  TEveGDoubleValuator* fValueInvMass;

private:
  AliEveHFListEditor(const AliEveHFListEditor&);            // Not implemented
  AliEveHFListEditor& operator=(const AliEveHFListEditor&); // Not implemented

  ClassDef(AliEveHFListEditor, 0); // GUI editor for AliEveV0List.
};

#endif
