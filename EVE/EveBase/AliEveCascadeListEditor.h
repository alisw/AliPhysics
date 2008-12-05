// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVECASCADELISTEDITOR_H
#define ALIEVECASCADELISTEDITOR_H

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TEveGDoubleValuator;

class AliEveCascadeList;

//______________________________________________________________________________
// Short description of AliEveCascadeListEditor
//

class AliEveCascadeListEditor : public TGedFrame
{
public:
  AliEveCascadeListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                     UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveCascadeListEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  void DoMinMaxRCut();
  void DoMinMaxDaughterDCA();
  void DoMinMaxPt();

protected:
  AliEveCascadeList            *fM; // Model object.

  // Declare widgets
  // TGSomeWidget*   fXYZZ;
  TEveGDoubleValuator* fMinMaxRCut;
  TEveGDoubleValuator* fMinMaxDaughterDCA;
  TEveGDoubleValuator* fMinMaxPt;

private:
  AliEveCascadeListEditor(const AliEveCascadeListEditor&);            // Not implemented
  AliEveCascadeListEditor& operator=(const AliEveCascadeListEditor&); // Not implemented

  ClassDef(AliEveCascadeListEditor, 0); // GUI editor for AliEveCascadeList.
};

#endif
