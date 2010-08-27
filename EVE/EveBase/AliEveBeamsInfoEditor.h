// $Id$
// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVEBEAMSINFOEDITOR_H
#define ALIEVEBEAMSINFOEDITOR_H

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGGroupFrame;
class TGComboBox;
class TGTextButton;

class AliEveBeamsInfo;

//______________________________________________________________________________
// Short description of AliEveBeamsInfoEditor
//

class AliEveBeamsInfoEditor : public TGedFrame
{
public:
   AliEveBeamsInfoEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
         UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
   virtual ~AliEveBeamsInfoEditor() {}

   virtual void SetModel(TObject* obj);

   // Declare callback/slot methods
   void ShowEventSelection();
   void SelectEventSelection(Int_t id);
   void ShowPrevEvent();
   void ShowNextEvent();

protected:
   AliEveBeamsInfo            *fM; // Model object.

private:
   TGGroupFrame  *fEventSelection;  // event selection group box
   TGCheckButton *fShowEvents;      // display information checkbox
   TGComboBox    *fSelect;          // combo box display
   TGTextButton  *fButtonPrev;      // previous event selection
   TGTextButton  *fButtonNext;      // next event selection

   AliEveBeamsInfoEditor(const AliEveBeamsInfoEditor&);            // Not implemented
   AliEveBeamsInfoEditor& operator=(const AliEveBeamsInfoEditor&); // Not implemented

   ClassDef(AliEveBeamsInfoEditor, 0); // GUI editor for AliEveBeamsInfo.
};

#endif
