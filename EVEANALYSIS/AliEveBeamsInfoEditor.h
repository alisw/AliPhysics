// $Id$
// Author: Stefano Carrazza 2010, CERN, stefano.carrazza@cern.ch

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVEBEAMSINFOEDITOR_H
#define ALIEVEBEAMSINFOEDITOR_H

#include "TGedFrame.h"

class TGButton;
class TGColorSelect;
class TGComboBox;
class TGCheckButton;
class TGGroupFrame;
class TGNumberEntry;
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

   // Set Methods
   void SetAlpha();

   // Slot methods
   void SelectEventSelection(Int_t id);
   void ShowEventSelection();
   void ShowPrevEvent();
   void ShowNextEvent();
   void SwitchDataType();

protected:
   AliEveBeamsInfo *fM;              // Model object.

private:
   TGCheckButton   *fIsMC;           // activating mc selection
   TGGroupFrame    *fEventSelection; // event selection group box
   TGCheckButton   *fShowEvents;     // display information checkbox
   TGComboBox      *fSelect;         // combo box display
   TGTextButton    *fButtonPrev;     // previous event selection
   TGTextButton    *fButtonNext;     // next event selection
   TGGroupFrame    *fSetAlpha;       // set alpha for overlay buttons
   TGNumberEntry   *fAlpha;          // alpha value

   AliEveBeamsInfoEditor(const AliEveBeamsInfoEditor&);            // Not implemented
   AliEveBeamsInfoEditor& operator=(const AliEveBeamsInfoEditor&); // Not implemented

   ClassDef(AliEveBeamsInfoEditor, 0); // GUI editor for AliEveBeamsInfo.
};

#endif
