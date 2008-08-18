// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTrackCounterEditor_H
#define AliEveTrackCounterEditor_H

#include "TGedFrame.h"

class TGComboBox;
class TGLabel;
class TGNumberEntry;

class AliEveTrackCounter;

//______________________________________________________________________________
// Short description of AliEveTrackCounterEditor
//

class AliEveTrackCounterEditor : public TGedFrame
{
public:
   AliEveTrackCounterEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                            UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
   virtual ~AliEveTrackCounterEditor();

   void UpdateModel();

   virtual void SetModel(TObject* obj);

   void DoPrev();
   void DoNext();
   void DoSetEvent();

   void DoPrintReport();
   void DoFileReport();
   void DoShowHistos();

   void DoClickAction(Int_t);

protected:
   AliEveTrackCounter *fM; // Model object.

   // Declare widgets
   TGComboBox       *fClickAction;
   TGLabel          *fInfoLabel;
   TGNumberEntry    *fEventId;

private:
   AliEveTrackCounterEditor(const AliEveTrackCounterEditor&);            // Not implemented
   AliEveTrackCounterEditor& operator=(const AliEveTrackCounterEditor&); // Not implemented

   ClassDef(AliEveTrackCounterEditor, 0); // GUI editor for AliEveTrackCounter.
};

#endif
