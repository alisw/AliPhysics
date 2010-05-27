// $Id$
// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveLegoEditor_H
#define AliEveLegoEditor_H

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGButtonGroup;
class TGRadioButton;
class TGLabel;
class TGComboBox;
class TGGroupFrame;

class AliEveLego;

//______________________________________________________________________________
// Short description of AliEveLegoEditor
//

class AliEveLegoEditor : public TGedFrame
{
public:
   AliEveLegoEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
         UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
   virtual ~AliEveLegoEditor() {}

   virtual void SetModel(TObject* obj);

   // Declare callback/slot methods
   void DoAllEvents();
   void ShowByCharge(Int_t id);
   void ShowByChargeAE(Int_t id);
   void SetThreshold();
   void SetThresholdAE();
   void SetMaxPt();
   void SetMaxPtAE();
   void ShowByTracks(Int_t id);
   void ShowByTracksAE(Int_t id);
   void ShowByEvents(Int_t id);
   void ShowEventSelection();
   void SelectEventSelection(Int_t id);
   void CreateAllEventsEditor();
   void ShowPrevEvent();
   void ShowNextEvent();

protected:
   AliEveLego            *fM; // Model object.

   // Single event GUI
   TGTextButton  *fAllEventsButton;
   TGButtonGroup *fParticlesBG;
   TGButtonGroup *fTrackSelection;
   TGGroupFrame *fEventSelection;
   TGCheckButton *fRevents;
   TGRadioButton *fRcharge[3];
   TGRadioButton *fRtracks[2];
   TGLabel       *fLabel;
   TGLabel       *fLabel1;
   TGNumberEntry *fThreshold;
   TGNumberEntry *fMaxPt;
   TGComboBox    *fSelect;
   TGTextButton  *fButtonPrev;
   TGTextButton  *fButtonNext;

   // All events GUI
   TGButtonGroup *fParticlesBGAE;
   TGButtonGroup *fTrackSelectionAE;
   TGGroupFrame *fEventSelectionAE;
   TGCheckButton *fReventsAE;
   TGRadioButton *fRchargeAE[3];
   TGRadioButton *fRtracksAE[2];
   TGLabel       *fLabelAE;
   TGLabel       *fLabel1AE;
   TGNumberEntry *fThresholdAE;
   TGNumberEntry *fMaxPtAE;
   TGComboBox    *fSelectAE;
   TGTextButton  *fButtonPrevAE;
   TGTextButton  *fButtonNextAE;


private:

   AliEveLegoEditor(const AliEveLegoEditor&);            // Not implemented
   AliEveLegoEditor& operator=(const AliEveLegoEditor&); // Not implemented

   ClassDef(AliEveLegoEditor, 0); // GUI editor for AliEveLego.
};

#endif
