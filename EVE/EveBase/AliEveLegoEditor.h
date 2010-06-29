// $Id$
// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVELEGOEDITOR_H
#define ALIEVELEGOEDITOR_H

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
   TGTextButton  *fAllEventsButton; // text button for all events
   TGButtonGroup *fParticlesBG;     // particle selection button
   TGButtonGroup *fTrackSelection;  // track selection button
   TGGroupFrame *fEventSelection;   // event selection button
   TGCheckButton *fRevents;         // check button for events
   TGRadioButton *fRcharge[3];      // radio button for charge selection
   TGRadioButton *fRtracks[2];      // radio button for track selection
   TGLabel       *fLabel;           // label for track selection
   TGLabel       *fLabel1;          // label for event selection
   TGNumberEntry *fThreshold;       // number entry to setup threshold
   TGNumberEntry *fMaxPt;           // number entry to setup max pt
   TGComboBox    *fSelect;          // combo box to filter events
   TGTextButton  *fButtonPrev;      // previous event selection button
   TGTextButton  *fButtonNext;      // next event selection button

   // All events GUI
   TGButtonGroup *fParticlesBGAE;    // particle selection button for all events
   TGButtonGroup *fTrackSelectionAE; // track selection for all events
   TGGroupFrame  *fEventSelectionAE; // event selection for all events
   TGCheckButton *fReventsAE;        // check button for event selection
   TGRadioButton *fRchargeAE[3];     // radio button event
   TGRadioButton *fRtracksAE[2];     // radio button track
   TGLabel       *fLabelAE;          // label for track selection
   TGLabel       *fLabel1AE;         // label for event selection
   TGNumberEntry *fThresholdAE;      // number entry to setup threshold
   TGNumberEntry *fMaxPtAE;          // number entry to setup max pt
   TGComboBox    *fSelectAE;         // combo box to filter events
   TGTextButton  *fButtonPrevAE;     // previous event selection button
   TGTextButton  *fButtonNextAE;     // next event selection button


private:

   AliEveLegoEditor(const AliEveLegoEditor&);            // Not implemented
   AliEveLegoEditor& operator=(const AliEveLegoEditor&); // Not implemented

   ClassDef(AliEveLegoEditor, 0); // GUI editor for AliEveLego.
};

#endif
