// $Id$
// Author: Stefano Carrazza 2010, CERN, stefano.carrazza@cern.ch

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVELEGOEDITOR_H
#define ALIEVELEGOEDITOR_H

#include "TGedFrame.h"

class TGButton;
class TGButtonGroup;
class TGCheckButton;
class TGColorSelect;
class TGComboBox;
class TGGroupFrame;
class TGLabel;
class TGNumberEntry;
class TGRadioButton;

class AliEveLego;

//______________________________________________________________________________
// AliEveLegoEditor is the class editor of AliEveLego
//

class AliEveLegoEditor : public TGedFrame
{
public:
   AliEveLegoEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
         UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
   virtual ~AliEveLegoEditor() {}

   virtual void SetModel(TObject* obj);

   // Slot methods
   void ApplyChanges();
   void CollisionCandidatesOnly();
   void CreateAllEventsEditor();
   void DataIsMC();
   void DoAllEvents();
   void ShowPosCharge();
   void ShowNegCharge();
   void ShowElectrons();
   void ShowMuons();
   void ShowPions();
   void ShowKaons();
   void ShowProtons();
   void ShowPosChargeAE();
   void ShowNegChargeAE();
   void ShowElectronsAE();
   void ShowMuonsAE();
   void ShowPionsAE();
   void ShowKaonsAE();
   void ShowProtonsAE();
   void SetThreshold();
   void SetThresholdAE();
   void SetMaxPt();
   void SetMaxPtAE();
   void ShowByTracks(Int_t id);
   void ShowByTracksAE(Int_t id);

protected:
   AliEveLego    *fM;                // Model object.

private:
   // Single event GUI
   TGTextButton  *fAllEventsButton;  // text button for all events
   TGGroupFrame  *fParticlesBG;      // particle selection button
   TGButtonGroup *fTrackSelection;   // track selection button
   TGCheckButton *fPosCharged;       // check button for positive only charged particles
   TGCheckButton *fNegCharged;       // check button for negative only charged particles
   TGCheckButton *fElectrons;        // check button for electrons
   TGCheckButton *fMuons;            // check button for muons
   TGCheckButton *fPions;            // check button for pions
   TGCheckButton *fKaons;            // check button for kaons
   TGCheckButton *fProtons;          // check button for protons
   TGRadioButton *fRtracks[2];       // radio button for track selection
   TGLabel       *fLabel;            // label for track selection
   TGLabel       *fLabel1;           // label for event selection
   TGNumberEntry *fThreshold;        // number entry to setup threshold
   TGNumberEntry *fMaxPt;            // number entry to setup max pt
   TGComboBox    *fSelect;           // combo box to filter events

   // All events GUI
   TGButtonGroup *fParticlesBGAE;    // particle selection button for all events
   TGButtonGroup *fTrackSelectionAE; // track selection for all events
   TGCheckButton *fPosChargedAE;     // check button for positive only charged particles
   TGCheckButton *fNegChargedAE;     // check button for negative only charged particles
   TGCheckButton *fElectronsAE;      // check button for electrons
   TGCheckButton *fMuonsAE;          // check button for muons
   TGCheckButton *fPionsAE;          // check button for pions
   TGCheckButton *fKaonsAE;          // check button for kaons
   TGCheckButton *fProtonsAE;        // check button for protons
   TGTextButton  *fApplyChanges;     // apply selections
   TGRadioButton *fRtracksAE[2];     // radio button track
   TGLabel       *fLabelAE;          // label for track selection
   TGLabel       *fLabel1AE;         // label for event selection
   TGNumberEntry *fThresholdAE;      // number entry to setup threshold
   TGNumberEntry *fMaxPtAE;          // number entry to setup max pt
   TGButtonGroup *fEventControl;     // event control panel
   TGCheckButton *fIsMC;             // check if data is from MC
   TGCheckButton *fCollisionCandidatesOnly; // fill all collision candidates events

   AliEveLegoEditor(const AliEveLegoEditor&);            // Not implemented
   AliEveLegoEditor& operator=(const AliEveLegoEditor&); // Not implemented

   ClassDef(AliEveLegoEditor, 0); // GUI editor for AliEveLego.
};

#endif
