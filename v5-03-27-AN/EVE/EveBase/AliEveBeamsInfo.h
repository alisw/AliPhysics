// $Id$
// Author: Stefano Carrazza 2010, CERN, stefano.carrazza@cern.ch

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVEBEAMSINFO_H
#define ALIEVEBEAMSINFO_H

#include "TEveElement.h"

class AliESDEvent;
class AliPhysicsSelection;
class AliEveEventSelector;
class AliEveMultiView;

class TGLOverlayButton;
class TEveViewer;

//______________________________________________________________________________
// Display beams and triggers information on AliEve viewers, from ESD.
//

class AliEveBeamsInfo : public TEveElementList
{
public:
  AliEveBeamsInfo(const char* name="AliEveBeamsInfo");
  virtual ~AliEveBeamsInfo();
  void ShowEventSelection(Bool_t status);

  // Set Methods
  void SetAlpha(Double_t val);

  // Get Methods
  TString * SepareTriggerClasses(Int_t &fNumberOfClasses, TString fTriggerSource);

  // Functions
  void AddOverlayButton(TGLOverlayButton *button);
  void AddTriggerClasses();
  void CreateEventPanel();
  void CreateRunPanel();
  void RemoveOverlayButton(TGLOverlayButton *button);
  void RemoveTriggerClasses();
  void SelectEventSelection(Int_t id);
  void ShowBeamsInfo(Bool_t show, Bool_t updateonly = kFALSE);
  void ShowPrevEvent();
  void ShowNextEvent();
  void SwitchDataType(Bool_t status);
  void Update();
  void UpdateTriggerClasses();  

private:
  Double_t            fAlpha;                // Alpha value
  Bool_t              fIsMC;                 // Check data type
  AliESDEvent         *fEsd;                 // Esd event
  Bool_t              fShowEventsInfo;       // Determine if show events info
  AliPhysicsSelection *fPhysicsSelection;    // Physics selection object

  TGLOverlayButton    *fEventNumber;         // Event number
  TGLOverlayButton    *fCollisionCandidate;  // AliPhysicsSelection button output
  TGLOverlayButton    *fCollisionBoolean;    // Collision boolean

  TGLOverlayButton    *fBeam1;               // Beam 1 information
  TGLOverlayButton    *fBeam1Boolean;        // Beam 1 boolean
  TGLOverlayButton    *fBeam2;               // Beam 2 information
  TGLOverlayButton    *fBeam2Boolean;        // Beam 2 boolean

  TGLOverlayButton    *fRunNumber;           // Show data run number
  TGLOverlayButton    *fEventType;           // Show event type
  TGLOverlayButton    *fEventTypeLabel;      // Show event type label
  TGLOverlayButton    *fPeriod;              // Show event period
  TGLOverlayButton    *fOrbit;               // Show orbit
  TGLOverlayButton    *fBC;                  // Show bc

  TGLOverlayButton    *fTimeStamp;           // Time stamp information
  TGLOverlayButton    *fMagnetField;         // Magnetic field
  TGLOverlayButton    *fTrigger;             // Trigger

  TGLOverlayButton    *fTriggerClassesPanel; // Active trigger classes panel
  Int_t               fNumberOfActiveTriggerClasses; // Number of active trigger classes
  TGLOverlayButton    **fTriggerClasses;      // Active trigger classes

  AliEveMultiView     *fAl;                   // Multiview instance
  TEveViewer          *fHisto2dv;             // 2D lego view
  AliEveEventSelector *fEventSelector;        // Current event selector

  AliEveBeamsInfo(const AliEveBeamsInfo&);            // Not implemented
  AliEveBeamsInfo& operator=(const AliEveBeamsInfo&); // Not implemented

  ClassDef(AliEveBeamsInfo, 0); // Short description.
};

#endif
