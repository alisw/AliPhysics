// $Id$
// Author: Stefano Carrazza 2010

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
// Display beams and triggers information on gl-viewers.
//

class AliEveBeamsInfo : public TEveElementList
{
public:
  AliEveBeamsInfo(const char* name="AliEveBeamsInfo");
  virtual ~AliEveBeamsInfo();

  // Set Methods
  void ShowEventSelection();

  // Functions
  void ShowBeamsInfo(Bool_t show, Bool_t updateonly = kFALSE);
  void Update();
  void SelectEventSelection(Int_t id);
  void ShowPrevEvent();
  void ShowNextEvent();

private:
  AliESDEvent *fEsd;                      // esd event
  Bool_t fShowEventsInfo;                 // determine if show events info
  AliPhysicsSelection *fPhysicsSelection; // physics selection object
  TGLOverlayButton *fCollisionCandidate;  // AliPhysicsSelection button output
  TGLOverlayButton *fBeam1;               // beam 1 information
  TGLOverlayButton *fBeam2;               // beam 2 information
  AliEveMultiView  *fAl;                  // multiview instance
  TEveViewer *fHisto2dv;                  // 2D lego view
  AliEveEventSelector *fEventSelector;    // current event selector

  AliEveBeamsInfo(const AliEveBeamsInfo&);            // Not implemented
  AliEveBeamsInfo& operator=(const AliEveBeamsInfo&); // Not implemented

  ClassDef(AliEveBeamsInfo, 0); // Short description.
};

#endif
