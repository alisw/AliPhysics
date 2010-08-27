// $Id$
// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveBeamsInfo.h"
#include "AliPhysicsSelection.h"
#include "AliESDEvent.h"
#include "AliEveEventManager.h"
#include "AliEveEventSelector.h"
#include "AliEveMultiView.h"

#include "TEveWindow.h"
#include "TEveManager.h"
#include "TEveBrowser.h"
#include "TEveViewer.h"
#include "TEveScene.h"
#include "TGLOverlayButton.h"


//______________________________________________________________________________
// This class provides the triggers information about beams.
//

ClassImp(AliEveBeamsInfo)

//______________________________________________________________________________
AliEveBeamsInfo::AliEveBeamsInfo(const char* name) :
  TEveElementList(name),
  fEsd(0),
  fShowEventsInfo(kTRUE),
  fPhysicsSelection(0),
  fCollisionCandidate(0),
  fBeam1(0),
  fBeam2(0),
  fAl(0),
  fHisto2dv(0),
  fEventSelector(0)
{
  // Constructor.
  gEve->AddToListTree(this,0);

  // get current ESD event
  fEsd = AliEveEventManager::AssertESD();
  fAl = AliEveMultiView::Instance();  
  fEventSelector = AliEveEventManager::GetMaster()->GetEventSelector();

  // AliPhysicsSelection
  fPhysicsSelection = new AliPhysicsSelection();
  fPhysicsSelection->SetAnalyzeMC(kFALSE);
  fPhysicsSelection->Initialize(fEsd->GetRunNumber());

  // loading physics selection and triggers buttons
  fCollisionCandidate = new TGLOverlayButton(0, "", 10.0, -10.0, 190.0, 20.0);
  fCollisionCandidate->SetAlphaValues(1.5,1.5);
  fBeam1 = new TGLOverlayButton(0, "", 10.0, -30.0, 190.0, 20.0);
  fBeam1->SetAlphaValues(1.5,1.5);
  fBeam2 = new TGLOverlayButton(0, "", 10.0, -50.0, 190.0, 20.0);
  fBeam2->SetAlphaValues(1.5,1.5);

  // show beams info
  ShowBeamsInfo(fShowEventsInfo, kFALSE);

}

//______________________________________________________________________________
AliEveBeamsInfo::~AliEveBeamsInfo()
{
  // deleting variables
  delete fEsd;
  delete fPhysicsSelection;
  delete fCollisionCandidate;
  delete fBeam1;
  delete fBeam2;
  delete fAl;
  delete fHisto2dv;

}

void AliEveBeamsInfo::ShowBeamsInfo(Bool_t show, Bool_t updateonly)
{
  // Collision candidate
  fHisto2dv = (TEveViewer*) gEve->GetViewers()->FindChild("2D Lego Viewer");

  if (show == 0)
  {

    gEve->GetDefaultGLViewer()->RemoveOverlayElement(fCollisionCandidate);
    fAl->Get3DView()->GetGLViewer()->RemoveOverlayElement(fCollisionCandidate);
    if(fHisto2dv)
      fHisto2dv->GetGLViewer()->RemoveOverlayElement(fCollisionCandidate);

    gEve->GetDefaultGLViewer()->RemoveOverlayElement(fBeam1);
    fAl->Get3DView()->GetGLViewer()->RemoveOverlayElement(fBeam1);
    if(fHisto2dv)
      fHisto2dv->GetGLViewer()->RemoveOverlayElement(fBeam1);

    gEve->GetDefaultGLViewer()->RemoveOverlayElement(fBeam2);
    fAl->Get3DView()->GetGLViewer()->RemoveOverlayElement(fBeam2);
    if(fHisto2dv)
      fHisto2dv->GetGLViewer()->RemoveOverlayElement(fBeam2);

  } else {

  if (updateonly == kFALSE) {
    gEve->GetDefaultGLViewer()->AddOverlayElement(fCollisionCandidate);
    fAl->Get3DView()->GetGLViewer()->AddOverlayElement(fCollisionCandidate);
    if(fHisto2dv)
      fHisto2dv->GetGLViewer()->AddOverlayElement(fCollisionCandidate);
  }

  Bool_t ev = fPhysicsSelection->IsCollisionCandidate(fEsd);

  if (ev == 1)
  {
     fCollisionCandidate->SetText("Collision candidate: YES");
  } else {
     fCollisionCandidate->SetText("Collision candidate: NO ");
  }

  // Beam 1 & 2 setup: method 1
  if (updateonly == kFALSE) {
     gEve->GetDefaultGLViewer()->AddOverlayElement(fBeam1);
     fAl->Get3DView()->GetGLViewer()->AddOverlayElement(fBeam1);
     if(fHisto2dv)
      fHisto2dv->GetGLViewer()->AddOverlayElement(fBeam1);

     gEve->GetDefaultGLViewer()->AddOverlayElement(fBeam2);
     fAl->Get3DView()->GetGLViewer()->AddOverlayElement(fBeam2);
     if(fHisto2dv)
      fHisto2dv->GetGLViewer()->AddOverlayElement(fBeam2);
  }

  Bool_t b1  = fEsd->IsTriggerClassFired("CINT1A-ABCE-NOPF-ALL");
  Bool_t b2  = fEsd->IsTriggerClassFired("CINT1C-ABCE-NOPF-ALL");
  Bool_t b12 = fEsd->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL");

  if (b1 == 1 || b12 == 1)
  {
     fBeam1->SetText("Beam 1: YES");
     fBeam1->SetBackColor(0x00ff00);
  } else {
     fBeam1->SetText("Beam 1: NO");
     fBeam1->SetBackColor(0xff0000);
  }

  if (b2 == 1 || b12 == 1)
  {
     fBeam2->SetText("Beam 2: YES");
     fBeam2->SetBackColor(0x00ff00);
  } else {
     fBeam2->SetText("Beam 2: NO");
     fBeam2->SetBackColor(0xff0000);
  }
}

gEve->Redraw3D(kTRUE);

}

//______________________________________________________________________________
void AliEveBeamsInfo::Update()
{
  // update beams information for current event
  ShowBeamsInfo(fShowEventsInfo, kTRUE);
}

//______________________________________________________________________________
void AliEveBeamsInfo::ShowEventSelection()
{
  // activate/deactivate info box
  if (fShowEventsInfo == 0)
  {
    fShowEventsInfo = 1;
  } else {
    fShowEventsInfo = 0;
  }

  ShowBeamsInfo(fShowEventsInfo);
}

//______________________________________________________________________________
void AliEveBeamsInfo::SelectEventSelection(Int_t id)
{
  // show trigger information
  if (id == 0)
  {
     fEventSelector->SetSelectOnTriggerType(kFALSE);
  } else {
     if (id == 1) fEventSelector->SetTriggerType("CINT1A-ABCE-NOPF-ALL");
     if (id == 2) fEventSelector->SetTriggerType("CINT1C-ABCE-NOPF-ALL");
     if (id == 3) fEventSelector->SetTriggerType("CINT1B-ABCE-NOPF-ALL");
     fEventSelector->SetSelectOnTriggerType(kTRUE);
  }
}

//______________________________________________________________________________
void AliEveBeamsInfo::ShowPrevEvent()
{
  // go to the previous event
  AliEveEventManager::GetMaster()->PrevEvent();
}

//______________________________________________________________________________
void AliEveBeamsInfo::ShowNextEvent()
{
  AliEveEventManager::GetMaster()->NextEvent();
}

/******************************************************************************/
