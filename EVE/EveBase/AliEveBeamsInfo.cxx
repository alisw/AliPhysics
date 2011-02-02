// $Id$
// Author: Stefano Carrazza 2010, CERN, stefano.carrazza@cern.ch

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
#include "AliRawEventHeaderBase.h"

#include "TEveWindow.h"
#include "TEveManager.h"
#include "TEveBrowser.h"
#include "TEveViewer.h"
#include "TEveScene.h"
#include "TGLOverlayButton.h"
#include "TTimeStamp.h"


//______________________________________________________________________________
// This class provides:
//
// 1) Run and event information from ESD as TGLOverlayButtons.
//
// 2) Determine for real and MC data if the event is collision candidate or not.
//
// 3) Possibility to filter events per type.
//

ClassImp(AliEveBeamsInfo)

//______________________________________________________________________________
AliEveBeamsInfo::AliEveBeamsInfo(const char* name) :
  TEveElementList(name),
  fAlpha(1.5),
  fIsMC(kFALSE),
  fEsd(0),
  fShowEventsInfo(kTRUE),
  fPhysicsSelection(0),
  fEventNumber(0),
  fCollisionCandidate(0),
  fCollisionBoolean(0),
  fBeam1(0),
  fBeam1Boolean(0),
  fBeam2(0),
  fBeam2Boolean(0),
  fRunNumber(0),
  fEventType(0),
  fEventTypeLabel(0),
  fPeriod(0),
  fOrbit(0),
  fBC(0),
  fTimeStamp(0),
  fMagnetField(0),
  fTrigger(0),
  fTriggerClassesPanel(0),
  fNumberOfActiveTriggerClasses(0),
  fTriggerClasses(0),
  fAl(0),
  fHisto2dv(0),
  fEventSelector(0)
{
  // Constructor.
  gEve->AddToListTree(this,0);

  // Get current ESD event
  fEsd = AliEveEventManager::AssertESD();
  fAl = AliEveMultiView::Instance();  
  fEventSelector = AliEveEventManager::GetMaster()->GetEventSelector();

  // AliPhysicsSelection
  fPhysicsSelection = new AliPhysicsSelection();
  fPhysicsSelection->SetAnalyzeMC(kFALSE);
  fPhysicsSelection->Initialize(fEsd);

  // Loading physics selection and triggers buttons
  CreateEventPanel();
  CreateRunPanel();

  // Show beams info
  ShowBeamsInfo(fShowEventsInfo, kFALSE);

}

//______________________________________________________________________________
AliEveBeamsInfo::~AliEveBeamsInfo()
{
  // Deleting variables
  RemoveTriggerClasses();
  delete fEsd;
  delete fPhysicsSelection;
  delete fEventNumber;
  delete fCollisionCandidate;
  delete fCollisionBoolean;
  delete fBeam1;
  delete fBeam1Boolean;
  delete fBeam2;
  delete fBeam2Boolean;
  delete fRunNumber;
  delete fPeriod;
  delete fOrbit;
  delete fBC;
  delete fTimeStamp;
  delete fMagnetField;
  delete fTrigger;
  delete fTriggerClassesPanel;  
  delete fEventType;
  delete fEventTypeLabel;
  delete fAl;
  delete fHisto2dv;
}

//______________________________________________________________________________
void AliEveBeamsInfo::CreateEventPanel()
{
  // Create vertical panel
  Double_t fPosx = 10.0;
  Double_t fPosy = -10.0;
  Double_t fXlengh = 220.0;
  Double_t fStep = 20.0; 

  fEventNumber = new TGLOverlayButton(0 , "", fPosx, fPosy, fXlengh, fStep);
  fPosy-=fStep;

  fCollisionCandidate = new TGLOverlayButton(0, "Collision:", fPosx, fPosy, fXlengh-fXlengh/2.0, fStep);
  fCollisionBoolean = new TGLOverlayButton(0, "", fPosx+fXlengh/2.0, fPosy, fXlengh/2.0, fStep);
  fPosy-=fStep;

  fBeam1 = new TGLOverlayButton(0, "Beam 1:", fPosx, fPosy, fXlengh-fXlengh/2.0, fStep);
  fBeam1Boolean = new TGLOverlayButton(0, "", fPosx+fXlengh/2.0, fPosy, fXlengh/2.0, fStep);
  fPosy-=fStep;

  fBeam2 = new TGLOverlayButton(0, "Beam 2:", fPosx, fPosy, fXlengh-fXlengh/2.0, fStep);
  fBeam2Boolean = new TGLOverlayButton(0, "", fPosx+fXlengh/2.0, fPosy, fXlengh/2.0, fStep);

  fPosy-=2*fStep;
  fTriggerClassesPanel = new TGLOverlayButton(0, "Active trigger classes:", fPosx, fPosy, fXlengh, fStep);
}

//______________________________________________________________________________
void AliEveBeamsInfo::CreateRunPanel()
{
  // Create horizontal panel
  Double_t fPosy = -10.0;
  Double_t fPosx = 250.0;
  Double_t fStep = 20.0;
  Double_t fXlengh = 120.0;
  TString fText;

  fText.Form("Run #: %d",fEsd->GetRunNumber());
  fRunNumber = new TGLOverlayButton(0, fText.Data(), fPosx, fPosy, fXlengh, fStep);
  fPosx+=fXlengh;

  fEventType = new TGLOverlayButton(0, "", fPosx, fPosy , fXlengh, fStep);
  fPosx+=fXlengh;

  fEventTypeLabel = new TGLOverlayButton(0, "", fPosx, fPosy, 2*fXlengh, fStep);

  // New line
  fPosx = 250.0;
  fPosy = -30.0;

  fPeriod = new TGLOverlayButton(0, "", fPosx, fPosy, fXlengh, fStep);
  fPosx+=fXlengh;

  fOrbit = new TGLOverlayButton(0, "", fPosx, fPosy, fXlengh, fStep);
  fPosx+=fXlengh;

  fBC = new TGLOverlayButton(0, "", fPosx, fPosy, fXlengh, fStep);
  fPosx+=fXlengh;

  fTrigger = new TGLOverlayButton(0, "", fPosx, fPosy, fXlengh, fStep);

  // New line
  fPosx = 250.0;
  fPosy = -50.0;

  fTimeStamp = new TGLOverlayButton(0, "", fPosx, fPosy, 2*fXlengh, fStep);
  fPosx+=2*fXlengh;

  fMagnetField = new TGLOverlayButton(0, "", fPosx, fPosy, 2*fXlengh, fStep);
  fPosx+=2*fXlengh;
}

//______________________________________________________________________________
void AliEveBeamsInfo::ShowBeamsInfo(Bool_t show, Bool_t updateonly)
{
  // Update & setup TGLOverlayButtons
  fHisto2dv = (TEveViewer*) gEve->GetViewers()->FindChild("2D Lego Viewer");

  if (!show)
  {
    RemoveOverlayButton(fEventNumber);
    RemoveOverlayButton(fCollisionCandidate);
    RemoveOverlayButton(fCollisionBoolean);
    RemoveOverlayButton(fBeam1);
    RemoveOverlayButton(fBeam1Boolean);
    RemoveOverlayButton(fBeam2);
    RemoveOverlayButton(fBeam2Boolean);

    RemoveOverlayButton(fRunNumber);
    RemoveOverlayButton(fEventType);
    RemoveOverlayButton(fEventTypeLabel);
    RemoveOverlayButton(fPeriod);
    RemoveOverlayButton(fOrbit);
    RemoveOverlayButton(fBC);
    RemoveOverlayButton(fTimeStamp);
    RemoveOverlayButton(fMagnetField);
    RemoveOverlayButton(fTrigger);

    RemoveOverlayButton(fTriggerClassesPanel);
    RemoveTriggerClasses();
  } else {

  if (!updateonly)
  {
    AddOverlayButton(fEventNumber);
    AddOverlayButton(fCollisionCandidate);
    AddOverlayButton(fCollisionBoolean);
    AddOverlayButton(fBeam1);
    AddOverlayButton(fBeam1Boolean);
    AddOverlayButton(fBeam2);
    AddOverlayButton(fBeam2Boolean);

    AddOverlayButton(fRunNumber);
    AddOverlayButton(fEventType);
    AddOverlayButton(fEventTypeLabel);
    AddOverlayButton(fPeriod);
    AddOverlayButton(fOrbit);
    AddOverlayButton(fBC);
    AddOverlayButton(fTimeStamp);
    AddOverlayButton(fMagnetField);
    AddOverlayButton(fTrigger);

    AddOverlayButton(fTriggerClassesPanel);
    AddTriggerClasses();
  }

  TString fText;
  fText.Form("Event #: %d", fEsd->GetEventNumberInFile());
  fEventNumber->SetText(fText.Data());

  fText.Form("Event type: %d",fEsd->GetEventType());
  fEventType->SetText(fText.Data());

  if(fEsd->GetEventType() == 0)
  {
    fText.Form("UNKNOW EVENT TYPE");
  } else {
    fText.Form("%s", AliRawEventHeaderBase::GetTypeName(fEsd->GetEventType()));
  }
  fEventTypeLabel->SetText(fText.Data());

  fText.Form("Period: %x", fEsd->GetPeriodNumber());
  fPeriod->SetText(fText.Data());

  fText.Form("Orbit: %x", fEsd->GetOrbitNumber());
  fOrbit->SetText(fText.Data());

  fText.Form("BC: %x", fEsd->GetBunchCrossNumber());
  fBC->SetText(fText.Data());

  TTimeStamp ts(fEsd->GetTimeStamp());
  fText.Form("Timestamp: %s",ts.AsString("s"));
  fTimeStamp->SetText(fText.Data());

  fText.Form("Magnetic field: %.2e kG", fEsd->GetMagneticField());
  fMagnetField->SetText(fText.Data());

  if (fEsd->GetTriggerMask() > (ULong64_t) 100 && fIsMC == kFALSE)
  {
    fText.Form("Trigger: #");
  } else {    
    fText.Form("Trigger: %llx", fEsd->GetTriggerMask());
  }
  fTrigger->SetText(fText.Data());

  UpdateTriggerClasses();

  Bool_t ev = fPhysicsSelection->IsCollisionCandidate(fEsd);

  if (ev == 1)
  {
     fCollisionBoolean->SetText("YES");
     fCollisionBoolean->SetBackColor(0x00ff00);
  } else {
     fCollisionBoolean->SetText("NO");
     fCollisionBoolean->SetBackColor(0xff0000);
  }

  Bool_t b1  = fEsd->IsTriggerClassFired("CINT1A-ABCE-NOPF-ALL");
  Bool_t b2  = fEsd->IsTriggerClassFired("CINT1C-ABCE-NOPF-ALL");
  Bool_t b12 = fEsd->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL");

  if (b1 == 1 || b12 == 1)
  {
     fBeam1Boolean->SetText("YES");
     fBeam1Boolean->SetBackColor(0x00ff00);
  } else {
     fBeam1Boolean->SetText("NO");
     fBeam1Boolean->SetBackColor(0xff0000);
  }

  if (b2 == 1 || b12 == 1)
  {
     fBeam2Boolean->SetText("YES");
     fBeam2Boolean->SetBackColor(0x00ff00);
  } else {
     fBeam2Boolean->SetText("NO");
     fBeam2Boolean->SetBackColor(0xff0000);
  }
}

gEve->Redraw3D(kTRUE);

}

//______________________________________________________________________________
void AliEveBeamsInfo::Update()
{
  // Update beams information for current event
  ShowBeamsInfo(fShowEventsInfo, kTRUE);
}

//______________________________________________________________________________
void AliEveBeamsInfo::AddOverlayButton(TGLOverlayButton *button)
{
  // Add buttons to viewers
  button->SetAlphaValues(fAlpha, fAlpha);

  gEve->GetDefaultGLViewer()->AddOverlayElement(button);
  fAl->Get3DView()->GetGLViewer()->AddOverlayElement(button);
  if(fHisto2dv)
    fHisto2dv->GetGLViewer()->AddOverlayElement(button);
}

//______________________________________________________________________________
void AliEveBeamsInfo::RemoveOverlayButton(TGLOverlayButton *button)
{
  // Remove buttons to viewers
  gEve->GetDefaultGLViewer()->RemoveOverlayElement(button);
  fAl->Get3DView()->GetGLViewer()->RemoveOverlayElement(button);
  if(fHisto2dv)
    fHisto2dv->GetGLViewer()->RemoveOverlayElement(button);
}

//______________________________________________________________________________
void AliEveBeamsInfo::ShowEventSelection(Bool_t status)
{
  // Activate/deactivate info box
  fShowEventsInfo = status;
  ShowBeamsInfo(fShowEventsInfo);
}

//______________________________________________________________________________
void AliEveBeamsInfo::SelectEventSelection(Int_t id)
{
  // Show trigger information
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
  // Go to the previous event
  AliEveEventManager::GetMaster()->PrevEvent();
}

//______________________________________________________________________________
void AliEveBeamsInfo::ShowNextEvent()
{
  // Go to the next event
  AliEveEventManager::GetMaster()->NextEvent();
}

//______________________________________________________________________________
void AliEveBeamsInfo::SwitchDataType(Bool_t status)
{
  // Activate/deactivate MC / real data type
  fIsMC = status;

  // Removing defaul physics selection
  delete fPhysicsSelection;
  fPhysicsSelection = NULL;

  fPhysicsSelection = new AliPhysicsSelection();
  fPhysicsSelection->SetAnalyzeMC(fIsMC);
  fPhysicsSelection->Initialize(fEsd);
  Update();
}

//______________________________________________________________________________
void AliEveBeamsInfo::AddTriggerClasses()
{
  // Add trigger classes
  if (!fTriggerClasses){

    TString fTriggerNameString = fEsd->GetESDRun()->GetActiveTriggerClasses();
    TString *fTriggerName = SepareTriggerClasses(fNumberOfActiveTriggerClasses, fTriggerNameString);

    fTriggerClasses = new TGLOverlayButton*[fNumberOfActiveTriggerClasses];

    Double_t fPosy = -130.0;
    Double_t fStep = 20.0;
    for (Int_t i = 0; i < fNumberOfActiveTriggerClasses; i++)
    {
      fTriggerClasses[i] = new TGLOverlayButton(0, fTriggerName[i].Data(), 10, fPosy, 220.0, 20);
      fTriggerClasses[i]->SetBackColor(0xce970a);
      fTriggerClasses[i]->SetAlphaValues(fAlpha, fAlpha);
      fPosy-=fStep;

      AddOverlayButton(fTriggerClasses[i]);
    }
  }
}

//______________________________________________________________________________
void AliEveBeamsInfo::RemoveTriggerClasses()
{
  // Remove overlay buttons
  for(Int_t i = 0; i < fNumberOfActiveTriggerClasses; i++)
  {
    RemoveOverlayButton(fTriggerClasses[i]);
  }

  delete[] fTriggerClasses;
  fTriggerClasses = 0;
}

//______________________________________________________________________________
void AliEveBeamsInfo::UpdateTriggerClasses()
{
  // Remove trigger information and update it
  RemoveTriggerClasses();
  AddTriggerClasses();
}

//______________________________________________________________________________
TString *AliEveBeamsInfo::SepareTriggerClasses(Int_t &fNumberOfClasses, TString fTriggerSource)
{
  // Get trigger string and separe triggers into TString's
  Int_t fStringLength = fTriggerSource.Length();
  fNumberOfClasses = 0;

  for (Int_t i = 1; i < fStringLength; i++)
  {
    TString fString = fTriggerSource(i,2);
    if (fString == "  "){
      fNumberOfClasses++;
    }
  }
  fNumberOfClasses++;

  TString *fTriggerResult = new TString[fNumberOfClasses];

  Int_t fIndex = 1;
  Int_t fLastIndex = 0;
  Int_t fFirstIndex = 1;
  Int_t fClassNumber = 0;
  for(;;)
  {
    if(fIndex >= fStringLength) break;
    TString fString = fTriggerSource(fIndex,1);
    fIndex++;
    fLastIndex++;
    if (fString == " "){
      fTriggerResult[fClassNumber] = fTriggerSource(fFirstIndex, fLastIndex-1);
      fFirstIndex = fIndex+1;
      fLastIndex = 0;
      fIndex++;
      fClassNumber++;
    }
  }
  return fTriggerResult;
}

//______________________________________________________________________________
void AliEveBeamsInfo::SetAlpha(Double_t val)
{
  // Set the new alpha value for TGLOverlayButton
  fAlpha = val;

  // First remove all buttons from viewers
  ShowBeamsInfo(kFALSE);

  // Then replot everything with the next alpha value
  ShowBeamsInfo(kTRUE);
}

/******************************************************************************/
