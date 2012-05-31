// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCLoaderEditor.h"

#include <EveDet/AliEveTPCLoader.h>
#include <EveDet/AliEveTPCData.h>

#include <TEveScene.h>
#include <TEveManager.h>
#include <TEveGValuators.h>
#include <TGDoubleSlider.h>

#include <TSystem.h>
#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGFileDialog.h>
#include <TGToolTip.h>


//______________________________________________________________________________
//
// Editor for AliEveTPCLoader.

ClassImp(AliEveTPCLoaderEditor)

AliEveTPCLoaderEditor::AliEveTPCLoaderEditor(const TGWindow *p,
					     Int_t width, Int_t height,
					     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),

  fM        (0),

  fFile     (0),
  fOpenFile (0),

  fEvent    (0),
  fDoubleSR (0),

  fDataLoadThreshold (0),
  fDataLoadPedestal  (0),
  fDataAutoPedestal  (0),

  fUpdateSectors   (0),
  fReloadSectors   (0),
  fCreateSectors3D (0),
  fDeleteSectors3D (0),

  gEtaRange(0),
  gCutOnEta(0)
  
{
  // Constructor.

  MakeTitle("AliEveTPCLoader");

  Int_t labelW;

  // File / event interface

  labelW = 42;
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGHorizontalFrame* g = new TGHorizontalFrame(f, labelW, 0, kFixedWidth);
    TGLabel* l = new TGLabel(g, "File: ");
    g->AddFrame(l, new TGLayoutHints(kLHintsLeft, 0,0,4,0));
    f->AddFrame(g);
    fFile = new TGTextEntry(f);
    fFile->SetWidth(140);
    f->AddFrame(fFile);
    fFile->Connect("DoubleClicked()",
		   "AliEveTPCLoaderEditor", this, "FileSelect()");
    fFile->Connect("TextChanged(const char *)",
		   "AliEveTPCLoaderEditor", this, "FileChanged()");
    fOpenFile = new TGTextButton(f, "Open");
    f->AddFrame(fOpenFile);
    fOpenFile->Connect("Clicked()",
		       "AliEveTPCLoaderEditor", this, "DoOpen()");
    AddFrame(f);
  }

  fEvent = new TEveGValuator(this, "Event:", 110, 0);
  fEvent->SetShowSlider(kFALSE);
  fEvent->SetLabelWidth(labelW);
  fEvent->SetNELength(6);
  fEvent->Build();
  fEvent->SetLimits(0, 999999);
  fEvent->SetToolTip("Current event number");
  fEvent->Connect("ValueSet(Double_t)",
		  "AliEveTPCLoaderEditor", this, "DoEvent()");
  // Reuse Event for DoubleSR button
  fDoubleSR = new TGCheckButton(fEvent, "Double SR");
  fDoubleSR->SetToolTipText("Double sampling rate");
  fEvent->AddFrame(fDoubleSR, new TGLayoutHints(kLHintsLeft, 12, 0, 2, 0));
  fDoubleSR->Connect("Toggled(Bool_t)",
                     "AliEveTPCLoaderEditor", this, "DoDoubleSR()");
  AddFrame(fEvent);


  // AliEveTPCData load settings

  labelW = 90;

  fDataLoadThreshold = new TEveGValuator(this, "Load threshold:", 110, 0);
  fDataLoadThreshold->SetShowSlider(kFALSE);
  fDataLoadThreshold->SetLabelWidth(labelW);
  fDataLoadThreshold->SetNELength(6);
  fDataLoadThreshold->Build();
  fDataLoadThreshold->SetLimits(0, 1000);
  fDataLoadThreshold->SetToolTip("Minimum signal that will be stored (pedestal subtracted first).");
  fDataLoadThreshold->Connect
    ("ValueSet(Double_t)",
     "AliEveTPCLoaderEditor", this, "DoDataLoadThreshold()");
  AddFrame(fDataLoadThreshold, new TGLayoutHints(kLHintsLeft, 0, 0, 6, 0));

  fDataLoadPedestal = new TEveGValuator(this, "Load pedestal:", 110, 0);
  fDataLoadPedestal->SetShowSlider(kFALSE);
  fDataLoadPedestal->SetLabelWidth(labelW);
  fDataLoadPedestal->SetNELength(6);
  fDataLoadPedestal->Build();
  fDataLoadPedestal->SetLimits(0, 1000);
  fDataLoadPedestal->SetToolTip("Pedestal that will be subtracted during data loading.");
  fDataLoadPedestal->Connect
    ("ValueSet(Double_t)",
     "AliEveTPCLoaderEditor", this, "DoDataLoadPedestal()");
  // Reuse DataLoadPedestal for DataAutoPedestal check button
  fDataAutoPedestal = new TGCheckButton(fDataLoadPedestal, "Automatic");
  fDataAutoPedestal->SetToolTipText("Determine per-pad pedestal during data loading.");
  fDataAutoPedestal->Connect
    ("Toggled(Bool_t)",
     "AliEveTPCLoaderEditor", this, "DoDataAutoPedestal()");
  fDataLoadPedestal->AddFrame(fDataAutoPedestal, new TGLayoutHints(kLHintsLeft, 12, 0, 2, 0));
  AddFrame(fDataLoadPedestal);

  // Steering buttons: update/reload sectors; show/hide 3d

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fUpdateSectors = new TGTextButton(f, "Update Sectors");
    f->AddFrame(fUpdateSectors, new TGLayoutHints(kLHintsExpandX, 0,4,0,0));
    fUpdateSectors->Connect("Clicked()",
			    "AliEveTPCLoaderEditor", this, "DoUpdateSectors()");
    fReloadSectors = new TGTextButton(f, "Reload Sectors");
    f->AddFrame(fReloadSectors, new TGLayoutHints(kLHintsExpandX, 4,0,0,0));
    fReloadSectors->Connect("Clicked()",
			    "AliEveTPCLoaderEditor", this, "DoReloadSectors()");
    AddFrame(f, new TGLayoutHints(kLHintsExpandX, 8,8,8,0));
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fCreateSectors3D = new TGTextButton(f, "Create 3D");
    f->AddFrame(fCreateSectors3D, new TGLayoutHints(kLHintsExpandX, 0,4,0,0));
    fCreateSectors3D->Connect("Clicked()",
			      "AliEveTPCLoaderEditor", this, "DoCreateSectors3D()");
    fDeleteSectors3D = new TGTextButton(f, "Delete 3D");
    f->AddFrame(fDeleteSectors3D, new TGLayoutHints(kLHintsExpandX, 4,0,0,0));
    fDeleteSectors3D->Connect("Clicked()",
			      "AliEveTPCLoaderEditor", this, "DoDeleteSectors3D()");
    AddFrame(f, new TGLayoutHints(kLHintsExpandX, 8,8,8,0));
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fCreateSectors3D = new TGTextButton(f, "Show 2D");
    f->AddFrame(fCreateSectors3D, new TGLayoutHints(kLHintsExpandX, 0,4,0,0));
    fCreateSectors3D->Connect("Clicked()",
			      "AliEveTPCLoaderEditor", this, "DoShowSectors2D()");
    fDeleteSectors3D = new TGTextButton(f, "Hide 2D");
    f->AddFrame(fDeleteSectors3D, new TGLayoutHints(kLHintsExpandX, 4,0,0,0));
    fDeleteSectors3D->Connect("Clicked()",
			      "AliEveTPCLoaderEditor", this, "DoHideSectors2D()");
    AddFrame(f, new TGLayoutHints(kLHintsExpandX, 8,8,8,0));
  }

   // Eta cuts slider

  {
   TGHorizontalFrame* f = new TGHorizontalFrame(this);

   gEtaRange = new TEveGDoubleValuator(f,"Eta range:", 40, 0);
   gEtaRange->SetNELength(6);
   gEtaRange->SetLabelWidth(50);
   gEtaRange->Build();
   gEtaRange->GetSlider()->SetWidth(180);
   gEtaRange->SetLimits(-1.5, 1.5, TGNumberFormat::kNESRealTwo);
   gEtaRange->SetValues(-1.5, 1.5, TGNumberFormat::kNESRealTwo);

   gCutOnEta = new TGCheckButton(f, "Set", 10);
   gCutOnEta->SetEnabled(kTRUE);

   f->AddFrame(gEtaRange, new TGLayoutHints(kLHintsExpandX, 10, 10, 10, 10));
   f->AddFrame(gCutOnEta, new TGLayoutHints(kLHintsNormal, 10, 10, 10, 10));

   AddFrame(f, new TGLayoutHints(kLHintsExpandX, 8,8,8,0));
  }

}

/******************************************************************************/

void AliEveTPCLoaderEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = static_cast<AliEveTPCLoader*>(obj);

  // !!!! order changed, need TGTextEntry::SetText NO BLOODY EMIT.
  fFile->SetText(fM->fFile);
  fEvent->SetValue(fM->fEvent);
  fEvent->SetEnabled(fM->fEvent >= 0);
  fDoubleSR->SetState(fM->fDoubleSR ? kButtonDown : kButtonUp);

  AliEveTPCData* tpcd = fM->GetData();
  Bool_t   tpcp = (tpcd != 0);
  fDataLoadThreshold->SetEnabled(tpcp);
  fDataLoadPedestal ->SetEnabled(tpcp && ! tpcd->GetAutoPedestal());
  fDataAutoPedestal ->SetEnabled(tpcp);
  if (tpcp) {
    fDataLoadThreshold->SetValue(tpcd->GetLoadThreshold());
    fDataLoadPedestal ->SetValue(tpcd->GetLoadPedestal());
    fDataAutoPedestal ->SetState(tpcd->GetAutoPedestal() ? kButtonDown : kButtonUp);
  }
}

/******************************************************************************/
/******************************************************************************/

void AliEveTPCLoaderEditor::FileSelect()
{
  // Slot for FileSelect.

  static const char *kTPCFileTypes[] = {
   "Root files",  "*.root",
   "All files",   "*.*",
    0,               0
};

  TGFileInfo fi;
  fi.fIniDir    = StrDup(gSystem->DirName (fM->fFile));
  fi.fFilename  = StrDup(gSystem->BaseName(fM->fFile));
  fi.fFileTypes = kTPCFileTypes;

  new TGFileDialog(fClient->GetRoot(), gEve->GetMainWindow(), kFDOpen, &fi);
  if (!fi.fFilename)
    return;

  fFile->SetText(fi.fFilename);
}

void AliEveTPCLoaderEditor::FileChanged()
{
  // Slot for FileChanged.

  fM->fFile = fFile->GetText();
}

void AliEveTPCLoaderEditor::DoOpen()
{
  // Slot for Open.

  fM->OpenFile();
  SetModel(fM);
}

/******************************************************************************/

void AliEveTPCLoaderEditor::DoEvent()
{
  // Slot for Event.

  fM->GotoEvent((Int_t) fEvent->GetValue());
  SetModel(fM);
}

void AliEveTPCLoaderEditor::DoDoubleSR()
{
  // Slot for DoubleSR.

  fM->SetDoubleSR(fDoubleSR->IsOn());
  Update();
}

/******************************************************************************/

void AliEveTPCLoaderEditor::DoDataLoadThreshold()
{
  // Slot for DataLoadThreshold.

  if (fM->GetData() == 0) return;
  fM->GetData()->SetLoadThreshold((Short_t) fDataLoadThreshold->GetValue());
}

void AliEveTPCLoaderEditor::DoDataLoadPedestal()
{
  // Slot for DataLoadPedestal.

  if (fM->GetData() == 0) return;
  fM->GetData()->SetLoadPedestal((Short_t) fDataLoadPedestal->GetValue());
}

void AliEveTPCLoaderEditor::DoDataAutoPedestal()
{
  // Slot for DataAutoPedestal.

  if (fM->GetData() == 0) return;
  fM->GetData()->SetAutoPedestal(fDataAutoPedestal->IsOn());
  fDataLoadPedestal->SetEnabled(!fDataAutoPedestal->IsOn());
}

/******************************************************************************/

void AliEveTPCLoaderEditor::DoUpdateSectors()
{
  // Slot for UpdateSectors.

  if(gCutOnEta)
    fM->SetCutOnEta(gCutOnEta->IsOn());

  if(gEtaRange)
  {
    fM->SetEtaMin(gEtaRange->GetMin());
    fM->SetEtaMax(gEtaRange->GetMax());
  }

  fM->UpdateSectors();
}

void AliEveTPCLoaderEditor::DoReloadSectors()
{
  // Slot for ReloadSectors.

  if(gCutOnEta)
    fM->SetCutOnEta(gCutOnEta->IsOn());

  if(gEtaRange)
  {
    fM->SetEtaMin(gEtaRange->GetMin());
    fM->SetEtaMax(gEtaRange->GetMax());
  }


  fM->ReloadSectors();
}

void AliEveTPCLoaderEditor::DoCreateSectors3D()
{
  // Slot for CreateSectors3D.

  if(gCutOnEta)
    fM->SetCutOnEta(gCutOnEta->IsOn());

  if(gEtaRange)
  {
    fM->SetEtaMin(gEtaRange->GetMin());
    fM->SetEtaMax(gEtaRange->GetMax());
  }

  fM->CreateSectors3D();
}


void AliEveTPCLoaderEditor::DoDeleteSectors3D()
{
  // Slot for DeleteSectors3D.

  fM->DeleteSectors3D();
}

void AliEveTPCLoaderEditor::DoShowSectors2D()
{

   for(Int_t i = 0; i< 36; i++)
   {
      if(gEve->GetEventScene()->FirstChild()->FindChild("AliEveTPCLoader")->FindChild(Form("Sector2D %d",i)))
      {
         gEve->GetEventScene()->FirstChild()->FindChild("AliEveTPCLoader")->FindChild(Form("Sector2D %d",i))->SetRnrSelf(kTRUE);
         gEve->GetEventScene()->FirstChild()->FindChild("AliEveTPCLoader")->FindChild(Form("Sector2D %d",i))->SetRnrChildren(kTRUE);
      }
   }

}

void AliEveTPCLoaderEditor::DoHideSectors2D()
{

   for(Int_t i = 0; i< 36; i++)
   {
      if(gEve->GetEventScene()->FirstChild()->FindChild("AliEveTPCLoader")->FindChild(Form("Sector2D %d",i)))
      {
         gEve->GetEventScene()->FirstChild()->FindChild("AliEveTPCLoader")->FindChild(Form("Sector2D %d",i))->SetRnrSelf(kFALSE);
         gEve->GetEventScene()->FirstChild()->FindChild("AliEveTPCLoader")->FindChild(Form("Sector2D %d",i))->SetRnrChildren(kFALSE);
      }
   }

}

