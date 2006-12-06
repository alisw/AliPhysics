// $Header$

#include "TPCLoaderEditor.h"
#include <Alieve/TPCLoader.h>
#include <Alieve/TPCData.h>
#include <Reve/RGTopFrame.h>
#include <Reve/RGValuators.h>

#include <TSystem.h>
#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGFileDialog.h>
#include <TGToolTip.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// TPCLoaderEditor
//

ClassImp(TPCLoaderEditor)

  TPCLoaderEditor::TPCLoaderEditor(const TGWindow *p,
				   Int_t width, Int_t height,
				   UInt_t options, Pixel_t back) :
    TGedFrame(p, width, height, options | kVerticalFrame, back),

    fM (0),

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
    fDeleteSectors3D (0)
{
  MakeTitle("TPCLoader");

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
		   "Alieve::TPCLoaderEditor", this, "FileSelect()");
    fFile->Connect("TextChanged(const char *)",
		   "Alieve::TPCLoaderEditor", this, "FileChanged()");
    fOpenFile = new TGTextButton(f, "Open");
    f->AddFrame(fOpenFile);
    fOpenFile->Connect("Clicked()",
		       "Alieve::TPCLoaderEditor", this, "DoOpen()");
    AddFrame(f);
  }

  fEvent = new RGValuator(this, "Event:", 110, 0);
  fEvent->SetShowSlider(kFALSE);
  fEvent->SetLabelWidth(labelW);
  fEvent->SetNELength(6);
  fEvent->Build();
  fEvent->SetLimits(0, 1000);
  fEvent->SetToolTip("Current event number");
  fEvent->Connect("ValueSet(Double_t)",
		  "Alieve::TPCLoaderEditor", this, "DoEvent()");
  // Reuse Event for DoubleSR button
  fDoubleSR = new TGCheckButton(fEvent, "Double SR");
  fDoubleSR->SetToolTipText("Double sampling rate");
  fEvent->AddFrame(fDoubleSR, new TGLayoutHints(kLHintsLeft, 12, 0, 2, 0));
  fDoubleSR->Connect("Toggled(Bool_t)",
                     "Alieve::TPCLoaderEditor", this, "DoDoubleSR()");
  AddFrame(fEvent);


  // TPCData load settings

  labelW = 90;

  fDataLoadThreshold = new RGValuator(this, "Load threshold:", 110, 0);
  fDataLoadThreshold->SetShowSlider(kFALSE);
  fDataLoadThreshold->SetLabelWidth(labelW);
  fDataLoadThreshold->SetNELength(6);
  fDataLoadThreshold->Build();
  fDataLoadThreshold->SetLimits(0, 1000);
  fDataLoadThreshold->SetToolTip("Minimum signal that will be stored (pedestal subtracted first).");
  fDataLoadThreshold->Connect
    ("ValueSet(Double_t)",
     "Alieve::TPCLoaderEditor", this, "DoDataLoadThreshold()");
  AddFrame(fDataLoadThreshold, new TGLayoutHints(kLHintsLeft, 0, 0, 6, 0));

  fDataLoadPedestal = new RGValuator(this, "Load pedestal:", 110, 0);
  fDataLoadPedestal->SetShowSlider(kFALSE);
  fDataLoadPedestal->SetLabelWidth(labelW);
  fDataLoadPedestal->SetNELength(6);
  fDataLoadPedestal->Build();
  fDataLoadPedestal->SetLimits(0, 1000);
  fDataLoadPedestal->SetToolTip("Pedestal that will be subtracted during data loading.");
  fDataLoadPedestal->Connect
    ("ValueSet(Double_t)",
     "Alieve::TPCLoaderEditor", this, "DoDataLoadPedestal()");
  // Reuse DataLoadPedestal for DataAutoPedestal check button
  fDataAutoPedestal = new TGCheckButton(fDataLoadPedestal, "Automatic");
  fDataAutoPedestal->SetToolTipText("Determine per-pad pedestal during data loading.");
  fDataAutoPedestal->Connect
    ("Toggled(Bool_t)",
     "Alieve::TPCLoaderEditor", this, "DoDataAutoPedestal()");
  fDataLoadPedestal->AddFrame(fDataAutoPedestal, new TGLayoutHints(kLHintsLeft, 12, 0, 2, 0));
  AddFrame(fDataLoadPedestal);

  // Steering buttons: update/reload sectors; show/hide 3d

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fUpdateSectors = new TGTextButton(f, "Update Sectors");
    f->AddFrame(fUpdateSectors, new TGLayoutHints(kLHintsExpandX, 0,4,0,0));
    fUpdateSectors->Connect("Clicked()",
			    "Alieve::TPCLoaderEditor", this, "DoUpdateSectors()");
    fReloadSectors = new TGTextButton(f, "Reload Sectors");
    f->AddFrame(fReloadSectors, new TGLayoutHints(kLHintsExpandX, 4,0,0,0));
    fReloadSectors->Connect("Clicked()",
			    "Alieve::TPCLoaderEditor", this, "DoReloadSectors()");
    AddFrame(f, new TGLayoutHints(kLHintsExpandX, 8,8,8,0));
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fCreateSectors3D = new TGTextButton(f, "Create 3D");
    f->AddFrame(fCreateSectors3D, new TGLayoutHints(kLHintsExpandX, 0,4,0,0));
    fCreateSectors3D->Connect("Clicked()",
			      "Alieve::TPCLoaderEditor", this, "DoCreateSectors3D()");
    fDeleteSectors3D = new TGTextButton(f, "Delete 3D");
    f->AddFrame(fDeleteSectors3D, new TGLayoutHints(kLHintsExpandX, 4,0,0,0));
    fDeleteSectors3D->Connect("Clicked()",
			      "Alieve::TPCLoaderEditor", this, "DoDeleteSectors3D()");
    AddFrame(f, new TGLayoutHints(kLHintsExpandX, 8,8,8,0));
  }
}

TPCLoaderEditor::~TPCLoaderEditor()
{}

/**************************************************************************/

void TPCLoaderEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<TPCLoader*>(obj);

  // !!!! order changed, need TGTextEntry::SetText NO BLOODY EMIT.
  fFile->SetToolTipText(gSystem->DirName(fM->fFile));
  fFile->SetText(gSystem->BaseName(fM->fFile));
  fEvent->SetValue(fM->fEvent);
  fEvent->SetEnabled(fM->fEvent >= 0);
  fDoubleSR->SetState(fM->fDoubleSR ? kButtonDown : kButtonUp);

  TPCData* tpcd = fM->GetData();
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

/**************************************************************************/
/**************************************************************************/

namespace {
const char *tpcfiletypes[] = {
   "Root files",  "*.root",
   "All files",   "*.*",
    0,               0
};
}

void TPCLoaderEditor::FileSelect()
{
  TGFileInfo fi;
  fi.fIniDir    = StrDup(gSystem->DirName (fM->fFile));
  fi.fFilename  = StrDup(gSystem->BaseName(fM->fFile));
  fi.fFileTypes = tpcfiletypes;

  new TGFileDialog(fClient->GetRoot(), gReve, kFDOpen, &fi);
  if (!fi.fFilename)
    return;

  fFile->SetToolTipText(gSystem->DirName (fi.fFilename));
  fFile->SetText       (gSystem->BaseName(fi.fFilename));
}

void TPCLoaderEditor::FileChanged()
{
  fM->fFile = Form("%s/%s", fFile->GetToolTip()->GetText()->Data(),
		   fFile->GetText());
}

void TPCLoaderEditor::DoOpen()
{
  fM->OpenFile();
  SetModel(fM);
}

/**************************************************************************/

void TPCLoaderEditor::DoEvent()
{
  fM->GotoEvent((Int_t) fEvent->GetValue());
  SetModel(fM);
}

void TPCLoaderEditor::DoDoubleSR()
{
  fM->SetDoubleSR(fDoubleSR->IsOn());
  Update();
}

/**************************************************************************/

void TPCLoaderEditor::DoDataLoadThreshold()
{
  if (fM->GetData() == 0) return;
  fM->GetData()->SetLoadThreshold((Short_t) fDataLoadThreshold->GetValue());
}

void TPCLoaderEditor::DoDataLoadPedestal()
{
  if (fM->GetData() == 0) return;
  fM->GetData()->SetLoadPedestal((Short_t) fDataLoadPedestal->GetValue());
}

void TPCLoaderEditor::DoDataAutoPedestal()
{
  if (fM->GetData() == 0) return;
  fM->GetData()->SetAutoPedestal(fDataAutoPedestal->IsOn());
  fDataLoadPedestal->SetEnabled(!fDataAutoPedestal->IsOn());
}

/**************************************************************************/

void TPCLoaderEditor::DoUpdateSectors()
{
  fM->UpdateSectors();
}

void TPCLoaderEditor::DoReloadSectors()
{
  fM->ReloadSectors();
}

void TPCLoaderEditor::DoCreateSectors3D()
{
  fM->CreateSectors3D();
}


void TPCLoaderEditor::DoDeleteSectors3D()
{
  fM->DeleteSectors3D();
}

