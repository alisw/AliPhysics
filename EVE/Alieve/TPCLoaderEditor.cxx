// $Header$

#include "TPCLoaderEditor.h"
#include <Alieve/TPCLoader.h>
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

  fUpdateSectors   (0),
  fCreateSectors3D (0),
  fDeleteSectors3D (0)
{
  MakeTitle("TPCLoader");

  Int_t labelW = 42;
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
  fEvent->AddFrame(fDoubleSR, new TGLayoutHints(kLHintsLeft, 12, 0, 1, 0));
  fDoubleSR->Connect("Toggled(Bool_t)",
                     "Alieve::TPCLoaderEditor", this, "DoDoubleSR()");
  AddFrame(fEvent);

  fUpdateSectors = new TGTextButton(this, "Update Sectors");
  AddFrame(fUpdateSectors, new TGLayoutHints(kLHintsExpandX, 8,8,8,0));
  fUpdateSectors->Connect("Clicked()",
		     "Alieve::TPCLoaderEditor", this, "DoUpdateSectors()");
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
  fDoubleSR->SetState(fM->fDoubleSR  ? kButtonDown : kButtonUp);
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

void TPCLoaderEditor::DoUpdateSectors()
{
  fM->UpdateSectors();
}

void TPCLoaderEditor::DoCreateSectors3D()
{
  fM->CreateSectors3D();
}


void TPCLoaderEditor::DoDeleteSectors3D()
{
  fM->DeleteSectors3D();
}

