// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTrackCounterEditor.h"
#include "AliEveTrackCounter.h"
#include "AliEveEventManager.h"

#include "TGedEditor.h"
#include "TVirtualPad.h"
#include "TColor.h"

// Cleanup these includes:
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TGComboBox.h"
#include "TGMsgBox.h"

#include "TH1F.h"

#include "TCanvas.h"
#include "TEveManager.h"

#include "TROOT.h"
#include "TSystem.h" // File input/output for track-count status.

//______________________________________________________________________________
// GUI editor for AliEveTrackCounter.
//

ClassImp(AliEveTrackCounterEditor)

//______________________________________________________________________________
AliEveTrackCounterEditor::AliEveTrackCounterEditor(const TGWindow *p, Int_t width, Int_t height,
                                               UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM(0),
   fClickAction (0),
   fInfoLabel   (0),
   fEventId     (0)
{
   // Constructor.

   MakeTitle("AliEveTrackCounter");

   Int_t labelW = 42;

   { // ClickAction
      TGHorizontalFrame* f = new TGHorizontalFrame(this);
      TGLabel* lab = new TGLabel(f, "Click:");
      f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 10, 1, 2));
      fClickAction = new TGComboBox(f);
      fClickAction->AddEntry("Print", 0);
      fClickAction->AddEntry("Toggle", 1);
      TGListBox* lb = fClickAction->GetListBox();
      lb->Resize(lb->GetWidth(), 2*16);
      fClickAction->Resize(70, 20);
      fClickAction->Connect("Selected(Int_t)", "AliEveTrackCounterEditor", this,
                            "DoClickAction(Int_t)");
      f->AddFrame(fClickAction, new TGLayoutHints(kLHintsLeft, 1, 2, 1, 1));

      AddFrame(f);
   }

   { // Status
      TGHorizontalFrame* f = new TGHorizontalFrame(this);
      TGLabel* lab = new TGLabel(f, "Status:");
      f->AddFrame(lab, new TGLayoutHints(kLHintsLeft, 1, 5, 1, 2));

      fInfoLabel = new TGLabel(f);
      f->AddFrame(fInfoLabel, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 9, 1, 2));

      AddFrame(f);
   }

   {
      TGHorizontalFrame* f = new TGHorizontalFrame(this, 210, 20, kFixedWidth);

      TGHorizontalFrame* g = new TGHorizontalFrame(f, labelW, 0, kFixedWidth);
      TGLabel* l = new TGLabel(g, "Event:");
      g->AddFrame(l, new TGLayoutHints(kLHintsLeft, 0,0,4,0));
      f->AddFrame(g);

      TGTextButton* b;

      b = new TGTextButton(f, "Prev");
      f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
      b->Connect("Clicked()", "AliEveTrackCounterEditor", this, "DoPrev()");

      fEventId = new TGNumberEntry(f, 0, 3, -1,TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
                                   TGNumberFormat::kNELLimitMinMax, 0, 10000);
      f->AddFrame(fEventId, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
      fEventId->Connect("ValueSet(Long_t)", "AliEveTrackCounterEditor", this, "DoSetEvent()");

      b = new TGTextButton(f, "Next");
      f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
      b->Connect("Clicked()", "AliEveTrackCounterEditor", this, "DoNext()");

      AddFrame(f);
   }

   {
      TGHorizontalFrame* f = new TGHorizontalFrame(this, 210, 20, kFixedWidth);

      TGHorizontalFrame* g = new TGHorizontalFrame(f, labelW, 0, kFixedWidth);
      TGLabel* l = new TGLabel(g, "Report:");
      g->AddFrame(l, new TGLayoutHints(kLHintsLeft, 0,0,4,0));
      f->AddFrame(g);

      TGTextButton* b;

      b = new TGTextButton(f, "Print");
      f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
      b->Connect("Clicked()", "AliEveTrackCounterEditor", this, "DoPrintReport()");

      b = new TGTextButton(f, "File");
      f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
      b->Connect("Clicked()", "AliEveTrackCounterEditor", this, "DoFileReport()");

      AddFrame(f, new TGLayoutHints(kLHintsLeft, 0, 0, 4, 0));
   }
   {
      TGHorizontalFrame* f = new TGHorizontalFrame(this, 210, 20, kFixedWidth);

      TGHorizontalFrame* g = new TGHorizontalFrame(f, labelW, 0, kFixedWidth);
      TGLabel* l = new TGLabel(g, "Histos:");
      g->AddFrame(l, new TGLayoutHints(kLHintsLeft, 0,0,4,0));
      f->AddFrame(g);

      TGTextButton* b;

      b = new TGTextButton(f, "Show");
      f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
      b->Connect("Clicked()", "AliEveTrackCounterEditor", this, "DoShowHistos()");

      AddFrame(f, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
   }

  AliEveEventManager::GetMaster()->Connect("NewEventLoaded()",
                        "AliEveTrackCounterEditor", this, "UpdateModel()");
}

AliEveTrackCounterEditor::~AliEveTrackCounterEditor()
{
  // Destructor.

  AliEveEventManager::GetMaster()->Disconnect("NewEventLoaded()", this);
}

/******************************************************************************/

void AliEveTrackCounterEditor::UpdateModel()
{
  if (fGedEditor && fM && fGedEditor->GetModel() == fM->GetEditorObject())
  {
    SetModel(fM->GetEditorObject());
  }
}

//______________________________________________________________________________
void AliEveTrackCounterEditor::SetModel(TObject* obj)
{
   // Set model object.

   fM = dynamic_cast<AliEveTrackCounter*>(obj);

   fClickAction->Select(fM->fClickAction, kFALSE);
   fInfoLabel->SetText(Form("All: %3d; Primaries: %3d", fM->fAllTracks, fM->fGoodTracks));
   fEventId->SetNumber(fM->GetEventId());
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrackCounterEditor::DoPrev()
{
   // Slot for Prev.

   AliEveEventManager::GetMaster()->PrevEvent();
}

//______________________________________________________________________________
void AliEveTrackCounterEditor::DoNext()
{
   // Slot for Next.

   AliEveEventManager::GetMaster()->NextEvent();
}

//______________________________________________________________________________
void AliEveTrackCounterEditor::DoSetEvent()
{
   // Slot for SetEvent.
   AliEveEventManager::GetMaster()->GotoEvent((Int_t) fEventId->GetNumber());
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrackCounterEditor::DoPrintReport()
{
   // Slot for PrintReport.

   fM->OutputEventTracks();
}

//______________________________________________________________________________
void AliEveTrackCounterEditor::DoFileReport()
{
   // Slot for FileReport.

   TString file(Form("ev-report-%03d.txt", fM->GetEventId()));
   if (gSystem->AccessPathName(file) == kFALSE)
   {
      Int_t ret;
      new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                   "File Exist",
                   Form("Event record for event %d already exist.\n Replace?", fM->GetEventId()),
                   kMBIconQuestion, kMBYes | kMBNo, &ret);
      if (ret == kMBNo)
         return;
   }
   FILE* out = fopen(file, "w");
   if (out) {
      fM->OutputEventTracks(out);
      fclose(out);
   } else {
      Error("AliEveTrackCounterEditor::DoFileReport",
            "Can not open file '%s' for writing.", file.Data());
   }
}

//______________________________________________________________________________
void AliEveTrackCounterEditor::DoShowHistos()
{
  // Slot for ShowHistos.

  TH1F* hcnt = new TH1F("cnt", "Primeries per event", 41, -0.5, 40.5);
  TH1F* hchg = new TH1F("chg", "Primary charge",       3, -1.5,  1.5);
  TH1F* hpt  = new TH1F("pt",  "pT distribution",     40,  0.0,  8.0);
  TH1F* heta = new TH1F("eta", "eta distribution",    40, -1.0,  1.0);

  Int_t nn; // fscanf return value

  for (Int_t i=0; i<1000; ++i)
  {
    TString file(Form("ev-report-%03d.txt", i));
    if (gSystem->AccessPathName(file) == kFALSE)
    {
      Int_t   ev, ntr;
      FILE* f = fopen(file, "read");
      nn = fscanf(f, "Event = %d  Ntracks = %d", &ev, &ntr);
      if (nn != 2) { printf("SAFR1 %d\n", nn); fclose(f); return;  }
      hcnt->Fill(ntr);
      for (Int_t t=0; t<ntr; ++t)
      {
        Int_t   id, chg;
        Float_t pt, eta;
        nn = fscanf(f, "%d: chg=%d pt=%f eta=%f", &id, &chg, &pt, &eta);
        if (nn != 4) { printf("SAFR2 %d\n", nn); fclose(f); return;  }
        hchg->Fill(chg);
        hpt ->Fill(pt);
        heta->Fill(eta);
      }
      fclose(f);
    }
  }

  TCanvas* c;
  if (gPad == 0 || gPad->GetCanvas()->IsEditable() == kFALSE) {
    c = new TCanvas("Scanwas", "Scanning Results", 800, 600);
  } else {
    c = gPad->GetCanvas();
    c->Clear();
  }
  c->Divide(2, 2);

  c->cd(1); hcnt->Draw();
  c->cd(2); hchg->Draw();
  c->cd(3); hpt ->Draw();
  c->cd(4); heta->Draw();

  c->Modified();
  c->Update();
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrackCounterEditor::DoClickAction(Int_t mode)
{
   // Slot for ClickAction.

   fM->SetClickAction(mode);
}
