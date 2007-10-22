// $Header$

#include "TrackEditors.h"
#include <Reve/Track.h>
#include <Reve/TrackRnrStyleEditor.h>
#include <Reve/ReveManager.h>

#include <Reve/RGValuators.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGedEditor.h>
#include <TAttMarkerEditor.h>
#include <TGLabel.h>
#include <TG3DLine.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>
#include <TGComboBox.h>

#include <TGMsgBox.h>
#include <TH1F.h>

#include <TSystem.h> // File input/output for track-count status.

using namespace Reve;


//______________________________________________________________________
// TrackEditor
//

ClassImp(TrackEditor)

TrackEditor::TrackEditor(const TGWindow *p, Int_t width, Int_t height,
			       UInt_t options, Pixel_t back) :
    TGedFrame(p, width, height, options | kVerticalFrame, back),
    fM (0),
    fRSEditor(0)
{
  MakeTitle("Track");

  TGHorizontalFrame* f = new TGHorizontalFrame(this);

  fRSEditor =  new TGTextButton(f, "Edit RnrStyle");
  fRSEditor->Connect("Clicked()", "Reve::TrackEditor", this, "DoEditRnrStyle()");
  f->AddFrame(fRSEditor, new TGLayoutHints(kLHintsLeft, 2, 1, 4, 4));
 
  AddFrame(f, new TGLayoutHints(kLHintsTop, 0,0,2,1));
}

void TrackEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<Track*>(obj); 
}

void TrackEditor::DoEditRnrStyle()
{
  fGedEditor->SetModel(fGedEditor->GetPad(), fM->GetRnrStyle(), kButton1Down);
}
//______________________________________________________________________
// TrackListEditor
//

ClassImp(TrackListEditor)

TrackListEditor::TrackListEditor(const TGWindow *p,
                                 Int_t width, Int_t height,
                                 UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),

  fTC         (0),
  fPtRange    (0),
  fPRange     (0),
  fRSSubEditor(0)
{
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
 
    fRnrPoints = new TGCheckButton(f, "Draw Marker");
    f->AddFrame(fRnrPoints, new TGLayoutHints(kLHintsLeft, 2,1,0,0));  
    fRnrPoints->Connect("Toggled(Bool_t)"," Reve::TrackListEditor", this, "DoRnrPoints()");

    fRnrLine  = new TGCheckButton(f, "Draw Line");
    f->AddFrame(fRnrLine, new TGLayoutHints(kLHintsLeft, 1,2,0,0));
    fRnrLine->Connect("Toggled(Bool_t)", "Reve::TrackListEditor", this, "DoRnrLine()");
  
    AddFrame(f, new TGLayoutHints(kLHintsTop, 0,0,2,1));
  }
  {  // --- Selectors
    Int_t labelW = 51;
    Int_t dbW    = 210;

    fPtRange = new RGDoubleValuator(this,"Pt rng:", 40, 0);
    fPtRange->SetNELength(6);
    fPtRange->SetLabelWidth(labelW);
    fPtRange->Build();
    fPtRange->GetSlider()->SetWidth(dbW);
    fPtRange->SetLimits(0, 10, TGNumberFormat::kNESRealTwo);
    fPtRange->Connect("ValueSet()",
                      "Reve::TrackListEditor", this, "DoPtRange()");
    AddFrame(fPtRange, new TGLayoutHints(kLHintsTop, 1, 1, 4, 1));

    fPRange = new RGDoubleValuator(this,"P rng:", 40, 0);
    fPRange->SetNELength(6);
    fPRange->SetLabelWidth(labelW);
    fPRange->Build();
    fPRange->GetSlider()->SetWidth(dbW);
    fPRange->SetLimits(0, 100, TGNumberFormat::kNESRealTwo);
    fPRange->Connect("ValueSet()",
                     "Reve::TrackListEditor", this, "DoPRange()");
    AddFrame(fPRange, new TGLayoutHints(kLHintsTop, 1, 1, 4, 1));
  }

  MakeTitle("RenderStyle");
  fRSSubEditor = new TrackRnrStyleSubEditor(this);
  fRSSubEditor->Connect("Changed()", "Reve::TrackListEditor", this, "Update()"); 
  AddFrame(fRSSubEditor, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0,0,0));
  CreateRefsTab();
}

TrackListEditor::~TrackListEditor()
{}

/**************************************************************************/

void TrackListEditor::CreateRefsTab()
{
  fRefs = CreateEditorTabSubFrame("Refs");

  TGCompositeFrame *title1 = new TGCompositeFrame(fRefs, 145, 10, 
						  kHorizontalFrame | 
						  kLHintsExpandX   | 
						  kFixedWidth      | 
						  kOwnBackground);
  title1->AddFrame(new TGLabel(title1, "PathMarks"), 
		   new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
  title1->AddFrame(new TGHorizontal3DLine(title1),
		   new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
  fRefs->AddFrame(title1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));
   
  // path marks 
  fRSSubEditor->CreateRefsContainer(fRefs);
  fRSSubEditor->fPMAtt->SetGedEditor((TGedEditor*)gReve->GetEditor()); 
  fRSSubEditor->fFVAtt->SetGedEditor((TGedEditor*)gReve->GetEditor());
}

/**************************************************************************/

void TrackListEditor::SetModel(TObject* obj)
{
  fTC = dynamic_cast<TrackList*>(obj); 
  fRnrLine  ->SetState(fTC->GetRnrLine()   ? kButtonDown : kButtonUp);
  fRnrPoints->SetState(fTC->GetRnrPoints() ? kButtonDown : kButtonUp);

  Float_t llim;
  fPtRange->SetValues(fTC->fMinPt, fTC->fMaxPt);
  llim = TMath::Log10(fTC->fLimPt);
  fPtRange->SetLimits(0, fTC->fLimPt, llim < 2 ? TGNumberFormat::kNESRealTwo : (llim < 3 ? TGNumberFormat::kNESRealOne : TGNumberFormat::kNESInteger));
  fPRange ->SetValues(fTC->fMinP, fTC->fMaxP);
  llim = TMath::Log10(fTC->fLimP);
  fPRange ->SetLimits(0, fTC->fLimP, llim < 2 ? TGNumberFormat::kNESRealTwo : (llim < 3 ? TGNumberFormat::kNESRealOne : TGNumberFormat::kNESInteger));

  fRSSubEditor->SetModel(fTC->GetRnrStyle());
}

/**************************************************************************/

void TrackListEditor::DoRnrLine()
{
  fTC->SetRnrLine(fRnrLine->IsOn());
  Update();
}

void TrackListEditor::DoRnrPoints()
{
  fTC->SetRnrPoints(fRnrPoints->IsOn());
  Update();
}

/**************************************************************************/

void TrackListEditor::DoPtRange()
{ 
  fTC->SelectByPt(fPtRange->GetMin(), fPtRange->GetMax());
  Update();
}

void TrackListEditor::DoPRange()
{
  fTC->SelectByP(fPRange->GetMin(), fPRange->GetMax());
  Update();
}


/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

#include <TCanvas.h>
#include <TGLViewer.h>
#include <Reve/ReveManager.h>

//______________________________________________________________________
// TrackCounterEditor
//

ClassImp(TrackCounterEditor)

TrackCounterEditor::TrackCounterEditor(const TGWindow *p, Int_t width, Int_t height,
				       UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fClickAction (0),
  fInfoLabel   (0),
  fEventId     (0)
{
  MakeTitle("TrackCounter");

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
    fClickAction->Connect("Selected(Int_t)", "Reve::TrackCounterEditor", this,
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
    TGLabel* l = new TGLabel(g, "View:");
    g->AddFrame(l, new TGLayoutHints(kLHintsLeft, 0,0,4,0));
    f->AddFrame(g);

    TGTextButton* b;

    b = new TGTextButton(f, "Orto XY");
    f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    b->Connect("Clicked()", "Reve::TrackCounterEditor", this, "DoOrtoXY()");

    b = new TGTextButton(f, "Orto ZY");
    f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    b->Connect("Clicked()", "Reve::TrackCounterEditor", this, "DoOrtoZY()");

    b = new TGTextButton(f, "Persp");
    f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    b->Connect("Clicked()", "Reve::TrackCounterEditor", this, "DoPersp()");

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
    b->Connect("Clicked()", "Reve::TrackCounterEditor", this, "DoPrev()");

    fEventId = new TGNumberEntry(f, 0, 3, -1,TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
				 TGNumberFormat::kNELLimitMinMax, 0, 1000);
    f->AddFrame(fEventId, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    fEventId->Connect("ValueSet(Long_t)", "Reve::TrackCounterEditor", this, "DoSetEvent()");

    b = new TGTextButton(f, "Next");
    f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    b->Connect("Clicked()", "Reve::TrackCounterEditor", this, "DoNext()");

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
    b->Connect("Clicked()", "Reve::TrackCounterEditor", this, "DoPrintReport()");

    b = new TGTextButton(f, "File");
    f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    b->Connect("Clicked()", "Reve::TrackCounterEditor", this, "DoFileReport()");

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
    b->Connect("Clicked()", "Reve::TrackCounterEditor", this, "DoShowHistos()");

    AddFrame(f, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
  }

}

TrackCounterEditor::~TrackCounterEditor()
{}

/**************************************************************************/

void TrackCounterEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<TrackCounter*>(obj);

  fClickAction->Select(fM->fClickAction, kFALSE);
  fInfoLabel->SetText(Form("All: %3d; Primaries: %3d", fM->fAllTracks, fM->fGoodTracks));
  fEventId->SetNumber(fM->GetEventId());
}

/**************************************************************************/

void TrackCounterEditor::DoOrtoXY()
{
  gReve->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY) ;
}

void TrackCounterEditor::DoOrtoZY()
{
  gReve->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoZOY) ;
}

void TrackCounterEditor::DoPersp()
{
  gReve->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraPerspXOZ) ;
}

/**************************************************************************/

void TrackCounterEditor::DoPrev()
{
  Reve::Macro("event_prev.C");
  gReve->EditRenderElement(fM);
}

void TrackCounterEditor::DoNext()
{
  Reve::Macro("event_next.C");
  gReve->EditRenderElement(fM);
}

void TrackCounterEditor::DoSetEvent()
{
  Reve::LoadMacro("event_goto.C");
  gROOT->ProcessLine(Form("event_goto(%d);", (Int_t) fEventId->GetNumber()));
  gReve->EditRenderElement(fM);
}

/**************************************************************************/

void TrackCounterEditor::DoPrintReport()
{
  fM->OutputEventTracks();
}

void TrackCounterEditor::DoFileReport()
{
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
  fM->OutputEventTracks(out);
  fclose(out);
}

void TrackCounterEditor::DoShowHistos()
{
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


/**************************************************************************/

void TrackCounterEditor::DoClickAction(Int_t mode)
{
  fM->SetClickAction(mode);
}
