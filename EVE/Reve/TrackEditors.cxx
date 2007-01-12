// $Header$

#include "TrackEditors.h"
#include <Reve/Track.h>

#include <Reve/RGValuators.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>
#include <TGComboBox.h>

#include <TGMsgBox.h>
#include <TH1F.h>

using namespace Reve;

//______________________________________________________________________
// TrackListEditor
//

ClassImp(TrackListEditor)

TrackListEditor::TrackListEditor(const TGWindow *p,
				 Int_t width, Int_t height,
				 UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),

  fTC (0),

  fMaxR(0),
  fMaxZ(0),
  fMaxOrbits(0),
  fMinAng(0),
  fDelta(0),

  fRnrTracks(0),
  fRnrMarkers(0),

  fFitDaughters(0),
  fFitDecay(0),

  fPtRange(0)
{
  MakeTitle("TrackList");
  Int_t labelW = 67;

  // --- Limits

  fMaxR = new RGValuator(this, "Max R:", 110, 0);
  fMaxR->SetLabelWidth(labelW);
  fMaxR->SetNELength(6);
  fMaxR->Build();
  fMaxR->SetLimits(0.1, 1000, 100, TGNumberFormat::kNESRealOne);
  fMaxR->SetToolTip("Maximum radius to which the tracks will be drawn.");
  fMaxR->Connect("ValueSet(Double_t)", "Reve::TrackListEditor", this, "DoMaxR()");
  AddFrame(fMaxR, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  fMaxZ = new RGValuator(this, "Max Z:", 110, 0);
  fMaxZ->SetLabelWidth(labelW);
  fMaxZ->SetNELength(6);
  fMaxZ->Build();
  fMaxZ->SetLimits(0.1, 2000, 100, TGNumberFormat::kNESRealOne);
  fMaxZ->SetToolTip("Maximum z-coordinate to which the tracks will be drawn.");
  fMaxZ->Connect("ValueSet(Double_t)", "Reve::TrackListEditor", this, "DoMaxZ()");
  AddFrame(fMaxZ, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Max Orbits:");
    f->AddFrame(l, new TGLayoutHints(kLHintsTop | kLHintsCenterY, 0, 2, 1, 1));
    fMaxOrbits = new TGNumberEntry(f, 0., 6, -1, 
				   TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive,
				   TGNumberFormat::kNELLimitMinMax, 0.1, 100.0);
    fMaxOrbits->GetNumberEntry()->SetToolTipText("Maximal angular path of tracks' orbits (1 ~ 2Pi).");
    f->AddFrame(fMaxOrbits, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fMaxOrbits->Associate(f);
    fMaxOrbits->Connect("ValueSet(Long_t)", "Reve::TrackListEditor", this, "DoMaxOrbits()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Min Angle:");
    f->AddFrame(l, new TGLayoutHints(kLHintsTop | kLHintsCenterY, 3, 2, 1, 1));
    fMinAng = new TGNumberEntry(f, 0., 6, -1, 
				TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive,
				TGNumberFormat::kNELLimitMinMax, 1, 180.0);
    fMinAng->GetNumberEntry()->SetToolTipText("Minimal angular step between two helix points.");
    f->AddFrame(fMinAng, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fMinAng->Associate(f);
    fMinAng->Connect("ValueSet(Long_t)", "Reve::TrackListEditor", this, "DoMinAng()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Delta:");
    f->AddFrame(l, new TGLayoutHints(kLHintsTop | kLHintsCenterY, 32, 2, 1, 1));
    fDelta = new TGNumberEntry(f, 0., 6, -1, 
			       TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAPositive,
			       TGNumberFormat::kNELLimitMinMax, 0.001, 100.0);
    fDelta->GetNumberEntry()->SetToolTipText("Maximal error at the mid-point of the line connecting to helix points.");
    f->AddFrame(fDelta, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fDelta->Associate(f);
    fDelta->Connect("ValueSet(Long_t)", "Reve::TrackListEditor", this, "DoDelta()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  // --- Rendering control

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fRnrTracks = new TGCheckButton(f, "Render tracks");
    f->AddFrame(fRnrTracks, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 0));
    fRnrTracks->Connect
      ("Toggled(Bool_t)", "Reve::TrackListEditor", this, "DoRnrTracks()");
    fWidthCombo = new TGLineWidthComboBox(f);
    fWidthCombo->Resize(80, 18);
    f->AddFrame(fWidthCombo, new TGLayoutHints(kLHintsLeft, 8, 1, 0, 0));

    fWidthCombo->Connect
      ("Selected(Int_t)", "Reve::TrackListEditor", this, "DoLineWidth(Int_t)"); 
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 3, 0));
  }

  fRnrMarkers = new TGCheckButton(this, "Render markers");
  AddFrame(fRnrMarkers, new TGLayoutHints(kLHintsTop, 3, 1, 2, 0));
  fRnrMarkers->Connect
    ("Toggled(Bool_t)",
     "Reve::TrackListEditor", this, "DoRnrMarkers()");  

  // --- Kinematics fitting

  fFitDaughters = new TGCheckButton(this, "Fit daughters");
  AddFrame(fFitDaughters, new TGLayoutHints(kLHintsTop, 3, 1, 2, 0));
  fFitDaughters->Connect("Toggled(Bool_t)","Reve::TrackListEditor", this, "DoFitDaughters()");

  fFitDecay = new TGCheckButton(this, "Fit decay");
  AddFrame(fFitDecay, new TGLayoutHints(kLHintsTop, 3, 1, 2, 0));
  fFitDecay->Connect("Toggled(Bool_t)","Reve::TrackListEditor", this, "DoFitDecay()");  

  // --- Selectors

  fPtRange = new RGDoubleValuator(this,"Pt Range", 200, 0);
  fPtRange->SetNELength(6);
  fPtRange->Build();
  fPtRange->GetSlider()->SetWidth(224);
  fPtRange->SetLimits(0.1, 10, TGNumberFormat::kNESRealTwo);
  fPtRange->Connect("ValueSet()",
                    "Reve::TrackListEditor", this, "DoPtRange()");
  AddFrame(fPtRange, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));
}

TrackListEditor::~TrackListEditor()
{}

/**************************************************************************/

void TrackListEditor::SetModel(TObject* obj)
{
  fTC = dynamic_cast<TrackList*>(obj);

  fMaxR->SetValue(fTC->GetMaxR());
  fMaxZ->SetValue(fTC->GetMaxZ());
  fMaxOrbits->SetNumber(fTC->GetMaxOrbs());
  fMinAng->SetNumber(fTC->GetMinAng());
  fDelta->SetNumber(fTC->GetDelta());

  fWidthCombo->Select(fTC->GetWidth());

  fRnrTracks->SetState(fTC->GetRnrTracks() ? kButtonDown : kButtonUp);
  fRnrMarkers->SetState(fTC->GetRnrMarkers() ? kButtonDown : kButtonUp);

  fFitDaughters->SetState(fTC->GetFitDaughters() ? kButtonDown : kButtonUp);
  fFitDecay->SetState(fTC->GetFitDecay() ? kButtonDown : kButtonUp);

  fPtRange->SetValues(0.1, 10);
}

/**************************************************************************/

void TrackListEditor::DoMaxR()
{
  fTC->SetMaxR(fMaxR->GetValue());
  Update();
}

void TrackListEditor::DoMaxZ()
{
  fTC->SetMaxZ(fMaxZ->GetValue());
  Update();
}

void TrackListEditor::DoMaxOrbits()
{
  fTC->SetMaxOrbs(fMaxOrbits->GetNumber());
  Update();
}

void TrackListEditor::DoMinAng()
{
  fTC->SetMinAng(fMinAng->GetNumber());
  Update();
}

void TrackListEditor::DoDelta()
{
  fTC->SetDelta(fDelta->GetNumber());
  Update();
}

/**************************************************************************/

void TrackListEditor::DoLineWidth(Int_t width)
{
  fTC->SetWidth(width);
  Update();
}

/**************************************************************************/

void TrackListEditor::DoRnrTracks()
{
  fTC->SetRnrTracks(fRnrTracks->IsOn());
  Update();
}

void TrackListEditor::DoRnrMarkers()
{
  fTC->SetRnrMarkers(fRnrMarkers->IsOn());
  Update();
}

/**************************************************************************/

void TrackListEditor::DoFitDaughters()
{
  fTC->SetFitDaughters(fFitDaughters->IsOn());
  Update();
}

void TrackListEditor::DoFitDecay()
{
  fTC->SetFitDecay(fFitDecay->IsOn());
  Update();
}

/**************************************************************************/

void TrackListEditor::DoPtRange()
{
  fTC->SelectByPt(fPtRange->GetMin(), fPtRange->GetMax());
  Update();
}



/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

#include <TCanvas.h>
#include <TGLViewer.h>
#include <Reve/RGTopFrame.h>

//______________________________________________________________________
// TrackCounterEditor
//

ClassImp(TrackCounterEditor)

TrackCounterEditor::TrackCounterEditor(const TGWindow *p, Int_t width, Int_t height,
				       UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fClickAction (0),
  fInfoLabel   (0)
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
}

/**************************************************************************/

// glv->SetOrthoCamera(TGLViewer::kCameraOrthoXOY, 2*left, 2*right, 2*top, bottom);

void TrackCounterEditor::DoOrtoXY()
{
  TGLViewer* glv = dynamic_cast<TGLViewer*>(gReve->GetGLCanvas()->GetViewer3D());
  glv->SetCurrentCamera(TGLViewer::kCameraOrthoXOY) ;
}

void TrackCounterEditor::DoOrtoZY()
{
  TGLViewer* glv = dynamic_cast<TGLViewer*>(gReve->GetGLCanvas()->GetViewer3D());
  glv->SetCurrentCamera(TGLViewer::kCameraOrthoZOY) ;
}

void TrackCounterEditor::DoPersp()
{
  TGLViewer* glv = dynamic_cast<TGLViewer*>(gReve->GetGLCanvas()->GetViewer3D());
  glv->SetCurrentCamera(TGLViewer::kCameraPerspXOZ) ;
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
