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

#include <TSystem.h> // File input/output for track-count status.

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
  fWidthCombo(0),
  fStyleCombo(0),

  fRnrMarkers(0),

  fPMFrame(0),
  fFitDaughters(0),
  fFitReferences(0),
  fFitDecay(0),
  fRnrDaughters(0),
  fRnrReferences(0),
  fRnrDecay(0),

  fPtRange(0),
  fPRange(0)
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

  // --- Selectors

  fPtRange = new RGDoubleValuator(this,"Pt Range", 200, 0);
  fPtRange->SetNELength(6);
  fPtRange->Build();
  fPtRange->GetSlider()->SetWidth(224);
  fPtRange->SetLimits(0, 10, TGNumberFormat::kNESRealTwo);
  fPtRange->Connect("ValueSet()",
		    "Reve::TrackListEditor", this, "DoPtRange()");
  AddFrame(fPtRange, new TGLayoutHints(kLHintsTop, 1, 1, 4, 1));

  fPRange = new RGDoubleValuator(this,"P Range", 200, 0);
  fPRange->SetNELength(6);
  fPRange->Build();
  fPRange->GetSlider()->SetWidth(224);
  fPRange->SetLimits(0, 100, TGNumberFormat::kNESRealTwo);
  fPRange->Connect("ValueSet()",
		   "Reve::TrackListEditor", this, "DoPRange()");
  AddFrame(fPRange, new TGLayoutHints(kLHintsTop, 1, 1, 4, 1));

  // --- Rendering control
  {
    TGHorizontalFrame* ft = new TGHorizontalFrame(this);
    fRnrTracks = new TGCheckButton(ft, "Render tracks");
    ft->AddFrame(fRnrTracks, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 0));
    fRnrTracks->Connect
      ("Toggled(Bool_t)", "Reve::TrackListEditor", this, "DoRnrTracks()");
    AddFrame(ft, new TGLayoutHints(kLHintsTop, 1, 1, 3, 0));
    
    TGHorizontalFrame* f = new TGHorizontalFrame(this); 
    TGLabel *wl = new TGLabel(f, "Width:");
    f->AddFrame(wl, new TGLayoutHints(kLHintsTop | kLHintsCenterY, 3, 2, 1, 1));
    fWidthCombo = new TGLineWidthComboBox(f);
    fWidthCombo->Resize(70, 18);
    f->AddFrame(fWidthCombo, new TGLayoutHints(kLHintsLeft, 8, 1, 0, 0));
    fWidthCombo->Connect
      ("Selected(Int_t)", "Reve::TrackListEditor", this, "DoLineWidth(Int_t)"); 

    TGLabel *sl = new TGLabel(f, "Style:");
    f->AddFrame(sl, new TGLayoutHints(kLHintsTop | kLHintsCenterY, 6, 2, 1, 1));
    fStyleCombo = new TGLineStyleComboBox(f);
    fStyleCombo->Resize(70, 18);
    f->AddFrame(fStyleCombo, new TGLayoutHints(kLHintsLeft, 8, 1, 0, 0));
    fStyleCombo->Connect
      ("Selected(Int_t)", "Reve::TrackListEditor", this, "DoLineStyle(Int_t)"); 


    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 3, 0));
  }

  // --- Kinematics fitting
  
  fPMFrame = new  TGHorizontalFrame(this);
  {
    TGGroupFrame* fitPM = new TGGroupFrame(fPMFrame, "PathMarks:", kLHintsTop | kLHintsCenterX);
    fitPM->SetTitlePos(TGGroupFrame::kLeft);
    fPMFrame->AddFrame( fitPM, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 3, 3, 3, 3));

    TGMatrixLayout *ml = new TGMatrixLayout(fitPM, 0,1,6);
    fitPM->SetLayoutManager(ml);

    fFitDaughters  = new TGCheckButton(fitPM, "Fit Daughters", PathMark::Daughter);
    fFitReferences = new TGCheckButton(fitPM, "Fit Refs",      PathMark::Reference);
    fFitDecay      = new TGCheckButton(fitPM, "Fit Decay",     PathMark::Decay);

    fitPM->AddFrame(fFitDaughters);
    fitPM->AddFrame(fFitReferences);
    fitPM->AddFrame(fFitDecay);

    fFitDecay->Connect("Clicked()","Reve::TrackListEditor", this, "DoFitPM()");  
    fFitReferences->Connect("Clicked()","Reve::TrackListEditor", this, "DoFitPM()");  
    fFitDaughters->Connect("Clicked()","Reve::TrackListEditor", this, "DoFitPM()");
  }
  {
    TGGroupFrame* rnrPM = new TGGroupFrame(fPMFrame, "PathMarks:", kLHintsTop | kLHintsCenterX);
    rnrPM->SetTitlePos(TGGroupFrame::kLeft);
    fPMFrame->AddFrame( rnrPM, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 3, 3, 3, 3));

    TGMatrixLayout *ml = new TGMatrixLayout(rnrPM, 0,1,6);
    rnrPM->SetLayoutManager(ml);

    fRnrDaughters  = new TGCheckButton(rnrPM, "Rnr Daughters", PathMark::Daughter);
    fRnrReferences = new TGCheckButton(rnrPM, "Rnr Refs",  PathMark::Reference);
    fRnrDecay      = new TGCheckButton(rnrPM, "Rnr Decay", PathMark::Decay);

    rnrPM->AddFrame(fRnrDaughters);
    rnrPM->AddFrame(fRnrReferences);
    rnrPM->AddFrame(fRnrDecay);

    fRnrDecay->Connect("Clicked()","Reve::TrackListEditor", this, "DoRnrPM()");  
    fRnrReferences->Connect("Clicked()","Reve::TrackListEditor", this, "DoRnrPM()");  
    fRnrDaughters->Connect("Clicked()","Reve::TrackListEditor", this, "DoRnrPM()");
  }
  AddFrame(fPMFrame, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
 
  fRnrMarkers = new TGCheckButton(this, "Render markers");
  AddFrame(fRnrMarkers, new TGLayoutHints(kLHintsTop, 3, 1, 2, 0));
  fRnrMarkers->Connect
    ("Toggled(Bool_t)",
     "Reve::TrackListEditor", this, "DoRnrMarkers()");  

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
  fStyleCombo->Select(fTC->GetStyle());

  fRnrTracks->SetState(fTC->GetRnrTracks() ? kButtonDown : kButtonUp);
 
  if(fTC->GetEditPathMarks()) 
  {
    HideFrame(fRnrMarkers);
    ShowFrame(fPMFrame);
    fRnrDaughters->SetState(fTC->GetRnrDaughters() ? kButtonDown : kButtonUp);
    fRnrReferences->SetState(fTC->GetRnrReferences() ? kButtonDown : kButtonUp);
    fRnrDecay->SetState(fTC->GetRnrDecay() ? kButtonDown : kButtonUp);

    fFitDaughters->SetState(fTC->GetFitDaughters() ? kButtonDown : kButtonUp);
    fFitReferences->SetState(fTC->GetFitReferences() ? kButtonDown : kButtonUp);
    fFitDecay->SetState(fTC->GetFitDecay() ? kButtonDown : kButtonUp);
  }
  else 
  {
    HideFrame(fPMFrame);
    ShowFrame(fRnrMarkers);
    fRnrMarkers->SetState(fTC->GetRnrMarkers() ? kButtonDown : kButtonUp);
  }
  fPtRange->SetValues(fTC->GetMinPt(), fTC->GetMaxPt());
  fPRange->SetValues(fTC->GetMinP(), fTC->GetMaxP());
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

void TrackListEditor::DoLineStyle(Int_t style)
{
  fTC->SetStyle(style);
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

void TrackListEditor::DoFitPM()
{
  TGButton* b = (TGButton *) gTQSender;
  PathMark::Type_e type = PathMark::Type_e(b->WidgetId());
  Bool_t on = b->IsOn();

  switch(type)
  {
    case PathMark::Daughter:
      fTC->SetFitDaughters(on);
      break; 
    case PathMark::Reference:
      fTC->SetFitReferences(on);
      break; 
    case PathMark::Decay:
      fTC->SetFitDecay(on);
      break;
     default:
      break;
  }
  Update();
}

void TrackListEditor::DoRnrPM()
{
  TGButton * b = (TGButton *) gTQSender;
  PathMark::Type_e type = PathMark::Type_e(b->WidgetId());
  Bool_t on = b->IsOn();
  switch(type){
    case  PathMark::Daughter:
      fTC->SetRnrDaughters(on);
      break; 
    case  PathMark::Reference:
      fTC->SetRnrReferences(on);
      break; 
    case  PathMark::Decay:
      fTC->SetRnrDecay(on);
      break;
 
    default:
      break;

  }
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
