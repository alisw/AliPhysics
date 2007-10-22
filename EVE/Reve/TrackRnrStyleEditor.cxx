// $Header$

#include "TrackRnrStyleEditor.h"
#include <Reve/Track.h>

#include <Reve/RGValuators.h>
#include <Reve/ReveManager.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TG3DLine.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>
#include <TGComboBox.h>
#include <TAttMarkerEditor.h>

using namespace Reve;

//______________________________________________________________________
// TrackRnrStyleSubEditor
//
//

ClassImp(TrackRnrStyleSubEditor)

TrackRnrStyleSubEditor::TrackRnrStyleSubEditor(const TGWindow *p):
    TGVerticalFrame(p),
    fM (0),

    fMaxR(0),
    fMaxZ(0),
    fMaxOrbits(0),
    fMinAng(0),
    fDelta(0),

    fRnrFV(0),

    fPMFrame(0),
    fFitDaughters(0),
    fFitReferences(0),
    fFitDecay(0),
    fRnrDaughters(0),
    fRnrReferences(0),
    fRnrDecay(0),

    fRefsCont(0),
    fPMAtt(0),
    fFVAtt(0)
{
  Int_t labelW = 51;

  // --- Limits
  fMaxR = new RGValuator(this, "Max R:", 90, 0);
  fMaxR->SetLabelWidth(labelW);
  fMaxR->SetNELength(6);
  fMaxR->Build();
  fMaxR->SetLimits(0.1, 1000, 101, TGNumberFormat::kNESRealOne);
  fMaxR->SetToolTip("Maximum radius to which the tracks will be drawn.");
  fMaxR->Connect("ValueSet(Double_t)", "Reve::TrackRnrStyleSubEditor", this, "DoMaxR()");
  AddFrame(fMaxR, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  fMaxZ = new RGValuator(this, "Max Z:", 90, 0);
  fMaxZ->SetLabelWidth(labelW);
  fMaxZ->SetNELength(6);
  fMaxZ->Build();
  fMaxZ->SetLimits(0.1, 2000, 101, TGNumberFormat::kNESRealOne);
  fMaxZ->SetToolTip("Maximum z-coordinate to which the tracks will be drawn.");
  fMaxZ->Connect("ValueSet(Double_t)", "Reve::TrackRnrStyleSubEditor", this, "DoMaxZ()");
  AddFrame(fMaxZ, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  fMaxOrbits = new RGValuator(this, "Orbits:", 90, 0);
  fMaxOrbits->SetLabelWidth(labelW);
  fMaxOrbits->SetNELength(6);
  fMaxOrbits->Build();
  fMaxOrbits->SetLimits(0.1, 10, 101, TGNumberFormat::kNESRealOne);
  fMaxOrbits->SetToolTip("Maximal angular path of tracks' orbits (1 ~ 2Pi).");
  fMaxOrbits->Connect("ValueSet(Double_t)", "Reve::TrackRnrStyleSubEditor", this, "DoMaxOrbits()");  
  AddFrame(fMaxOrbits, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  fMinAng = new RGValuator(this, "Angle:", 90, 0);
  fMinAng->SetLabelWidth(labelW);
  fMinAng->SetNELength(6);
  fMinAng->Build();
  fMinAng->SetLimits(1, 160, 81, TGNumberFormat::kNESRealOne);
  fMinAng->SetToolTip("Minimal angular step between two helix points.");
  fMinAng->Connect("ValueSet(Double_t)", "Reve::TrackRnrStyleSubEditor", this, "DoMinAng()");
  AddFrame(fMinAng, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  fDelta = new RGValuator(this, "Delta:", 90, 0);
  fDelta->SetLabelWidth(labelW);
  fDelta->SetNELength(6);
  fDelta->Build();
  fDelta->SetLimits(0.001, 10, 101, TGNumberFormat::kNESRealThree);
  fDelta->SetToolTip("Maximal error at the mid-point of the line connecting to helix points.");
  fDelta->Connect("ValueSet(Double_t)", "Reve::TrackRnrStyleSubEditor", this, "DoDelta()");
  AddFrame(fDelta, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
}

//______________________________________________________________________
void TrackRnrStyleSubEditor::CreateRefsContainer(TGVerticalFrame* p)
{
  fRefsCont = new TGCompositeFrame(p, 80, 20, kVerticalFrame);
  fPMFrame  = new TGVerticalFrame(fRefsCont);
  {
    // --- Rendering control
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

    fFitDecay->Connect("Clicked()","Reve::TrackRnrStyleSubEditor", this, "DoFitPM()");  
    fFitReferences->Connect("Clicked()","Reve::TrackRnrStyleSubEditor", this, "DoFitPM()");  
    fFitDaughters->Connect("Clicked()","Reve::TrackRnrStyleSubEditor", this, "DoFitPM()");
  }
  { 
    // --- Kinematics fitting
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

    fRnrDecay->Connect("Clicked()","Reve::TrackRnrStyleSubEditor", this, "DoRnrPM()");  
    fRnrReferences->Connect("Clicked()","Reve::TrackRnrStyleSubEditor", this, "DoRnrPM()");  
    fRnrDaughters->Connect("Clicked()","Reve::TrackRnrStyleSubEditor", this, "DoRnrPM()");
  
    fRefsCont->AddFrame(fPMFrame, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }
  // marker attributes
  {
    fPMAtt = new TAttMarkerEditor(fRefsCont);
    TGFrameElement *el = (TGFrameElement*) fPMAtt->GetList()->First();
    TGFrame *f = el->fFrame; fPMAtt->RemoveFrame(f); delete f;
    fRefsCont->AddFrame(fPMAtt, new TGLayoutHints(kLHintsTop, 1, 1, 3, 1));
  }

  // FIRST VERTEX
  TGCompositeFrame *title1 = new TGCompositeFrame(fRefsCont, 145, 10, 
						  kHorizontalFrame | 
						  kLHintsExpandX   | 
						  kFixedWidth      | 
						  kOwnBackground);
  title1->AddFrame(new TGLabel(title1, "FirstVertex"), 
		   new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
  title1->AddFrame(new TGHorizontal3DLine(title1),
		   new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 5));
  fRefsCont->AddFrame(title1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));

  fRnrFV = new TGCheckButton(fRefsCont, "Rnr");
  fRnrFV->Connect("Clicked()","Reve::TrackRnrStyleSubEditor", this, "DoRnrFV()");
  fRefsCont->AddFrame(fRnrFV, new TGLayoutHints(kLHintsTop, 5, 1, 2, 0));
  {
    fFVAtt = new TAttMarkerEditor(fRefsCont);
    TGFrameElement *el = (TGFrameElement*) fFVAtt->GetList()->First();
    TGFrame *f = el->fFrame; fFVAtt->RemoveFrame(f); delete f;
    fRefsCont->AddFrame(fFVAtt, new TGLayoutHints(kLHintsTop, 1, 1, 3, 1));
  }
  p->AddFrame(fRefsCont,new TGLayoutHints(kLHintsTop| kLHintsExpandX));
}

//______________________________________________________________________
void TrackRnrStyleSubEditor::SetModel(TrackRnrStyle* m)
{
  fM = m;

  fMaxR->SetValue(fM->fMaxR);
  fMaxZ->SetValue(fM->fMaxZ);
  fMaxOrbits->SetValue(fM->fMaxOrbs);
  fMinAng->SetValue(fM->fMinAng);
  fDelta->SetValue(fM->fDelta);

  if(fM->fEditPathMarks) 
  {
    ShowFrame(fPMFrame);
    fRnrDaughters->SetState(fM->fRnrDaughters ? kButtonDown : kButtonUp);
    fRnrReferences->SetState(fM->fRnrReferences ? kButtonDown : kButtonUp);
    fRnrDecay->SetState(fM->fRnrDecay ? kButtonDown : kButtonUp);

    fFitDaughters->SetState(fM->fFitDaughters ? kButtonDown : kButtonUp);
    fFitReferences->SetState(fM->fFitReferences ? kButtonDown : kButtonUp);
    fFitDecay->SetState(fM->fFitDecay ? kButtonDown : kButtonUp);

    fPMAtt->SetModel(&fM->fPMAtt);
  }
  else 
  {
   fRefsCont->HideFrame(fPMFrame);
  }

  fRnrFV->SetState(fM->fRnrFV ? kButtonDown : kButtonUp);
  fFVAtt->SetModel(&fM->fFVAtt);
}

/**************************************************************************/
void TrackRnrStyleSubEditor::Changed()
{
  fM->UpdateBackPtrItems();
  Emit("Changed()");
}

/**************************************************************************/
void TrackRnrStyleSubEditor::DoMaxR()
{
  fM->SetMaxR(fMaxR->GetValue());
  Changed();
}

void TrackRnrStyleSubEditor::DoMaxZ()
{
  fM->SetMaxZ(fMaxZ->GetValue());
  Changed();
}

void TrackRnrStyleSubEditor::DoMaxOrbits()
{
  fM->SetMaxOrbs(fMaxOrbits->GetValue());
  Changed();
}

void TrackRnrStyleSubEditor::DoMinAng()
{
  fM->SetMinAng(fMinAng->GetValue());
  Changed();
}

void TrackRnrStyleSubEditor::DoDelta()
{
  fM->SetDelta(fDelta->GetValue());
  Changed();
}

/**************************************************************************/
void TrackRnrStyleSubEditor::DoFitPM()
{
  TGButton* b = (TGButton *) gTQSender;
  PathMark::Type_e type = PathMark::Type_e(b->WidgetId());
  Bool_t on = b->IsOn();

  switch(type)
  {
    case PathMark::Daughter:
      fM->SetFitDaughters(on);
      break; 
    case PathMark::Reference:
      fM->SetFitReferences(on);
      break; 
    case PathMark::Decay:
      fM->SetFitDecay(on);
      break;
    default:
      break;
  }
  Changed();
}

void TrackRnrStyleSubEditor::DoRnrPM()
{
  TGButton * b = (TGButton *) gTQSender;
  PathMark::Type_e type = PathMark::Type_e(b->WidgetId());
  Bool_t on = b->IsOn();
  switch(type){
    case  PathMark::Daughter:
      fM->SetRnrDaughters(on);
      break; 
    case  PathMark::Reference:
      fM->SetRnrReferences(on);
      break; 
    case  PathMark::Decay:
      fM->SetRnrDecay(on);
      break;
 
    default:
      break;
  }
  Changed();
}

void TrackRnrStyleSubEditor::DoRnrFV()
{
  fM->SetRnrFV(fRnrFV->IsOn());
  Changed();
}


//______________________________________________________________________
// TrackRnrStyleEditor
//

ClassImp(TrackRnrStyleEditor)

TrackRnrStyleEditor::TrackRnrStyleEditor(const TGWindow *p,
                                         Int_t width, Int_t height,
                                         UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),

  fM (0),
  fRSSubEditor(0)
{
  MakeTitle("RenderStyle");

  fRSSubEditor = new TrackRnrStyleSubEditor(this);
  fRSSubEditor->Connect("Changed()", "Reve::TrackRnrStyleEditor", this, "Update()"); 
  AddFrame(fRSSubEditor, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0,0,0));

  TGVerticalFrame* refsFrame = CreateEditorTabSubFrame("Refs");
  TGCompositeFrame *title1 = new TGCompositeFrame(refsFrame, 145, 10, 
						  kHorizontalFrame | 
						  kLHintsExpandX   | 
						  kFixedWidth      | 
						  kOwnBackground);
  title1->AddFrame(new TGLabel(title1, "PathMarks"), 
		   new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
  title1->AddFrame(new TGHorizontal3DLine(title1),
		   new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
  refsFrame->AddFrame(title1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));
   
  // path marks 
  fRSSubEditor->CreateRefsContainer(refsFrame);
  fRSSubEditor->fPMAtt->SetGedEditor((TGedEditor*)gReve->GetEditor()); 
  fRSSubEditor->fFVAtt->SetGedEditor((TGedEditor*)gReve->GetEditor());

}

TrackRnrStyleEditor::~TrackRnrStyleEditor()
{
}

/**************************************************************************/

void TrackRnrStyleEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<TrackRnrStyle*>(obj); 
  fRSSubEditor->SetModel(fM);
}
