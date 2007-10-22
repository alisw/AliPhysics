// $Header$

#include "GeoNodeEditor.h"
#include "RGValuators.h"

#include <Reve/GeoNode.h>
#include <TGeoNode.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// GeoNodeRnrElEditor
//

ClassImp(GeoNodeRnrElEditor)

GeoNodeRnrElEditor::GeoNodeRnrElEditor(const TGWindow *p,
                                       Int_t width, Int_t height,
                                       UInt_t options, Pixel_t back) :
  TGedFrame(p,width, height, options | kVerticalFrame, back),

  fNodeRE (0),

  fVizNode(0),
  fVizNodeDaughters(0),
  fVizVolume(0),
  fVizVolumeDaughters(0),

  fTransparency(0)
{
  MakeTitle("GeoNode");

  // --- Visibility control

  fVizNode = new TGCheckButton(this, "VizNode");
  AddFrame(fVizNode, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fVizNode->Connect
    ("Toggled(Bool_t)",
     "Reve::GeoNodeRnrElEditor", this, "DoVizNode()");

  fVizNodeDaughters = new TGCheckButton(this, "VizNodeDaughters");
  AddFrame(fVizNodeDaughters, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fVizNodeDaughters->Connect
    ("Toggled(Bool_t)",
     "Reve::GeoNodeRnrElEditor", this, "DoVizNodeDaughters()");

  fVizVolume = new TGCheckButton(this, "VizVolume");
  AddFrame(fVizVolume, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fVizVolume->Connect
    ("Toggled(Bool_t)",
     "Reve::GeoNodeRnrElEditor", this, "DoVizVolume()");

  fVizVolumeDaughters = new TGCheckButton(this, "VizVolumeDaughters");
  AddFrame(fVizVolumeDaughters, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fVizVolumeDaughters->Connect
    ("Toggled(Bool_t)",
     "Reve::GeoNodeRnrElEditor", this, "DoVizVolumeDaughters()");


  // --- Color props

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Transparency:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 1, 2, 1, 1));
    fTransparency = new TGNumberEntry(f, 0., 6, -1, 
			      TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
			      TGNumberFormat::kNELLimitMinMax, 0, 100);
    fTransparency->GetNumberEntry()->SetToolTipText("0 is opaque, 100 fully transparent.");
    f->AddFrame(fTransparency, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fTransparency->Connect("ValueSet(Long_t)",
			   "Reve::GeoNodeRnrElEditor", this, "DoTransparency()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }
}

/**************************************************************************/

void GeoNodeRnrElEditor::SetModel(TObject* obj)
{
  fNodeRE = dynamic_cast<GeoNodeRnrEl*>(obj);
  TGeoNode*  node = fNodeRE->fNode;
  TGeoVolume* vol = node->GetVolume();

  fVizNode->SetState(node->TGeoAtt::IsVisible() ? kButtonDown : kButtonUp);
  fVizNodeDaughters->SetState(node->TGeoAtt::IsVisDaughters() ? kButtonDown : kButtonUp);
  fVizVolume->SetState(vol->IsVisible() ? kButtonDown : kButtonUp);
  fVizVolumeDaughters->SetState(vol->IsVisDaughters() ? kButtonDown : kButtonUp);

  fTransparency->SetNumber(vol->GetTransparency());
}

/**************************************************************************/

void GeoNodeRnrElEditor::DoVizNode()
{
  fNodeRE->SetRnrSelf(fVizNode->IsOn());
  Update();
}

void GeoNodeRnrElEditor::DoVizNodeDaughters()
{
  fNodeRE->SetRnrChildren(fVizNodeDaughters->IsOn());
  Update();
}

void GeoNodeRnrElEditor::DoVizVolume()
{
  fNodeRE->fNode->GetVolume()->SetVisibility(fVizVolume->IsOn());
  Update();
}

void GeoNodeRnrElEditor::DoVizVolumeDaughters()
{
  fNodeRE->fNode->GetVolume()->VisibleDaughters(fVizVolumeDaughters->IsOn());
  Update();
}

/**************************************************************************/

void GeoNodeRnrElEditor::DoTransparency()
{
  fNodeRE->fNode->GetVolume()->SetTransparency(char(fTransparency->GetNumber()));
  Update();
}


/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// GeoTopNodeRnrElEditor
//

ClassImp(GeoTopNodeRnrElEditor)

GeoTopNodeRnrElEditor::GeoTopNodeRnrElEditor(const TGWindow *p,
                                             Int_t width, Int_t height,
					     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),

  fTopNodeRE (0),
  fVisOption (0),
  fVisLevel  (0)
{
  MakeTitle("GeoTopNode");

  Int_t labelW = 58;

  fVisOption = new RGValuator(this, "VisOption:", 90, 0);
  fVisOption->SetLabelWidth(labelW);
  fVisOption->SetShowSlider(kFALSE);
  fVisOption->SetNELength(4);
  fVisOption->Build();
  fVisOption->SetLimits(0, 2, 10, TGNumberFormat::kNESInteger);
  fVisOption->SetToolTip("Visualization option passed to TGeoPainter.");
  fVisOption->Connect("ValueSet(Double_t)", "Reve::GeoTopNodeRnrElEditor", this, "DoVisOption()");
  AddFrame(fVisOption, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  fVisLevel = new RGValuator(this, "VisLevel:", 90, 0);
  fVisLevel->SetLabelWidth(labelW);
  fVisLevel->SetShowSlider(kFALSE);
  fVisLevel->SetNELength(4);
  fVisLevel->Build();
  fVisLevel->SetLimits(0, 10, 10, TGNumberFormat::kNESInteger);
  fVisLevel->SetToolTip("Level (depth) to which the geometry is traversed.");
  fVisLevel->Connect("ValueSet(Double_t)", "Reve::GeoTopNodeRnrElEditor", this, "DoVisLevel()");
  AddFrame(fVisLevel, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
}

/**************************************************************************/

void GeoTopNodeRnrElEditor::SetModel(TObject* obj)
{
  fTopNodeRE = dynamic_cast<GeoTopNodeRnrEl*>(obj);

  fVisOption->SetValue(fTopNodeRE->GetVisOption());
  fVisLevel ->SetValue(fTopNodeRE->GetVisLevel());
}

/**************************************************************************/

void GeoTopNodeRnrElEditor::DoVisOption()
{
  fTopNodeRE->SetVisOption(Int_t(fVisOption->GetValue()));
  Update();
}

void GeoTopNodeRnrElEditor::DoVisLevel()
{
  fTopNodeRE->SetVisLevel(Int_t(fVisLevel->GetValue()));
  Update();
}
