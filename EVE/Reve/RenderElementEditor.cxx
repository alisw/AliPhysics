// $Header$

#include "RenderElementEditor.h"
#include <Reve/RenderElement.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// RenderElementEditor
//

ClassImp(RenderElementEditor)

  RenderElementEditor::RenderElementEditor(const TGWindow *p,
					   Int_t width, Int_t height,
					   UInt_t options, Pixel_t back) :
    TGedFrame(p, width, height, options | kVerticalFrame, back),

    fRE           (0),
    fHFrame       (0),
    fRnrSelf   (0),
    fRnrChildren (0),
    fMainColor    (0)
{
  MakeTitle("RenderElement");
  fPriority = 0;

  fHFrame = new TGHorizontalFrame(this);

  fRnrSelf = new TGCheckButton(fHFrame, "RnrSelf");
  fHFrame->AddFrame(fRnrSelf, new TGLayoutHints(kLHintsLeft, 1, 2, 1, 1));
  fRnrSelf->Connect
    ("Toggled(Bool_t)",
     "Reve::RenderElementEditor", this, "DoRnrSelf()");

  fRnrChildren = new TGCheckButton(fHFrame, "RnrChildren");
  fHFrame->AddFrame(fRnrChildren, new TGLayoutHints(kLHintsLeft, 2, 1, 1, 1));
  fRnrChildren->Connect
    ("Toggled(Bool_t)",
     "Reve::RenderElementEditor", this, "DoRnrChildren()");

  fMainColor = new TGColorSelect(fHFrame, 0, -1);
  fHFrame->AddFrame(fMainColor, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
  fMainColor->Connect
    ("ColorSelected(Pixel_t)",
     "Reve::RenderElementEditor", this, "DoMainColor(Pixel_t)");

  AddFrame(fHFrame, new TGLayoutHints(kLHintsTop, 0, 0, 1, 1));    
}

RenderElementEditor::~RenderElementEditor()
{}

/**************************************************************************/

void RenderElementEditor::SetModel(TObject* obj)
{
  fRE = dynamic_cast<RenderElement*>(obj);

  if (fRE->CanEditRnrElement()) {
    fRnrSelf->SetState(fRE->GetRnrSelf() ? kButtonDown : kButtonUp);
    fRnrChildren->SetState(fRE->GetRnrChildren() ? kButtonDown : kButtonUp);
    fRnrSelf->MapWindow();
    fRnrChildren->MapWindow();
  } else {
    fRnrSelf->UnmapWindow();
    fRnrChildren->UnmapWindow();
  }

  if (fRE->CanEditMainColor()) {
    fMainColor->SetColor(TColor::Number2Pixel(fRE->GetMainColor()), kFALSE);
    fMainColor->MapWindow();
  } else {
    fMainColor->UnmapWindow();
  }

  fHFrame->Layout();
}

/**************************************************************************/


void RenderElementEditor::DoRnrSelf()
{
  fRE->SetRnrSelf(fRnrSelf->IsOn());
  Update();
}


void RenderElementEditor::DoRnrChildren()
{
  fRE->SetRnrChildren(fRnrChildren->IsOn());
  Update();
}

void RenderElementEditor::DoMainColor(Pixel_t color)
{
  fRE->SetMainColor(color);
  Update();
}
