// $Header$

#include "LineEditor.h"
#include <Reve/Line.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// LineEditor
//

ClassImp(LineEditor)

  LineEditor::LineEditor(const TGWindow *p, Int_t width, Int_t height,
			 UInt_t options, Pixel_t back) :
    TGedFrame(p, width, height, options | kVerticalFrame, back),
    fM(0),
    fRnrLine   (0),
    fRnrPoints (0)
{
  // MakeTitle("Line");

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);

    fRnrLine  = new TGCheckButton(f, "Draw Line");
    f->AddFrame(fRnrLine, new TGLayoutHints(kLHintsLeft, 1,2,0,0));
    fRnrLine->Connect("Toggled(Bool_t)", "Reve::LineEditor", this, "DoRnrLine()");
    fRnrPoints = new TGCheckButton(f, "Draw Marker");
    f->AddFrame(fRnrPoints, new TGLayoutHints(kLHintsLeft, 2,1,0,0));  
    fRnrPoints->Connect("Toggled(Bool_t)"," Reve::LineEditor", this, "DoRnrPoints()");

    AddFrame(f, new TGLayoutHints(kLHintsTop, 0,0,2,1));
  }
}

LineEditor::~LineEditor()
{}

/**************************************************************************/

void LineEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<Line*>(obj);

  fRnrLine  ->SetState(fM->fRnrLine  ? kButtonDown : kButtonUp);
  fRnrPoints->SetState(fM->fRnrPoints ? kButtonDown : kButtonUp);
}

/**************************************************************************/

void LineEditor::DoRnrLine()
{
  fM->SetRnrLine(fRnrLine->IsOn());
  Update();
}

void LineEditor::DoRnrPoints()
{
  fM->SetRnrPoints(fRnrPoints->IsOn());
  Update();
}
