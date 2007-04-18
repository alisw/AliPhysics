// $Header$

#include "QuadSetEditor.h"
#include <Reve/QuadSet.h>

#include <Reve/RGValuators.h>
#include <Reve/ZTransEditor.h>
#include <Reve/RGBAPaletteEditor.h>

#include <TVirtualPad.h>
#include <TColor.h>
#include <TH1F.h>
#include <TStyle.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// QuadSetEditor
//

ClassImp(QuadSetEditor)

QuadSetEditor::QuadSetEditor(const TGWindow *p, Int_t width, Int_t height,
			     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fHMTrans   (0),
  fPalette   (0)
  // Initialize widget pointers to 0
{
  MakeTitle("Transformation matrix");

  fHMTrans = new ZTransSubEditor(this);
  fHMTrans->Connect("UseTrans()",     "Reve::QuadSetEditor", this, "Update()");
  fHMTrans->Connect("TransChanged()", "Reve::QuadSetEditor", this, "Update()");
  AddFrame(fHMTrans, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));


  MakeTitle("Palette controls");

  fPalette = new RGBAPaletteSubEditor(this);
  fPalette->Connect("Changed", "Reve::QuadSetEditor", this, "Update()");
  AddFrame(fPalette, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));


  MakeTitle("QuadSet");

  fHistoButtFrame = new TGHorizontalFrame(this);
  {
    TGTextButton* b = 0;

    b = new TGTextButton(fHistoButtFrame, "Histo");
    b->SetToolTipText("Show histogram over full range.");
    fHistoButtFrame->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    b->Connect("Clicked()", "Reve::QuadSetEditor", this, "DoHisto()");

    b = new TGTextButton(fHistoButtFrame, "Range Histo");
    b->SetToolTipText("Show histogram over selected range.");
    fHistoButtFrame->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    b->Connect("Clicked()", "Reve::QuadSetEditor", this, "DoRangeHisto()");
  }
  AddFrame(fHistoButtFrame, new TGLayoutHints(kLHintsExpandX, 2, 0, 0, 0));
}

QuadSetEditor::~QuadSetEditor()
{}

/**************************************************************************/

void QuadSetEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<QuadSet*>(obj);

  fHMTrans->SetDataFromTrans(&fM->fHMTrans);

  if (fM->fValueIsColor || fM->fPalette == 0) {
    fPalette->UnmapWindow();
  } else {
    fPalette->SetModel(fM->fPalette);
    fPalette->MapWindow();
  }


  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());

  if (fM->fHistoButtons)
    fHistoButtFrame->MapWindow();
  else
    fHistoButtFrame->UnmapWindow();
}

/**************************************************************************/

void QuadSetEditor::DoHisto()
{
  Int_t min, max;
  if (fM->fPalette) {
    min = fM->fPalette->GetLowLimit();
    max = fM->fPalette->GetHighLimit();
  } else {
    fM->ScanMinMaxValues(min, max);
  }
  PlotHisto(min, max);
}

void QuadSetEditor::DoRangeHisto()
{
  Int_t min, max;
  if (fM->fPalette) {
    min = fM->fPalette->GetMinVal();
    max = fM->fPalette->GetMaxVal();
  } else {
    fM->ScanMinMaxValues(min, max);
  }
  PlotHisto(min, max);
}

void QuadSetEditor::PlotHisto(Int_t min, Int_t max)
{
  Int_t nbins = max-min+1;
  while (nbins > 200)
    nbins /= 2;

  TH1F* h = new TH1F(fM->GetName(), fM->GetTitle(), nbins, min-0.5, max+0.5);
  h->SetDirectory(0);
  h->SetBit(kCanDelete);
  VoidCPlex::iterator qi(fM->fPlex);
  while (qi.next())
    h->Fill(((QuadSet::QuadBase*)qi())->fValue);

  gStyle->SetOptStat(1111111);
  h->Draw();
  gPad->Modified();
  gPad->Update();
}
