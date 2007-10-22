// $Header$

#include "DigitSetEditor.h"
#include <Reve/DigitSet.h>

#include <Reve/RGValuators.h>
#include <Reve/ZTransEditor.h>
#include <Reve/RGBAPaletteEditor.h>
#include <Reve/RGEditor.h>

#include <TVirtualPad.h>
#include <TColor.h>
#include <TH1F.h>
#include <TStyle.h>

#include <TGLabel.h>
#include <TG3DLine.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// DigitSetEditor
//

ClassImp(DigitSetEditor)

DigitSetEditor::DigitSetEditor(const TGWindow *p, Int_t width, Int_t height,
			     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fHMTrans   (0),
  fPalette   (0),
  
  fHistoButtFrame(0),
  fInfoFrame(0)
{
  MakeTitle("Transformation matrix");

  fHMTrans = new ZTransSubEditor(this);
  fHMTrans->Connect("UseTrans()",     "Reve::DigitSetEditor", this, "Update()");
  fHMTrans->Connect("TransChanged()", "Reve::DigitSetEditor", this, "Update()");
  AddFrame(fHMTrans, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));


  MakeTitle("Palette controls");

  fPalette = new RGBAPaletteSubEditor(this);
  fPalette->Connect("Changed", "Reve::DigitSetEditor", this, "Update()");
  AddFrame(fPalette, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));

  CreateInfoTab();
}

DigitSetEditor::~DigitSetEditor()
{}

/*************************************************************************/
void DigitSetEditor::CreateInfoTab()
{

  fInfoFrame = CreateEditorTabSubFrame("Info");

  TGCompositeFrame *title1 = new TGCompositeFrame(fInfoFrame, 180, 10, 
						  kHorizontalFrame | 
						  kLHintsExpandX   | 
						  kFixedWidth      | 
						  kOwnBackground);

  title1->AddFrame(new TGLabel(title1, "DigitSet Info"), 
		   new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
  title1->AddFrame(new TGHorizontal3DLine(title1),
		   new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
  fInfoFrame->AddFrame(title1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));

  
  fHistoButtFrame = new TGHorizontalFrame(fInfoFrame);
  TGTextButton* b = 0;
  b = new TGTextButton(fHistoButtFrame, "Histo");
  b->SetToolTipText("Show histogram over full range.");
  fHistoButtFrame->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
  b->Connect("Clicked()", "Reve::DigitSetEditor", this, "DoHisto()");

  b = new TGTextButton(fHistoButtFrame, "Range Histo");
  b->SetToolTipText("Show histogram over selected range.");
  fHistoButtFrame->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
  b->Connect("Clicked()", "Reve::DigitSetEditor", this, "DoRangeHisto()");
  fInfoFrame->AddFrame(fHistoButtFrame, new TGLayoutHints(kLHintsExpandX, 2, 0, 0, 0));
}

/**************************************************************************/

void DigitSetEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<DigitSet*>(obj);

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

void DigitSetEditor::DoHisto()
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

void DigitSetEditor::DoRangeHisto()
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

void DigitSetEditor::PlotHisto(Int_t min, Int_t max)
{
  Int_t nbins = max-min+1;
  while (nbins > 200)
    nbins /= 2;

  TH1F* h = new TH1F(fM->GetName(), fM->GetTitle(), nbins, min-0.5, max+0.5);
  h->SetDirectory(0);
  h->SetBit(kCanDelete);
  VoidCPlex::iterator qi(fM->fPlex);
  while (qi.next())
    h->Fill(((DigitSet::DigitBase*)qi())->fValue);

  gStyle->SetOptStat(1111111);
  h->Draw();
  gPad->Modified();
  gPad->Update();
}
