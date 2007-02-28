// $Header$

#include "PMDModuleEditor.h"
#include <Alieve/PMDModule.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

#include <TCanvas.h>
#include <TGLViewer.h>
#include <Reve/RGTopFrame.h>


#include <TH1F.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// PMDModuleEditor
//

ClassImp(PMDModuleEditor)

PMDModuleEditor::PMDModuleEditor(const TGWindow *p, Int_t width, Int_t height,
				 UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fInfoLabel0(0),
  fInfoLabel1(0),
  fInfoLabel2(0),
  fInfoLabel3(0),
  fInfoLabel4(0),
  fInfoLabel5(0)
  // Initialize widget pointers to 0
{
  MakeTitle("PMDModule");

  Int_t labelW = 67;

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Alieve::PMDModuleEditor", this, "DoXYZZ()");

  fInfoLabel0 = new TGLabel(this);
  fInfoLabel0->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel0, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));

  fInfoLabel1 = new TGLabel(this);
  fInfoLabel1->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel1, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));

  fInfoLabel2 = new TGLabel(this);
  fInfoLabel2->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel2, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));

  fInfoLabel3 = new TGLabel(this);
  fInfoLabel3->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel3, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));

  fInfoLabel4 = new TGLabel(this);
  fInfoLabel4->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel4, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));

  fInfoLabel5 = new TGLabel(this);
  fInfoLabel5->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel5, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));


  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this, 210, 20, kFixedWidth);

    TGHorizontalFrame* g = new TGHorizontalFrame(f, labelW, 0, kFixedWidth);
    TGLabel* l = new TGLabel(g, "Histos:");
    g->AddFrame(l, new TGLayoutHints(kLHintsLeft, 0,0,4,0));
    f->AddFrame(g);

    TGTextButton* b;

    b = new TGTextButton(f, "Show");
    f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    b->Connect("Clicked()", "Alieve::PMDModuleEditor", this, "DisplayHistos()");

    AddFrame(f, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
  }


}

PMDModuleEditor::~PMDModuleEditor()
{}

/**************************************************************************/

void PMDModuleEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<PMDModule*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());

  fInfoLabel0->SetText(Form("Cells hit per Module : %d", fM->GetNPads()));
  fInfoLabel1->SetText(Form("ADC       per Module : %d", fM->GetAdc()));
  fInfoLabel2->SetText(Form("Tot Cells for PRE    : %d", fM->GetPRETotPads()));
  fInfoLabel3->SetText(Form("Tot ADC   for PRE    : %d", fM->GetPRETotAdc()));
  fInfoLabel4->SetText(Form("Tot Cells for CPV    : %d", fM->GetCPVTotPads()));
  fInfoLabel5->SetText(Form("Tot ADC   for CPV    : %d", fM->GetCPVTotAdc()));
}

void PMDModuleEditor::DisplayHistos()
{
  fM->GetHisto()->Draw();
  gPad->Modified();
  gPad->Update();
}




/**************************************************************************/

// Implements callback/slot methods

// void PMDModuleEditor::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }
