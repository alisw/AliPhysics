// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEvePMDModuleEditor.h"
#include <EveDet/AliEvePMDModule.h>

#include <TVirtualPad.h>
#include <TH1F.h>

#include <TGLabel.h>
#include <TG3DLine.h>
#include <TGButton.h>

//==============================================================================
//==============================================================================
// AliEvePMDModuleEditor
//==============================================================================

//______________________________________________________________________________
//
// GUI editor for AliEvePMDModule.

ClassImp(AliEvePMDModuleEditor)

AliEvePMDModuleEditor::AliEvePMDModuleEditor(const TGWindow *p, Int_t width, Int_t height,
                                             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fInfoFrame(0),
  fInfoLabel0(0),
  fInfoLabel1(0),
  fInfoLabel2(0),
  fInfoLabel3(0),
  fInfoLabel4(0),
  fInfoLabel5(0)
{
  // Constructor.

  MakeTitle("AliEvePMDModule");

  CreateInfoFrame();
}

void AliEvePMDModuleEditor::CreateInfoFrame()
{
  // Create frame holding information.

  fInfoFrame = CreateEditorTabSubFrame("Info");

  TGCompositeFrame *title1 = new TGCompositeFrame(fInfoFrame, 145, 10,
						  kHorizontalFrame |
						  kLHintsExpandX   |
						  kFixedWidth      |
						  kOwnBackground);

  title1->AddFrame(new TGLabel(title1, "AliEvePMDModule Info"),
		   new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
  title1->AddFrame(new TGHorizontal3DLine(title1),
		   new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
  fInfoFrame->AddFrame(title1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));

  Int_t labelW = 67;

  fInfoLabel0 = new TGLabel(fInfoFrame);
  fInfoLabel0->SetTextJustify(kTextLeft);
  fInfoFrame->AddFrame(fInfoLabel0, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));

  fInfoLabel1 = new TGLabel(fInfoFrame);
  fInfoLabel1->SetTextJustify(kTextLeft);
  fInfoFrame->AddFrame(fInfoLabel1, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));

  fInfoLabel2 = new TGLabel(fInfoFrame);
  fInfoLabel2->SetTextJustify(kTextLeft);
  fInfoFrame->AddFrame(fInfoLabel2, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));

  fInfoLabel3 = new TGLabel(fInfoFrame);
  fInfoLabel3->SetTextJustify(kTextLeft);
  fInfoFrame->AddFrame(fInfoLabel3, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));

  fInfoLabel4 = new TGLabel(fInfoFrame);
  fInfoLabel4->SetTextJustify(kTextLeft);
  fInfoFrame->AddFrame(fInfoLabel4, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));

  fInfoLabel5 = new TGLabel(fInfoFrame);
  fInfoLabel5->SetTextJustify(kTextLeft);
  fInfoFrame->AddFrame(fInfoLabel5, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));


  {
    TGHorizontalFrame* f = new TGHorizontalFrame(fInfoFrame, 210, 20, kFixedWidth);

    TGHorizontalFrame* g = new TGHorizontalFrame(f, labelW, 0, kFixedWidth);
    TGLabel* l = new TGLabel(g, "Histos:");
    g->AddFrame(l, new TGLayoutHints(kLHintsLeft, 0,0,4,0));
    f->AddFrame(g);

    TGTextButton* b;

    b = new TGTextButton(f, "Show");
    f->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
    b->Connect("Clicked()", "AliEvePMDModuleEditor", this, "DisplayHistos()");

    fInfoFrame->AddFrame(f, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
  }
}

/******************************************************************************/

void AliEvePMDModuleEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = static_cast<AliEvePMDModule*>(obj);

  fInfoLabel0->SetText(Form("Cells hit per Module : %d", fM->GetNPads()));
  fInfoLabel1->SetText(Form("ADC       per Module : %d", fM->GetAdc()));
  fInfoLabel2->SetText(Form("Tot Cells for PRE    : %d", fM->GetPRETotPads()));
  fInfoLabel3->SetText(Form("Tot ADC   for PRE    : %d", fM->GetPRETotAdc()));
  fInfoLabel4->SetText(Form("Tot Cells for CPV    : %d", fM->GetCPVTotPads()));
  fInfoLabel5->SetText(Form("Tot ADC   for CPV    : %d", fM->GetCPVTotAdc()));
}

void AliEvePMDModuleEditor::DisplayHistos()
{
  // Slot for displaying histograms with module data.

  fM->GetHisto()->Draw();
  gPad->Modified();
  gPad->Update();
}
