// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveITSScaledModuleEditor.h"
#include <EveDet/AliEveITSScaledModule.h>

#include <TVirtualPad.h>
#include <TGedEditor.h>
#include <TGLabel.h>
#include <TG3DLine.h>
#include <TGNumberEntry.h>
#include <TGComboBox.h>

//==============================================================================
//==============================================================================
// AliEveITSScaledModuleEditor
//==============================================================================

//______________________________________________________________________________
//
// Editor for AliEveITSScaledModule.

ClassImp(AliEveITSScaledModuleEditor)

AliEveITSScaledModuleEditor::AliEveITSScaledModuleEditor(const TGWindow *p, Int_t width, Int_t height,
                                                         UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fModule(0),
  fScale(0),
  fStatistic(0),
  fInfoFrame(0),
  fInfoLabel0(0),
  fInfoLabel1(0)
{
  // Constructor.

  MakeTitle("AliEveITSScaledModule");
  // Create widgets
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Scale:");
    f->AddFrame(l, new TGLayoutHints(kLHintsTop | kLHintsCenterY, 0, 5, 1, 1));
    fScale = new TGNumberEntry(f, 0, 2, -1,
			       TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive,
			       TGNumberFormat::kNELLimitMinMax, 1, 5);
    fScale->GetNumberEntry()->SetToolTipText("Set cell size.");
    f->AddFrame(fScale, new TGLayoutHints(kLHintsLeft, 1, 7, 1, 1));
    fScale->Associate(f);
    fScale->Connect("ValueSet(Long_t)", "AliEveITSScaledModuleEditor", this, "DoScale()");

    TGLabel* lab = new TGLabel(f, "Statistic:");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 2, 1, 2));
    fStatistic = new TGComboBox(f);
    fStatistic->AddEntry("Occup", 0);
    fStatistic->AddEntry("Average", 1);
    fStatistic->AddEntry("RMS", 2);
    TGListBox* lb = fStatistic->GetListBox();
    lb->Resize(lb->GetWidth(), 3*16);
    fStatistic->Resize(74, 20);
    fStatistic->Connect("Selected(Int_t)", "AliEveITSScaledModuleEditor", this, "DoStatType(Int_t)");
    f->AddFrame(fStatistic, new TGLayoutHints(kLHintsLeft, 1, 2, 1, 1));
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  CreateInfoFrame();
}

/******************************************************************************/

void AliEveITSScaledModuleEditor::CreateInfoFrame()
{
  // Create a frame under tab "Info".

  fInfoFrame = CreateEditorTabSubFrame("Info");
  TGCompositeFrame *title1 = new TGCompositeFrame(fInfoFrame, 145, 10,
						  kHorizontalFrame |
						  kLHintsExpandX   |
						  kFixedWidth      |
						  kOwnBackground);

  title1->AddFrame(new TGLabel(title1, "ScaledDigits Info"),
		   new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
  title1->AddFrame(new TGHorizontal3DLine(title1),
		   new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
  fInfoFrame->AddFrame(title1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));


  Int_t lp = 2;
  fInfoLabel0 = new TGLabel(fInfoFrame);
  fInfoLabel0->SetTextJustify(kTextLeft);
  fInfoFrame->AddFrame(fInfoLabel0, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                                      lp, 0, 8, 0));

  fInfoLabel1 = new TGLabel(fInfoFrame);
  fInfoLabel1->SetTextJustify(kTextLeft);
  fInfoFrame->AddFrame(fInfoLabel1, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                                      lp, 0, 2, 8));
}

/******************************************************************************/

void AliEveITSScaledModuleEditor::SetModel(TObject* obj)
{
  // Set model object.

  fModule = static_cast<AliEveITSScaledModule*>(obj);

  // widgets
  fScale->SetIntNumber(fModule->GetScaleInfo()->GetScale());
  fStatistic->Select(fModule->GetScaleInfo()->GetStatType(), kFALSE);

  // text info
  Int_t cnx, cnz, total;
  fModule->GetScaleData(cnx, cnz, total);
  fInfoLabel0->SetText(Form("Cell size:  Nx=%d Nz=%d", cnx, cnz));
  fInfoLabel1->SetText(Form("Num cells:  %d", total));
}


/******************************************************************************/

void AliEveITSScaledModuleEditor::DoScale()
{
  // Slot for Scale.

  fModule->GetScaleInfo()->ScaleChanged(fScale->GetIntNumber());

  Int_t cnx, cnz, total;
  fModule->GetScaleData(cnx, cnz, total);
  fInfoLabel0->SetText(Form("Cell size:  Nx=%d Nz=%d", cnx, cnz));
  Update();
  fGedEditor->SetModel(fGedEditor->GetPad(), fGedEditor->GetModel(), kButton1Down);
}

/******************************************************************************/

void AliEveITSScaledModuleEditor::DoStatType(Int_t v)
{
  // Slot for StatType.

  fModule->GetScaleInfo()->StatTypeChanged(v);
  Update();
  fGedEditor->SetModel(fGedEditor->GetPad(), fGedEditor->GetModel(), kButton1Down);
}
