// $Header$

#include "ITSScaledModuleEditor.h"
#include <Alieve/ITSScaledModule.h>
#include <Reve/ZTransEditor.h>
#include <Reve/RGValuators.h>

#include <TVirtualPad.h>
#include <TColor.h>
#include <TMath.h>

#include <TGedEditor.h>
#include <TGLabel.h>
#include <TG3DLine.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGComboBox.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// ITSScaledModuleEditor 
//

ClassImp(ITSScaledModuleEditor)

  ITSScaledModuleEditor::ITSScaledModuleEditor(const TGWindow *p, Int_t width, Int_t height,
					       UInt_t options, Pixel_t back) :
    TGedFrame(p, width, height, options | kVerticalFrame, back),

    fInfoFrame(0),

    fModule(0), 

    fScale(0),
    fStatistic(0),
    fInfoLabel0(0),
    fInfoLabel1(0)
{
  MakeTitle("ITSScaledModule");
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
    fScale->Connect("ValueSet(Long_t)", "Alieve::ITSScaledModuleEditor", this, "DoScale()");
    
    TGLabel* lab = new TGLabel(f, "Statistic:");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 2, 1, 2));
    fStatistic = new TGComboBox(f);
    fStatistic->AddEntry("Occup", 0);
    fStatistic->AddEntry("Average", 1);
    fStatistic->AddEntry("RMS", 2);
    TGListBox* lb = fStatistic->GetListBox();
    lb->Resize(lb->GetWidth(), 3*16);
    fStatistic->Resize(74, 20);
    fStatistic->Connect("Selected(Int_t)", "Alieve::ITSScaledModuleEditor", this, "DoStatType(Int_t)");
    f->AddFrame(fStatistic, new TGLayoutHints(kLHintsLeft, 1, 2, 1, 1));
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  CreateInfoFrame();
}

/*************************************************************************/
ITSScaledModuleEditor::~ITSScaledModuleEditor()
{}

/*************************************************************************/
void ITSScaledModuleEditor::CreateInfoFrame()
{
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

/**************************************************************************/

void ITSScaledModuleEditor::SetModel(TObject* obj)
{
  fModule = dynamic_cast<ITSScaledModule*>(obj); 

  // widgets
  fScale->SetIntNumber(fModule->GetScaleInfo()->GetScale());
  fStatistic->Select(fModule->GetScaleInfo()->GetStatType(), kFALSE);

  // text info  
  Int_t cnx, cnz, total;
  fModule->GetScaleData(cnx, cnz, total);
  fInfoLabel0->SetText(Form("Cell size:  Nx=%d Nz=%d", cnx, cnz));
  fInfoLabel1->SetText(Form("Num cells:  %d", total));
}


/**************************************************************************/

void ITSScaledModuleEditor::DoScale()
{
  fModule->GetScaleInfo()->ScaleChanged(fScale->GetIntNumber());

  Int_t cnx, cnz, total;
  fModule->GetScaleData(cnx, cnz, total);
  fInfoLabel0->SetText(Form("Cell size:  Nx=%d Nz=%d", cnx, cnz));
  Update();
  fGedEditor->SetModel(fGedEditor->GetPad(), fGedEditor->GetModel(), kButton1Down);
}

/**************************************************************************/

void ITSScaledModuleEditor::DoStatType(Int_t v)
{
  fModule->GetScaleInfo()->StatTypeChanged(v);
  Update();
  fGedEditor->SetModel(fGedEditor->GetPad(), fGedEditor->GetModel(), kButton1Down);
}
