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
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGComboBox.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// SDPaletteSubEditor
//
ITSSDSubEditor::ITSSDSubEditor(const TGWindow* p) :
  RGBAPaletteSubEditor(p),
  fModule(0),
  fScale(0),
  fStatistic(0),
  fInfoLabel0(0),
  fInfoLabel1(0),
  fInfoLabel2(0)
{
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
    fScale->Connect("ValueSet(Long_t)", "Alieve::ITSSDSubEditor", this, "DoScale()");
    
    TGLabel* lab = new TGLabel(f, "Statistic:");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 2, 1, 2));
    fStatistic = new TGComboBox(f);
    fStatistic->AddEntry("Occup", 0);
    fStatistic->AddEntry("Average", 1);
    fStatistic->AddEntry("RMS", 2);
    TGListBox* lb = fStatistic->GetListBox();
    lb->Resize(lb->GetWidth(), 3*16);
    fStatistic->Resize(74, 20);
    fStatistic->Connect("Selected(Int_t)", "Alieve::ITSSDSubEditor", this, "DoStatType(Int_t)");
    f->AddFrame(fStatistic, new TGLayoutHints(kLHintsLeft, 1, 2, 1, 1));
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  Int_t lp = 2;
  fInfoLabel0 = new TGLabel(this);
  fInfoLabel0->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel0, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 lp, 0, 8, 0));

  fInfoLabel1 = new TGLabel(this);
  fInfoLabel1->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel1, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 lp, 0, 2, 0));

  fInfoLabel2 = new TGLabel(this);
  fInfoLabel2->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel2, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 lp, 0, 2, 0));
}

/**************************************************************************/

void ITSSDSubEditor::SetModel(ITSScaledModule* mod)
{
  fModule = mod;
  RGBAPaletteSubEditor::SetModel(fModule->GetPalette());
  
  fScale->SetIntNumber(ITSScaledModule::fgDigitScaleInfo->GetScale());
  fStatistic->Select(ITSScaledModule::fgDigitScaleInfo->GetStatType(), kFALSE);
  
  Int_t cnx, cnz, total;
  Float_t  maxoc;
  GetSubDetScaleData(cnx, cnz, total, maxoc);
  fInfoLabel0->SetText(Form("Cell size:  Nx=%d Nz=%d", cnx, cnz));
  fInfoLabel1->SetText(Form("Num cells:  %d", total));
  fInfoLabel2->SetText(Form("Max occupancy:  %5.3f%%", maxoc));

  SetPaletteFromDigitInfo();
}

/**************************************************************************/

void ITSSDSubEditor::GetSubDetScaleData(Int_t& cnx, Int_t& cnz, Int_t& total, Float_t& maxoc)
{
  Int_t scale = fScale->GetIntNumber() -1;
  Alieve::ITSDigitsInfo* di = fModule->GetDigitsInfo();
  switch(fModule->fDetID)
  {
    case 0: 
      cnx =   di->fSPDScaleX[scale], cnz = di->fSPDScaleZ[scale];
      total =   di->fSegSPD->Npx()*di->fSegSPD->Npz();
      maxoc = di->fSPDMaxOcc;;
      break;
    case 1:
      cnx =   di->fSDDScaleX[scale], cnz = di->fSDDScaleZ[scale];
      total = di->fSegSDD->Npx()*di->fSegSDD->Npz();
      maxoc = di->fSDDMaxOcc;;
      break;
    case 2:
      cnx =   di->fSSDScale[scale], cnz = 1;
      total = di->fSegSSD->Npx()*di->fSegSSD->Npz();
      maxoc = di->fSSDMaxOcc;;
      break;
  }
}

/**************************************************************************/

void ITSSDSubEditor::SetPaletteFromDigitInfo()
{  
  // apply values for color palette
  if(ITSScaledModule::fgDigitScaleInfo->fAutoUpdatePalette) 
  {   
    Int_t cnx, cnz, total;
    Float_t  maxoc;
    GetSubDetScaleData(cnx, cnz, total, maxoc);
    Alieve::ITSDigitsInfo* di = fModule->GetDigitsInfo();
    if(fStatistic->GetSelected() == DigitScaleInfo::ST_Occup) 
    {
      Int_t scale = fScale->GetIntNumber() -1;
      fMinMax->SetValues(0, TMath::Max(Int_t(cnx*cnz*maxoc),1), kFALSE);
      ITSModule::fgSPDPalette->SetMinMax(0,TMath::Max(Int_t(di->fSPDScaleZ[scale]*di->fSPDScaleX[scale]*di->fSPDMaxOcc), 1));
      ITSModule::fgSDDPalette->SetMinMax(0,TMath::Max(Int_t(di->fSDDScaleZ[scale]*di->fSDDScaleX[scale]*di->fSPDMaxOcc), 1));
      ITSModule::fgSSDPalette->SetMinMax(0,TMath::Max(Int_t(di->fSSDScale[scale]*di->fSPDMaxOcc), 1));
    }
    else 
    {
      ITSModule::fgSPDPalette->SetMinMax(di->fSPDMinVal,di->fSPDMaxVal);
      ITSModule::fgSDDPalette->SetMinMax(di->fSDDMinVal,di->fSDDMaxVal);
      ITSModule::fgSSDPalette->SetMinMax(di->fSSDMinVal,di->fSSDMaxVal);
    }
  }
}

/**************************************************************************/

void ITSSDSubEditor::DoScale()
{
  Int_t cnx, cnz, total; Float_t  maxoc;
  GetSubDetScaleData(cnx, cnz, total, maxoc);
  fInfoLabel0->SetText(Form("Cell size:  Nx=%d Nz=%d", cnx, cnz));

  SetPaletteFromDigitInfo();

  ITSScaledModule::fgDigitScaleInfo->ScaleChanged(fScale->GetIntNumber());
  Changed();
}

/**************************************************************************/

void ITSSDSubEditor::DoStatType(Int_t v)
{
  // update palette
  SetPaletteFromDigitInfo();

  ITSScaledModule::fgDigitScaleInfo->StatTypeChanged(v);
  
  Changed(); 
}

//______________________________________________________________________
// ITSScaledModuleEditor 
//

ClassImp(ITSScaledModuleEditor)

ITSScaledModuleEditor::ITSScaledModuleEditor(const TGWindow *p, Int_t width, Int_t height,
					       UInt_t options, Pixel_t back) :
    TGedFrame(p, width, height, options | kVerticalFrame, back),
    fM(0), 
    fHMTrans   (0),
    fSDPalette   (0)
{
  MakeTitle("ITSScaledModule");
  
  fHMTrans = new ZTransSubEditor(this);
  fHMTrans->Connect("UseTrans()",     "Alieve::ITSScaledModuleEditor", this, "Update()");
  fHMTrans->Connect("TransChanged()", "Alieve::ITSScaledModuleEditor", this, "Update()");
  AddFrame(fHMTrans, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));

  MakeTitle("Palette controls");

  fSDPalette = new ITSSDSubEditor(this);
  fSDPalette->Connect("Changed", "Alieve::ITSScaledModuleEditor", this, "Update()");
  AddFrame(fSDPalette, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));
}

/*************************************************************************/
ITSScaledModuleEditor::~ITSScaledModuleEditor()
{}

/**************************************************************************/

void ITSScaledModuleEditor::ActivateBaseClassEditors(TClass* cl)
{
  // exclude QuadSet editor
  fGedEditor->ExcludeClassEditor(QuadSet::Class());
  TGedFrame::ActivateBaseClassEditors(cl);
}

/**************************************************************************/

void ITSScaledModuleEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<ITSScaledModule*>(obj);

  fHMTrans->SetDataFromTrans(&fM->RefHMTrans());

  if (fM->GetValueIsColor() || fM->GetPalette() == 0) {
    fSDPalette->UnmapWindow();
  } else {
    fSDPalette->SetModel(fM);
    fSDPalette->MapWindow();
  } 
}
