
/***********************************************************************
  This editor appears in the Reve window when v0 are visualize.
It allows to select the v0 as a function of some useful parameters.

Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/




#include "CascadeEditors.h"
#include <Reve/Cascade.h>

#include <Reve/RGValuators.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>


using namespace Reve;

//______________________________________________________________________
// CascadeListEditor
//

ClassImp(CascadeListEditor)

CascadeListEditor::CascadeListEditor(const TGWindow *p,
                                 Int_t width, Int_t height,
                                 UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fMList(0),
  fRnrV0Daughters(0),
  fRnrV0path(0),
  fRnrVtx(0),
  fRnrBach(0),
  fRnrCasPath(0),
  fMainTabA(0),
  fMainTabB(0)
{
  MakeTitle("CascadeList");
 
  //TGHorizontalFrame* frame = new TGHorizontalFrame(this);
 
  // --- Rendering control

  fRnrVtx = new TGCheckButton(this, "Render v0 and cascade vertices");
  AddFrame(fRnrVtx, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrVtx->Connect("Toggled(Bool_t)",
		     "Reve::CascadeListEditor", this, "DoRnrVtx()");

  fRnrV0path = new TGCheckButton(this, "Render v0 path");
  AddFrame(fRnrV0path, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrV0path->Connect("Toggled(Bool_t)",
		     "Reve::CascadeListEditor", this, "DoRnrV0path()");  

  fRnrCasPath = new TGCheckButton(this, "Render cascade path");
  AddFrame(fRnrCasPath, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrCasPath->Connect("Toggled(Bool_t)",
		       "Reve::CascadeListEditor", this, "DoRnrCasPath()");  

  fRnrV0Daughters = new TGCheckButton(this, "Render v0 daughter tracks");
  AddFrame(fRnrV0Daughters, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrV0Daughters->Connect("Toggled(Bool_t)",
			"Reve::CascadeListEditor", this, "DoRnrV0Daughters()");

  fRnrBach = new TGCheckButton(this, "Render bachelor track");
  AddFrame(fRnrBach, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrBach->Connect("Toggled(Bool_t)",
			"Reve::CascadeListEditor", this, "DoRnrBach()");
    
  for (Int_t i=0; i<fgkNRange; i++) fRange[i]=0;
  for (Int_t i=0; i<fgkNCanvas; i++) fCanvasA[i]=0;
  for (Int_t i=0; i<fgkNCanvas; i++) fCanvasB[i]=0;

  AddSelectTab();
  AddSeeTab();

  TGTextButton* resetCutsButton = new TGTextButton(this, "Reset all cuts", 40);
  resetCutsButton->Connect("Clicked()", "Reve::CascadeListEditor", this, "ResetCuts()");
  AddFrame(resetCutsButton, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
}

CascadeListEditor::~CascadeListEditor()
{}


//_________________________________________________________________________
void CascadeListEditor::SetModel(TObject* obj)
{
  fMList = dynamic_cast<CascadeList*>(obj);

  for (Int_t i=0; i<fgkNRange; i++)
    if (fRange[i])
      fRange[i]->SetValues( fMList->GetMin(i), fMList->GetMax(i) );

  fRnrV0Daughters->SetState(fMList->GetRnrV0Daughters() ? kButtonDown : kButtonUp);
  fRnrV0path->SetState(fMList->GetRnrV0path() ? kButtonDown : kButtonUp);
  fRnrVtx->SetState(fMList->GetRnrCasVtx() ? kButtonDown : kButtonUp);
  fRnrBach->SetState(fMList->GetRnrBachelor() ? kButtonDown : kButtonUp);
  fRnrCasPath->SetState(fMList->GetRnrCasPath() ? kButtonDown : kButtonUp);

  FillCanvas();
}


//_________________________________________________________________________
TGCompositeFrame* CascadeListEditor::AddTab(TGTab* tab, Int_t i, Int_t can,
				       char *name) {

  TGCompositeFrame* frameTab = tab->AddTab(name);
  frameTab->SetLayoutManager(new TGVerticalLayout( frameTab ));

  TRootEmbeddedCanvas** embeddedCanvas = 0;

  if (can==1) embeddedCanvas = &fCanvasA[i];
  else if (can==2) embeddedCanvas = &fCanvasB[i];

  *embeddedCanvas = new TRootEmbeddedCanvas("EmbeddedCanvas", frameTab, 200, 200);
  TCanvas *c1 = (*embeddedCanvas)->GetCanvas();
  c1->SetBorderMode(0);
  c1->SetTicks(1,0);
  c1->SetGrid();
  c1->SetBorderMode(0);
 
  frameTab->AddFrame(*embeddedCanvas, new TGLayoutHints(kLHintsTop|kLHintsExpandX,
						       0, 0, 0, 0));
  return frameTab;
}


//_________________________________________________________________________
TGCompositeFrame** CascadeListEditor::CreateTab(TGTab **pMainTab, TGTab **ptab, Int_t can) {

  //------
  // tab widget
  pMainTab[0] = new TGTab(this,0,0);
  pMainTab[0]->Connect("Selected(Int_t)", "Reve::CascadeListEditor", this, "UpdateSelectedTab()");
  this->AddFrame(pMainTab[0], new TGLayoutHints( kLHintsTop | kLHintsExpandX,2,2,2,2));

  //------
  // container of "Tab1"
  TGCompositeFrame *frameTab1 = pMainTab[0]->AddTab("ident.");
  frameTab1->SetLayoutManager(new TGVerticalLayout(frameTab1));
  
  // tab widget
  ptab[0] = new TGTab(frameTab1,2,2);
  ptab[0]->Resize(ptab[0]->GetDefaultSize());
  // The following is for updating the canvas of a tab if this one is selected
  // (it updates every canvas)
  ptab[0]->Connect("Selected(Int_t)", "Reve::CascadeListEditor", this, "UpdateSelectedTab()");
  frameTab1->AddFrame(ptab[0], new TGLayoutHints(kLHintsLeft| kLHintsExpandX,0,0,0,0));
  
  //------
  // container of "Tab2"
  TGCompositeFrame *frameTab2 = pMainTab[0]->AddTab("cascade");
  frameTab2->SetLayoutManager(new TGVerticalLayout(frameTab2));
  
  // tab widget
  ptab[1] = new TGTab(frameTab2,440,299);
  ptab[1]->Resize(ptab[1]->GetDefaultSize());
  ptab[1]->Connect("Selected(Int_t)", "Reve::CascadeListEditor", this, "UpdateSelectedTab()");
  frameTab2->AddFrame(ptab[1], new TGLayoutHints(kLHintsLeft| kLHintsExpandX ,0,0,0,0));
  
  //------
  // container of "Tab3"
  TGCompositeFrame *frameTab3 = pMainTab[0]->AddTab("V0daugh/bach.");
  frameTab3->SetLayoutManager(new TGVerticalLayout(frameTab3));
  
  // tab widget
  ptab[2] = new TGTab(frameTab3,440,299);
  ptab[2]->Resize(ptab[2]->GetDefaultSize());
  ptab[2]->Connect("Selected(Int_t)", "Reve::CascadeListEditor", this, "UpdateSelectedTab()");
  frameTab3->AddFrame(ptab[2], new TGLayoutHints(kLHintsLeft| kLHintsExpandX ,0,0,0,0));

  //------
  TGCompositeFrame **frameTab = new TGCompositeFrame*[fgkNCanvas];
  
  frameTab[0] = AddTab(ptab[0], 0, can, "Xi");
  frameTab[1] = AddTab(ptab[0], 1, can, "Omega");
  frameTab[2] = AddTab(ptab[0], 2, can, "Arm.Podo.");
  frameTab[3] = AddTab(ptab[0], 3, can, "Index");

  frameTab[4] = AddTab(ptab[1], 4, can, "cosPointing");
  frameTab[5] = AddTab(ptab[1], 5, can, "bachV0DCA");
  frameTab[6] = AddTab(ptab[1], 6, can, "radius");
  frameTab[7] = AddTab(ptab[1], 7, can, "Pt");
  frameTab[8] = AddTab(ptab[1], 8, can, "eta");

  frameTab[9] = AddTab(ptab[2], 9, can, "neg Pt");
  frameTab[10] = AddTab(ptab[2], 10, can, "neg eta");
  frameTab[11] = AddTab(ptab[2], 11, can, "pos Pt");
  frameTab[12] = AddTab(ptab[2], 12, can, "pos eta");
  frameTab[13] = AddTab(ptab[2], 13, can, "bach Pt");
  frameTab[14] = AddTab(ptab[2], 14, can, "bach eta");

  pMainTab[0]->SetTab(0);
  ptab[0]->SetTab(0);
  ptab[1]->SetTab(0);
  ptab[2]->SetTab(0);

  Pixel_t darkgrey;
  gClient->GetColorByName("grey50", darkgrey);
  ptab[0]->SetBackgroundColor(darkgrey);
  ptab[1]->SetBackgroundColor(darkgrey);
  ptab[2]->SetBackgroundColor(darkgrey);

  return frameTab;
}


//_________________________________________________________________________
void CascadeListEditor::AddValuator(TGCompositeFrame* frame, char *name,
			       Float_t min, Float_t max, Int_t pres,
			       char *func, Int_t iHist) {

  TGCompositeFrame* downFrame = new TGCompositeFrame(frame,
				    60, 60, kHorizontalFrame);

   // --- Selectors
  fRange[iHist] = new RGDoubleValuator(downFrame, name, 200, 0);
  fRange[iHist]->SetNELength(6);
  fRange[iHist]->Build();
  fRange[iHist]->GetSlider()->SetWidth(200);

  fRange[iHist]->SetLimits(min, max, TGNumberFormat::kNESRealTwo);
  if (pres==0)
    fRange[iHist]->SetLimits(min, max, TGNumberFormat::kNESInteger);
  else if (pres==3)
    fRange[iHist]->SetLimits(min, max, TGNumberFormat::kNESRealThree);
  else if (pres==4)
    fRange[iHist]->SetLimits(min, max, TGNumberFormat::kNESRealFour);
  else if (pres==5)
    fRange[iHist]->SetLimits(min, max, TGNumberFormat::kNESReal);

  fRange[iHist]->Connect("ValueSet()",
			 "Reve::CascadeListEditor", this, func);
  downFrame->AddFrame(fRange[iHist], new TGLayoutHints(kLHintsLeft,
						      0, 0, 0, 0));

  TGTextButton* adjustButton = new TGTextButton(downFrame, "Adjust Hist", 40);

  char ch[40];
  sprintf(ch,"AdjustHist(=%i)",iHist);
  adjustButton->Connect("Clicked()", "Reve::CascadeListEditor", this, ch);
  downFrame->AddFrame(adjustButton, new TGLayoutHints(kLHintsTop, 0, 0, 0, 0));

  frame->AddFrame(downFrame, new TGLayoutHints(kLHintsTop| kLHintsExpandY,
						0, 0, 0, 0));
}


//_________________________________________________________________________
void CascadeListEditor::AddSelectTab() {
 
  TGCompositeFrame** tab = CreateTab(&fMainTabA, fTabA, 1);

  AddValuator(tab[0],  "mass Xi",       0,  5, 3, "MassXiRange()",     0);
  AddValuator(tab[1],  "mass Omega",    0,  5, 3, "MassOmegaRange()",  1);

  AddValuator(tab[3],  "Index",                  0,   1e5, 0, "IndexRange()",       2);
  AddValuator(tab[4],  "cos pointing angle",     0.8,   1, 5, "CosPointingRange()", 3);
  AddValuator(tab[5],  "bach-V0 DCA",            0,     5, 3, "BachV0DCARange()",   4);
  AddValuator(tab[6],  "radius",                 0,   100, 2, "RadiusRange()",      5);
  AddValuator(tab[7],  "Pt",                     0,    10, 2, "PtRange()",          6);
  AddValuator(tab[8],  "Pseudo-rapidity",       -2,     2, 2, "PseudoRapRange()",   7);

  AddValuator(tab[9],  "neg Pt",                 0,    10, 2, "NegPtRange()",       8);
  AddValuator(tab[10], "neg pseudo-rapidity",   -2,     2, 2, "NegEtaRange()",      9);
  AddValuator(tab[11], "pos Pt",                 0,    10, 2, "PosPtRange()",      10);
  AddValuator(tab[12], "pos pseudo-rapidity",   -2,     2, 2, "PosEtaRange()",     11);
  AddValuator(tab[13], "bach. Pt",               0,    10, 2, "BachPtRange()",     12);
  AddValuator(tab[14], "bach. pseudo-rapidity", -2,     2, 2, "BachEtaRange()",    13);

  delete[] tab;
}


//_________________________________________________________________________
void CascadeListEditor::AddSeeTab() {

  TGCompositeFrame** tab = CreateTab(&fMainTabB, fTabB, 2);
  delete[] tab;
}


//_________________________________________________________________________
void CascadeListEditor::AdjustHist(Int_t iHist) {

  if (! fMList) return;
  fMList->AdjustHist(iHist);

  FillCanvas();
}

//_________________________________________________________________________
void CascadeListEditor::ResetCuts() {

  if (! fMList) return;

  Float_t min,max;

  for (Int_t i=0; i<fgkNRange;i++) {
    
    if (i==2) continue;
    min = fRange[i]->GetLimitMin();
    max = fRange[i]->GetLimitMax();
    fMList->SetMin(i, min);
    fMList->SetMax(i, max);
    fRange[i]->SetValues(min, max);
    fMList->AdjustHist(i);
  }

  // for the Index we scan its actual range
  Int_t imin, imax;
  fMList->GetCasIndexRange(imin, imax);
  if (imin<imax) {
    Int_t minH = imin-(imax-imin)/20;
    Int_t maxH = imax+(imax-imin)/20;
    fMList->SetMin(2, minH);
    fMList->SetMax(2, maxH);
    fRange[2]->SetLimits(minH, maxH, TGNumberFormat::kNESInteger);
    fRange[2]->SetValues(minH, maxH);
    fMList->AdjustHist(2);

  }
  FillCanvas();
  Update();
}


//_________________________________________________________________________
void CascadeListEditor::FillCanvas() {

  fMList->FilterAll();

  TCanvas *c1, *c1b;
  TH1F *hist = 0;
  TH2F *hist2D = 0;
  Bool_t is2D;

  Int_t canvasMap[fgkNCanvas]={0,1,1000,2,3,4,5,6,7,8,9,10,11,12,13};

  for (Int_t i=0; i<fgkNCanvas; i++)
    if (fCanvasA[i]) {

      is2D = canvasMap[i]>999;

      if (is2D) hist2D = fMList->GetHist2D(canvasMap[i]-1000);
      else hist = fMList->GetHist(canvasMap[i]);

      c1 = fCanvasA[i]->GetCanvas();
      c1->cd();
      if (is2D) hist2D->Draw("colz"); else hist->Draw();
      c1->Modified();
      c1->Update();
      
      c1b = fCanvasB[i]->GetCanvas();
      c1b->cd();
      if (is2D) hist2D->Draw("colz"); else hist->Draw();
      c1b->Modified();
      c1b->Update();
  }
  UpdateSelectedTab();
}
//_________________________________________________________________________
void CascadeListEditor::UpdateSelectedTab() {

  Pixel_t yellow;
  Pixel_t grey;
  gClient->GetColorByName("yellow", yellow);
  gClient->GetColorByName("grey", grey);

  TGTabElement *tabElem; 
  for (Int_t i=0; i<fMainTabA->GetNumberOfTabs(); i++) {

    tabElem = fMainTabA->GetTabTab(i);
    tabElem->ChangeBackground(grey);

    for (Int_t j=0; j<fTabA[i]->GetNumberOfTabs();j++) {
      tabElem = fTabA[i]->GetTabTab(j);
      tabElem->ChangeBackground(grey);
    }
  }

  Int_t currentTab = fMainTabA->GetCurrent();
  Int_t currentSubTab = fTabA[currentTab]->GetCurrent();
  tabElem = fMainTabA->GetTabTab(currentTab);
  tabElem->ChangeBackground(yellow);
  tabElem = fTabA[currentTab]->GetTabTab(currentSubTab);
  tabElem->ChangeBackground(yellow);

  TCanvas *c1;
  Int_t iCan = 0;
  Int_t i=0;

  while (currentTab>i) {
    iCan += fTabA[i]->GetNumberOfTabs();
    i++;
  }
  iCan += currentSubTab;
  c1 = fCanvasA[iCan]->GetCanvas();
  c1->GetCanvas()->Modified();
  c1->GetCanvas()->Update();

  //---

  for (Int_t i=0; i<fMainTabB->GetNumberOfTabs(); i++) {

    tabElem = fMainTabB->GetTabTab(i);
    tabElem->ChangeBackground(grey);

    for (Int_t j=0; j<fTabB[i]->GetNumberOfTabs();j++) {
      tabElem = fTabB[i]->GetTabTab(j);
      tabElem->ChangeBackground(grey);
    }
  }

  currentTab = fMainTabB->GetCurrent();
  currentSubTab = fTabB[currentTab]->GetCurrent();
  tabElem = fMainTabB->GetTabTab(currentTab);
  tabElem->ChangeBackground(yellow);
  tabElem = fTabB[currentTab]->GetTabTab(currentSubTab);
  tabElem->ChangeBackground(yellow);

  iCan = 0;
  i=0;

  while (currentTab>i) {
    iCan += fTabB[i]->GetNumberOfTabs();
    i++;
  }
  iCan += currentSubTab;
  c1 = fCanvasB[iCan]->GetCanvas();
  c1->GetCanvas()->Modified();
  c1->GetCanvas()->Update();
}


//_________________________________________________________________________
void CascadeListEditor::UpdateAll(Int_t iCanA) {

  TCanvas *c1 = fCanvasA[iCanA]->GetCanvas();
  c1->Modified();
  c1->Update();

  static Int_t iCan, i;
  iCan=0;
  i=0;
  while (fMainTabB->GetCurrent()>i) {
    iCan += fTabB[i]->GetNumberOfTabs();
    i++;
  }
  iCan += fTabB[i]->GetCurrent();
  c1 = fCanvasB[iCan]->GetCanvas();
  c1->GetCanvas()->Modified();
  c1->GetCanvas()->Update();

  Update();
}


//_________________________________________________________________________

void CascadeListEditor::DoRnrVtx()
{
  fMList->SetRnrCasVtx(fRnrVtx->IsOn());
  Update();
}

void CascadeListEditor::DoRnrV0path()
{
  fMList->SetRnrV0path(fRnrV0path->IsOn());
  Update();
}

void CascadeListEditor::DoRnrV0Daughters()
{
  fMList->SetRnrV0Daughters(fRnrV0Daughters->IsOn());
  Update();
}

void CascadeListEditor::DoRnrBach()
{
  fMList->SetRnrBachelor(fRnrBach->IsOn());
  Update();
}

void CascadeListEditor::DoRnrCasPath()
{
  fMList->SetRnrCasPath(fRnrCasPath->IsOn());
  Update();
}


//_________________________________________________________________________
void CascadeListEditor::MassXiRange() {

  fMList->XiMassFilter(fRange[0]->GetMin(), fRange[0]->GetMax());
  UpdateAll(0);
}

//_________________________________________________________________________
  void CascadeListEditor::MassOmegaRange() {
  fMList->OmegaMassFilter(fRange[1]->GetMin(), fRange[1]->GetMax());
  UpdateAll(1);
}

//_________________________________________________________________________
  void CascadeListEditor::IndexRange() {
  fMList->IndexFilter(fRange[2]->GetMin(), fRange[2]->GetMax());
  UpdateAll(3);
}
//_________________________________________________________________________
  void CascadeListEditor::CosPointingRange() {
  fMList->CosPointingFilter(fRange[3]->GetMin(), fRange[3]->GetMax());
  UpdateAll(4);
}

//_________________________________________________________________________
  void CascadeListEditor::BachV0DCARange() {
  fMList->BachV0DCAFilter(fRange[4]->GetMin(), fRange[4]->GetMax());
  UpdateAll(5);
}

//_________________________________________________________________________
  void CascadeListEditor::RadiusRange() {
  fMList->RadiusFilter(fRange[5]->GetMin(), fRange[5]->GetMax());
  UpdateAll(6);
}

//_________________________________________________________________________
  void CascadeListEditor::PtRange() {
  fMList->PtFilter(fRange[6]->GetMin(), fRange[6]->GetMax());
  UpdateAll(7);
}

//_________________________________________________________________________
  void CascadeListEditor::PseudoRapRange() {
  fMList->PseudoRapFilter(fRange[7]->GetMin(), fRange[7]->GetMax());
  UpdateAll(8);
}

//_________________________________________________________________________
  void CascadeListEditor::NegPtRange() {
  fMList->NegPtFilter(fRange[8]->GetMin(), fRange[8]->GetMax());
  UpdateAll(9);
}

//_________________________________________________________________________
  void CascadeListEditor::NegEtaRange() {
  fMList->NegEtaFilter(fRange[9]->GetMin(), fRange[9]->GetMax());
  UpdateAll(10);
}

//_________________________________________________________________________
  void CascadeListEditor::PosPtRange() {
  fMList->PosPtFilter(fRange[10]->GetMin(), fRange[10]->GetMax());
  UpdateAll(11);
}

//_________________________________________________________________________
  void CascadeListEditor::PosEtaRange() {
  fMList->PosEtaFilter(fRange[11]->GetMin(), fRange[11]->GetMax());
  UpdateAll(12);
}

//_________________________________________________________________________
  void CascadeListEditor::BachPtRange() {
  fMList->BachPtFilter(fRange[12]->GetMin(), fRange[12]->GetMax());
  UpdateAll(13);
}

//_________________________________________________________________________
  void CascadeListEditor::BachEtaRange() {
  fMList->BachEtaFilter(fRange[13]->GetMin(), fRange[13]->GetMax());
  UpdateAll(14);
}


