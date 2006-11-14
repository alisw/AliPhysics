/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


/***********************************************************************
  This editor appears in the Reve window when v0 are visualize.
It allows to select the v0 as a function of some useful parameters.

Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/


#include "V0Editors.h"
#include <Reve/V0.h>

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
// V0ListEditor
//

ClassImp(V0ListEditor)

V0ListEditor::V0ListEditor(const TGWindow *p,
                                 Int_t width, Int_t height,
                                 UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fMList(0),
  fRnrV0sDaugh(0),
  fRnrV0vtx(0),
  fRnrV0path(0),
  fMainTabA(0),
  fMainTabB(0)
{
  MakeTitle("V0List");
 
  //TGHorizontalFrame* frame = new TGHorizontalFrame(this);
 
  // --- Rendering control

  fRnrV0path = new TGCheckButton(this, "Render v0 path");
  AddFrame(fRnrV0path, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrV0path->Connect("Toggled(Bool_t)",
		     "Reve::V0ListEditor", this, "DoRnrV0path()");  

  fRnrV0vtx = new TGCheckButton(this, "Render v0 vertices");
  AddFrame(fRnrV0vtx, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrV0vtx->Connect("Toggled(Bool_t)",
		     "Reve::V0ListEditor", this, "DoRnrV0vtx()");

  fRnrV0sDaugh = new TGCheckButton(this, "Render v0 daughters");
  AddFrame(fRnrV0sDaugh, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrV0sDaugh->Connect("Toggled(Bool_t)",
			"Reve::V0ListEditor", this, "DoRnrDaughters()");
    
  for (Int_t i=0; i<fgkNRange; i++) fRange[i]=0;
  for (Int_t i=0; i<fgkNCanvas; i++) fCanvasA[i]=0;

  AddSelectTab();
  AddSeeTab();

  TGTextButton* resetCutsButton = new TGTextButton(this, "Reset index cut", 40);
  resetCutsButton->Connect("Clicked()", "Reve::V0ListEditor", this, "ResetCuts()");
  AddFrame(resetCutsButton, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
}

V0ListEditor::~V0ListEditor() {}


//_________________________________________________________________________
void V0ListEditor::SetModel(TObject* obj)
{
  fMList = dynamic_cast<V0List*>(obj);

  for (Int_t i=0; i<fgkNRange; i++)
    if (fRange[i]) {
      fRange[i]->SetValues( fMList->GetMin(i), fMList->GetMax(i) );
    }
  fRnrV0sDaugh->SetState(fMList->GetRnrDaughters() ? kButtonDown : kButtonUp);
  fRnrV0vtx->SetState(fMList->GetRnrV0vtx() ? kButtonDown : kButtonUp);
  fRnrV0path->SetState(fMList->GetRnrV0path() ? kButtonDown : kButtonUp);

  FillCanvas();
}


//_________________________________________________________________________
TGCompositeFrame* V0ListEditor::AddTab(TGTab* tab, Int_t i, Int_t can,
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
TGCompositeFrame** V0ListEditor::CreateTab(TGTab **pMainTab, TGTab **ptab, Int_t can) {

  //------
  // tab widget
  pMainTab[0] = new TGTab(this,0,0);
  pMainTab[0]->Connect("Selected(Int_t)", "Reve::V0ListEditor", this, "UpdateSelectedTab()");
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
  ptab[0]->Connect("Selected(Int_t)", "Reve::V0ListEditor", this, "UpdateSelectedTab()");
  frameTab1->AddFrame(ptab[0], new TGLayoutHints(kLHintsLeft| kLHintsExpandX,0,0,0,0));
  
  //------
  // container of "Tab2"
  TGCompositeFrame *frameTab2 = pMainTab[0]->AddTab("v0");
  frameTab2->SetLayoutManager(new TGVerticalLayout(frameTab2));
  
  // tab widget
  ptab[1] = new TGTab(frameTab2,440,299);
  ptab[1]->Resize(ptab[1]->GetDefaultSize());
  ptab[1]->Connect("Selected(Int_t)", "Reve::V0ListEditor", this, "UpdateSelectedTab()");
  frameTab2->AddFrame(ptab[1], new TGLayoutHints(kLHintsLeft| kLHintsExpandX ,0,0,0,0));
  
  //------
  // container of "Tab3"
  TGCompositeFrame *frameTab3 = pMainTab[0]->AddTab("daughters");
  frameTab3->SetLayoutManager(new TGVerticalLayout(frameTab3));
  
  // tab widget
  ptab[2] = new TGTab(frameTab3,440,299);
  ptab[2]->Resize(ptab[2]->GetDefaultSize());
  ptab[2]->Connect("Selected(Int_t)", "Reve::V0ListEditor", this, "UpdateSelectedTab()");
  frameTab3->AddFrame(ptab[2], new TGLayoutHints(kLHintsLeft| kLHintsExpandX ,0,0,0,0));

  //------
  TGCompositeFrame **frameTab = new TGCompositeFrame*[fgkNCanvas];
  
  frameTab[0] = AddTab(ptab[0], 0, can, "K0s");
  frameTab[1] = AddTab(ptab[0], 1, can, "Lambda");
  frameTab[2] = AddTab(ptab[0], 2, can, "Anti-Lambda");
  frameTab[3] = AddTab(ptab[0], 3, can, "Arm-Podo");
  frameTab[4] = AddTab(ptab[0], 4, can, "Index");
  frameTab[5] = AddTab(ptab[1], 5, can, "cosPointing");
  frameTab[6] = AddTab(ptab[1], 6, can, "daughterDCA");
  frameTab[7] = AddTab(ptab[1], 7, can, "radius");
  frameTab[8] = AddTab(ptab[1], 8, can, "Pt");
  frameTab[9] = AddTab(ptab[1], 9, can, "eta");
  frameTab[10] = AddTab(ptab[2], 10, can, "neg Pt");
  frameTab[11] = AddTab(ptab[2], 11, can, "neg eta");
  frameTab[12] = AddTab(ptab[2], 12, can, "pos Pt");
  frameTab[13] = AddTab(ptab[2], 13, can, "pos eta");

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
void V0ListEditor::AddValuator(TGCompositeFrame* frame, char *name,
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
			 "Reve::V0ListEditor", this, func);
  downFrame->AddFrame(fRange[iHist], new TGLayoutHints(kLHintsLeft,
						      0, 0, 0, 0));

  TGTextButton* adjustButton = new TGTextButton(downFrame, "Adjust Hist", 40);

  char ch[40];
  sprintf(ch,"AdjustHist(=%i)",iHist);
  adjustButton->Connect("Clicked()", "Reve::V0ListEditor", this, ch);
  downFrame->AddFrame(adjustButton, new TGLayoutHints(kLHintsTop, 0, 0, 0, 0));

  frame->AddFrame(downFrame, new TGLayoutHints(kLHintsTop| kLHintsExpandY,
						0, 0, 0, 0));
}


//_________________________________________________________________________
void V0ListEditor::AddSelectTab() {
 
  TGCompositeFrame** tab = CreateTab(&fMainTabA, fTabA, 1);

  AddValuator(tab[0],  "mass K0s",           0,  5, 3, "MassK0sRange()",     1);
  AddValuator(tab[1],  "mass Lambda",        0,  5, 3, "MassLamRange()",     2);
  AddValuator(tab[2],  "mass Anti-Lambda",   0,  5, 3, "MassAntiLamRange()", 3);
  AddValuator(tab[4],  "ESD v0 index",       0,1e5, 0, "ESDv0IndexRange()", 12);
  AddValuator(tab[5],  "cos pointing angle", 0.8,1, 5, "CosPointingRange()", 5);
  AddValuator(tab[6],  "daughter DCA",       0, 10, 2, "DaughterDcaRange()", 4);
  AddValuator(tab[7],  "radius",             0,100, 2, "RadiusRange()",      6);
  AddValuator(tab[8],  "Pt range",           0, 10, 2, "PtRange()",          0);
  AddValuator(tab[9],  "eta",               -2,  2, 2, "EtaRange()",         7);
  AddValuator(tab[10],  "neg Pt",            0, 10, 2, "NegPtRange()",       8);
  AddValuator(tab[11], "neg eta",           -2,  2, 2, "NegEtaRange()",      9);
  AddValuator(tab[12], "pos Pt",             0, 10, 2, "PosPtRange()",       10);
  AddValuator(tab[13], "pos eta",           -2,  2, 2, "PosEtaRange()",      11);

  delete[] tab;
}


//_________________________________________________________________________
void V0ListEditor::AddSeeTab() {

  TGCompositeFrame** tab = CreateTab(&fMainTabB, fTabB, 2);
  delete[] tab;
}


//_________________________________________________________________________
void V0ListEditor::AdjustHist(Int_t iHist) {

  if (! fMList) return;
  fMList->AdjustHist(iHist);

  FillCanvas();
}

//_________________________________________________________________________
void V0ListEditor::ResetCuts() {

  if (! fMList) return;

  Int_t imin, imax;
  fMList->GetV0IndexRange(imin, imax);
  if (imin<imax) {
    Int_t minH = imin-(imax-imin)/20;
    Int_t maxH = imax+(imax-imin)/20;
    fMList->SetMin(12, minH);
    fMList->SetMax(12, maxH);
    fRange[12]->SetLimits(minH, maxH, TGNumberFormat::kNESInteger);
    fRange[12]->SetValues(minH, maxH);
    AdjustHist(12);

  }
  FillCanvas();
}



//_________________________________________________________________________
void V0ListEditor::FillCanvas() {

  fMList->FilterAll();

  TCanvas *c1, *c1b;
  TH1F *hist = 0;
  TH2F *hist2D = 0;
  Bool_t is2D;

  Int_t canvasMap[fgkNCanvas]={1,2,3,1000,12,5,4,6,0, 7,8,9,10,11};

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
void V0ListEditor::UpdateSelectedTab() {

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
void V0ListEditor::UpdateAll(Int_t iCanA) {

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

void V0ListEditor::DoRnrV0vtx()
{
  fMList->SetRnrV0vtx(fRnrV0vtx->IsOn());
  Update();
}

void V0ListEditor::DoRnrV0path()
{
  fMList->SetRnrV0path(fRnrV0path->IsOn());
  Update();
}

void V0ListEditor::DoRnrDaughters()
{
  fMList->SetRnrDaughters(fRnrV0sDaugh->IsOn());
  Update();
}


//_________________________________________________________________________
void V0ListEditor::MassK0sRange() {

  fMList->K0sMFilter(fRange[1]->GetMin(), fRange[1]->GetMax());
  UpdateAll(0);
}

//_________________________________________________________________________
  void V0ListEditor::MassLamRange() {
  fMList->LamMFilter(fRange[2]->GetMin(), fRange[2]->GetMax());
  UpdateAll(1);
}

//_________________________________________________________________________
  void V0ListEditor::MassAntiLamRange() {
  fMList->ALamMFilter(fRange[3]->GetMin(), fRange[3]->GetMax());
  UpdateAll(2);
}

//_________________________________________________________________________
void V0ListEditor::ESDv0IndexRange() {
  fMList->IndexFilter(fRange[12]->GetMin(), fRange[12]->GetMax());
  UpdateAll(4);
}

//_________________________________________________________________________
  void V0ListEditor::CosPointingRange() {
  fMList->CosPointingFilter(fRange[5]->GetMin(), fRange[5]->GetMax());
  UpdateAll(5);
}

//_________________________________________________________________________
  void V0ListEditor::DaughterDcaRange() {
  fMList->DaughterDCAFilter(fRange[4]->GetMin(), fRange[4]->GetMax());
  UpdateAll(6);
}

//_________________________________________________________________________
  void V0ListEditor::RadiusRange() {
  fMList->RadiusFilter(fRange[6]->GetMin(), fRange[6]->GetMax());
  UpdateAll(7);
}

//_________________________________________________________________________
void V0ListEditor::PtRange()
{
  fMList->PtFilter(fRange[0]->GetMin(), fRange[0]->GetMax());
  UpdateAll(8);
}

//_________________________________________________________________________
  void V0ListEditor::EtaRange() {
  fMList->EtaFilter(fRange[7]->GetMin(), fRange[7]->GetMax());
  UpdateAll(9);
}

//_________________________________________________________________________
  void V0ListEditor::NegPtRange() {
  fMList->NegPtFilter(fRange[8]->GetMin(), fRange[8]->GetMax());
  UpdateAll(10);
}

//_________________________________________________________________________
  void V0ListEditor::NegEtaRange() {
  fMList->NegEtaFilter(fRange[9]->GetMin(), fRange[9]->GetMax());
  UpdateAll(11);
}

//_________________________________________________________________________
  void V0ListEditor::PosPtRange() {
  fMList->PosPtFilter(fRange[10]->GetMin(), fRange[10]->GetMax());
  UpdateAll(12);
}

//_________________________________________________________________________
  void V0ListEditor::PosEtaRange() {
  fMList->PosEtaFilter(fRange[11]->GetMin(), fRange[11]->GetMax());
  UpdateAll(12);
}

