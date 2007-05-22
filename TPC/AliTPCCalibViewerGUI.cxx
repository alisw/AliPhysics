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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  GUI for the AliTPCCalibViewer                                            //
//  used for the calibration monitor   
//  Example usage: 
/*
  aliroot

  AliTPCCalibViewerGUI v(gClient->GetRoot(), 1000, 600, "/u/sgaertne/calibration/localFit/AliTPCCalibViewer/allInOne6.root")

 - Resize windows - (BUG to BE FIXED)

*/                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTPCCalibViewerGUI.h"

#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPad.h>

#include <TObjArray.h>
#include <TObjString.h>
#include <TVector.h>

ClassImp(AliTPCCalibViewerGUI)

AliTPCCalibViewerGUI::AliTPCCalibViewerGUI(const TGWindow *p, UInt_t w, UInt_t h, char* fileName)
  : TGMainFrame(p, w, h),
    fViewer(0),
    fContAll(0),
    fContLeft(0),
    fContRight(0),
    fContCenter(0),
    fContPlotOpt(0),
    fContDrawOpt(0),
    fContNormalized(0),
    fContCustom(0),
    fContCuts(0),
    fContSector(0),
    fContAddCuts(0),
    fListVariables(0),
    fBtnDraw(0),
    fCanvMain(0),
    fRadioRaw(0),
    fRadioNormalized(0),
    fRadioCustom(0),
    fRadio1D(0),
    fRadio2D(0),
    fRadioTPC(0),
    fRadioSideA(0),
    fRadioSideC(0),
    fRadioSector(0),
    fChkAuto(0),
    fComboMethod(0),
    fListNormalization(0),
    fTxtCustom(0),
    fNmbSector(0),
    fChkAddCuts(0),
    fTxtAddCuts(0)
//
// AliTPCCalibViewerGUI constructor; fileName specifies the ROOT tree used for drawing
//
{
   SetCleanup(kDeepCleanup);
   
   // ************************* content of this MainFrame *************************
   // top level container with horizontal layout
   fContAll = new TGCompositeFrame(this, w, h, kHorizontalFrame | kFixedWidth | kFixedHeight);
   AddFrame(fContAll, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

   // ************************* content of fContAll *************************
   // left container
   fContLeft = new TGCompositeFrame(fContAll, 200, 200, kVerticalFrame | kFixedWidth | kFitHeight);
   fContAll->AddFrame(fContLeft, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY, 5, 3, 3, 3));
   
   // left vertical splitter
   TGVSplitter *splitLeft = new TGVSplitter(fContAll);
   splitLeft->SetFrame(fContLeft, kTRUE);
   fContAll->AddFrame(splitLeft, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 0, 0, 0, 0));

   // right container
   fContRight = new TGCompositeFrame(fContAll, 150, 200, kVerticalFrame | kFixedWidth | kFitHeight);
   fContAll->AddFrame(fContRight, new TGLayoutHints(kLHintsTop | kLHintsRight | kLHintsExpandY, 3, 5, 3, 3));
   
   // center container
   fContCenter = new TGCompositeFrame(fContAll, 200, 200, kVerticalFrame | kFixedWidth | kFitHeight);
   fContAll->AddFrame(fContCenter, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

   // right vertical splitter
   TGVSplitter *splitRight = new TGVSplitter(fContAll);
   splitRight->SetFrame(fContRight, kFALSE);
   fContAll->AddFrame(splitRight, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 0, 0, 0, 0));
   
   // ************************* content of fContLeft *************************
   // list of variables
   fListVariables = new TGListBox(fContLeft);
   fContLeft->AddFrame(fListVariables, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   fListVariables->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");

   // plot options container
   //fContPlotOpt = new TGCompositeFrame(fContLeft, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContPlotOpt = new TGGroupFrame(fContLeft, "Plot options", kVerticalFrame | kFitWidth | kFitHeight);
   fContLeft->AddFrame(fContPlotOpt, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

   // draw options container
   fContDrawOpt = new TGCompositeFrame(fContLeft, 200, 20, kHorizontalFrame | kFitWidth | kFixedHeight);
   fContLeft->AddFrame(fContDrawOpt, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
   
   // draw button
   fBtnDraw = new TGTextButton(fContLeft, "&Draw");
   fContLeft->AddFrame(fBtnDraw, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   //fBtnDraw->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoTest(=\"fBtnDraw clicked\")");
   fBtnDraw->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoDraw()");
   
   // ************************* content of fContRight *************************
   // cut options container
   //fContCuts = new TGCompositeFrame(fContRight, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContCuts = new TGGroupFrame(fContRight, "Cuts", kVerticalFrame | kFitWidth | kFitHeight);
   fContRight->AddFrame(fContCuts, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   // ************************* content of fContCenter *************************
   // main drawing canvas
   fCanvMain = new TRootEmbeddedCanvas("Main Canvas", fContCenter, 200, 200, kFitWidth | kFitHeight);
   fContCenter->AddFrame(fCanvMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   

   // ************************* content of fContPlotOpt *************************
   //TGButtonGroup *fBtngrpPlotOpt = new TGButtonGroup(fContPlotOpt, "Plot options", 
   // raw radio button
   fRadioRaw = new TGRadioButton(fContPlotOpt, "Raw", 10);
   fContPlotOpt->AddFrame(fRadioRaw, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fRadioRaw->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   // normalized radio button
   fRadioNormalized = new TGRadioButton(fContPlotOpt, "Normalized", 11);
   fContPlotOpt->AddFrame(fRadioNormalized, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fRadioNormalized->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   //fContPlotOpt->Show();

   // normalized options container
   fContNormalized = new TGCompositeFrame(fContPlotOpt, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContPlotOpt->AddFrame(fContNormalized, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 15, 0, 0, 0));

   // custom radio button
   fRadioCustom = new TGRadioButton(fContPlotOpt, "Custom", 12);
   fContPlotOpt->AddFrame(fRadioCustom, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fRadioCustom->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   // custom options container
   fContCustom = new TGCompositeFrame(fContPlotOpt, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContPlotOpt->AddFrame(fContCustom, new TGLayoutHints(kLHintsExpandX, 15, 0, 0, 0));

   // ************************* content of fContDrawOpt *************************
   // 1D radio button
   fRadio1D = new TGRadioButton(fContDrawOpt, "1D", 30);
   fContDrawOpt->AddFrame(fRadio1D, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
   fRadio1D->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");
   
   // 2D radio button
   fRadio2D = new TGRadioButton(fContDrawOpt, "2D", 31);
   fContDrawOpt->AddFrame(fRadio2D, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
   fRadio2D->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   // automatic redraw check button
   fChkAuto = new TGCheckButton(fContDrawOpt, "auto redraw");
   fContDrawOpt->AddFrame(fChkAuto, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));

   // ************************* content of fContCuts *************************
   // TPC radio button
   fRadioTPC = new TGRadioButton(fContCuts, "whole TPC", 20);
   fContCuts->AddFrame(fRadioTPC, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fRadioTPC->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   // side A radio button
   fRadioSideA = new TGRadioButton(fContCuts, "side A", 21);
   fContCuts->AddFrame(fRadioSideA, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fRadioSideA->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   // side C radio button
   fRadioSideC = new TGRadioButton(fContCuts, "side C", 22);
   fContCuts->AddFrame(fRadioSideC, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fRadioSideC->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   // sector radio button
   fRadioSector = new TGRadioButton(fContCuts, "sector", 23);
   fContCuts->AddFrame(fRadioSector, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fRadioSector->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   // sector options container
   fContSector = new TGCompositeFrame(fContCuts, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContCuts->AddFrame(fContSector, new TGLayoutHints(kLHintsExpandX, 15, 0, 0, 0));
   
   // additional cuts check button
   fChkAddCuts = new TGCheckButton(fContCuts, "additional cuts");
   fContCuts->AddFrame(fChkAddCuts, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
   fChkAddCuts->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");

   // additional cuts container
   fContAddCuts = new TGCompositeFrame(fContCuts, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContCuts->AddFrame(fContAddCuts, new TGLayoutHints(kLHintsExpandX, 15, 0, 0, 0));
   
   // ************************* content of fContNormalized *************************
   // method drop down combo box
   fComboMethod = new TGComboBox(fContNormalized);
   fComboMethod->Resize(0, fBtnDraw->GetDefaultHeight());
   fContNormalized->AddFrame(fComboMethod, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   fComboMethod->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");

   // list of normalization variables
   fListNormalization = new TGListBox(fContNormalized);
   fContNormalized->AddFrame(fListNormalization, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   fListNormalization->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");

   // ************************* content of fContCustom *************************
   // text field for custom draw command
   fTxtCustom = new TGTextEntry(fContCustom);
   fContCustom->AddFrame(fTxtCustom, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   fTxtCustom->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
   
   // ************************* content of fContSector *************************
   // sector number entry
   fNmbSector = new TGNumberEntry(fContSector, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 71);
   fContSector->AddFrame(fNmbSector, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   fNmbSector->Connect("ValueSet(Long_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
   
   // ************************* content of fContAddCuts *************************
   // text field for additional cuts
   fTxtAddCuts = new TGTextEntry(fContAddCuts);
   fContAddCuts->AddFrame(fTxtAddCuts, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   fTxtAddCuts->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");

   // Display everything
   Initialize(fileName);
   SetWindowName("AliTPCCalibViewer GUI");
   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
}

AliTPCCalibViewerGUI::AliTPCCalibViewerGUI(const AliTPCCalibViewerGUI &c)
   : TGMainFrame(c.fParent, c.fWidth, c.fHeight),
    fViewer(0),
    fContAll(0),
    fContLeft(0),
    fContRight(0),
    fContCenter(0),
    fContPlotOpt(0),
    fContDrawOpt(0),
    fContNormalized(0),
    fContCustom(0),
    fContCuts(0),
    fContSector(0),
    fContAddCuts(0),
    fListVariables(0),
    fBtnDraw(0),
    fCanvMain(0),
    fRadioRaw(0),
    fRadioNormalized(0),
    fRadioCustom(0),
    fRadio1D(0),
    fRadio2D(0),
    fRadioTPC(0),
    fRadioSideA(0),
    fRadioSideC(0),
    fRadioSector(0),
    fChkAuto(0),
    fComboMethod(0),
    fListNormalization(0),
    fTxtCustom(0),
    fNmbSector(0),
    fChkAddCuts(0),
    fTxtAddCuts(0)
{
  //
  // dummy AliTPCCalibViewerGUI copy constructor
  //
}

AliTPCCalibViewerGUI & AliTPCCalibViewerGUI::operator =(const AliTPCCalibViewerGUI & param) {
   //
   // dummy assignment operator
   //
   return (*this);
}

AliTPCCalibViewerGUI::~AliTPCCalibViewerGUI() {
   Cleanup();
}

void AliTPCCalibViewerGUI::CloseWindow() {
   DeleteWindow();
}

void AliTPCCalibViewerGUI::Initialize(char* fileName) {
   //
   // initializes the GUI with default settings and opens tree for drawing
   //
   
   // create AliTPCCalibViewer object, which will be used for generating all drawings
   fViewer = new AliTPCCalibViewer(fileName);

   // fill fListVariables
   TObjArray* arr = fViewer->GetListOfVariables();
   TIterator* iter = arr->MakeIterator();
   iter->Reset();
   TObjString* currentStr = 0;
   Int_t id = 0;
   while (currentStr = (TObjString*)(iter->Next())) {
      fListVariables->AddEntry(currentStr->GetString().Data(), id);
      id++;
   }
   delete iter;
   arr->Delete();
   delete arr;

   // fill fComboMethod
   fComboMethod->AddEntry("subtract", 0);
   fComboMethod->AddEntry("divide by", 1);

   // fill fListNorm
   arr = fViewer->GetListOfNormalizationVariables();
   iter = arr->MakeIterator();
   iter->Reset();
   currentStr = 0;
   id = 0;
   while (currentStr = (TObjString*)(iter->Next())) {
      fListNormalization->AddEntry(currentStr->GetString().Data(), id);
      id++;
   }
   delete iter;
   arr->Delete();
   delete arr;

   // set default button states
   fRadioRaw->SetState(kButtonDown);
   fRadioTPC->SetState(kButtonDown);
   //fRadioCustom->SetState(kButtonDisabled);
   fRadio1D->SetState(kButtonDown);
   fChkAuto->SetState(kButtonDown);
   fChkAddCuts->SetState(kButtonUp);
   fListVariables->Select(0);
   fListNormalization->Select(0);
   fComboMethod->Select(0);
}

void AliTPCCalibViewerGUI::HandleButtons(Int_t id) {
   //
   // handles mutual radio button exclusions
   //
   if (id == -1) {
      TGButton *btn = (TGButton *) gTQSender;
      id = btn->WidgetId();
   }
   
   switch (id) {
      case 10:             // fRadioRaw
         fRadioNormalized->SetState(kButtonUp);
         fRadioCustom->SetState(kButtonUp);
         //fComboMethod->UnmapWindow();
         //fListNormalization->UnmapWindow();
         break;
      case 11:             // fRadioNormalized
         fRadioRaw->SetState(kButtonUp);
         fRadioCustom->SetState(kButtonUp);
         break;
      case 12:             // fRadioCustom
         fRadioRaw->SetState(kButtonUp);
         fRadioNormalized->SetState(kButtonUp);
         break;
      //--------
      case 20:             // fRadioTPC
         fRadioSideA->SetState(kButtonUp);
         fRadioSideC->SetState(kButtonUp);
         fRadioSector->SetState(kButtonUp);
         break;
      case 21:             // fRadioSideA
         fRadioTPC->SetState(kButtonUp);
         fRadioSideC->SetState(kButtonUp);
         fRadioSector->SetState(kButtonUp);
         break;
      case 22:             // fRadioSideC
         fRadioTPC->SetState(kButtonUp);
         fRadioSideA->SetState(kButtonUp);
         fRadioSector->SetState(kButtonUp);
         break;
      case 23:             // fRadioSector
         fRadioTPC->SetState(kButtonUp);
         fRadioSideA->SetState(kButtonUp);
         fRadioSideC->SetState(kButtonUp);
         break;
      //--------
      case 30:             // fRadio1D
         fRadio2D->SetState(kButtonUp);
         break;
      case 31:             // fRadio2D
         fRadio1D->SetState(kButtonUp);
         break;
   }
   //fRadioCustom->SetState(kButtonDisabled);

   if (fChkAuto->GetState() == kButtonDown) DoDraw();
}

void AliTPCCalibViewerGUI::DoNewSelection() {
   //
   // decides whether to redraw if user makes another selection
   //
   
   if (fChkAuto->GetState() == kButtonDown) DoDraw();
}

void AliTPCCalibViewerGUI::DoDraw() {
   //
   // main method for drawing according to user selection
   //
   
   // specify data to plot
   TString desiredData("");
   desiredData += ((TGTextLBEntry*)(fListVariables->GetSelectedEntry()))->GetTitle();
   desiredData += ".fElements";

   // specify normalization
   if (fRadioNormalized->GetState() == kButtonDown) {
      TString op("");
      switch (fComboMethod->GetSelected()) {
         case 0:        // subtraction
            op += "-";
            break;
         case 1:        // division
            op += "/";
            break;
      }
      TString normalizationData("");
      normalizationData += ((TGTextLBEntry*)(fListNormalization->GetSelectedEntry()))->GetTitle();
      desiredData += op;
      desiredData += ((TGTextLBEntry*)(fListVariables->GetSelectedEntry()))->GetTitle();
      //desiredData += "_";
      desiredData += normalizationData;
   }
   else if (fRadioCustom->GetState() == kButtonDown) {
      desiredData = fTxtCustom->GetText();
      if (desiredData == "") return;
   }

   // specify cuts
   TString sectorStr("");
   if (fRadioTPC->GetState() == kButtonDown)
      sectorStr += "ALL";
   if (fRadioSideA->GetState() == kButtonDown)
      sectorStr += "A"; //cuts += "(sector/18)%2==0";
   if (fRadioSideC->GetState() == kButtonDown)
      sectorStr+= "C"; //cuts += "(sector/18)%2==1";
   if (fRadioSector->GetState() == kButtonDown) {
      Int_t sector = (Int_t)(fNmbSector->GetNumber());
      sectorStr += sector; //cuts += "sector==";
   }
   TString cutsStr("");
   if (fChkAddCuts->GetState() == kButtonDown)
      cutsStr += fTxtAddCuts->GetText();
   
   // draw finally
   fCanvMain->GetCanvas()->cd();
   //fViewer->Draw(desiredData.Data(), cuts.Data());
   if (fRadio1D->GetState() == kButtonDown)
      fViewer->EasyDraw1D(desiredData.Data(), sectorStr.Data(), cutsStr.Data());
   else if (fRadio2D->GetState() == kButtonDown)
      fViewer->EasyDraw(desiredData.Data(), sectorStr.Data(), cutsStr.Data());
   
   fCanvMain->GetCanvas()->Update();
}
