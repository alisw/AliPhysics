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
//  used for the calibration monitor                                         //
//  Example usage:                                                           //
/*
  aliroot
  AliTPCCalibViewerGUI::showGUI("allInOne22.root")
*/
// - Resize windows - (BUG to BE FIXED -> ROOT bug)                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTPCCalibViewerGUI.h"

#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPad.h>

#include <TObjArray.h>
#include <TObjString.h>
#include <TVector.h>
#include <string.h>
#include <TH1.h>


ClassImp(AliTPCCalibViewerGUI)

AliTPCCalibViewerGUI::AliTPCCalibViewerGUI(const TGWindow *p, UInt_t w, UInt_t h, char* fileName)
  : TGCompositeFrame(p, w, h),
    fViewer(0),
    fContTopBottom(0),
    fContLCR(0),
    fContLeft(0),
    fContRight(0),
    fContCenter(0),
    fContPlotOpt(0),
    fContDrawOpt(0),
    fContDrawOptSub1D2D(0),
    fContNormalized(0),
    fContCustom(0),
    fContCuts(0),
    fContSector(0),
    fContAddCuts(0),
    fContFit(0),
    fContAddFit(0),
    fContScaling(0),
    fContSetMax(0),
    fContSetMin(0),
    fListVariables(0),
    fBtnDraw(0),
    fBtnFit(0),
    fBtnAddFitFunction(0),
    fBtnGetMinMax(0),
    fCanvMain(0),
    fRadioRaw(0),
    fRadioNormalized(0),
    fRadioPredefined(0),
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
    fComboCustom(0),
    fNmbSector(0),
    fLblSector(0),
    fChkAddCuts(0),
    fComboAddCuts(0), 
    fComboCustomFit(0),
    fChkSetMax(0),
    fChkSetMin(0),
    fChkGetMinMaxAuto(0),
    fTxtSetMax(0),
    fTxtSetMin(0)
//
// AliTPCCalibViewerGUI constructor; fileName specifies the ROOT tree used for drawing
//
{
   SetCleanup(kDeepCleanup);
   
   // ************************* content of this MainFrame *************************
   // top level container with horizontal layout
   fContTopBottom = new TGCompositeFrame(this, w, h, kVerticalFrame | kFixedWidth | kFixedHeight);
   AddFrame(fContTopBottom, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   
   fContLCR = new TGCompositeFrame(fContTopBottom, w, h, kHorizontalFrame | kFixedWidth | kFixedHeight);
   fContTopBottom->AddFrame(fContLCR, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   

   // ************************* content of fContLCR *************************
   // left container
   fContLeft = new TGCompositeFrame(fContLCR, 200, 200, kVerticalFrame | kFixedWidth | kFitHeight);
   fContLCR->AddFrame(fContLeft, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY, 5, 3, 3, 3));
   
   // left vertical splitter
   TGVSplitter *splitLeft = new TGVSplitter(fContLCR);
   splitLeft->SetFrame(fContLeft, kTRUE);
   fContLCR->AddFrame(splitLeft, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 0, 0, 0, 0));

   // right container
   fContRight = new TGCompositeFrame(fContLCR, 150, 200, kVerticalFrame | kFixedWidth | kFitHeight);
   fContLCR->AddFrame(fContRight, new TGLayoutHints(kLHintsTop | kLHintsRight | kLHintsExpandY, 3, 5, 3, 3));
   
   // center container
   fContCenter = new TGCompositeFrame(fContLCR, 200, 200, kVerticalFrame | kFixedWidth | kFitHeight);
   fContLCR->AddFrame(fContCenter, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

   // right vertical splitter
   TGVSplitter *splitRight = new TGVSplitter(fContLCR);
   splitRight->SetFrame(fContRight, kFALSE);
   fContLCR->AddFrame(splitRight, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 0, 0, 0, 0));
   
   // ************************* content of fContLeft *************************
   // draw button
   fBtnDraw = new TGTextButton(fContLeft, "&Draw");
   fContLeft->AddFrame(fBtnDraw, new TGLayoutHints(kLHintsExpandX, 10, 10, 0, 0));
   //fBtnDraw->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoTest(=\"fBtnDraw clicked\")");
   fBtnDraw->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoDraw()");
   
   // draw options container
   fContDrawOpt = new TGGroupFrame(fContLeft, "Plot options", kVerticalFrame | kFitWidth | kFitHeight);
   fContLeft->AddFrame(fContDrawOpt, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
   fContDrawOptSub1D2D = new TGCompositeFrame(fContDrawOpt, 200, 20, kHorizontalFrame | kFitWidth | kFixedHeight);
   fContDrawOpt->AddFrame(fContDrawOptSub1D2D, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   
   
   // predefined radio button
   fRadioPredefined = new TGRadioButton(fContLeft, "Predefined: ", 13);
   fContLeft->AddFrame(fRadioPredefined, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fRadioPredefined->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");
   
   // list of variables
   fListVariables = new TGListBox(fContLeft);
   fContLeft->AddFrame(fListVariables, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY, 10, 0, 0, 0));
   fListVariables->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");

   // plot options container
   //fContPlotOpt = new TGCompositeFrame(fContLeft, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContPlotOpt = new TGGroupFrame(fContLeft, "Normalization options", kVerticalFrame | kFitWidth | kFitHeight);
   fContLeft->AddFrame(fContPlotOpt, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 0, 0, 0));
   
   // custom radio button
   fRadioCustom = new TGRadioButton(fContLeft, "Custom: ", 12);
   fContLeft->AddFrame(fRadioCustom, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fRadioCustom->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   // custom options container
   fContCustom = new TGCompositeFrame(fContTopBottom, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContTopBottom->AddFrame(fContCustom, new TGLayoutHints(kLHintsExpandX, 10, 0, 0, 0));

   // ************************* content of fContRight *************************
   // cut options container
   //fContCuts = new TGCompositeFrame(fContRight, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContCuts = new TGGroupFrame(fContRight, "Cuts", kVerticalFrame | kFitWidth | kFitHeight);
   fContRight->AddFrame(fContCuts, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   // Fit options container
   fContFit = new TGGroupFrame(fContRight, "Custom Fit", kVerticalFrame | kFitWidth | kFitHeight);
   fContRight->AddFrame(fContFit, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   
   // Scaling options container
   fContScaling = new TGGroupFrame(fContRight, "Scaling", kVerticalFrame | kFitWidth | kFitHeight);
   fContRight->AddFrame(fContScaling, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   
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

   // ************************* content of fContDrawOpt *************************
   // 1D radio button
   fRadio1D = new TGRadioButton(fContDrawOptSub1D2D, "1D", 30);
//   fContDrawOpt->AddFrame(fRadio1D, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
   fContDrawOptSub1D2D->AddFrame(fRadio1D, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
   fRadio1D->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");
   
   // 2D radio button
   fRadio2D = new TGRadioButton(fContDrawOptSub1D2D, "2D", 31);
   fContDrawOptSub1D2D->AddFrame(fRadio2D, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
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
   fContSector = new TGCompositeFrame(fContCuts, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
   fContCuts->AddFrame(fContSector, new TGLayoutHints(kLHintsExpandX, 5, 0, 0, 0));
   
   // additional cuts check button
   fChkAddCuts = new TGCheckButton(fContCuts, "additional cuts");
   fContCuts->AddFrame(fChkAddCuts, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
   fChkAddCuts->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");

   // additional cuts container
   fContAddCuts = new TGCompositeFrame(fContCuts, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContCuts->AddFrame(fContAddCuts, new TGLayoutHints(kLHintsExpandX, -5, -5, 0, 0));
   
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
   fComboCustom = new TGComboBox(fContCustom);
   fComboCustom->Resize(0, fBtnDraw->GetDefaultHeight());
   fComboCustom->EnableTextInput(kTRUE);
   fContCustom->AddFrame(fComboCustom, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   fComboCustom->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtons(=42)");
   fComboCustom->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
   
   // ************************* content of fContSector *************************
   // sector number entry
   fNmbSector = new TGNumberEntry(fContSector, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 71);
   fContSector->AddFrame(fNmbSector, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   fNmbSector->Connect("ValueSet(Long_t)", "AliTPCCalibViewerGUI", this, "ChangeSector()");
   
   // sector number label
   fLblSector = new TGLabel(fContSector, "IROC, A");
   fContSector->AddFrame(fLblSector, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
   
   // ************************* content of fContAddCuts *************************
   // combo text field for additional cuts
   fComboAddCuts = new TGComboBox(fContAddCuts);
   fComboAddCuts->Resize(0, fBtnDraw->GetDefaultHeight());
   fComboAddCuts->EnableTextInput(kTRUE);
   fContAddCuts->AddFrame(fComboAddCuts, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   fComboAddCuts->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
   fComboAddCuts->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
   
   // ************************* content of fContFit *************************
   // container for additional fits
   fContAddFit = new TGCompositeFrame(fContFit, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContFit->AddFrame(fContAddFit, new TGLayoutHints(kLHintsExpandX, -5, -5, 0, 0));
   
   // ************************* content of fContAddFit *************************
   // text field for custom fit
   fComboCustomFit = new TGComboBox(fContAddFit);
   fComboCustomFit->Resize(0, fBtnDraw->GetDefaultHeight());
   fComboCustomFit->EnableTextInput(kTRUE);
   fContAddFit->AddFrame(fComboCustomFit, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   fComboCustomFit->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "DoFit()");
   fComboCustomFit->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoFit()");
   
   // fit button
   fBtnFit = new TGTextButton(fContAddFit, "&Fit");
   fContAddFit->AddFrame(fBtnFit, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fBtnFit->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoFit()");

   // add fit function button
   //fBtnAddFitFunction = new TGTextButton(fContAddFit, "&Add fit function to normalization");
   //fContAddFit->AddFrame(fBtnAddFitFunction, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   //fBtnAddFitFunction->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "AddFitFunction()");

   // ************************* content of fContScaling *************************
   // SetMaximum container
   fContSetMax = new TGCompositeFrame(fContScaling, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContScaling->AddFrame(fContSetMax, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   // SetMinimum container
   fContSetMin = new TGCompositeFrame(fContScaling, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContScaling->AddFrame(fContSetMin, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   
   // get Min & Max from Plot - button
   fBtnGetMinMax = new TGTextButton(fContScaling, "&Get scale from plot");
   fContScaling->AddFrame(fBtnGetMinMax, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   fBtnGetMinMax->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "GetMinMax()");
   
   // GetMinMaxAuto - checkbox
   fChkGetMinMaxAuto = new TGCheckButton(fContScaling, "Get Min + Max auto.");
   fContScaling->AddFrame(fChkGetMinMaxAuto, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
   fChkGetMinMaxAuto->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");

   
   // ************************* content of fContSetMax *************************
   // SetMaximum - checkbox
   fChkSetMax = new TGCheckButton(fContSetMax, "Set fixed max.");
   fContSetMax->AddFrame(fChkSetMax, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
   fChkSetMax->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons()");
   
   // text field for maximum value
   fTxtSetMax = new TGTextEntry(fContSetMax, "", 41);
   fContSetMax->AddFrame(fTxtSetMax, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   fTxtSetMax->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   // ************************* content of fContSetMin *************************
   // SetMinimum - checkbox
   fChkSetMin = new TGCheckButton(fContSetMin, "Set fixed min.");
   fContSetMin->AddFrame(fChkSetMin, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
   fChkSetMin->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
   
   // text field for minimum value
   fTxtSetMin = new TGTextEntry(fContSetMin, "", 40);
   fContSetMin->AddFrame(fTxtSetMin, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   fTxtSetMin->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtons()");

   
   // Display everything
   Initialize(fileName);
   SetWindowName("AliTPCCalibViewer GUI");
   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
}

AliTPCCalibViewerGUI::AliTPCCalibViewerGUI(const AliTPCCalibViewerGUI &c)
   : TGCompositeFrame(c.fParent, c.fWidth, c.fHeight),
    fViewer(0),
    fContTopBottom(0),
    fContLCR(0),
    fContLeft(0),
    fContRight(0),
    fContCenter(0),
    fContPlotOpt(0),
    fContDrawOpt(0),
    fContDrawOptSub1D2D(0),
    fContNormalized(0),
    fContCustom(0),
    fContCuts(0),
    fContSector(0),
    fContAddCuts(0),
    fContFit(0),
    fContAddFit(0),
    fContScaling(0),
    fContSetMax(0),
    fContSetMin(0),
    fListVariables(0),
    fBtnDraw(0),
    fBtnFit(0),
    fBtnAddFitFunction(0),
    fBtnGetMinMax(0),
    fCanvMain(0),
    fRadioRaw(0),
    fRadioNormalized(0),
    fRadioPredefined(0),
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
    fComboCustom(0),
    fNmbSector(0),
    fLblSector(0),
    fChkAddCuts(0),
    fComboAddCuts(0), 
    fComboCustomFit(0),
    fChkSetMax(0),
    fChkSetMin(0),
    fChkGetMinMaxAuto(0),
    fTxtSetMax(0),
    fTxtSetMin(0)
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
   if (fCanvMain && fCanvMain->GetCanvas()) {
      for (Int_t i = 0; i < fCanvMain->GetCanvas()->GetListOfPrimitives()->GetEntries(); i++) {
         if (strcmp(fCanvMain->GetCanvas()->GetListOfPrimitives()->At(i)->ClassName(), "TFrame") != 0)
            fCanvMain->GetCanvas()->GetListOfPrimitives()->At(i)->Delete();
      }
   }
   Cleanup();
   if (fViewer) fViewer->Delete();
}

/*
void AliTPCCalibViewerGUI::CloseWindow() {
   DeleteWindow();
}
*/

void AliTPCCalibViewerGUI::Initialize(char* fileName) {
   //
   // initializes the GUI with default settings and opens tree for drawing
   //
   
   // create AliTPCCalibViewer object, which will be used for generating all drawings
   if (fViewer) delete fViewer;
   fViewer = new AliTPCCalibViewer(fileName);

   // fill fListVariables
   TObjArray* arr = fViewer->GetListOfVariables();
   TIterator* iter = arr->MakeIterator();
   iter->Reset();
   TObjString* currentStr = 0;
   Int_t id = 0;
   while ((currentStr = (TObjString*)(iter->Next()))) {
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
   while ((currentStr = (TObjString*)(iter->Next()))) {
      fListNormalization->AddEntry(currentStr->GetString().Data(), id);
      id++;
   }
   delete iter;
   arr->Delete();
   delete arr;

   // set default button states
   fRadioPredefined->SetState(kButtonDown);
   fRadioRaw->SetState(kButtonDown);
   fRadioTPC->SetState(kButtonDown);
   fRadio1D->SetState(kButtonDown);
   fChkAuto->SetState(kButtonDown);
   fChkAddCuts->SetState(kButtonUp);
   fListVariables->Select(0);
   fListNormalization->Select(0);
   fComboMethod->Select(0);
   fChkGetMinMaxAuto->SetState(kButtonDown);
   fChkSetMin->SetState(kButtonUp);
   fChkSetMax->SetState(kButtonUp);

   //fCanvMain->GetCanvas()->ToggleEventStatus(); // klappt nicht
   //fCanvMain->GetCanvas()->GetCanvasImp()->ShowStatusBar(kTRUE); // klappt auch nicht
   fListVariables->IntegralHeight(kFALSE);         // naja
   fListNormalization->IntegralHeight(kFALSE);     // naja
   DoDraw();
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
         fRadioPredefined->SetState(kButtonDown);
         fRadioCustom->SetState(kButtonUp);
         //fComboMethod->UnmapWindow();
         //fListNormalization->UnmapWindow();
         break;
      case 11:             // fRadioNormalized
         fRadioRaw->SetState(kButtonUp);
         fRadioPredefined->SetState(kButtonDown);
         fRadioCustom->SetState(kButtonUp);
         break;
      case 12:             // fRadioCustom
         fRadioPredefined->SetState(kButtonUp);
         //fRadioNormalized->SetState(kButtonUp);
         break;
      case 13:             // fRadioPredefined
         fRadioCustom->SetState(kButtonUp);
         //fRadioNormalized->SetState(kButtonUp);
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
      //--------
      case 40:             // fTxtSetMin
         fChkSetMin->SetState(kButtonDown);
         break;
      case 41:             // fTxtSetMax
         fChkSetMax->SetState(kButtonDown);
         break;
      case 42:             // fComboCustom
         fRadioCustom->SetState(kButtonDown);
         fRadioPredefined->SetState(kButtonUp);
         break;
   }
   DoNewSelection();
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
   if (!fListVariables->GetSelectedEntry()) return;
   desiredData += ((TGTextLBEntry*)(fListVariables->GetSelectedEntry()))->GetTitle();
   desiredData += ".fElements";

   // specify normalization
   if (fRadioPredefined->GetState() == kButtonDown && fRadioNormalized->GetState() == kButtonDown) {
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
      if (!fListNormalization->GetSelectedEntry()) return;
      normalizationData += ((TGTextLBEntry*)(fListNormalization->GetSelectedEntry()))->GetTitle();
      
      if ( normalizationData.BeginsWith("Fit")) {
         // create fit formula, evaluate it an replace normalizationData-String
         // ********** create cut string **********
         TString cutStr("");
         if (fRadioTPC->GetState() == kButtonDown)
            cutStr += ""; // whole TPC is used for fitting
         if (fRadioSideA->GetState() == kButtonDown)
            cutStr += "(sector/18)%2==0"; // side A
         if (fRadioSideC->GetState() == kButtonDown)
            cutStr+= "(sector/18)%2==1"; // side C
         if (fRadioSector->GetState() == kButtonDown) {
            Int_t sector = (Int_t)(fNmbSector->GetNumber());
            cutStr += "sector==";
            cutStr += sector; 
         }
         if (fChkAddCuts->GetState() == kButtonDown && strcmp(fComboAddCuts->GetTextEntry()->GetText(), "") != 0){
            if (fRadioTPC->GetState() != kButtonDown) cutStr += " && ";
            cutStr += fComboAddCuts->GetTextEntry()->GetText();  
         }
         Double_t chi2 = 0;
         TVectorD fitParam(0);
         TMatrixD covMatrix(0,0);
         TString formulaStr("");
         if (normalizationData.CompareTo("FitLinLocal") == 0)
            formulaStr = "lx~ ++ ly~";
         if (normalizationData.CompareTo("FitLinGlobal") == 0) 
            formulaStr = "gx~ ++ gy~";
         normalizationData = *fViewer->Fit(desiredData.Data(), formulaStr.Data(), cutStr.Data(), chi2, fitParam, covMatrix);
      }

      desiredData += op;
      if (! (TString(((TGTextLBEntry*)(fListNormalization->GetSelectedEntry()))->GetTitle())).BeginsWith("Fit"))
         desiredData += ((TGTextLBEntry*)(fListVariables->GetSelectedEntry()))->GetTitle();
      desiredData += normalizationData;
   }
   else if (fRadioCustom->GetState() == kButtonDown) {
      desiredData = fComboCustom->GetTextEntry()->GetText();
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
      cutsStr += fComboAddCuts->GetTextEntry()->GetText();
   
   // draw finally
   for (Int_t i = 0; i < fCanvMain->GetCanvas()->GetListOfPrimitives()->GetEntries(); i++) {
      if (strcmp(fCanvMain->GetCanvas()->GetListOfPrimitives()->At(i)->ClassName(), "TFrame") != 0)
         fCanvMain->GetCanvas()->GetListOfPrimitives()->At(i)->Delete();
   }
   //fCanvMain->GetCanvas()->Clear();
   fCanvMain->GetCanvas()->cd();
   Int_t entries = -1;
   if (fRadio1D->GetState() == kButtonDown)
      entries = fViewer->EasyDraw1D(desiredData.Data(), sectorStr.Data(), cutsStr.Data());
   else if (fRadio2D->GetState() == kButtonDown)
      entries = fViewer->EasyDraw(desiredData.Data(), sectorStr.Data(), cutsStr.Data());
   if (entries == -1) return;
   
   TList* listOfPrimitives = fCanvMain->GetCanvas()->GetListOfPrimitives();
   TObject* ptr = 0;
   for (Int_t i = 0; i < listOfPrimitives->GetEntries(); i++) {
      ptr = listOfPrimitives->At(i);
      if ( ptr->InheritsFrom("TH1") ) break;
   }
   if ( ptr != 0 && !ptr->InheritsFrom("TH1") ) return;      // if the loop did not find a TH1
   TH1 *hist = (TH1*)ptr; 
   TString minTxt(fTxtSetMin->GetText());
   TString maxTxt(fTxtSetMax->GetText());
   if (fChkSetMax->GetState() == kButtonDown && (maxTxt.IsDigit() || maxTxt.IsFloat()) )
      hist->SetMaximum(maxTxt.Atof());
   if (fChkSetMin->GetState() == kButtonDown && (minTxt.IsDigit() || minTxt.IsFloat()) )
      hist->SetMinimum(minTxt.Atof());
      
   if (fChkGetMinMaxAuto->GetState() == kButtonDown) {
      if (fChkSetMax->GetState() == kButtonUp)
         fTxtSetMax->SetText(Form("%f", hist->GetMaximum()));
      if (fChkSetMin->GetState() == kButtonUp)
         fTxtSetMin->SetText(Form("%f", hist->GetMinimum()));
   }
   
   fCanvMain->GetCanvas()->Update();
}


void AliTPCCalibViewerGUI::DoFit() {
   //
   // main method for fitting
   //
   
   Double_t chi2 = 0;
   TVectorD fitParam(0);
   TMatrixD covMatrix(0,0);
   TString drawStr("");
   TString cutStr("");
   TString formulaStr("");
   TString *returnStr = new TString("");

   
   // ******** create draw string *********
   if (fRadioCustom->GetState() == kButtonDown) {
   // take custom text as draw string
      drawStr = fComboCustom->GetTextEntry()->GetText();
      if (drawStr == "") return;
   }
   else if (fRadioPredefined->GetState() == kButtonDown) {
   // create drawStr out of selection
      drawStr += ((TGTextLBEntry*)(fListVariables->GetSelectedEntry()))->GetTitle();
      drawStr += ".fElements";
      if (fRadioNormalized->GetState() == kButtonDown) {
      // normalize data by selection
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
         drawStr += op;
         drawStr += ((TGTextLBEntry*)(fListVariables->GetSelectedEntry()))->GetTitle();
         drawStr += normalizationData;
      }
   }

   // ********** create cut string **********
   if (fRadioTPC->GetState() == kButtonDown)
      cutStr += ""; // whole TPC is used for fitting
   if (fRadioSideA->GetState() == kButtonDown)
      cutStr += "(sector/18)%2==0"; // side A
   if (fRadioSideC->GetState() == kButtonDown)
      cutStr+= "(sector/18)%2==1"; // side C
   if (fRadioSector->GetState() == kButtonDown) {
      Int_t sector = (Int_t)(fNmbSector->GetNumber());
      cutStr += "sector==";
      cutStr += sector; 
   }
   if (fChkAddCuts->GetState() == kButtonDown && strcmp(fComboAddCuts->GetTextEntry()->GetText(), "") != 0){
      if (fRadioTPC->GetState() != kButtonDown) cutStr += " && ";
      cutStr += fComboAddCuts->GetTextEntry()->GetText();  
   }
   
   // ********** get formula string **********
   formulaStr += fComboCustomFit->GetTextEntry()->GetText();

   // ********** call AliTPCCalibViewer's fit-function
   returnStr = fViewer->Fit(drawStr.Data(), formulaStr.Data(), cutStr.Data(), chi2, fitParam, covMatrix);
   
   std::cout << std::endl;
   std::cout << "Your fit formula reads as follows:" << std::endl;
   std::cout << returnStr->Data() << std::endl;
   std::cout << "chi2 = " << chi2 << std::endl;
}

void AliTPCCalibViewerGUI::GetMinMax() {
   //
   // Read current Min & Max from the plot and set it to fTxtSetMin & fTxtSetMax
   //
   TList* listOfPrimitives = fCanvMain->GetCanvas()->GetListOfPrimitives();
   TObject* ptr = 0;
   for (Int_t i = 0; i < listOfPrimitives->GetEntries(); i++) {
      ptr = listOfPrimitives->At(i);
      if ( ptr->InheritsFrom("TH1") ) break;
   }
   if ( ptr != 0 && !ptr->InheritsFrom("TH1") ) return;      // if the loop did not find a TH1
   TH1 *hist = (TH1*)ptr;
   Double_t histMax = hist->GetMaximum();
   Double_t histMin = hist->GetMinimum();
   fTxtSetMax->SetText(Form("%f",histMax));
   fTxtSetMin->SetText(Form("%f",histMin));
}

void AliTPCCalibViewerGUI::ChangeSector(){
   // 
   // function that is called, when the number of the sector is changed
   // to change the sector label
   // 
   Int_t sector = (Int_t)(fNmbSector->GetNumber());
   char* secLabel = "";
   if (sector >= 0 && sector <= 17) // IROC, Side A
      secLabel = "IROC, A";
   if (sector >= 18 && sector <= 35) // IROC, Side C
      secLabel = "IROC, C";
   if (sector >= 36 && sector <= 53) // OROC, Side A
      secLabel = "OROC, A";
   if (sector >= 54 && sector <= 71) // OROC, Side C
      secLabel = "OROC, C";
   fLblSector->SetText(secLabel);
   DoNewSelection();
}

void AliTPCCalibViewerGUI::AddFitFunction(){ 
   //
   // adds the last fit function to the normalization list
   // 
   std::cout << "Not yet implemented." << std::endl;
}
   

void AliTPCCalibViewerGUI::ShowGUI(const char* fileName) {
   //
   // initialize and show GUI for presentation
   // 
   TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
   frmMain->SetWindowName("AliTPCCalibViewer GUI");
   frmMain->SetCleanup(kDeepCleanup);
   
   TGTab* tabMain = new TGTab(frmMain, 1000, 600);
   frmMain->AddFrame(tabMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

   TGCompositeFrame* tabCont1 = tabMain->AddTab("Viewer 1");
   TGCompositeFrame* tabCont2 = tabMain->AddTab("Viewer 2");
   
   AliTPCCalibViewerGUI* calibViewer1 = new AliTPCCalibViewerGUI(tabCont1, 1000, 600, (char*)fileName);
   tabCont1->AddFrame(calibViewer1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

   AliTPCCalibViewerGUI* calibViewer2 = new AliTPCCalibViewerGUI(tabCont2, 1000, 600, (char*)fileName);
   tabCont2->AddFrame(calibViewer2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   
   frmMain->MapSubwindows();
   frmMain->Resize();
   frmMain->MapWindow();

}

