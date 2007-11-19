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
//  All functionalities of the AliTPCCalibViewer are here available
//  
//  Example usage:                                                           //
/*
  aliroot
  AliTPCCalibViewerGUI::ShowGUI("CalibTree.root")
  Begin_macro(source,gui)
  {   
      char* fileName = "CalibTreeEmpty.root";
      AliTPCCalibViewer::MakeTree(fileName, 0, "$ALICE_ROOT/TPC/Calib/MapCalibrationObjects.root");
      gROOT->SetStyle("Plain");
      // content of AliTPCCalibViewerGUI::ShowGUI(...)
      TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
      frmMain->SetWindowName("AliTPCCalibViewer GUI");
      frmMain->SetCleanup(kDeepCleanup);
      
      TGTab* tabMain = new TGTab(frmMain, 1000, 600);
      frmMain->AddFrame(tabMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   
      TGCompositeFrame* tabCont1 = tabMain->AddTab("Viewer 1");
      TGCompositeFrame* tabCont2 = tabMain->AddTab("Viewer 2");
   
      AliTPCCalibViewerGUI* calibViewer1 = new AliTPCCalibViewerGUI(tabCont1, 1000, 600, fileName);
      tabCont1->AddFrame(calibViewer1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   
      AliTPCCalibViewerGUI* calibViewer2 = new AliTPCCalibViewerGUI(tabCont2, 1000, 600, fileName);
      tabCont2->AddFrame(calibViewer2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
      
      frmMain->MapSubwindows();
      frmMain->Resize();
      frmMain->MapWindow();
      
      return frmMain;
     
  }
  End_macro
*/
//                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTPCCalibViewerGUI.h"
#include <iostream>

#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPad.h>

#include <TObjArray.h>
#include <TObjString.h>
#include <TVector.h>
#include <string.h>
#include <TH1.h>
#include "TStyle.h"
#include "AliTPCCalibViewer.h"

// #include "TGListBox.h"
// #include "TGNumberEntry"
// #include "TGSplitter"
// #include "TGTab"
// #include "TGLabel"
// #include "TGButtonGroup"
// #include "TGComboBox"
// #include "TRootEmbeddedCanvas"
// #include "TGButton"
// #include "TGRadioButton"
// #include "GTCheckButton"
// #include "TGTextEntry"




ClassImp(AliTPCCalibViewerGUI)

AliTPCCalibViewerGUI::AliTPCCalibViewerGUI(const TGWindow *p, UInt_t w, UInt_t h, char* fileName)
  : TGCompositeFrame(p, w, h),
    fViewer(0),
    fContTopBottom(0),
    fContLCR(0),
    fContLeft(0),
    ftabLeft(0),
    ftabLeft0(0),
    ftabLeft1(0),
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
    fContAddDrawOpt(0),
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
    fComboAddDrawOpt(0),
    fChkAuto(0),
    fComboMethod(0),
    fListNormalization(0),
    fComboCustom(0),
    fChkAddDrawOpt(0),
    fNmbSector(0),
    fLblSector(0),
    fChkCutZero(0),
    fChkAddCuts(0),
    fComboAddCuts(0), 
    fComboCustomFit(0),
    fChkSetMax(0),
    fChkSetMin(0),
    fChkGetMinMaxAuto(0),
    fTxtSetMax(0),
    fTxtSetMin(0) ,
    fContDrawOpt1D(0), 
    fcontDrawOpt1DSubLR(0),
    fContDrawOpt1DSubNSC(0), 
    fRadioNorm(0),
    fRadioSigma(0),
    fTxtSigmas(0),
    fContCumuLR(0),
    fContCumLeft(0),
    fContCumRight(0),
    fLblSigmaMax(0),
    fTxtSigmaMax(0),
    fRadioCumulative(0),
    fCheckCumulativePM(0),
    fRadioIntegrate(0),
    fContDrawOpt1DSubMML(0),
    fChkMean(0),
    fChkMedian(0),
    fChkLTM(0),
    fContStatOpt(0),
    fChkStatName(0),
    fChkStatEntries(0),
    fContStatMean(0),
    fChkStatMean(0),
    fChkStatMeanPM(0),
    fContStatRMS(0),
    fChkStatRMS(0),
    fChkStatRMSPM(0),
    fChkStatUnderflow(0),
    fChkStatOverflow(0),
    fChkStatIntegral(0),
    fContStatSkew(0),
    fChkStatSkewness(0),
    fChkStatSkewnessPM(0),
    fContStatKurt(0),
    fChkStatKurtosis(0),
    fChkStatKurtosisPM(0)
{
   //
   // AliTPCCalibViewerGUI constructor; fileName specifies the ROOT tree used for drawing 
   //
   SetCleanup(kDeepCleanup);
   
   // *****************************************************************************
   // ************************* content of this MainFrame *************************
   // *****************************************************************************
   // top level container with horizontal layout
   fContTopBottom = new TGCompositeFrame(this, w, h, kVerticalFrame | kFixedWidth | kFixedHeight);
   AddFrame(fContTopBottom, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   
   fContLCR = new TGCompositeFrame(fContTopBottom, w, h, kHorizontalFrame | kFixedWidth | kFixedHeight);
   fContTopBottom->AddFrame(fContLCR, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   
   // ***********************************************************************
   // ************************* content of fContLCR *************************
   // ***********************************************************************
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
   
   
   // ========================================================================
   // ************************* content of fContLeft *************************
   // ========================================================================
   // --- draw button and tabLeft ---
   // draw button
   fBtnDraw = new TGTextButton(fContLeft, "&Draw");
   fContLeft->AddFrame(fBtnDraw, new TGLayoutHints(kLHintsExpandX, 10, 10, 0, 0));
   //fBtnDraw->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoTest(=\"fBtnDraw clicked\")");
   fBtnDraw->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoDraw()");
   
   // three tabs on the left side:
   ftabLeft = new TGTab(fContLeft);
   fContLeft->AddFrame(ftabLeft, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 8, 0));
   ftabLeft0 = ftabLeft->AddTab("General");
   ftabLeft1 = ftabLeft->AddTab("More plot options");

   
      // **************************** content of tabLeft0 *******************************
      
      // draw options container *** fcontDrawOpt ***  " Plot options "
      fContDrawOpt = new TGGroupFrame(ftabLeft0, "Plot options", kVerticalFrame | kFitWidth | kFitHeight);
      ftabLeft0->AddFrame(fContDrawOpt, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
      fContDrawOptSub1D2D = new TGCompositeFrame(fContDrawOpt, 200, 23, kHorizontalFrame | kFitWidth | kFixedHeight);
      fContDrawOpt->AddFrame(fContDrawOptSub1D2D, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
         // ------------------------- content of fContDrawOpt -------------------------
         // -- radio1D, radio2D, chkAuto
         // 1D radio button
         fRadio1D = new TGRadioButton(fContDrawOptSub1D2D, "1D", 30);
         fContDrawOptSub1D2D->AddFrame(fRadio1D, new TGLayoutHints(kLHintsNormal, 0, 2, 0, 0));
         fRadio1D->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral()");
         
         // 2D radio button
         fRadio2D = new TGRadioButton(fContDrawOptSub1D2D, "2D", 31);
         fContDrawOptSub1D2D->AddFrame(fRadio2D, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
         fRadio2D->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral()");
         
         // additional draw options container
         fContAddDrawOpt = new TGCompositeFrame(fContDrawOpt, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
         fContDrawOpt->AddFrame(fContAddDrawOpt, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

            //  content of --- fContAddDrawOpt ---
            // addition draw options label
            fChkAddDrawOpt = new TGCheckButton(fContAddDrawOpt, "Draw options:");
            //fChkAddDrawOpt->SetTextJustify(kTextLeft);
            fContAddDrawOpt->AddFrame(fChkAddDrawOpt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
            fChkAddDrawOpt->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
            
            // additional draw options combo box
            fComboAddDrawOpt = new TGComboBox(fContAddDrawOpt);
            fComboAddDrawOpt->Resize(0, fBtnDraw->GetDefaultHeight());
            fComboAddDrawOpt->EnableTextInput(kTRUE);
            fContAddDrawOpt->AddFrame(fComboAddDrawOpt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
            fComboAddDrawOpt->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
            fComboAddDrawOpt->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
                  
         // automatic redraw check button
         fChkAuto = new TGCheckButton(fContDrawOpt, "auto redraw");
         fContDrawOpt->AddFrame(fChkAuto, new TGLayoutHints(kLHintsNormal, 0, 2, 0, 0));
               
      
      // *** predefined radio button ***  " Predefined "
      fRadioPredefined = new TGRadioButton(ftabLeft0, "Predefined: ", 13);
      ftabLeft0->AddFrame(fRadioPredefined, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      fRadioPredefined->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral()");
      
      // list of variables
      fListVariables = new TGListBox(ftabLeft0);
      ftabLeft0->AddFrame(fListVariables, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY, 10, 0, 0, 0));
      fListVariables->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
   
      
      // normalization options container *** fContPlotOpt ***
      //fContPlotOpt = new TGCompositeFrame(fContLeft, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
      fContPlotOpt = new TGGroupFrame(ftabLeft0, "Normalization options", kVerticalFrame | kFitWidth | kFitHeight);
      ftabLeft0->AddFrame(fContPlotOpt, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 0, 0, 0));

         // ------------------------- content of fContPlotOpt -------------------------
         // raw radio button
         fRadioRaw = new TGRadioButton(fContPlotOpt, "Raw", 10);
         fContPlotOpt->AddFrame(fRadioRaw, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fRadioRaw->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral()");
      
         // normalized radio button
         fRadioNormalized = new TGRadioButton(fContPlotOpt, "Normalized", 11);
         fContPlotOpt->AddFrame(fRadioNormalized, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fRadioNormalized->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral()");
      
         // normalized options container *** fContNormalized ***
         fContNormalized = new TGCompositeFrame(fContPlotOpt, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
         fContPlotOpt->AddFrame(fContNormalized, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 15, 0, 0, 0));
      
            // --- content of fContNormalized ---
            // --- combo box to select 'subtract' or 'divide', list of normalization variables
            // method drop down combo box
            fComboMethod = new TGComboBox(fContNormalized);
            fComboMethod->Resize(0, fBtnDraw->GetDefaultHeight());
            fContNormalized->AddFrame(fComboMethod, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
            fComboMethod->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
         
            // list of normalization variables
            fListNormalization = new TGListBox(fContNormalized);
            fContNormalized->AddFrame(fListNormalization, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
            fListNormalization->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");

      // custom radio button
      fRadioCustom = new TGRadioButton(ftabLeft0, "Custom: ", 12);
      ftabLeft0->AddFrame(fRadioCustom, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      fRadioCustom->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral()");
      
      // custom options container
      // --- fComboCustom --- the custom draw line
      fContCustom = new TGCompositeFrame(fContTopBottom, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
      fContTopBottom->AddFrame(fContCustom, new TGLayoutHints(kLHintsExpandX, 10, 0, 0, 0));
   
         // ------------------------- content of fContCustom -------------------------
         // text field for custom draw command
         fComboCustom = new TGComboBox(fContCustom);
         fComboCustom->Resize(0, fBtnDraw->GetDefaultHeight());
         fComboCustom->EnableTextInput(kTRUE);
         fContCustom->AddFrame(fComboCustom, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
         fComboCustom->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral(=42)");
         fComboCustom->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
      

      
      // **************************** content of tabLeft1 *******************************
      
      // draw options container *** fcontDrawOpt1D ***  " Plot options "
      fContDrawOpt1D = new TGGroupFrame(ftabLeft1, "1D Plot options", kVerticalFrame | kFitWidth | kFitHeight);
      ftabLeft1->AddFrame(fContDrawOpt1D, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
      
      fcontDrawOpt1DSubLR = new TGCompositeFrame(fContDrawOpt1D, 1, 1, kVerticalFrame | kFitWidth | kFitHeight);
      fContDrawOpt1D->AddFrame(fcontDrawOpt1DSubLR, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
         // ***** content of fContDrawOpt1DSubLR *****
         fContDrawOpt1DSubNSC = new TGCompositeFrame(fcontDrawOpt1DSubLR, 200, 23, kVerticalFrame | kFitWidth | kFitHeight);
         fcontDrawOpt1DSubLR->AddFrame(fContDrawOpt1DSubNSC, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         
            // --------------------------- content of fContDrawOpt1DSubNSC -----------------
            fRadioNorm = new TGRadioButton(fContDrawOpt1DSubNSC, "Normal", 110);
            fContDrawOpt1DSubNSC->AddFrame(fRadioNorm, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fRadioNorm->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
               
            fRadioSigma = new TGRadioButton(fContDrawOpt1DSubNSC, "Sigma", 111);
            fContDrawOpt1DSubNSC->AddFrame(fRadioSigma, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 0));
            fRadioSigma->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");

            fTxtSigmas = new TGTextEntry(fContDrawOpt1DSubNSC, "2; 4; 6", 111);
            fContDrawOpt1DSubNSC->AddFrame(fTxtSigmas, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 10, 15, 0, 0));
            fTxtSigmas->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtons1D(=111)");
               
            fContCumuLR = new TGCompositeFrame(fContDrawOpt1DSubNSC, 200, 23, kHorizontalFrame | kFitWidth | kFitHeight);
            fContDrawOpt1DSubNSC->AddFrame(fContCumuLR, new TGLayoutHints(kLHintsExpandX, 0, 0, 5, 0));
            
               fContCumLeft = new TGCompositeFrame(fContCumuLR, 200, 23, kVerticalFrame | kFitWidth | kFitHeight);
               fContCumuLR->AddFrame(fContCumLeft, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
                           
                  fRadioCumulative = new TGRadioButton(fContCumLeft, "Cumulative", 112);
                  fContCumLeft->AddFrame(fRadioCumulative, new TGLayoutHints(kLHintsNormal, 0, 10, 0, 0));
                  fRadioCumulative->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
                  
                  fCheckCumulativePM = new TGCheckButton(fContCumLeft, "Plus/Minus");
                  fContCumLeft->AddFrame(fCheckCumulativePM, new TGLayoutHints(kLHintsNormal, 10, 15, 0, 0));
                  fCheckCumulativePM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
                  
                  fRadioIntegrate = new TGRadioButton(fContCumLeft, "Integrate", 113);
                  fContCumLeft->AddFrame(fRadioIntegrate, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 0));
                  fRadioIntegrate->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
                  
               fContCumRight = new TGCompositeFrame(fContCumuLR, 200, 23, kVerticalFrame | kFitWidth | kFitHeight);
               fContCumuLR->AddFrame(fContCumRight, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
               
                  fLblSigmaMax = new TGLabel(fContCumRight, "SigmaMax:");
                  fLblSigmaMax->SetTextJustify(kTextLeft);
                  fContCumRight->AddFrame(fLblSigmaMax, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 5, 0, 0, 0));

                  fTxtSigmaMax = new TGTextEntry(fContCumRight, "5", 112);
                  fContCumRight->AddFrame(fTxtSigmaMax, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 10, 15, 0, 0));
                  fTxtSigmaMax->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtons1D(=112)");
             
            
         fContDrawOpt1DSubMML = new TGCompositeFrame(fcontDrawOpt1DSubLR, 200, 23, kHorizontalFrame | kFitWidth | kFitHeight);
         fcontDrawOpt1DSubLR->AddFrame(fContDrawOpt1DSubMML, new TGLayoutHints(kLHintsExpandX, 0, 0, 5, 0));
         
            // -------------- content of fcontDrawOpt1DSubLR
            fChkMean = new TGCheckButton(fContDrawOpt1DSubMML, "Mean");
            fContDrawOpt1DSubMML->AddFrame(fChkMean, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkMean->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");

            fChkMedian = new TGCheckButton(fContDrawOpt1DSubMML, "Median");
            fContDrawOpt1DSubMML->AddFrame(fChkMedian, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkMedian->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");

            fChkLTM = new TGCheckButton(fContDrawOpt1DSubMML, "LTM");
            fContDrawOpt1DSubMML->AddFrame(fChkLTM, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkLTM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
            
      
      // statistic options container *** fcontStatOpt1D ***  " Statistic options "      
      fContStatOpt = new TGGroupFrame(ftabLeft1, "Statistic options", kVerticalFrame | kFitWidth | kFitHeight);
      ftabLeft1->AddFrame(fContStatOpt, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
      
         fChkStatName = new TGCheckButton(fContStatOpt, "Name");
         fContStatOpt->AddFrame(fChkStatName, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkStatName->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
      
         fChkStatEntries = new TGCheckButton(fContStatOpt, "Entries");
         fContStatOpt->AddFrame(fChkStatEntries, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkStatEntries->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
      
         fContStatMean = new TGCompositeFrame(fContStatOpt, 1, 1, kHorizontalFrame | kFitWidth | kFitHeight);
         fContStatOpt->AddFrame(fContStatMean, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
            fChkStatMean = new TGCheckButton(fContStatMean, "Mean");
            fContStatMean->AddFrame(fChkStatMean, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkStatMean->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            
            fChkStatMeanPM = new TGCheckButton(fContStatMean, "+- Error");
            fContStatMean->AddFrame(fChkStatMeanPM, new TGLayoutHints(kLHintsNormal, 10, 0, 0, 0));
            fChkStatMeanPM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");

         fContStatRMS = new TGCompositeFrame(fContStatOpt, 1, 1, kHorizontalFrame | kFitWidth | kFitHeight);
         fContStatOpt->AddFrame(fContStatRMS, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
            fChkStatRMS = new TGCheckButton(fContStatRMS, "RMS");
            fContStatRMS->AddFrame(fChkStatRMS, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkStatRMS->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            
            fChkStatRMSPM = new TGCheckButton(fContStatRMS, "+- Error");
            fContStatRMS->AddFrame(fChkStatRMSPM, new TGLayoutHints(kLHintsNormal, 10, 0, 0, 0));
            fChkStatRMSPM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");

         fChkStatUnderflow = new TGCheckButton(fContStatOpt, "Underflow");
         fContStatOpt->AddFrame(fChkStatUnderflow, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkStatUnderflow->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
      
         fChkStatOverflow = new TGCheckButton(fContStatOpt, "Overflow");
         fContStatOpt->AddFrame(fChkStatOverflow, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkStatOverflow->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
      
         fChkStatIntegral = new TGCheckButton(fContStatOpt, "Integral");
         fContStatOpt->AddFrame(fChkStatIntegral, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkStatIntegral->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
      
         fContStatSkew = new TGCompositeFrame(fContStatOpt, 1, 1, kHorizontalFrame | kFitWidth | kFitHeight);
         fContStatOpt->AddFrame(fContStatSkew, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
            fChkStatSkewness = new TGCheckButton(fContStatSkew, "Skewness");
            fContStatSkew->AddFrame(fChkStatSkewness, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkStatSkewness->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            
            fChkStatSkewnessPM = new TGCheckButton(fContStatSkew, "+- Error");
            fContStatSkew->AddFrame(fChkStatSkewnessPM, new TGLayoutHints(kLHintsNormal, 10, 0, 0, 0));
            fChkStatSkewnessPM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");

         fContStatKurt = new TGCompositeFrame(fContStatOpt, 1, 1, kHorizontalFrame | kFitWidth | kFitHeight);
         fContStatOpt->AddFrame(fContStatKurt, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
            fChkStatKurtosis = new TGCheckButton(fContStatKurt, "Kurtosis");
            fContStatKurt->AddFrame(fChkStatKurtosis, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkStatKurtosis->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            
            fChkStatKurtosisPM = new TGCheckButton(fContStatKurt, "+- Error");
            fContStatKurt->AddFrame(fChkStatKurtosisPM, new TGLayoutHints(kLHintsNormal, 10, 0, 0, 0));
            fChkStatKurtosisPM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");


         
         

   // ==========================================================================
   // ************************* content of fContCenter *************************
   // ========================================================================
   // main drawing canvas
   fCanvMain = new TRootEmbeddedCanvas("Main Canvas", fContCenter, 200, 200, kFitWidth | kFitHeight);
   fContCenter->AddFrame(fCanvMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
      
   
   
   
   // =========================================================================   
   // ************************* content of fContRight *************************
   // ========================================================================
   // cut options container
   //fContCuts = new TGCompositeFrame(fContRight, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
   fContCuts = new TGGroupFrame(fContRight, "Cuts", kVerticalFrame | kFitWidth | kFitHeight);
   fContRight->AddFrame(fContCuts, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   
      // ************************* content of fContCuts *************************
      // TPC radio button
      fRadioTPC = new TGRadioButton(fContCuts, "whole TPC", 20);
      fContCuts->AddFrame(fRadioTPC, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      fRadioTPC->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsRight()");
   
      // side A radio button
      fRadioSideA = new TGRadioButton(fContCuts, "side A", 21);
      fContCuts->AddFrame(fRadioSideA, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      fRadioSideA->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsRight()");
   
      // side C radio button
      fRadioSideC = new TGRadioButton(fContCuts, "side C", 22);
      fContCuts->AddFrame(fRadioSideC, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      fRadioSideC->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsRight()");
   
      // sector radio button
      fRadioSector = new TGRadioButton(fContCuts, "sector", 23);
      fContCuts->AddFrame(fRadioSector, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      fRadioSector->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsRight()");
   
      // sector options container
      fContSector = new TGCompositeFrame(fContCuts, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
      fContCuts->AddFrame(fContSector, new TGLayoutHints(kLHintsExpandX, 5, 0, 0, 0));
      
         // ------------------------- content of fContSector -------------------------
         // sector number entry
         fNmbSector = new TGNumberEntry(fContSector, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 71);
         fContSector->AddFrame(fNmbSector, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
         fNmbSector->Connect("ValueSet(Long_t)", "AliTPCCalibViewerGUI", this, "ChangeSector()");
         
         // sector number label
         fLblSector = new TGLabel(fContSector, "IROC, A");
         fLblSector->SetTextJustify(kTextLeft);
         fContSector->AddFrame(fLblSector, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 5, 0, 0, 0));
      
      // additional cuts check button
      fChkCutZero = new TGCheckButton(fContCuts, "Cut zeros");
      fContCuts->AddFrame(fChkCutZero, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
      fChkCutZero->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
   
      // additional cuts check button
      fChkAddCuts = new TGCheckButton(fContCuts, "additional cuts");
      fContCuts->AddFrame(fChkAddCuts, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
      fChkAddCuts->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
   
      // additional cuts container
      fContAddCuts = new TGCompositeFrame(fContCuts, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
      fContCuts->AddFrame(fContAddCuts, new TGLayoutHints(kLHintsExpandX, -5, -5, 0, 0));
      
         // ------------------------- content of fContAddCuts -------------------------
         // combo text field for additional cuts
         fComboAddCuts = new TGComboBox(fContAddCuts);
         fComboAddCuts->Resize(0, fBtnDraw->GetDefaultHeight());
         fComboAddCuts->EnableTextInput(kTRUE);
         fContAddCuts->AddFrame(fComboAddCuts, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
         fComboAddCuts->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
         fComboAddCuts->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
   
   
   // Scaling options container
   fContScaling = new TGGroupFrame(fContRight, "Scaling", kVerticalFrame | kFitWidth | kFitHeight);
   fContRight->AddFrame(fContScaling, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

      // ************************* content of fContScaling *************************
      // SetMaximum container
      fContSetMax = new TGCompositeFrame(fContScaling, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
      fContScaling->AddFrame(fContSetMax, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
         // ------------------------- content of fContSetMax -------------------------
         // SetMaximum - checkbox
         fChkSetMax = new TGCheckButton(fContSetMax, "Set fixed max.");
         fContSetMax->AddFrame(fChkSetMax, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkSetMax->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsRight()");
         
         // text field for maximum value
         fTxtSetMax = new TGTextEntry(fContSetMax, "", 41);
         fContSetMax->AddFrame(fTxtSetMax, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
         fTxtSetMax->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsRight()");
   
      // SetMinimum container
      fContSetMin = new TGCompositeFrame(fContScaling, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
      fContScaling->AddFrame(fContSetMin, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
         // ------------------------- content of fContSetMin -------------------------
         // SetMinimum - checkbox
         fChkSetMin = new TGCheckButton(fContSetMin, "Set fixed min.");
         fContSetMin->AddFrame(fChkSetMin, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkSetMin->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
         
         // text field for minimum value
         fTxtSetMin = new TGTextEntry(fContSetMin, "", 40);
         fContSetMin->AddFrame(fTxtSetMin, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
         fTxtSetMin->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsRight()");
      
      // get Min & Max from Plot - button
      fBtnGetMinMax = new TGTextButton(fContScaling, "&Get scale from plot");
      fContScaling->AddFrame(fBtnGetMinMax, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      fBtnGetMinMax->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "GetMinMax()");
      
      // GetMinMaxAuto - checkbox
      fChkGetMinMaxAuto = new TGCheckButton(fContScaling, "Get Min + Max auto.");
      fContScaling->AddFrame(fChkGetMinMaxAuto, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
      fChkGetMinMaxAuto->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");

      
   // Fit options container
   fContFit = new TGGroupFrame(fContRight, "Custom Fit", kVerticalFrame | kFitWidth | kFitHeight);
   fContRight->AddFrame(fContFit, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

      // ------------------------- content of fContFit -------------------------
      // container for additional fits
      fContAddFit = new TGCompositeFrame(fContFit, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
      fContFit->AddFrame(fContAddFit, new TGLayoutHints(kLHintsExpandX, -5, -5, 0, 0));
   
         // --- content of fContAddFit ---
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
         
   // set default button states
   fRadioPredefined->SetState(kButtonDown);
   fRadioRaw->SetState(kButtonDown);
   fRadioTPC->SetState(kButtonDown);
   fRadio1D->SetState(kButtonDown);
   fChkAuto->SetState(kButtonDown);
   fChkAddCuts->SetState(kButtonUp);
   fChkGetMinMaxAuto->SetState(kButtonDown);
   fChkSetMin->SetState(kButtonUp);
   fChkSetMax->SetState(kButtonUp);
   fRadioNorm->SetState(kButtonDown);
   fRadioSigma->SetState(kButtonUp);
   fRadioCumulative->SetState(kButtonUp);
   fChkMean->SetState(kButtonDown);
   fCheckCumulativePM->SetState(kButtonUp);
   
   fChkStatName->SetState(kButtonDown);
   fChkStatEntries->SetState(kButtonDown);
   fChkStatMean->SetState(kButtonDown);
   fChkStatRMS->SetState(kButtonDown);
//    fChkStatMeanPM->SetState(kButtonUp);
//    fChkStatRMSPM->SetState(kButtonUp);
//    fChkStatUnderflow->SetState(kButtonUp);
//    fChkStatOverflow->SetState(kButtonUp);
//    fChkStatIntegral->SetState(kButtonUp);
//    fChkStatSkewness->SetState(kButtonUp);
//    fChkStatSkewnessPM->SetState(kButtonUp);
//    fChkStatKurtosis->SetState(kButtonUp);
//    fChkStatKurtosisPM->SetState(kButtonUp);
      
   // ======================================================================   
   // ************************* Display everything *************************
   // ======================================================================

   if (fileName) Initialize(fileName);
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
    ftabLeft(0),
    ftabLeft0(0),
    ftabLeft1(0),
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
    fContAddDrawOpt(0),
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
    fComboAddDrawOpt(0),
    fChkAuto(0),
    fComboMethod(0),
    fListNormalization(0),
    fComboCustom(0),
    fChkAddDrawOpt(0),
    fNmbSector(0),
    fLblSector(0),
    fChkCutZero(0),
    fChkAddCuts(0),
    fComboAddCuts(0), 
    fComboCustomFit(0),
    fChkSetMax(0),
    fChkSetMin(0),
    fChkGetMinMaxAuto(0),
    fTxtSetMax(0),
    fTxtSetMin(0), 
    fContDrawOpt1D(0),
    fcontDrawOpt1DSubLR(0),
    fContDrawOpt1DSubNSC(0), 
    fRadioNorm(0),
    fRadioSigma(0),
    fTxtSigmas(0),
    fContCumuLR(0),
    fContCumLeft(0),
    fContCumRight(0),
    fLblSigmaMax(0),
    fTxtSigmaMax(0),
    fRadioCumulative(0),
    fCheckCumulativePM(0),
    fRadioIntegrate(0),
    fContDrawOpt1DSubMML(0),
    fChkMean(0),
    fChkMedian(0),
    fChkLTM(0), 
    fContStatOpt(0),
    fChkStatName(0),
    fChkStatEntries(0),
    fContStatMean(0),
    fChkStatMean(0),
    fChkStatMeanPM(0),
    fContStatRMS(0),
    fChkStatRMS(0),
    fChkStatRMSPM(0),
    fChkStatUnderflow(0),
    fChkStatOverflow(0),
    fChkStatIntegral(0),
    fContStatSkew(0),
    fChkStatSkewness(0),
    fChkStatSkewnessPM(0),
    fContStatKurt(0),
    fChkStatKurtosis(0),
    fChkStatKurtosisPM(0)
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
   // 
   // Destructor
   // 
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

   // fill fComboAddDrawOpt with some additional drawing options
   fComboAddDrawOpt->AddEntry("profbox", 0);
   fComboAddDrawOpt->AddEntry("profcolz", 1);
   fComboAddDrawOpt->AddEntry("profcont0", 2);
   fComboAddDrawOpt->AddEntry("proflego", 3);
   fComboAddDrawOpt->AddEntry("proflego2", 4);
   fComboAddDrawOpt->AddEntry("profsurf", 5);
   fComboAddDrawOpt->AddEntry("profsurf1", 6);
   fComboAddDrawOpt->AddEntry("profsurf2", 7);
   fComboAddDrawOpt->AddEntry("box", 8);
   fComboAddDrawOpt->AddEntry("colz", 9);
   fComboAddDrawOpt->AddEntry("cont0", 10);
   fComboAddDrawOpt->AddEntry("lego", 11);
   fComboAddDrawOpt->AddEntry("lego2", 12);
   fComboAddDrawOpt->AddEntry("surf", 13);
   fComboAddDrawOpt->AddEntry("surf1", 14);
   fComboAddDrawOpt->AddEntry("surf2", 15);


   fListVariables->Select(0);
   fListNormalization->Select(0);
   fComboMethod->Select(0);

   //fCanvMain->GetCanvas()->ToggleEventStatus(); // klappt nicht
   //fCanvMain->GetCanvas()->GetCanvasImp()->ShowStatusBar(kTRUE); // klappt auch nicht
   fListVariables->IntegralHeight(kFALSE);         // naja
   fListNormalization->IntegralHeight(kFALSE);     // naja
   DoDraw();
}



void AliTPCCalibViewerGUI::HandleButtonsGeneral(Int_t id) {
   //
   // handles mutual radio button exclusions
   // for general Tab
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
         // fComboMethod->UnmapWindow();
         // fListNormalization->UnmapWindow();
         break;
      case 11:             // fRadioNormalized
         fRadioRaw->SetState(kButtonUp);
         fRadioPredefined->SetState(kButtonDown);
         fRadioCustom->SetState(kButtonUp);
         break;
      case 12:             // fRadioCustom
         fRadioPredefined->SetState(kButtonUp);
         // fComboCustom->SetEnabled(kFALSE);
         // fRadioNormalized->SetState(kButtonUp);
         break;
      case 42:             // fComboCustom
         fRadioCustom->SetState(kButtonDown);
         fRadioPredefined->SetState(kButtonUp);
         break;
      case 13:             // fRadioPredefined
         fRadioCustom->SetState(kButtonUp);
         // fComboCustom->SetEnabled(kTRUE);
         //f RadioNormalized->SetState(kButtonUp);
         break;
      //--------
      case 30:             // fRadio1D
         fRadio2D->SetState(kButtonUp);
         break;
      case 31:             // fRadio2D
         fRadio1D->SetState(kButtonUp);
         break;
   }
   DoNewSelection();
}


void AliTPCCalibViewerGUI::HandleButtons1D(Int_t id) {
   //
   // handles mutual radio button exclusions
   // 1D-Tab buttons
   //
   
   if (id == -1) {
      TGButton *btn = (TGButton *) gTQSender;
      id = btn->WidgetId();
   }
   switch (id) {
      case 110:            // 1D draw normal
         fRadioNorm->SetState(kButtonDown);
         fRadioSigma->SetState(kButtonUp);
         fRadioCumulative->SetState(kButtonUp);
         fRadioIntegrate->SetState(kButtonUp);
         break;
      case 111:            // 1D draw sigma
         fRadioNorm->SetState(kButtonUp);
         fRadioSigma->SetState(kButtonDown);
         fRadioCumulative->SetState(kButtonUp);
         fRadioIntegrate->SetState(kButtonUp);
         break;
      case 112:            // 1D draw cumulative
         fRadioNorm->SetState(kButtonUp);
         fRadioSigma->SetState(kButtonUp);
         fRadioCumulative->SetState(kButtonDown);
         fRadioIntegrate->SetState(kButtonUp);
         break;
      case 113:            // 1D draw integral
         fRadioNorm->SetState(kButtonUp);
         fRadioSigma->SetState(kButtonUp);
         fRadioCumulative->SetState(kButtonUp);
         fRadioIntegrate->SetState(kButtonDown);
         break;
   }
   DoNewSelection();
}


void AliTPCCalibViewerGUI::HandleButtons2D(Int_t id) {
   //
   // handles mutual radio button exclusions
   // 2D-Tab buttons
   //
   if (id == -1) {
      TGButton *btn = (TGButton *) gTQSender;
      id = btn->WidgetId();
   }

   switch (id) {
      case 211:
      break;
   }
   DoNewSelection();
}


void AliTPCCalibViewerGUI::HandleButtonsStat(Int_t id) {
   // 
   // handles statistic check boxes 
   // checks each checkbox if checked
   // if the checkbox is checked, appends 'n' for name, 'e' for entries, ...
   // to a TString, passes this TString to gStyle->SetOptStat(...)
   // 
   id = id; // to avoid compiler warnings 
   TString statOpt("");
   if (fChkStatName->GetState() == kButtonDown) statOpt.Append("n");
   if (fChkStatEntries->GetState() == kButtonDown) statOpt.Append("e");
   if (fChkStatMean->GetState() == kButtonDown && fChkStatMeanPM->GetState() == kButtonUp) statOpt.Append("m");
   if (fChkStatMeanPM->GetState() == kButtonDown) statOpt.Append("M");
   if (fChkStatRMS->GetState() == kButtonDown && fChkStatRMSPM->GetState() == kButtonUp) statOpt.Append("r");
   if (fChkStatRMSPM->GetState() == kButtonDown) statOpt.Append("R");
   if (fChkStatUnderflow->GetState() == kButtonDown) statOpt.Append("u");
   if (fChkStatOverflow->GetState() == kButtonDown) statOpt.Append("o");
   if (fChkStatIntegral->GetState() == kButtonDown) statOpt.Append("i");
   if (fChkStatSkewness->GetState() == kButtonDown && fChkStatSkewnessPM->GetState() == kButtonUp) statOpt.Append("s");
   if (fChkStatSkewnessPM->GetState() == kButtonDown) statOpt.Append("S");
   if (fChkStatKurtosis->GetState() == kButtonDown && fChkStatKurtosisPM->GetState() == kButtonUp) statOpt.Append("k");
   if (fChkStatKurtosisPM->GetState() == kButtonDown) statOpt.Append("K");

   gStyle->SetOptStat(statOpt);
   DoNewSelection();
}


void AliTPCCalibViewerGUI::HandleButtonsRight(Int_t id) {
   //
   // handles mutual radio button exclusions
   // right side buttons
   //
    if (id == -1) {
      TGButton *btn = (TGButton *) gTQSender;
      id = btn->WidgetId();
   }

   switch (id) {
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
      case 40:             // fTxtSetMin
         fChkSetMin->SetState(kButtonDown);
         break;
      case 41:             // fTxtSetMax
         fChkSetMax->SetState(kButtonDown);
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
   desiredData += "~";

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
         if (normalizationData.CompareTo("FitParLocal") == 0)
            formulaStr = "lx~ ++ ly~ ++ lx~^2 ++ ly~^2 ++ lx~*ly~";
         if (normalizationData.CompareTo("FitParGlobal") == 0)
            formulaStr = "gx~ ++ gy~ ++ gx~^2 ++ gy~^2 ++ gx~*gy~";
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
   
         
   if (fChkCutZero->GetState() == kButtonDown) {
      cutsStr += desiredData.Data();
      cutsStr += "!=0";
      if (fChkAddCuts->GetState() == kButtonDown) cutsStr += " && ";
   }
   if (fChkAddCuts->GetState() == kButtonDown)
      cutsStr += fComboAddCuts->GetTextEntry()->GetText();

   TString addDrawOpt("");
   if (fChkAddDrawOpt->GetState() == kButtonDown)
      addDrawOpt += fComboAddDrawOpt->GetTextEntry()->GetText();
   
   // draw finally
   for (Int_t i = 0; i < fCanvMain->GetCanvas()->GetListOfPrimitives()->GetEntries(); i++) {
      if (strcmp(fCanvMain->GetCanvas()->GetListOfPrimitives()->At(i)->ClassName(), "TFrame") != 0)
         fCanvMain->GetCanvas()->GetListOfPrimitives()->At(i)->Delete();
   }
   //fCanvMain->GetCanvas()->Clear();
   fCanvMain->GetCanvas()->cd();
   Int_t entries = -1;
   if (fRadio1D->GetState() == kButtonDown){
      // 1D-Drawing
      TString strSigmaMax(fTxtSigmaMax->GetText());  // get sigmaMax from text enty
      Double_t sigmaMax = (strSigmaMax.IsFloat()) ? strSigmaMax.Atof() : 5; // convert to double, if not convertable, set to 5
      Bool_t plotMean   = fChkMean->GetState() == kButtonDown;
      Bool_t plotMedian = fChkMedian->GetState() == kButtonDown;
      Bool_t plotLTM    = fChkLTM->GetState() == kButtonDown;
      if (fRadioNorm->GetState() == kButtonDown)  // normal 1D drawing
         entries = fViewer->EasyDraw1D(desiredData.Data(), sectorStr.Data(), cutsStr.Data(), addDrawOpt.Data());
      if (fRadioSigma->GetState() == kButtonDown) // sigma 1D drawing
         entries = fViewer->DrawHisto1D(desiredData.Data(), sectorStr.Data(), cutsStr.Data(), // 
            fTxtSigmas->GetText(), plotMean, plotMedian, plotLTM);
      if (fRadioCumulative->GetState() == kButtonDown)  // cumulative 1D drawing
         entries = fViewer->SigmaCut(desiredData.Data(), sectorStr.Data(), cutsStr.Data(), //
            sigmaMax, plotMean, plotMedian, plotLTM, // 
            fCheckCumulativePM->GetState() == kButtonDown, fTxtSigmas->GetText(), /* Float_t sigmaStep =*/ -1);
      if (fRadioIntegrate->GetState() == kButtonDown)  // integral 1D drawing  
         entries = fViewer->Integrate(desiredData.Data(), sectorStr.Data(), cutsStr.Data(), //
            sigmaMax, plotMean, plotMedian, plotLTM, // 
            fTxtSigmas->GetText(), /* Float_t sigmaStep =*/ -1);            
   }
   else if (fRadio2D->GetState() == kButtonDown) {
      // 2D-Drawing
      entries = fViewer->EasyDraw(desiredData.Data(), sectorStr.Data(), cutsStr.Data(), addDrawOpt.Data());
   }
   if (entries == -1) return; // nothing was drawn, there is no histogram to get min and max
   
   
   // get or set Min & Max 
   // 
   // search for histogram
   TList* listOfPrimitives = fCanvMain->GetCanvas()->GetListOfPrimitives();
   TObject* ptr = 0;
   for (Int_t i = 0; i < listOfPrimitives->GetEntries(); i++) {
      ptr = listOfPrimitives->At(i);
      if ( ptr->InheritsFrom("TH1") ) break;
   }
   if ( ptr == 0 || !ptr->InheritsFrom("TH1") ) {  // if the loop did not find a TH1
      fCanvMain->GetCanvas()->Update();
      return;
      // unable to find histogram, no min and max wil be read out
   }
   TH1 *hist = (TH1*)ptr; 
   TString minTxt(fTxtSetMin->GetText());
   TString maxTxt(fTxtSetMax->GetText());
   // set min and max according to specified values, if checkbox is checked
   if (fChkSetMax->GetState() == kButtonDown && (maxTxt.IsDigit() || maxTxt.IsFloat()) )
      hist->SetMaximum(maxTxt.Atof());
   if (fChkSetMin->GetState() == kButtonDown && (minTxt.IsDigit() || minTxt.IsFloat()) )
      hist->SetMinimum(minTxt.Atof());
   // get min and max from plot       
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

void AliTPCCalibViewerGUI::AddFitFunction() const { 
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

