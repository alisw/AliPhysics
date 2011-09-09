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
  
*/
//                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTPCCalibViewerGUI.h"
#include "AliTPCCalPad.h"
#include <iostream>

#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPad.h>

#include <TROOT.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TVector.h>
#include <string.h>
#include <TH1.h>
#include "TStyle.h"
#include "TGFileDialog.h"
#include "TGInputDialog.h"
#include "AliTPCCalibViewer.h"
#include "AliTPCPreprocessorOnline.h"

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
    fPreprocessor(0),
    fContTopBottom(0),
    fContLCR(0),
    fContLeft(0),
    ftabLeft(0),
    ftabLeft0(0),
    ftabLeft1(0),
    ftabRight(0),
    fTabRight0(0),
    fTabRight1(0),
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
    fRadioROC(0),
    fRadioSector(0),
    fComboAddDrawOpt(0),
    fChkAuto(0),
    fChkAutoAppend(0),
    fComboMethod(0),
    fListNormalization(0),
    fComboCustom(0),
    fLblCustomDraw(0),
    fChkAddDrawOpt(0),
    fNmbSector(0),
    fLblSector(0),
    fChkCutZero(0),
    fChkAddCuts(0),
    fLblAddCuts(0),
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
    fChkStatKurtosisPM(0),
    fBtnUnchekAll(0),
    fContLabeling(0),
    fChkLabelTitle(0),
    fTxtLabelTitle(0),
    fChkLabelXaxis(0),
    fTxtLabelXaxis(0),
    fChkLabelYaxis(0),
    fTxtLabelYaxis(0),
    fChkLabelGetAuto(0),
    fContSave(0),
    fBtnSave(0),
    fContAddSaveOpt(0),
    fChkAddSaveOpt(0),
    fComboAddSaveOpt(0),
    fContExport(0),
    fContAddExport(0),
    fComboExportName(0),
    fBtnExport(0),
    fBtnAddNorm(0), 
    fContTree(0),
    fBtnDumpToFile(0),
    fBtnLoadTree(0),
    fChkAddAsReference(0),
    fTxtRefName(0), 
    fInitialized(0)
{
   //
   // AliTPCCalibViewerGUI constructor; fileName specifies the ROOT tree used for drawing 
   //

   // draw the GUI:
   fPreprocessor = new AliTPCPreprocessorOnline();
   DrawGUI(p, w, h);
   // initialize the AliTPCCalibViewer:
   if (fileName) Initialize(fileName);
   // set default button states:
   SetInitialValues();
   // do first drawing: 
   if (fileName) DoDraw();
}


void AliTPCCalibViewerGUI::DrawGUI(const TGWindow */*p*/, UInt_t w, UInt_t h) {
   // 
   // draw the GUI
   // 
   // ======================================================================   
   // ************************* Display everything *************************
   // ======================================================================
   
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
   fBtnDraw->SetToolTipText("Press here to draw according to selections.");
   
   // tabs on the left side:
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
         fRadio1D->SetToolTipText("1D drawing \nSelect this if you want to have the full control for the custom draw.");
         
         // 2D radio button
         fRadio2D = new TGRadioButton(fContDrawOptSub1D2D, "2D", 31);
         fContDrawOptSub1D2D->AddFrame(fRadio2D, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
         fRadio2D->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral()");
         fRadio2D->SetToolTipText("2D drawing");
         
         // additional draw options container
         fContAddDrawOpt = new TGCompositeFrame(fContDrawOpt, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
         fContDrawOpt->AddFrame(fContAddDrawOpt, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

            //  content of --- fContAddDrawOpt ---
            // addition draw options label
            fChkAddDrawOpt = new TGCheckButton(fContAddDrawOpt, "Draw options:");
            //fChkAddDrawOpt->SetTextJustify(kTextLeft);
            fContAddDrawOpt->AddFrame(fChkAddDrawOpt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
            fChkAddDrawOpt->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
            fChkAddDrawOpt->SetToolTipText("Enter additional draw options like 'prof' or 'colz' here.\nBe careful with the option 'same' for 2D drawings as it will crash (ROOT feature).");
            
            // additional draw options combo box
            fComboAddDrawOpt = new TGComboBox(fContAddDrawOpt);
            fComboAddDrawOpt->Resize(0, fBtnDraw->GetDefaultHeight());
            fComboAddDrawOpt->EnableTextInput(kTRUE);
            fContAddDrawOpt->AddFrame(fComboAddDrawOpt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
            fComboAddDrawOpt->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral(=14)");
            fComboAddDrawOpt->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
                  
         // automatic redraw check button
         fChkAuto = new TGCheckButton(fContDrawOpt, "Auto redraw");
         fContDrawOpt->AddFrame(fChkAuto, new TGLayoutHints(kLHintsNormal, 0, 2, 0, 0));
         fChkAuto->SetToolTipText("Decide if you want an automatic redraw on each new selection.\nNot recommended on a slow machine, during remote connection or if your draw option is 'same'.");
         
         // automatic append ending check button
         fChkAutoAppend = new TGCheckButton(fContDrawOpt, "Auto add appending");
         fContDrawOpt->AddFrame(fChkAutoAppend, new TGLayoutHints(kLHintsNormal, 0, 2, 0, 0));
         fChkAutoAppend->SetToolTipText("Tries to repair your custom draw string or custom cut string, if you forgot '~' or '.fElements' \nThis function may be buggy!");
               
      
      // *** predefined radio button ***  " Predefined "
      fRadioPredefined = new TGRadioButton(ftabLeft0, "Predefined: ", 13);
      ftabLeft0->AddFrame(fRadioPredefined, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      fRadioPredefined->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral()");
      fRadioPredefined->SetToolTipText("Draw predefined variables according to selection.");
      
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
         fRadioRaw->SetToolTipText("Plot without normalization");
      
         // normalized radio button
         fRadioNormalized = new TGRadioButton(fContPlotOpt, "Normalized", 11);
         fContPlotOpt->AddFrame(fRadioNormalized, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fRadioNormalized->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral()");
         fRadioNormalized->SetToolTipText("Normalize data");
      
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
      fRadioCustom->SetToolTipText("Draw data according to user specific text entry in the 'Custom Draw' line. Remember '~' (= '.fElements')!");
      // custom options container is located further down
      
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
            fRadioNorm->SetToolTipText("Produce a normal 1D plot, a histogram of the selected data.");
               
            fRadioSigma = new TGRadioButton(fContDrawOpt1DSubNSC, "Sigma", 111);
            fContDrawOpt1DSubNSC->AddFrame(fRadioSigma, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 0));
            fRadioSigma->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
            fRadioSigma->SetToolTipText("Draw a normal histogram, but also lines that indicate the mean/median/LTM \nand sigmas of the selected data.");

            fTxtSigmas = new TGTextEntry(fContDrawOpt1DSubNSC, "2; 4; 6", 111);
            fContDrawOpt1DSubNSC->AddFrame(fTxtSigmas, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 10, 15, 0, 0));
            fTxtSigmas->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtons1D(=111)");
            fTxtSigmas->SetToolTipText("Enter sigma intervals you would like to be indicated by lines. \nExample: '2; 4; 6'");
               
            fContCumuLR = new TGCompositeFrame(fContDrawOpt1DSubNSC, 200, 23, kHorizontalFrame | kFitWidth | kFitHeight);
            fContDrawOpt1DSubNSC->AddFrame(fContCumuLR, new TGLayoutHints(kLHintsExpandX, 0, 0, 5, 0));
            
               fContCumLeft = new TGCompositeFrame(fContCumuLR, 200, 23, kVerticalFrame | kFitWidth | kFitHeight);
               fContCumuLR->AddFrame(fContCumLeft, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
                           
                  fRadioCumulative = new TGRadioButton(fContCumLeft, "Cumulative", 112);
                  fContCumLeft->AddFrame(fRadioCumulative, new TGLayoutHints(kLHintsNormal, 0, 10, 0, 0));
                  fRadioCumulative->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
                  fRadioCumulative->SetToolTipText("Draw the cumulative (SigmaCut) of the given selection. \nThe data distribution is integrated, starting from the mean/median/LTM.");
                  
                  fCheckCumulativePM = new TGCheckButton(fContCumLeft, "Plus/Minus");
                  fContCumLeft->AddFrame(fCheckCumulativePM, new TGLayoutHints(kLHintsNormal, 10, 15, 0, 0));
                  fCheckCumulativePM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
                  fCheckCumulativePM->SetToolTipText("Decide whether you want the cumulative integration for each direction (+/-) \nor only for the absolute distance to the mean/median/LTM value.");
                  
                  fRadioIntegrate = new TGRadioButton(fContCumLeft, "Integrate", 113);
                  fContCumLeft->AddFrame(fRadioIntegrate, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 0));
                  fRadioIntegrate->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
                  fRadioIntegrate->SetToolTipText("Draw the integral of the given selection.");
                  
               fContCumRight = new TGCompositeFrame(fContCumuLR, 200, 23, kVerticalFrame | kFitWidth | kFitHeight);
               fContCumuLR->AddFrame(fContCumRight, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
               
                  fLblSigmaMax = new TGLabel(fContCumRight, "SigmaMax:");
                  fLblSigmaMax->SetTextJustify(kTextLeft);
                  fContCumRight->AddFrame(fLblSigmaMax, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 5, 0, 0, 0));

                  fTxtSigmaMax = new TGTextEntry(fContCumRight, "5", 112);
                  fContCumRight->AddFrame(fTxtSigmaMax, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 10, 15, 0, 0));
                  fTxtSigmaMax->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtons1D(=112)");
                  fTxtSigmaMax->SetToolTipText("Enter up to which multiple of sigma you want to integrate.");
             
            
         fContDrawOpt1DSubMML = new TGCompositeFrame(fcontDrawOpt1DSubLR, 200, 23, kHorizontalFrame | kFitWidth | kFitHeight);
         fcontDrawOpt1DSubLR->AddFrame(fContDrawOpt1DSubMML, new TGLayoutHints(kLHintsExpandX, 0, 0, 5, 0));
         
            // -------------- content of fcontDrawOpt1DSubLR
            fChkMean = new TGCheckButton(fContDrawOpt1DSubMML, "Mean");
            fContDrawOpt1DSubMML->AddFrame(fChkMean, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkMean->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
            fChkMean->SetToolTipText("Activate Mean for Sigma/Cumulative/Integrate");

            fChkMedian = new TGCheckButton(fContDrawOpt1DSubMML, "Median");
            fContDrawOpt1DSubMML->AddFrame(fChkMedian, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkMedian->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
            fChkMedian->SetToolTipText("Activate Median for Sigma/Cumulative/Integrate");

            fChkLTM = new TGCheckButton(fContDrawOpt1DSubMML, "LTM");
            fContDrawOpt1DSubMML->AddFrame(fChkLTM, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkLTM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtons1D()");
            fChkLTM->SetToolTipText("Activate LTM for Sigma/Cumulative/Integrate");
            
      
      // statistic options container *** fcontStatOpt1D ***  " Statistic options "      
      fContStatOpt = new TGGroupFrame(ftabLeft1, "Statistic options", kVerticalFrame | kFitWidth | kFitHeight);
      ftabLeft1->AddFrame(fContStatOpt, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
      
         fChkStatName = new TGCheckButton(fContStatOpt, "Name");
         fContStatOpt->AddFrame(fChkStatName, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkStatName->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
         fChkStatName->SetToolTipText("Display the name in the statistics legend.");
      
         fChkStatEntries = new TGCheckButton(fContStatOpt, "Entries");
         fContStatOpt->AddFrame(fChkStatEntries, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkStatEntries->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
         fChkStatEntries->SetToolTipText("Display the number of entries in the statistics legend.");
      
         fContStatMean = new TGCompositeFrame(fContStatOpt, 1, 1, kHorizontalFrame | kFitWidth | kFitHeight);
         fContStatOpt->AddFrame(fContStatMean, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
            fChkStatMean = new TGCheckButton(fContStatMean, "Mean");
            fContStatMean->AddFrame(fChkStatMean, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkStatMean->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            fChkStatMean->SetToolTipText("Display the mean value of the data in the statistics legend.");
            
            fChkStatMeanPM = new TGCheckButton(fContStatMean, "+- Error");
            fContStatMean->AddFrame(fChkStatMeanPM, new TGLayoutHints(kLHintsNormal, 10, 0, 0, 0));
            fChkStatMeanPM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            fChkStatMeanPM->SetToolTipText("Display the mean value's error in the statistics legend.");

         fContStatRMS = new TGCompositeFrame(fContStatOpt, 1, 1, kHorizontalFrame | kFitWidth | kFitHeight);
         fContStatOpt->AddFrame(fContStatRMS, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
            fChkStatRMS = new TGCheckButton(fContStatRMS, "RMS");
            fContStatRMS->AddFrame(fChkStatRMS, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkStatRMS->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            fChkStatRMS->SetToolTipText("Display the RMS value of the data in the statistics legend.");
            
            fChkStatRMSPM = new TGCheckButton(fContStatRMS, "+- Error");
            fContStatRMS->AddFrame(fChkStatRMSPM, new TGLayoutHints(kLHintsNormal, 10, 0, 0, 0));
            fChkStatRMSPM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            fChkStatRMSPM->SetToolTipText("Display the RMS value's error in the statistics legend.");

         fChkStatUnderflow = new TGCheckButton(fContStatOpt, "Underflow");
         fContStatOpt->AddFrame(fChkStatUnderflow, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkStatUnderflow->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
         fChkStatUnderflow->SetToolTipText("Display the number of entries in the underflow bin.");
      
         fChkStatOverflow = new TGCheckButton(fContStatOpt, "Overflow");
         fContStatOpt->AddFrame(fChkStatOverflow, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkStatOverflow->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
         fChkStatOverflow->SetToolTipText("Display the number of entries in the overflow bin.");
      
         fChkStatIntegral = new TGCheckButton(fContStatOpt, "Integral");
         fContStatOpt->AddFrame(fChkStatIntegral, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkStatIntegral->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
         fChkStatIntegral->SetToolTipText("Display the integral of the data in the statistics legend.");
      
         fContStatSkew = new TGCompositeFrame(fContStatOpt, 1, 1, kHorizontalFrame | kFitWidth | kFitHeight);
         fContStatOpt->AddFrame(fContStatSkew, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
            fChkStatSkewness = new TGCheckButton(fContStatSkew, "Skewness");
            fContStatSkew->AddFrame(fChkStatSkewness, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkStatSkewness->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            fChkStatSkewness->SetToolTipText("Display the skewness of the data in the statistics legend. \nBe careful! Sometimes the skewness causes a floating point exception that hangs the GUI!");
            
            fChkStatSkewnessPM = new TGCheckButton(fContStatSkew, "+- Error");
            fContStatSkew->AddFrame(fChkStatSkewnessPM, new TGLayoutHints(kLHintsNormal, 10, 0, 0, 0));
            fChkStatSkewnessPM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            fChkStatSkewnessPM->SetToolTipText("Display the skewness' error in the statistics legend.");

         fContStatKurt = new TGCompositeFrame(fContStatOpt, 1, 1, kHorizontalFrame | kFitWidth | kFitHeight);
         fContStatOpt->AddFrame(fContStatKurt, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      
            fChkStatKurtosis = new TGCheckButton(fContStatKurt, "Kurtosis");
            fContStatKurt->AddFrame(fChkStatKurtosis, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkStatKurtosis->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            fChkStatKurtosis->SetToolTipText("Display the kurtosis of the data in the statistics legend.");
            
            fChkStatKurtosisPM = new TGCheckButton(fContStatKurt, "+- Error");
            fContStatKurt->AddFrame(fChkStatKurtosisPM, new TGLayoutHints(kLHintsNormal, 10, 0, 0, 0));
            fChkStatKurtosisPM->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsStat()");
            fChkStatKurtosisPM->SetToolTipText("Display the kurtosis' error in the statistics legend.");
       
      fBtnUnchekAll = new TGTextButton(fContStatOpt, "&Uncheck all");
      fContStatOpt->AddFrame(fBtnUnchekAll, new TGLayoutHints(kLHintsExpandX, 10, 10, 0, 0));
      //fBtnDraw->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoTest(=\"fBtnDraw clicked\")");
      fBtnUnchekAll->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "UnchekAllStat()");
      fBtnUnchekAll->SetToolTipText("Disable all statistics legend entries, \nno statistics legend.");

      
      // custom options container
      // --- fComboCustom --- the custom draw line on the very low
      fContCustom = new TGCompositeFrame(fContTopBottom, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
      fContTopBottom->AddFrame(fContCustom, new TGLayoutHints(kLHintsExpandX, 10, 0, 0, 0));
   
         // ------------------------- content of fContCustom -------------------------
         fLblCustomDraw = new TGLabel(fContCustom, "Custom draw: ");
         fLblCustomDraw->SetTextJustify(kTextLeft);
         fContCustom->AddFrame(fLblCustomDraw, new TGLayoutHints(kLHintsNormal, 5, 0, 0, 0));
         // text field for custom draw command
         fComboCustom = new TGComboBox(fContCustom);
         fComboCustom->Resize(0, fBtnDraw->GetDefaultHeight());
         fComboCustom->EnableTextInput(kTRUE);
         fContCustom->AddFrame(fComboCustom, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
         fComboCustom->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral(=42)");
         fComboCustom->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
      
   
      // additional cuts container
      fContAddCuts = new TGCompositeFrame(fContTopBottom, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
      fContTopBottom->AddFrame(fContAddCuts, new TGLayoutHints(kLHintsExpandX, 10, 0, 0, 0));
      
         // ------------------------- content of fContAddCuts -------------------------
         fLblAddCuts = new TGLabel(fContAddCuts, "Custom cuts:  ");
         fLblAddCuts->SetTextJustify(kTextLeft);
         fContAddCuts->AddFrame(fLblAddCuts, new TGLayoutHints(kLHintsNormal, 5, 0, 0, 0));
         // combo text field for additional cuts
         fComboAddCuts = new TGComboBox(fContAddCuts);
         fComboAddCuts->Resize(0, fBtnDraw->GetDefaultHeight());
         fComboAddCuts->EnableTextInput(kTRUE);
         fContAddCuts->AddFrame(fComboAddCuts, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
         fComboAddCuts->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsCuts(=31)");
         fComboAddCuts->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
         
   // ==========================================================================
   // ************************* content of fContCenter *************************
   // ========================================================================
   // main drawing canvas
   fCanvMain = new TRootEmbeddedCanvas("Main_Canvas", fContCenter, 200, 200, kFitWidth | kFitHeight);
   fContCenter->AddFrame(fCanvMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   fCanvMain->GetCanvas()->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)", "AliTPCCalibViewerGUI", this, "MouseMove(Int_t, Int_t, Int_t, TObject*)");
   fCanvMain->GetCanvas()->Connect("RangeAxisChanged()", "AliTPCCalibViewerGUI", this, "GetMinMax()");
   fCanvMain->GetCanvas()->SetToolTipText("The Main_Canvas, here your plots are displayed.");
   
   
   // =========================================================================   
   // ************************* content of fContRight *************************
   // ========================================================================
   
   // tabs on the right side:
   ftabRight = new TGTab(fContRight);
   fContRight->AddFrame(ftabRight, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 8, 0));
   fTabRight0 = ftabRight->AddTab("Basic");
   fTabRight1 = ftabRight->AddTab("Advanced");

   
      // **************************** content of tabLeft0 *******************************
      // cut options container
      fContCuts = new TGGroupFrame(fTabRight0, "Cuts", kVerticalFrame | kFitWidth | kFitHeight);
      fTabRight0->AddFrame(fContCuts, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   
      
         // ************************* content of fContCuts *************************
         // TPC radio button
         fRadioTPC = new TGRadioButton(fContCuts, "Whole TPC", 20);
         fContCuts->AddFrame(fRadioTPC, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fRadioTPC->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsCuts()");
         fRadioTPC->SetToolTipText("No cuts, use the whole TPC. \nIn 2D mode, A side and C side are superimposed!");
      
         // side A radio button
         fRadioSideA = new TGRadioButton(fContCuts, "Side A", 21);
         fContCuts->AddFrame(fRadioSideA, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fRadioSideA->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsCuts()");
         fRadioSideA->SetToolTipText("Use only side A.");
      
         // side C radio button
         fRadioSideC = new TGRadioButton(fContCuts, "Side C", 22);
         fContCuts->AddFrame(fRadioSideC, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fRadioSideC->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsCuts()");
         fRadioSideC->SetToolTipText("Use only side C.");
      
         // roc radio button
         fRadioROC = new TGRadioButton(fContCuts, "ROC", 23);
         fContCuts->AddFrame(fRadioROC, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fRadioROC->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsCuts()");
         fRadioROC->SetToolTipText("Use only one ROC (readout chamber).");

	 // sector radio button
         fRadioSector = new TGRadioButton(fContCuts, "Sector", 24);
         fContCuts->AddFrame(fRadioSector, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fRadioSector->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsCuts()");
         fRadioSector->SetToolTipText("Use only one Sector (0-17 A-Side, 18-35 C-Side).");
      
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
         
         // cut zeros check button
         fChkCutZero = new TGCheckButton(fContCuts, "Cut zeros");
         fContCuts->AddFrame(fChkCutZero, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkCutZero->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
         fChkCutZero->SetToolTipText("Use only values that are not zero. \nWhen doing a custom draw this cut is composed of the draw string\nup to the first ':' or '>>' and the condition that this should not be zero.");
      
         // additional cuts check button
         fChkAddCuts = new TGCheckButton(fContCuts, "Custom cuts");
         fContCuts->AddFrame(fChkAddCuts, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkAddCuts->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
         fChkAddCuts->SetToolTipText("Activate the custom cuts on the bottom. Remember '~' (= '.fElements')!");
         // fContAddCuts' content is locaed further op
      
      
      // Scaling options container
      fContScaling = new TGGroupFrame(fTabRight0, "Scaling", kVerticalFrame | kFitWidth | kFitHeight);
      fTabRight0->AddFrame(fContScaling, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   
         // ************************* content of fContScaling *************************
         // SetMaximum container
         fContSetMax = new TGCompositeFrame(fContScaling, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
         fContScaling->AddFrame(fContSetMax, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         
            // ------------------------- content of fContSetMax -------------------------
            // SetMaximum - checkbox
            fChkSetMax = new TGCheckButton(fContSetMax, "Set fixed max.");
            fContSetMax->AddFrame(fChkSetMax, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkSetMax->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw()");
            fChkSetMax->SetToolTipText("Set the maximum fixed to the value specified here.");
            
            // text field for maximum value
            fTxtSetMax = new TGTextEntry(fContSetMax, "", 41);
            fContSetMax->AddFrame(fTxtSetMax, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
            fTxtSetMax->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw()");
            fTxtSetMax->SetToolTipText("maximum value for the drawing");
      
         // SetMinimum container
         fContSetMin = new TGCompositeFrame(fContScaling, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
         fContScaling->AddFrame(fContSetMin, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         
            // ------------------------- content of fContSetMin -------------------------
            // SetMinimum - checkbox
            fChkSetMin = new TGCheckButton(fContSetMin, "Set fixed min.");
            fContSetMin->AddFrame(fChkSetMin, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkSetMin->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw()");
            fChkSetMin->SetToolTipText("Set the minimum fixed to the value specified here.");
            
            // text field for minimum value
            fTxtSetMin = new TGTextEntry(fContSetMin, "", 40);
            fContSetMin->AddFrame(fTxtSetMin, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
            fTxtSetMin->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw()");
            fTxtSetMin->SetToolTipText("minimum value for the drawing");
         
         // get Min & Max from Plot - button
         fBtnGetMinMax = new TGTextButton(fContScaling, "&Get scale from plot");
         fContScaling->AddFrame(fBtnGetMinMax, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fBtnGetMinMax->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "GetMinMax()");
         fBtnGetMinMax->SetToolTipText("Get min and max from plot, e.g. after rescaling by dragging the palette. \nObsolete! The button's function will change to 'Unzoom all'.");
         
         // GetMinMaxAuto - checkbox
         fChkGetMinMaxAuto = new TGCheckButton(fContScaling, "Get Min + Max auto.");
         fContScaling->AddFrame(fChkGetMinMaxAuto, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
         fChkGetMinMaxAuto->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw()");
         fChkGetMinMaxAuto->SetToolTipText("Get minimum and maximum automatically from each new plot. \nDeactivate this, if you want to 'save' your specified minimum and maximum.");
         
      // labeling container *** fContLabeling ***  " Labeling "      
         fContLabeling = new TGGroupFrame(fTabRight0, "Labeling", kVerticalFrame | kFitWidth | kFitHeight);
         fTabRight0->AddFrame(fContLabeling, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         
            fChkLabelTitle = new TGCheckButton(fContLabeling, "Set title:");
            fContLabeling->AddFrame(fChkLabelTitle, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkLabelTitle->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw()");
            fChkLabelTitle->SetToolTipText("Set the plot title.");
               
            fTxtLabelTitle = new TGTextEntry(fContLabeling, "Title", 500);
            fContLabeling->AddFrame(fTxtLabelTitle, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
            fTxtLabelTitle->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw(=50)");
            fTxtLabelTitle->SetToolTipText("plot title");
   
            fChkLabelXaxis = new TGCheckButton(fContLabeling, "Set X-axis label:");
            fContLabeling->AddFrame(fChkLabelXaxis, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkLabelXaxis->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw()");
            fChkLabelXaxis->SetToolTipText("Set the X-axis label.");
               
            fTxtLabelXaxis = new TGTextEntry(fContLabeling, "XaxisLabel", 500);
            fContLabeling->AddFrame(fTxtLabelXaxis, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
            fTxtLabelXaxis->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw(=51)");
            fTxtLabelXaxis->SetToolTipText("X-axis label");
   
            fChkLabelYaxis = new TGCheckButton(fContLabeling, "Set Y-axis label:");
            fContLabeling->AddFrame(fChkLabelYaxis, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkLabelYaxis->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw()");
            fChkLabelYaxis->SetToolTipText("Set the Y-axis label.");
               
            fTxtLabelYaxis = new TGTextEntry(fContLabeling, "YaxisLabel", 500);
            fContLabeling->AddFrame(fTxtLabelYaxis, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
            fTxtLabelYaxis->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw(=52)");
            fTxtLabelYaxis->SetToolTipText("Y-axis label");
   
            fChkLabelGetAuto = new TGCheckButton(fContLabeling, "Get labels auto.");
            fContLabeling->AddFrame(fChkLabelGetAuto, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
            fChkLabelGetAuto->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw()");
            fChkLabelGetAuto->SetToolTipText("Get labels automatically from each new plot \nDeactivate this, if you want to 'save' your specified labels.");

      
      // **************************** content of ftabRight1 *******************************
      // Save container
      fContSave = new TGGroupFrame(fTabRight1, "Save", kVerticalFrame | kFitWidth | kFitHeight);
      fTabRight1->AddFrame(fContSave, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         // save button
         fBtnSave = new TGTextButton(fContSave, "&Save picture");
         fContSave->AddFrame(fBtnSave, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fBtnSave->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "SavePicture()");
         fBtnSave->SetToolTipText("Open a 'Save as...' dialog to save the current plot as picture or macro.");

         // additional save options container
         fContAddSaveOpt = new TGCompositeFrame(fContSave, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
         fContSave->AddFrame(fContAddSaveOpt, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   
            //  content of --- fContAddSaveOpt ---
            // addition save options label
            fChkAddSaveOpt = new TGCheckButton(fContAddSaveOpt, "Save options:");
            fContAddSaveOpt->AddFrame(fChkAddSaveOpt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
            fChkAddSaveOpt->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
            fChkAddSaveOpt->SetToolTipText("Additional save options (see documentation for TPad::Print()).");
            
            // additional save options combo box
            fComboAddSaveOpt = new TGComboBox(fContAddSaveOpt);
            fContAddSaveOpt->AddFrame(fComboAddSaveOpt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
            fComboAddSaveOpt->Resize(0, fBtnDraw->GetDefaultHeight());
            fComboAddSaveOpt->EnableTextInput(kTRUE);
            fComboAddSaveOpt->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "SavePicture()");
            // fComboAddSaveOpt->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "SavePicture()");
     
      // calPad export container
      fContExport = new TGGroupFrame(fTabRight1, "Export AliTPCCalPad", kVerticalFrame | kFitWidth | kFitHeight);
      fTabRight1->AddFrame(fContExport, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   
         // ------------------------- content of fContExport -------------------------
         // container for export name
         fContAddExport = new TGCompositeFrame(fContExport, 200, 200, kVerticalFrame | kFitWidth | kFitHeight);
         fContExport->AddFrame(fContAddExport, new TGLayoutHints(kLHintsExpandX, -5, -5, 0, 0));
         
         fComboExportName = new TGComboBox(fContAddExport);
         fComboExportName->Resize(0, fBtnDraw->GetDefaultHeight());
         fContAddExport->AddFrame(fComboExportName, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
         fComboExportName->AddEntry("calPad",  0);  // first default value
         fComboExportName->Select(0);               // select default value before connecting
         fComboExportName->EnableTextInput(kTRUE);
         fComboExportName->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "DoExport()");
         fComboExportName->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoExport()");
  
         // export button
         fBtnExport = new TGTextButton(fContExport, "&Export to CINT");
         fContExport->AddFrame(fBtnExport, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fBtnExport->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoExport()");
         fBtnExport->SetToolTipText("Export the current 2D view as AliTPCCalPad to the CINT command line interpreter, use the specified name. \nThis works only in 2D mode.");
      
         // add to normalisation button
         fBtnAddNorm = new TGTextButton(fContExport, "&Add to normalization");
         fContExport->AddFrame(fBtnAddNorm, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fBtnAddNorm->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoExportNorm()");
         fBtnAddNorm->SetToolTipText("Use the current 2D view as normalization variable, use the specified name. \nNot yet working!");

      // Tree container
      fContTree = new TGGroupFrame(fTabRight1, "Tree", kVerticalFrame | kFitWidth | kFitHeight);
      fTabRight1->AddFrame(fContTree, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

         // dump tree to file button
         fBtnDumpToFile = new TGTextButton(fContTree, "&Dump to File");
         fContTree->AddFrame(fBtnDumpToFile, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fBtnDumpToFile->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoDumpToFile()");
         fBtnDumpToFile->SetToolTipText("Write the exported CalPads to a new CalibTree");

         // dump tree to file button
         fBtnLoadTree = new TGTextButton(fContTree, "&Load Tree");
         fContTree->AddFrame(fBtnLoadTree, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         fBtnLoadTree->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoLoadTree()");
         fBtnLoadTree->SetToolTipText("Load and initialize a new calibration tree. ");

         fChkAddAsReference = new TGCheckButton(fContTree, "as reference:");
         fContTree->AddFrame(fChkAddAsReference, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
         fChkAddAsReference->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoLoadTree()");
         fChkAddAsReference->SetToolTipText("To add a new tree as reference tree.");
            
         fTxtRefName = new TGTextEntry(fContTree, "R", 500);
         fContTree->AddFrame(fTxtRefName, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 15, 0, 0, 0));
         // fTxtRefName->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsNoRedraw(=50)");
         fTxtRefName->SetToolTipText("Reference Name");
        
            
      // Fit options container
      fContFit = new TGGroupFrame(fTabRight1, "Custom fit", kVerticalFrame | kFitWidth | kFitHeight);
      fTabRight1->AddFrame(fContFit, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   
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
         fBtnFit->SetToolTipText("Fit a whole TPC side, e.g. with gx~ ++ gy~, the result is printed to the console \nNot yet final status.");
      
         // add fit function button
         //fBtnAddFitFunction = new TGTextButton(fContAddFit, "&Add fit function to normalization");
         //fContAddFit->AddFrame(fBtnAddFitFunction, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
         //fBtnAddFitFunction->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "AddFitFunction()");

   SetWindowName("AliTPCCalibViewer GUI");
   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
}


void AliTPCCalibViewerGUI::SetInitialValues() {
   // 
   // Set the default button states
   // 
   fChkAuto->SetState(kButtonUp);
   fRadioPredefined->SetState(kButtonDown);
   fRadioRaw->SetState(kButtonDown);
   fRadioTPC->SetState(kButtonDown);
   fRadio1D->SetState(kButtonDown);
   fChkAddCuts->SetState(kButtonUp);
   fChkGetMinMaxAuto->SetState(kButtonDown);
   fChkSetMin->SetState(kButtonUp);
   fChkSetMax->SetState(kButtonUp);
   fRadioNorm->SetState(kButtonDown);
   fRadioSigma->SetState(kButtonUp);
   fRadioCumulative->SetState(kButtonUp);
   fChkMean->SetState(kButtonDown);
   fCheckCumulativePM->SetState(kButtonUp);
   
   fChkLabelGetAuto->SetState(kButtonDown);

   // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
   // TODO Set the checkboxes state as it is really     TODO 
   // TODO in gStyle                                    TODO 
   // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
//          k = 1;  kurtosis printed
//          k = 2;  kurtosis and kurtosis error printed
//          s = 1;  skewness printed
//          s = 2;  skewness and skewness error printed
//          i = 1;  integral of bins printed
//          o = 1;  number of overflows printed
//          u = 1;  number of underflows printed
//          r = 1;  rms printed
//          r = 2;  rms and rms error printed
//          m = 1;  mean value printed
//          m = 2;  mean and mean error values printed
//          e = 1;  number of entries printed
//          n = 1;  name of histogram is printed    
//          (default = 000001111)
         
   Int_t statOpt = gStyle->GetOptStat();
   if (statOpt == 1) statOpt = 1111;
   if (statOpt / 200000000 >= 1) {
      fChkStatKurtosis->SetState(kButtonDown);
      fChkStatKurtosisPM->SetState(kButtonDown);
      statOpt -= 200000000;
   }
   if (statOpt / 100000000 >= 1) {
      fChkStatKurtosis->SetState(kButtonDown);
      statOpt -= 100000000;
   }
   if (statOpt / 20000000 >= 1) {
      fChkStatSkewness->SetState(kButtonDown);
      fChkStatSkewnessPM->SetState(kButtonDown);
      statOpt -= 20000000;
   }
   if (statOpt / 10000000 >= 1) {
      fChkStatSkewness->SetState(kButtonDown);
      statOpt -= 10000000;
   }
   if (statOpt / 1000000 >= 1) {
      fChkStatIntegral->SetState(kButtonDown);
      statOpt -= 1000000;
   }
   if (statOpt / 100000 >= 1) {
      fChkStatOverflow->SetState(kButtonDown);
      statOpt -= 100000;
   }
   if (statOpt / 10000 >= 1) {
      fChkStatUnderflow->SetState(kButtonDown);
      statOpt -= 10000;
   }
   if (statOpt / 2000 >= 1) {
      fChkStatRMS->SetState(kButtonDown);
      fChkStatRMSPM->SetState(kButtonDown);
      statOpt -= 2000;
   }
   if (statOpt / 1000 >= 1) {
      fChkStatRMS->SetState(kButtonDown);
      statOpt -= 1000;
   }
   if (statOpt / 200 >= 1) {
      fChkStatMean->SetState(kButtonDown);
      fChkStatMeanPM->SetState(kButtonDown);
      statOpt -= 200;
   }
   if (statOpt / 100 >= 1) {
      fChkStatMean->SetState(kButtonDown);
      statOpt -= 100;
   }
   if (statOpt / 10 >= 1) {
      fChkStatEntries->SetState(kButtonDown);
      statOpt -= 10;
   }
   if (statOpt / 1 >= 1) {
      fChkStatName->SetState(kButtonDown);
      statOpt -= 1;
   }
   
  /* 
   fChkStatName->SetState(kButtonDown);
   fChkStatEntries->SetState(kButtonDown);
   fChkStatMean->SetState(kButtonDown);
   fChkStatRMS->SetState(kButtonDown);
   */
//    fChkStatMeanPM->SetState(kButtonUp);
//    fChkStatRMSPM->SetState(kButtonUp);
//    fChkStatUnderflow->SetState(kButtonUp);
//    fChkStatOverflow->SetState(kButtonUp);
//    fChkStatIntegral->SetState(kButtonUp);
//    fChkStatSkewness->SetState(kButtonUp);
//    fChkStatSkewnessPM->SetState(kButtonUp);
//    fChkStatKurtosis->SetState(kButtonDown);
//    fChkStatKurtosisPM->SetState(kButtonDown);
   
   // fill fComboAddDrawOpt with some additional drawing options
   fComboAddDrawOpt->AddEntry("same",      0);
   fComboAddDrawOpt->AddEntry("profbox",   1);
   fComboAddDrawOpt->AddEntry("profcolz",  2);
   fComboAddDrawOpt->AddEntry("profcont0", 3);
   fComboAddDrawOpt->AddEntry("proflego",  4);
   fComboAddDrawOpt->AddEntry("proflego2", 5);
   fComboAddDrawOpt->AddEntry("profsurf",  6);
   fComboAddDrawOpt->AddEntry("profsurf1", 7);
   fComboAddDrawOpt->AddEntry("profsurf2", 8);
   fComboAddDrawOpt->AddEntry("box",    9);
   fComboAddDrawOpt->AddEntry("colz",  10);
   fComboAddDrawOpt->AddEntry("cont0", 11);
   fComboAddDrawOpt->AddEntry("lego",  12);
   fComboAddDrawOpt->AddEntry("lego2", 13);
   fComboAddDrawOpt->AddEntry("surf",  14);
   fComboAddDrawOpt->AddEntry("surf1", 15);
   fComboAddDrawOpt->AddEntry("surf2", 16);

   // fill fComboAddSaveOpt with some additional drawing options
   fComboAddSaveOpt->AddEntry("Portrait",  0);
   fComboAddSaveOpt->AddEntry("Landscape", 1);
   fComboAddSaveOpt->AddEntry("Preview",   2);
   fComboAddSaveOpt->AddEntry("+50",       3);

   // fill fComboMethod
   fComboMethod->AddEntry("subtract",  0);
   fComboMethod->AddEntry("divide by", 1);
   
   // fill fComboExportName
//    fComboExportName->AddEntry("calPad",  0);
   // fComboExportName->AddEntry("calPad2", 1);
   fBtnExport->SetEnabled(kFALSE);
   fBtnAddNorm->SetEnabled(kFALSE);

   // select initial variables
   fListVariables->Select(0);
   fListNormalization->Select(0);
   fComboMethod->Select(0);
//    fComboExportName->Select(0);
//    fComboExportName->EnableTextInput(kTRUE);

   //fCanvMain->GetCanvas()->ToggleEventStatus(); // klappt nicht
   //fCanvMain->GetCanvas()->GetCanvasImp()->ShowStatusBar(kTRUE); // klappt auch nicht
   fListVariables->IntegralHeight(kFALSE);         // naja
   fListNormalization->IntegralHeight(kFALSE);     // naja
   fChkAuto->SetState(kButtonDown);
  
   // Make first drawing:
   // DoDraw();
}


AliTPCCalibViewerGUI::AliTPCCalibViewerGUI(const AliTPCCalibViewerGUI &c)
   : TGCompositeFrame(c.fParent, c.fWidth, c.fHeight),
    fViewer(0),
    fPreprocessor(0),
    fContTopBottom(0),
    fContLCR(0),
    fContLeft(0),
    ftabLeft(0),
    ftabLeft0(0),
    ftabLeft1(0),
    ftabRight(0),
    fTabRight0(0),
    fTabRight1(0),
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
    fRadioROC(0),
    fRadioSector(0),
    fComboAddDrawOpt(0),
    fChkAuto(0),
    fChkAutoAppend(0),
    fComboMethod(0),
    fListNormalization(0),
    fComboCustom(0),
    fLblCustomDraw(0),
    fChkAddDrawOpt(0),
    fNmbSector(0),
    fLblSector(0),
    fChkCutZero(0),
    fChkAddCuts(0),
    fLblAddCuts(0),
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
    fChkStatKurtosisPM(0),
    fBtnUnchekAll(0),
    fContLabeling(0),
    fChkLabelTitle(0),
    fTxtLabelTitle(0),
    fChkLabelXaxis(0),
    fTxtLabelXaxis(0),
    fChkLabelYaxis(0),
    fTxtLabelYaxis(0),
    fChkLabelGetAuto(0),
    fContSave(0),
    fBtnSave(0),
    fContAddSaveOpt(0),
    fChkAddSaveOpt(0),
    fComboAddSaveOpt(0),
    fContExport(0),
    fContAddExport(0),
    fComboExportName(0),
    fBtnExport(0),
    fBtnAddNorm(0), 
    fContTree(0),
    fBtnDumpToFile(0),
    fBtnLoadTree(0),
    fChkAddAsReference(0),
    fTxtRefName(0),
    fInitialized(0)
{
  //
  // dummy AliTPCCalibViewerGUI copy constructor
  //
}


AliTPCCalibViewerGUI & AliTPCCalibViewerGUI::operator =(const AliTPCCalibViewerGUI & /*param*/) {
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
   delete fPreprocessor;
}


/*
void AliTPCCalibViewerGUI::CloseWindow() {
   DeleteWindow();
}
*/


void AliTPCCalibViewerGUI::Initialize(const char* fileName, const char* treeName) {
   // 
   // initialize the GUI with a calibrationTree from fileName
   // 
   
   // create AliTPCCalibViewer object, which will be used for generating all drawings
   if (fViewer) delete fViewer;
   fViewer = new AliTPCCalibViewer(fileName, treeName);
   Initialize(fViewer);   
}


void AliTPCCalibViewerGUI::Initialize(AliTPCCalibViewer *viewer) {
   //
   // initializes the GUI with default settings and opens tree for drawing
   //
   
   fViewer = viewer;
   TString selectedVariable("");
   TString selectedNormalization("");
   Int_t variableId = -1;
   Int_t normalizationId = -1;
   if (fInitialized) {
      // remember the selected entry
      if (fListVariables->GetSelectedEntry()) selectedVariable = fListVariables->GetSelectedEntry()->GetTitle();
      if (fListNormalization->GetSelectedEntry()) selectedNormalization = fListNormalization->GetSelectedEntry()->GetTitle();
   }
   
   // fill fListVariables, list of drawable variables:
   TObjArray* arr = fViewer->GetListOfVariables();
   if (!arr) return;
   TIterator* iter = arr->MakeIterator();
   iter->Reset();
   TObjString* currentStr = 0;
   Int_t id = 0;
   fListVariables->RemoveAll();
   while ((currentStr = (TObjString*)(iter->Next()))) {
      fListVariables->AddEntry(currentStr->GetString().Data(), id);
      if (fInitialized && currentStr->GetString() == selectedVariable) variableId = id;
      id++;
   }

   // fill fListNorm, list of normalization variables:
   TObjArray *arrNorm = fViewer->GetListOfNormalizationVariables();
   TIterator *iterNorm = arrNorm->MakeIterator();
   iterNorm->Reset();
   currentStr = 0;
   id = 0;
   fListNormalization->RemoveAll();
   while ((currentStr = (TObjString*)(iterNorm->Next()))) {
      fListNormalization->AddEntry(currentStr->GetString().Data(), id);
      if (fInitialized && currentStr->GetString() == selectedNormalization) normalizationId = id;
      id++;
   }
   currentStr = 0;
   iter->Reset();
   //Add draw variables to the list of normalisation
   while ((currentStr = (TObjString*)(iter->Next()))) {
      if (currentStr->GetString().BeginsWith("Map")) continue; //don't add mapping information
      fListNormalization->AddEntry(currentStr->GetString().Data(), id);
      if (fInitialized && currentStr->GetString() == selectedNormalization) normalizationId = id;
      id++;
   }

   delete iterNorm;
   arrNorm->Delete();
   delete arrNorm;

   delete iter;
   arr->Delete();
   delete arr;
   
   // trick do display the entries corectly after reinitialization
   // otherwise all the entries would appear as one kryptic entry
   // resizing the listbox somehow fixes the problem...
   if (fInitialized) fListVariables->Resize(fListVariables->GetWidth()-1, fListVariables->GetHeight());
   if (fInitialized) fListVariables->Resize(fListVariables->GetWidth()+1, fListVariables->GetHeight());
   if (fInitialized) fListNormalization->Resize(fListNormalization->GetWidth()-1, fListNormalization->GetHeight());
   if (fInitialized) fListNormalization->Resize(fListNormalization->GetWidth()+1, fListNormalization->GetHeight());
   
   // select the last selected variable and normalization
   if (fInitialized && variableId != -1)     fListVariables->Select(variableId);
   if (fInitialized && normalizationId != -1)fListVariables->Select(normalizationId);
   
   if (fInitialized) Info("Initialize", "AliTPCCalibViewerGUI new initialized.");
   fInitialized = kTRUE;

}

void AliTPCCalibViewerGUI::Reset(){
  //
  // reset variables, delete calib viewer
  //
  if (fViewer) delete fViewer;
  fListVariables->RemoveAll();
  fListNormalization->RemoveAll();
  fInitialized = kFALSE;
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
      case 14:             // select Draw options fComboAddDrawOpt
         fChkAddDrawOpt->SetState(kButtonDown);
         break;
      case 13:             // fRadioPredefined
         fRadioCustom->SetState(kButtonUp);
         // fComboCustom->SetEnabled(kTRUE);
         //f RadioNormalized->SetState(kButtonUp);
         break;
      //--------
      case 30:             // fRadio1D
         fRadio2D->SetState(kButtonUp);
         fBtnExport->SetEnabled(kFALSE);
         fBtnAddNorm->SetEnabled(kFALSE);
         break;
      case 31:             // fRadio2D
         fRadio1D->SetState(kButtonUp);
         fBtnExport->SetEnabled(kTRUE);
         fBtnAddNorm->SetEnabled(kTRUE);
         break;
      case 42:             // fComboCustom
         fRadioCustom->SetState(kButtonDown);
         fRadioPredefined->SetState(kButtonUp);
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


void AliTPCCalibViewerGUI::HandleButtonsStat(Int_t id) {
   // 
   // handles statistic check boxes 
   // checks each checkbox if checked
   // if the checkbox is checked, appends 'n' for name, 'e' for entries, ...
   // to a TString, passes this TString to gStyle->SetOptStat(...)
   // 
   if (id == -1) {
      TGButton *btn = (TGButton *) gTQSender;
      id = btn->WidgetId();
   }
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


void AliTPCCalibViewerGUI::HandleButtonsCuts(Int_t id) {
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
         fRadioROC->SetState(kButtonUp);
         fRadioSector->SetState(kButtonUp);
         break;
      case 21:             // fRadioSideA
         fRadioTPC->SetState(kButtonUp);
         fRadioSideC->SetState(kButtonUp);
         fRadioROC->SetState(kButtonUp);
         fRadioSector->SetState(kButtonUp);
         break;
      case 22:             // fRadioSideC
         fRadioTPC->SetState(kButtonUp);
         fRadioSideA->SetState(kButtonUp);
         fRadioROC->SetState(kButtonUp);
         fRadioSector->SetState(kButtonUp);
         break;
      case 23:             // fRadioROC
         fRadioTPC->SetState(kButtonUp);
         fRadioSideA->SetState(kButtonUp);
         fRadioSideC->SetState(kButtonUp);
         fRadioSector->SetState(kButtonUp);
         break;
      case 24:             // fRadioSector
         fRadioTPC->SetState(kButtonUp);
         fRadioSideA->SetState(kButtonUp);
         fRadioSideC->SetState(kButtonUp);
         fRadioROC->SetState(kButtonUp);
         break;
      case 31:            // fComboAddCuts
         fChkAddCuts->SetState(kButtonDown);
         break;
   }
   DoNewSelection();
}


void AliTPCCalibViewerGUI::HandleButtonsNoRedraw(Int_t id) {
   //
   // handles label & scaling checkboxes 
   // without redrawing (not necessary, faster like this)
   //
    if (id == -1) {
      TGButton *btn = (TGButton *) gTQSender;
      id = btn->WidgetId();
   }

   switch (id) {
      case 40:             // fTxtSetMin
         fChkSetMin->SetState(kButtonDown);
         break;
      case 41:             // fTxtSetMax
         fChkSetMax->SetState(kButtonDown);
         break;
      case 50:             // fTxtLabelTitle
         fChkLabelTitle->SetState(kButtonDown);
         break;
      case 51:             // fTxtLabelXaxis
         fChkLabelXaxis->SetState(kButtonDown);
         break;
      case 52:             // fTxtLabelXaxis
         fChkLabelYaxis->SetState(kButtonDown);
         break;
   }
   SetMinMaxLabel();
}

void AliTPCCalibViewerGUI::ReplacePlaceHolders(TString &str)
{
    //
    // replace the defined placeholders in the custom draw string and cut string
    //
    TString drawPlaceHolder("#draw#");
    TString normPlaceHolder("#norm#");

    //current draw variable
    TString desiredData("");
    if (fListVariables->GetSelectedEntry()){
      desiredData += ((TGTextLBEntry*)(fListVariables->GetSelectedEntry()))->GetTitle();
      str.ReplaceAll(drawPlaceHolder,desiredData);
    }
//    desiredData += fViewer->GetAbbreviation();

    //current normalisation
    TString normalizationData("");
    if (fListNormalization->GetSelectedEntry()){
      normalizationData += ((TGTextLBEntry*)(fListNormalization->GetSelectedEntry()))->GetTitle();
      if (! (TString(((TGTextLBEntry*)(fListNormalization->GetSelectedEntry()))->GetTitle())).BeginsWith("Fit"))
        if ( normalizationData.BeginsWith("_") ) normalizationData = desiredData+normalizationData;
      if ( fListVariables->FindEntry(normalizationData.Data()) )
        normalizationData+="~";
      str.ReplaceAll(normPlaceHolder,normalizationData);
    }
}

void AliTPCCalibViewerGUI::DoNewSelection() {
   //
   // decides whether to redraw if user makes another selection
   //
   if (fChkAuto->GetState() == kButtonDown) DoDraw();
}


TString* AliTPCCalibViewerGUI::GetDrawString() {
   // 
   // create the draw string out of selection
   // 
   
   // specify data to plot
   TString desiredData("");
   if (!fListVariables->GetSelectedEntry()) return 0;
   desiredData += ((TGTextLBEntry*)(fListVariables->GetSelectedEntry()))->GetTitle();
   desiredData += fViewer->GetAbbreviation();

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
      if (!fListNormalization->GetSelectedEntry()) return 0;
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
         if (fRadioROC->GetState() == kButtonDown) {
            Int_t sector = (Int_t)(fNmbSector->GetNumber());
            cutStr += "sector==";
            cutStr += sector; 
         }
         if (fRadioSector->GetState() == kButtonDown) {
            Int_t sector = ((Int_t)(fNmbSector->GetNumber()))%36;
            cutStr += "sector%36==";
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
         formulaStr.ReplaceAll("~", fViewer->GetAbbreviation());
         normalizationData = *fViewer->Fit(desiredData.Data(), formulaStr.Data(), cutStr.Data(), chi2, fitParam, covMatrix);
      }  // if ( normalizationData.BeginsWith("Fit")

      desiredData += op;
      if (! (TString(((TGTextLBEntry*)(fListNormalization->GetSelectedEntry()))->GetTitle())).BeginsWith("Fit"))
         if ( normalizationData.BeginsWith("_") ) desiredData += ((TGTextLBEntry*)(fListVariables->GetSelectedEntry()))->GetTitle();
      if ( fListVariables->FindEntry(normalizationData.Data()) )
          normalizationData+="~";
      desiredData += normalizationData;
   }
   else if (fRadioCustom->GetState() == kButtonDown) {
      desiredData = fComboCustom->GetTextEntry()->GetText();
      if (desiredData == "") return 0;
      ReplacePlaceHolders(desiredData);
   }
   
   // try to add forgotten '~'
   if (fChkAutoAppend->GetState() == kButtonDown) 
      desiredData = TString(fViewer->AddAbbreviations((char*)desiredData.Data()));
   return new TString(desiredData.Data());
}   

  
TString* AliTPCCalibViewerGUI::GetSectorString() {
   // 
   // create the sector string out of selection
   // 

   TString sectorStr("");
   if (fRadioTPC->GetState() == kButtonDown)
      sectorStr += "ALL";
   if (fRadioSideA->GetState() == kButtonDown)
      sectorStr += "A"; //cuts += "(sector/18)%2==0";
   if (fRadioSideC->GetState() == kButtonDown)
      sectorStr+= "C"; //cuts += "(sector/18)%2==1";
   if (fRadioROC->GetState() == kButtonDown) {
      Int_t sector = (Int_t)(fNmbSector->GetNumber());
      sectorStr += sector; //cuts += "sector==";
   }
   if (fRadioSector->GetState() == kButtonDown) {
      Int_t sector = ((Int_t)(fNmbSector->GetNumber()))%36;
      sectorStr += "S";
      sectorStr += sector; //cuts += "sector%36==";
   }
   return new TString(sectorStr.Data());
}   
  
 
 TString* AliTPCCalibViewerGUI::GetCutString() {
   // 
   // create the cut string out of selection
   // 
  
   TString cutsStr("");
   if (fChkCutZero->GetState() == kButtonDown) {
      TString cutZerosStr(GetDrawString()->Data());
      if (cutZerosStr.Contains(">>")) {
         cutZerosStr.Remove(cutZerosStr.First(">>"));
      }
      if (cutZerosStr.Contains(":")) {
         cutZerosStr.Remove(cutZerosStr.First(":"));
      }
      cutsStr += cutZerosStr.Data();
      cutsStr += "!=0";
      if (fChkAddCuts->GetState() == kButtonDown) cutsStr += " && ";
   }
   if (fChkAddCuts->GetState() == kButtonDown){
       cutsStr += fComboAddCuts->GetTextEntry()->GetText();
       ReplacePlaceHolders(cutsStr);
   }

   
   // try to add forgotten '~'
   if (fChkAutoAppend->GetState() == kButtonDown) 
      cutsStr = TString(fViewer->AddAbbreviations((char*)cutsStr.Data()));
   return new TString(cutsStr.Data());
}


void AliTPCCalibViewerGUI::DoDraw() {
   //
   // main method for drawing according to user selection
   //
   
   // specify data to plot:
  if (!GetDrawString()) return;
   TString desiredData(GetDrawString()->Data());
   // specify sector:
   TString sectorStr(GetSectorString()->Data());
   // specify cuts:
   TString cutsStr(GetCutString()->Data());

   TString addDrawOpt("");
   if (fChkAddDrawOpt->GetState() == kButtonDown)
      addDrawOpt += fComboAddDrawOpt->GetTextEntry()->GetText();
   
   // remove last picture
   if (!addDrawOpt.Contains("same"))
      for (Int_t i = 0; i < fCanvMain->GetCanvas()->GetListOfPrimitives()->GetEntries(); i++) {
         if (strcmp(fCanvMain->GetCanvas()->GetListOfPrimitives()->At(i)->ClassName(), "TFrame") != 0)
            fCanvMain->GetCanvas()->GetListOfPrimitives()->At(i)->Delete();
      }
   //fCanvMain->GetCanvas()->Clear();
   fCanvMain->GetCanvas()->cd();
   Int_t entries = -1;
   // draw finally
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
   
   SetMinMaxLabel();
   fCanvMain->GetCanvas()->Update();
}


void AliTPCCalibViewerGUI::DoFit() {
   //
   // main method for fitting
   //
   
   Double_t chi2 = 0;
   TVectorD fitParam(0);
   TMatrixD covMatrix(0,0);
   TString cutStr("");
   TString formulaStr("");

   // specify data to plot:
   TString drawStr(GetDrawString()->Data());
   
   // ********** create cut string **********
   if (fRadioTPC->GetState() == kButtonDown)
      cutStr += ""; // whole TPC is used for fitting
   if (fRadioSideA->GetState() == kButtonDown)
      cutStr += "(sector/18)%2==0"; // side A
   if (fRadioSideC->GetState() == kButtonDown)
      cutStr+= "(sector/18)%2==1"; // side C
   if (fRadioROC->GetState() == kButtonDown) {
      Int_t sector = (Int_t)(fNmbSector->GetNumber());
      cutStr += "sector==";
      cutStr += sector; 
   }
   if (fRadioSector->GetState() == kButtonDown) {
      Int_t sector = (Int_t)(fNmbSector->GetNumber())%36;
      cutStr += "sector%36==";
      cutStr += sector; 
   }
   if (fChkAddCuts->GetState() == kButtonDown && strcmp(fComboAddCuts->GetTextEntry()->GetText(), "") != 0){
      if (fRadioTPC->GetState() != kButtonDown) cutStr += " && ";
      cutStr += fComboAddCuts->GetTextEntry()->GetText();  
   }
   // try to add forgotten '~'
   if (fChkAutoAppend->GetState() == kButtonDown) 
      cutStr = TString(fViewer->AddAbbreviations((char*)cutStr.Data()));
   
   // ********** get formula string **********
   formulaStr += fComboCustomFit->GetTextEntry()->GetText();
   if (fChkAutoAppend->GetState() == kButtonDown) 
      formulaStr = TString(fViewer->AddAbbreviations((char*)formulaStr.Data()));

   // ********** call AliTPCCalibViewer's fit-function
  ReplacePlaceHolders(drawStr);
  ReplacePlaceHolders(cutStr);
  TString *returnStr = fViewer->Fit(drawStr.Data(), formulaStr.Data(), cutStr.Data(), chi2, fitParam, covMatrix);
   
   std::cout << std::endl;
   std::cout << "Your fit formula reads as follows:" << std::endl;
   std::cout << returnStr->Data() << std::endl;
   std::cout << "chi2 = " << chi2 << std::endl;

  delete returnStr;
}


void AliTPCCalibViewerGUI::DoExport() {
   //
   // function to export a CalPad to Cint
   //
   if ( fRadio2D->GetState() != kButtonDown){
      Error("ExportCalPad", "Export of AliTPCCalPad to CINT works only in 2D mode.");
      return;
   }
   // specify data to plot:
   TString desiredData(GetDrawString()->Data());
   // specify cuts:
   TString cutsStr(GetCutString()->Data());
   // get name for the calPad   
   const char* calPadName = fComboExportName->GetTextEntry()->GetText();
   // create calPad according to drawCommand and cuts
   AliTPCCalPad *calPad = fViewer->GetCalPad(desiredData.Data(), (char*)cutsStr.Data(), (char*)calPadName);
   // finally export calPad to Cint:
   gROOT->ProcessLine(Form("AliTPCCalPad* %s = (AliTPCCalPad*)0x%lx;", calPadName, (ULong_t)calPad));
   fPreprocessor->AddComponent(calPad);
   Info("ExportCalPad", "Current 2D view has been exported to an AliTPCCalPad* with name '%s'", calPadName);
}


void AliTPCCalibViewerGUI::DoExportNorm() {
   //
   // function to export a CalPad to Cint
   //
   
   if ( fRadio2D->GetState() != kButtonDown){
      Error("ExportCalPad", "Adding an AliTPCCalPad to the normalization works only in 2D mode.");
      return;
   }
   
   Error("DoExportNorm", "Not yet implemented.");
/*  
   return;
   
   // specify data to plot:
   TString desiredData(GetDrawString()->Data());
   // specify sector:
   TString sectorStr(GetSectorString()->Data());
   // specify cuts:
   TString cutsStr(GetCutString()->Data());
   
   // get name for the calPad   
   const char* calPadName = fComboExportName->GetTextEntry()->GetText();
   // create calPad according to drawCommand and cuts
   AliTPCCalPad *calPad = fViewer->GetCalPad(desiredData.Data(), (char*)cutsStr.Data(), (char*)fComboExportName->GetTextEntry()->GetText());
   // finally export calPad to Cint:
   gROOT->ProcessLine(Form("AliTPCCalPad* %s = (AliTPCCalPad*)0x%lx;", calPadName, (ULong_t) calPad));
   Info("ExportCalPad", "Current 2D view has been exported to an AliTPCCalPad* with name '%s'", calPadName);
  */
}


void AliTPCCalibViewerGUI::GetMinMax() {
   //
   // Read current Min & Max from the plot and set it to fTxtSetMin & fTxtSetMax
   //
   if (fChkGetMinMaxAuto->GetState() == kButtonUp) return;
   TList* listOfPrimitives = fCanvMain->GetCanvas()->GetListOfPrimitives();
   TObject* ptr = 0;
   for (Int_t i = 0; i < listOfPrimitives->GetEntries(); i++) {
      ptr = listOfPrimitives->At(i);
      if ( ptr->InheritsFrom("TH1") ) break;
   }
   if ( !ptr || !ptr->InheritsFrom("TH1") ) return;      // if the loop did not find a TH1
   TH1 *hist = (TH1*)ptr;

//    Double_t histMax = hist->GetMaximum();
//    Double_t histMin = hist->GetMinimum();
//    fTxtSetMax->SetText(Form("%f",histMax));
//    fTxtSetMin->SetText(Form("%f",histMin));

   if (fRadio2D->GetState() == kButtonDown) {
      if (fChkSetMax->GetState() == kButtonUp)
         fTxtSetMax->SetText(Form("%f", hist->GetMaximum()));
      if (fChkSetMin->GetState() == kButtonUp)
         fTxtSetMin->SetText(Form("%f", hist->GetMinimum()));
   }
   else if (fRadio1D->GetState() == kButtonDown) {
      if (fChkSetMax->GetState() == kButtonUp)
         fTxtSetMax->SetText( Form("%f", hist->GetXaxis()->GetXmax()) );
      if (fChkSetMin->GetState() == kButtonUp)
         fTxtSetMin->SetText( Form("%f", hist->GetXaxis()->GetXmin()) );
   }
}


void AliTPCCalibViewerGUI::SetMinMaxLabel() {
   // 
   // Set Minimum, Maximum and labels without redrawing the plot
   // (faster)
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
      Warning("SetMinMaxLabel","No Histogram found!");
      return;
      // unable to find histogram, no min and max wil be read out
   }
   
   TH1 *hist = (TH1*)ptr; 
   TString minTxt(fTxtSetMin->GetText());
   TString maxTxt(fTxtSetMax->GetText());
   
   // set min and max according to specified values, if checkbox is checked
   if (fRadio2D->GetState() == kButtonDown) {
      if (fChkSetMax->GetState() == kButtonDown && fChkSetMax->GetState() == kButtonDown &&(maxTxt.IsDigit() || maxTxt.IsFloat()) )
         hist->SetMaximum(maxTxt.Atof());
      if (fChkSetMax->GetState() == kButtonUp)
         hist->SetMaximum(-1111);  // default value, to unzoom
      if (fChkSetMin->GetState() == kButtonDown && (minTxt.IsDigit() || minTxt.IsFloat()) )
         hist->SetMinimum(minTxt.Atof());
      if (fChkSetMin->GetState() == kButtonUp)
         hist->SetMinimum(-1111);  // default value, to unzoom
   }
   else if (fRadio2D->GetState() == kButtonDown) {
      if (fChkSetMin->GetState() == kButtonDown && 
          fChkSetMax->GetState() == kButtonDown && hist->GetXaxis())
         hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
      else if (fChkSetMax->GetState() == kButtonDown && hist->GetXaxis())
         hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(), maxTxt.Atof());
      else if (fChkSetMin->GetState() == kButtonDown && hist->GetXaxis())
         hist->GetXaxis()->SetRangeUser(minTxt.Atof(), hist->GetXaxis()->GetXmax());
      hist->SetTitle(hist->GetTitle());  // trick to update the histogram
   }
   
   // get min and max from plot       
   GetMinMax();
   
   // set labels according to specification, if cehckboxes are checked
   if (fChkLabelTitle->GetState() == kButtonDown) 
      hist->SetTitle(fTxtLabelTitle->GetText());
   if (fChkLabelXaxis->GetState() == kButtonDown)
      hist->GetXaxis()->SetTitle(fTxtLabelXaxis->GetText());
   if (fChkLabelYaxis->GetState() == kButtonDown)
      hist->GetYaxis()->SetTitle(fTxtLabelYaxis->GetText());
   // get and/or set labels and title
   if (fChkLabelGetAuto->GetState() == kButtonDown) {
      fTxtLabelTitle->SetText(hist->GetTitle());
      fTxtLabelXaxis->SetTitle(hist->GetXaxis()->GetTitle());
      fTxtLabelYaxis->SetTitle(hist->GetYaxis()->GetTitle());
   }
   hist->SetTitle(hist->GetTitle());  // trick to update the histogram
   fCanvMain->GetCanvas()->Update();
}
   
   
void AliTPCCalibViewerGUI::ChangeSector(){
   // 
   // function that is called, when the number of the sector is changed
   // to change the sector label
   //
   if ( fRadioROC->GetState()!=kButtonDown && fRadioSector->GetState()!=kButtonDown ){
       fLblSector->SetText("not used");
       return;
   }

   Int_t sector = (Int_t)(fNmbSector->GetNumber());
   TString secLabel = "";
   if ( sector < 36 )
       secLabel = "IROC";
   else
       secLabel = "OROC";

   if (fRadioSector->GetState()==kButtonDown)
       secLabel="Sector";

   if ( sector%36<18 ) //A-Side
       secLabel += ", A";
   else
       secLabel += ", C";

   secLabel += Form("%02d",sector%18);

   fLblSector->SetText(secLabel);
   DoNewSelection();
}


void AliTPCCalibViewerGUI::AddFitFunction() const { 
   //
   // adds the last fit function to the normalization list
   // 
   std::cout << "Not yet implemented." << std::endl;
}


void AliTPCCalibViewerGUI::UnchekAllStat() {
   // 
   // Disable all statistical legend entries, no statistical legend.
   // 
   fChkStatName->SetState(kButtonUp);
   fChkStatEntries->SetState(kButtonUp);
   fChkStatMean->SetState(kButtonUp);
   fChkStatMeanPM->SetState(kButtonUp);
   fChkStatRMS->SetState(kButtonUp);
   fChkStatRMSPM->SetState(kButtonUp);
   fChkStatUnderflow->SetState(kButtonUp);
   fChkStatOverflow->SetState(kButtonUp);
   fChkStatIntegral->SetState(kButtonUp);
   fChkStatSkewness->SetState(kButtonUp);
   fChkStatSkewnessPM->SetState(kButtonUp);
   fChkStatKurtosis->SetState(kButtonUp);
   fChkStatKurtosisPM->SetState(kButtonUp);
   
   HandleButtonsStat(0);
}


void AliTPCCalibViewerGUI::MouseMove(Int_t event, Int_t x, Int_t y, TObject *selectedObject) { 
   //
   // mouse move
   // zoom to sector works ONLY in 2D mode, if one side is specified
   // 
   Double_t pi = TMath::Pi();
   if (event != kButton1Double )
      return;
   if (!selectedObject->InheritsFrom("TH2")) return;
   // zoom to sector works ONLY in 2D mode, if one side is specified
   if (fRadio2D->GetState() == kButtonUp) return;
   if (fRadioROC->GetState() == kButtonDown) { // return to full side view
      // return to full side view
      Int_t sector = (Int_t)(fNmbSector->GetNumber());
      if ( (sector >= 0 && sector <= 17) || (sector >= 36 &&  sector <= 53) ) {
         // A-Side
         fRadioSideA->Clicked();
         fRadioSideA->SetState(kButtonDown);
         DoNewSelection();
      }
      if ( (sector >= 18 && sector <= 35) || (sector >= 54 &&  sector <= 71) ) {
         // C-Side
         fRadioSideC->Clicked();
         fRadioSideC->SetState(kButtonDown);
         DoNewSelection();
      }
      return;
   }
   if (!(fRadioSideA->GetState() == kButtonDown || fRadioSideC->GetState() == kButtonDown)) return;
   
   // Int_t    px        = gPad->GetEventX();
   // Int_t    py        = gPad->GetEventY();
   Float_t  upy       = gPad->AbsPixeltoY(y);
   Float_t  upx       = gPad->AbsPixeltoX(x);
   Float_t  gy         = gPad->PadtoY(upy);
   Float_t  gx         = gPad->PadtoX(upx);
   Int_t quadrant = -1;
   if (gx >= 0 && gy >= 0) quadrant = 1;
   if (gx <  0 && gy >= 0) quadrant = 2;
   if (gx <  0 && gy <  0) quadrant = 3;
   if (gx >= 0 && gy <  0) quadrant = 4;
   gx = TMath::Abs(gx);
   gy = TMath::Abs(gy);
   Double_t phi = TMath::ATan(gy/gx);        // angle phi is in Pi- units
   Double_t r   = TMath::Sqrt(gx*gx + gy*gy);
   if (quadrant == 2) phi = pi - phi;
   if (quadrant == 3) phi = pi + phi;
   if (quadrant == 4) phi = 2 * pi - phi;
   Double_t phiGrad = phi / pi * 180;
   Int_t sector = (Int_t) phiGrad / 20;  // one sector coresponds to 20�
   // IROC starts at 84.5 cm
   // IROC ends at 135.5 cm, OROC begins
   // OROC ends at 250 cm
   if (r < 84.5 || r > 250) return; // outside TPC
   if (r < 135.5) { // IROC 
      if (fRadioSideC->GetState() == kButtonDown) sector += 18;
   }
   else {// OROC
      sector += 36;
      if (fRadioSideC->GetState() == kButtonDown) sector += 18;
   }
   // printf("r: %f, phi: %f, phiGrad: %f, gy/gx: %f, quadrant: %i, sector: %i \n", r, phi, phiGrad, gy/gx, quadrant, sector);
   fNmbSector->SetNumber(sector);
   fRadioROC->Clicked();
   fRadioROC->SetState(kButtonDown);
   ChangeSector();   
}


void AliTPCCalibViewerGUI::SavePicture() {
   // 
   // saves the current picture
   // 
   // use the following combination of file type and save options:
   // (see also TCanvas::Print)
   // 
   //       "ps"  - Postscript file is produced (see special cases below)
   //    "Portrait" - Postscript file is produced (Portrait)
   // "Landscape" - Postscript file is produced (Landscape)
   //       "eps" - an Encapsulated Postscript file is produced
   //    "Preview" - an Encapsulated Postscript file with preview is produced.
   //       "pdf" - a PDF file is produced
   //       "svg" - a SVG file is produced
   //       "gif" - a GIF file is produced
   //       "gif+NN" - an animated GIF file is produced, where NN is delay in 10ms units
   //       "xpm" - a XPM file is produced
   //       "png" - a PNG file is produced
   //       "jpg" - a JPEG file is produced
   //       "tiff" - a TIFF file is produced
   //       "cxx" - a C++ macro file is produced
   //       "xml" - a XML file
   //       "root" - a ROOT binary file
   
   const char *kSaveAsTypes[] = {
      "Postscript",  "*.ps",
      "Encapsulated Postscript",   "*.eps",
      "PDF",   "*.pdf",
      "JPEG",   "*.jpg",
      "PNG",   "*.png",
      "TIFF",   "*.tiff",
      "GIF",   "*.gif",
      "XPM",   "*.xpm",
      "SVG",   "*.svg",
      "XML",   "*.xml",
      "C++ macro",   "*.cxx",
      "Macro file",  "*.C",
      "ROOT file",   "*.root",
      "All file",    "*",
       0,              0
   };
   TString addSaveOpt("");
   if (fChkAddSaveOpt->GetState() == kButtonDown)
   addSaveOpt += fComboAddSaveOpt->GetTextEntry()->GetText();
   TString dir(".");
   TGFileInfo fi;
   fi.fFileTypes = kSaveAsTypes;
   // fi.fIniDir    = StrDup(dir);
   fi.fOverwrite = kFALSE;
   new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDSave, &fi);
   if (fi.fFilename && strlen(fi.fFilename)) {
      if (addSaveOpt != "")
         fCanvMain->GetCanvas()->Print(fi.fFilename, addSaveOpt.Data());
      else 
         fCanvMain->GetCanvas()->Print(fi.fFilename);
   }

// TList*	fFileNamesList	list of selected file names
// Int_t	fFileTypeIdx	selected file type, index in fFileTypes
// const char**	fFileTypes	file types used to filter selectable files
// char*	fFilename	selected file name
// char*	fIniDir	on input: initial directory, on output: new directory
// Bool_t	fMultipleSelection	if true, allow multiple file selection
// Bool_t	fOverwrite	if true overwrite the file with existing name on save
   
}
   

void AliTPCCalibViewerGUI::DoDumpToFile() {
   // 
   // This function is called, when the "Dump to File" button is pressed. 
   // All the exported CalPads will be written into an new CalibTree, 
   // a Save File dialog will appear to specify the filename
   // 
   const char *kSaveAsTypes[] = {
      "ROOT file",   "*.root",
       0,              0
   };
   TString dir(".");
   TGFileInfo fi;
   fi.fFileTypes = kSaveAsTypes;
   // fi.fIniDir    = StrDup(dir);
   fi.fOverwrite = kFALSE;
   new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDSave, &fi);
   if (fi.fFilename && strlen(fi.fFilename)) {
      fPreprocessor->DumpToFile(fi.fFilename);
      Info("DumpToFile", Form("New CalibTree has been writen to file '%s'", fi.fFilename));
   }
}


void AliTPCCalibViewerGUI::DoLoadTree() {
   // function to load a new calib tree
   // 
   // 
   const char *kFileTypes[] = {
      "ROOT file",   "*.root",
       0,              0
   };
   TString dir(".");
   TGFileInfo fi;
   fi.fFileTypes = kFileTypes;
   // fi.fIniDir    = StrDup(dir);
   fi.fOverwrite = kFALSE;
   new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen, &fi);
   if (fi.fFilename && strlen(fi.fFilename) && fChkAddAsReference->GetState() == kButtonUp) {
      Initialize(fi.fFilename);
      Reload();
   }
   else if (fi.fFilename && strlen(fi.fFilename) && fChkAddAsReference->GetState() == kButtonDown) {
      fViewer->AddReferenceTree(fi.fFilename, "calPads", fTxtRefName->GetText());
      Reload();
   }
}


TObjArray* AliTPCCalibViewerGUI::ShowGUI(const char* fileName) {
   //
   // Initialize and show GUI for presentation for demonstration purposes
   // or for fast standalone use
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
   
   TObjArray *guiArray = new TObjArray();
   guiArray->Add(calibViewer1);
   guiArray->Add(calibViewer2);
   
   frmMain->MapSubwindows();
   frmMain->Resize();
   frmMain->MapWindow();
   
   return guiArray;
}

