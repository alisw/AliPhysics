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
//  GUI for displaying calibration entries over time                         //
//                                                                           //
//  Authors:     Marian Ivanov (Marian.Ivanov@cern.ch)                       //
//               Jens Wiechula (Jens.Wiechula@cern.ch)                       //
//               Ionut Arsene  (iarsene@cern.ch)                             //
//                                                                           //
//  Usage:   AliCalibViewerGUItime::ShowGUI(TChain* chain)                   //
//           AliCalibViewerGUItime::ShowGUI()                                //
///////////////////////////////////////////////////////////////////////////////


#include <iostream>
//Root includes
#include <TROOT.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TVector.h>
#include <TH1.h>
#include <TCut.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TIterator.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TTimeStamp.h>
#include <TMath.h>
#include <TMap.h>
//
#include <TGFileDialog.h>
#include <TGInputDialog.h>
//
#include <TGButton.h>
#include <TGListBox.h>
#include <TGComboBox.h>
#include <TGNumberEntry.h>
#include <TGLayout.h>
#include <TRootEmbeddedCanvas.h>
#include <TGSplitter.h>
#include <TGButtonGroup.h>
#include <TGLabel.h>
#include <TGTab.h>
#include <TGString.h>

//AliRoot includes
#include <AliLog.h>
#include "AliBaseCalibViewer.h"
#include "AliBaseCalibViewerGUI.h"
#include "AliCalibViewerGUItime.h"


ClassImp(AliCalibViewerGUItime)

//_______________________________________________________________________________________________________________
AliCalibViewerGUItime::AliCalibViewerGUItime(const TGWindow *p, UInt_t w, UInt_t h, const Char_t* det) :
TGCompositeFrame(p,w,h),
  fDetector(det),
  fFile(0x0),
  fTree(0x0),
  fCalibViewerGUI(0x0),
  fCalibViewerGUItab(0x0),
  fCurrentHist(0x0),
  fCurrentGraph(0x0),
  fCurrentRunDetails(-1),
  fOutputCacheDir("/tmp"),
  fDrawString(""),
  fIsCustomDraw(kFALSE),
  fRunNumbers(10),
  fTimeStamps(10),
  fValuesX(10),
  fValuesY(10),
  fNoGraph(kFALSE),
  fGraphLimitEntries(10000),
  fMapRefTrees(new TMap),
  //GUI elements
  //main canvas Top part, bottom part
  fContTopBottom(0x0),
  //top left, centre, right
  fContLCR(0x0),
  //content left
  fContLeft(0x0),
  fContDrawOpt(0x0),
  fChkDrawOptSame(0x0),
  fComboAddDrawOpt(0x0),
  fContDrawSel(0x0),
  fContDrawSelSubRunTime(0x0),
  fRadioXhist(0x0),
  fRadioXrun(0x0),
  fRadioXtime(0x0),
  fListVariables(0x0),
  fComboRunType(0x0),
  fLblRunType(0x0),
  fNmbPar(0x0),
  fLblPar(0x0),
  fListCalibType(0x0),
  fContCalibType(0x0),
  //content centre
  fContCenter(0x0),
  fCanvMain(0x0),
  //content right
  fContRight(0x0),
  fContValues(0x0),
  fLblRunNumber(0x0),
  fLblRunTime(0x0),
  fLblValueX(0x0),
  fLblValueY(0x0),
  fLblRunNumberVal(0x0),
  fLblRunTimeVal(0x0),
  fLblValueXVal(0x0),
  fLblValueYVal(0x0),
  fBtnDumpRuns(0x0),
  fBtnSave(0x0),
  //load a file
  fContLoad(0x0),
  fContFilename(0x0),
  fContConfigFile(0x0),
  fContTreeName(0x0),
  fLblTreeName(0x0),
  fTxtFilename(0x0),
  fTxtConfigFile(0x0),
  fTxtTreeName(0x0),
  fBtnLoadFile(0x0),
  // extract information from OCDB
  //content bottom
  fContCustom(0x0),
  fContCustomCuts(0x0),
  fLblCustomDraw(0x0),
  fLblCustomCuts(0x0),
  fComboCustomDraw(0x0),
  fComboCustomCuts(0x0),
  fTrashBox(new TObjArray)
{
  //
  // ctor
  //
  fDetector.ToUpper();
  fMapRefTrees->SetOwnerKeyValue();
  fTrashBox->SetOwner();
  DrawGUI(p,w,h);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.5);
  SetInitialValues();
}

//______________________________________________________________________________
AliCalibViewerGUItime::~AliCalibViewerGUItime(){
  //
  // dtor
  //
  //  delete fConfigParser;
  delete fTrashBox;
  delete fMapRefTrees;
}

//______________________________________________________________________________
void AliCalibViewerGUItime::DrawGUI(const TGWindow */*p*/, UInt_t w, UInt_t h) {
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
  fContRight = new TGCompositeFrame(fContLCR, 200, 200, kVerticalFrame | kFixedWidth | kFitHeight);
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
  // draw options
  fContDrawOpt = new TGGroupFrame(fContLeft, "Draw options", kVerticalFrame | kFitWidth | kFitHeight);
  fContLeft->AddFrame(fContDrawOpt, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
  fChkDrawOptSame = new TGCheckButton(fContDrawOpt, "Same");
  fContDrawOpt->AddFrame(fChkDrawOptSame, new TGLayoutHints(kLHintsNormal, 0, 2, 0, 0));
  fChkDrawOptSame->SetToolTipText("Add draw option 'same'");
  // additional draw options combo box
  fComboAddDrawOpt = new TGComboBox(fContDrawOpt);
  fComboAddDrawOpt->Resize(0, 22);
  fComboAddDrawOpt->EnableTextInput(kTRUE);
  fContDrawOpt->AddFrame(fComboAddDrawOpt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
    
  // draw selection group
  fContDrawSel = new TGGroupFrame(fContLeft, "Draw selection", kVerticalFrame | kFitWidth | kFitHeight);
  fContLeft->AddFrame(fContDrawSel, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 10, 0));
  //x-axis variables selection, Run of Time
  fContDrawSelSubRunTime = new TGCompositeFrame(fContDrawSel, 200, 23, kHorizontalFrame | kFitWidth | kFixedHeight);
  fContDrawSel->AddFrame(fContDrawSelSubRunTime, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX , 0, 0, 0, 0));
  
  // ------------------------- content of fContDrawOpt -------------------------
  // 
  // Run radio button
    // Time radio button
  fRadioXhist = new TGRadioButton(fContDrawSelSubRunTime, "hist", kRadioXhist);
  fContDrawSelSubRunTime->AddFrame(fRadioXhist, new TGLayoutHints(kLHintsNormal, 0, 2, 0, 0));
  fRadioXhist->Connect("Clicked()", "AliCalibViewerGUItime", this, "HandleButtonsDrawSel()");
  fRadioXhist->SetToolTipText("Draw the distribution of the variable");
  
  fRadioXrun = new TGRadioButton(fContDrawSelSubRunTime, ":Run", kRadioXrun);
  fContDrawSelSubRunTime->AddFrame(fRadioXrun, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
  fRadioXrun->Connect("Clicked()", "AliCalibViewerGUItime", this, "HandleButtonsDrawSel()");
  fRadioXrun->SetToolTipText("Use run number as x-value");
  
  // Time radio button
  fRadioXtime = new TGRadioButton(fContDrawSelSubRunTime, ":Time", kRadioXtime);
  fContDrawSelSubRunTime->AddFrame(fRadioXtime, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
  fRadioXtime->Connect("Clicked()", "AliCalibViewerGUItime", this, "HandleButtonsDrawSel()");
  fRadioXtime->SetToolTipText("Use time stamp number as x-value");
  
  
  // list of variables
  fListVariables = new TGListBox(fContDrawSel);
  fContDrawSel->AddFrame(fListVariables, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  fListVariables->Connect("Selected(Int_t)", "AliCalibViewerGUItime", this, "DoNewSelection()");

  
//-------------------- run type selection ------------------------
  // Parameter label
  fLblRunType = new TGLabel(fContDrawSel, "Run Type:");
  fLblRunType->SetTextJustify(kTextLeft);
  fContDrawSel->AddFrame(fLblRunType, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  
  fComboRunType = new TGComboBox(fContDrawSel);
  fComboRunType->EnableTextInput(kFALSE);
  fContDrawSel->AddFrame(fComboRunType, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  fComboRunType->Connect("Selected(Int_t)", "AliCalibViewerGUItime", this, "DoDraw()");
  fComboRunType->Resize(0, 22);
  
  //-------------------- parameter selection ------------------------
  // Parameter label
  fLblPar = new TGLabel(fContDrawSel, "Parameter:");
  fLblPar->SetTextJustify(kTextLeft);
  fContDrawSel->AddFrame(fLblPar, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  
  fNmbPar = new TGNumberEntry(fContDrawSel, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 71);
  fContDrawSel->AddFrame(fNmbPar, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  fNmbPar->Connect("ValueSet(Long_t)", "AliCalibViewerGUItime", this, "DoParLimitChange()");
  fNmbPar->SetState(kFALSE);
  
  //-------------------- calibration type selection ------------------------
  // label
  // draw selection group
  fContCalibType = new TGGroupFrame(fContLeft, "Calib type selection", kVerticalFrame | kFitWidth | kFitHeight);
  fContLeft->AddFrame(fContCalibType, new TGLayoutHints(kLHintsExpandX , 0, 0, 10, 0));
    
    // list of variables
  fListCalibType = new TGListBox(fContCalibType);
  fContCalibType->AddFrame(fListCalibType, new TGLayoutHints(kLHintsNormal | kLHintsExpandX , 0, 0, 0, 0));
  fListCalibType->Connect("Selected(Int_t)", "AliCalibViewerGUItime", this, "DoChangeSelectionList()");
  fListCalibType->Resize(0,88);
  fListCalibType->SetMultipleSelections();
  
  
     // ==========================================================================
   // ************************* content of fContCenter *************************
   // ========================================================================
   // main drawing canvas
  fCanvMain = new TRootEmbeddedCanvas("GUItime_Canvas", fContCenter, 200, 200, kFitWidth | kFitHeight);
  fContCenter->AddFrame(fCanvMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  fCanvMain->GetCanvas()->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)", "AliCalibViewerGUItime", this, "MouseMove(Int_t, Int_t, Int_t, TObject*)");
  fCanvMain->GetCanvas()->SetToolTipText("The Main_Canvas, here your plots are displayed.");
  fCanvMain->GetCanvas()->SetRightMargin(0.062);
  fCanvMain->GetCanvas()->SetLeftMargin(0.15);
  
   // =========================================================================
   // ************************* content of fContRight *************************
   // ========================================================================
  //group frame for value information
  fContValues = new TGGroupFrame(fContRight, "Data point info", kVerticalFrame | kFitWidth | kFitHeight);
  fContRight->AddFrame(fContValues, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
  //set layout manager
  fContValues->SetLayoutManager(new TGMatrixLayout(fContValues, 0, 2, 4));
  //value information labels

  //run number label
  fLblRunNumber = new TGLabel(fContValues, "Run:");
  fLblRunNumber->SetTextJustify(kTextLeft);
  fContValues->AddFrame(fLblRunNumber, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
  //run number
  fLblRunNumberVal = new TGLabel(fContValues, "0000000");
  fLblRunNumberVal->SetTextJustify(kTextLeft);
  fContValues->AddFrame(fLblRunNumberVal, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 0, 0, 0, 0));
  //time stamp label
  fLblRunTime = new TGLabel(fContValues, "Time:");
  fLblRunTime->SetTextJustify(kTextLeft);
  fContValues->AddFrame(fLblRunTime, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
  //run number
  fLblRunTimeVal = new TGLabel(fContValues, "00.00.0000\n00:00:00");
  fLblRunTimeVal->SetTextJustify(kTextLeft);
  fContValues->AddFrame(fLblRunTimeVal, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 0, 0, 0, 0));
  //value label x
  fLblValueX = new TGLabel(fContValues, "x-Value:");
  fLblValueX->SetTextJustify(kTextLeft);
  fContValues->AddFrame(fLblValueX, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
  //value x
  fLblValueXVal = new TGLabel(fContValues, "00.000e+00");
  fLblValueXVal->SetTextJustify(kTextRight);
  fContValues->AddFrame(fLblValueXVal, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  //value label y
  fLblValueY = new TGLabel(fContValues, "y-Value:");
  fLblValueY->SetTextJustify(kTextLeft);
  fContValues->AddFrame(fLblValueY, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
  //value y
  fLblValueYVal = new TGLabel(fContValues, "00.000e+00");
  fLblValueYVal->SetTextJustify(kTextRight);
  fContValues->AddFrame(fLblValueYVal, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   // draw button
  fBtnDumpRuns = new TGTextButton(fContRight, "&Dump runs");
  fContRight->AddFrame(fBtnDumpRuns, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
  fBtnDumpRuns->Connect("Clicked()", "AliCalibViewerGUItime", this, "DoDumpRuns()");
  fBtnDumpRuns->SetToolTipText("Press to dump the run numbers of the current selection.");
  // save button
  fBtnSave = new TGTextButton(fContRight, "&Save picture");
  fContRight->AddFrame(fBtnSave, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
  fBtnSave->SetToolTipText("Press to save the current picture in the canvas");
  fBtnSave->Connect("Clicked()", "AliCalibViewerGUItime", this, "SavePicture()");

  //group frame for loading a file
  fContLoad = new TGGroupFrame(fContRight, "Load file", kVerticalFrame | kFitWidth | kFitHeight);
  fContRight->AddFrame(fContLoad, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));

    fContFilename = new TGCompositeFrame(fContLoad, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
    fContLoad->AddFrame(fContFilename, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 0));

      fTxtFilename = new TGTextEntry(fContFilename, "Input file", 100);
      fContFilename->AddFrame(fTxtFilename, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      fTxtFilename->Connect("DoubleClicked()", "AliCalibViewerGUItime", this, "HandleLoadRunTextEntry()");
     
    fContTreeName = new TGCompositeFrame(fContLoad, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
    fContLoad->AddFrame(fContTreeName, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 0));

      fLblTreeName = new TGLabel(fContTreeName, "Tree Name:");
      fLblTreeName->SetTextJustify(kTextLeft);
      fContTreeName->AddFrame(fLblTreeName, new TGLayoutHints(kLHintsNormal | kLHintsCenterY, 0, 1, 0, 0));

      fTxtTreeName = new TGTextEntry(fContTreeName, "trdTree");
      if(fDetector.Contains("TPC")) fTxtTreeName->SetText("dcs");
      fContTreeName->AddFrame(fTxtTreeName, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));

    fBtnLoadFile = new TGTextButton(fContLoad, "Load file", 100);
    fContLoad->AddFrame(fBtnLoadFile, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
    fBtnLoadFile->SetToolTipText("Load a file into viewer");
    fBtnLoadFile->Connect("Clicked()", "AliCalibViewerGUItime", this, "HandleLoadRunButtons()");

   // =========================================================================
   // ****************** bottom content of fContTopBottom *********************
   // =========================================================================
  
  // custom options container
  // --- fComboCustom --- the custom draw line on the very low
  fContCustom = new TGCompositeFrame(fContTopBottom, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
  fContTopBottom->AddFrame(fContCustom, new TGLayoutHints(kLHintsExpandX, 10, 0, 0, 0));

         // ------------------------- content of fContCustom -------------------------
  fLblCustomDraw = new TGLabel(fContCustom, "Custom draw: ");
  fLblCustomDraw->SetTextJustify(kTextLeft);
  fContCustom->AddFrame(fLblCustomDraw, new TGLayoutHints(kLHintsNormal, 5, 0, 0, 0));
         // text field for custom draw command
  fComboCustomDraw = new TGComboBox(fContCustom);
  fComboCustomDraw->Resize(0, 22);
  fComboCustomDraw->EnableTextInput(kTRUE);
  fContCustom->AddFrame(fComboCustomDraw, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  fComboCustomDraw->Connect("ReturnPressed()", "AliCalibViewerGUItime", this, "DoCustomDraw()");
  fComboCustomDraw->Connect("Selected(Int_t)", "AliCalibViewerGUItime", this, "DoCustomDraw()");
  
  
      // additional cuts container
  fContCustomCuts = new TGCompositeFrame(fContTopBottom, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
  fContTopBottom->AddFrame(fContCustomCuts, new TGLayoutHints(kLHintsExpandX, 10, 0, 0, 0));
  
         // ------------------------- content of fContCustomCuts -------------------------
  fLblCustomCuts = new TGLabel(fContCustomCuts, "Custom cuts:  ");
  fLblCustomCuts->SetTextJustify(kTextLeft);
  fContCustomCuts->AddFrame(fLblCustomCuts, new TGLayoutHints(kLHintsNormal, 5, 0, 0, 0));
         // combo text field for additional cuts
  fComboCustomCuts = new TGComboBox(fContCustomCuts);
//   fComboCustomCuts->Resize(0, fLblValueY->GetDefaultHeight());
  fComboCustomCuts->Resize(0, 22);
  fComboCustomCuts->EnableTextInput(kTRUE);
  fContCustomCuts->AddFrame(fComboCustomCuts, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  fComboCustomCuts->Connect("ReturnPressed()", "AliCalibViewerGUItime", this, "DoCustomCutsDraw()");
  fComboCustomCuts->Connect("Selected(Int_t)", "AliCalibViewerGUItime", this, "DoCustomCutsDraw()");

  if(fDetector.Contains("TPC")) SetWindowName("AliTPCCalibViewer GUI - Time");
  else if (fDetector.Contains("TRD")) SetWindowName("AliTRDCalibViewer GUI - Time");
  else SetWindowName("Unknown detector");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
}

//______________________________________________________________________________
void AliCalibViewerGUItime::SetInitialValues(){
  //
  // Set inital selections of the gui
  //
  fRadioXrun->SetState(kButtonDown);
  fRadioXtime->SetState(kButtonUp);
}

//______________________________________________________________________________
void AliCalibViewerGUItime::UseFile(const char* fileName, const char* treeName) {
  //
  // retrieve tree from file
  //
  TString s(fileName);
  TObjArray *arr=s.Tokenize(" ");
  TIter next(arr);
  TObject *o=0;
  if (fTree) delete fTree;
  fTree=new TChain(treeName);
  while ( (o=next()) ){
    TString tmpString(o->GetName());
    std::cout << "Adding " << tmpString.Data() << " to the chain" << std::endl;
    if(tmpString.Contains(".root"))
      fTree->AddFile(tmpString.Data());
    else {
      tmpString += "/*.root";
      fTree->Add(tmpString.Data());
    }
  }
  delete arr;
  if (!CheckChain())
    return;
  //  UseConfigFile(fConfigFile.Data());
  FillCalibTypes();
  Reload();
}

//______________________________________________________________________________
void AliCalibViewerGUItime::UseChain(TChain* chain)
{
  //
  // load a chain
  //
  fTree=chain;
  if (!CheckChain()) {
    std::cout << "AliCalibViewerGUItime::UseChain !checkchain OUT" << std::endl;
    return;
  }
  //set configuration file
  //  UseConfigFile(fConfigFile.Data());
  FillCalibTypes();
  Reload();
}

//______________________________________________________________________________
Bool_t AliCalibViewerGUItime::CheckChain()
{
  //
  // check whether cahin has entries
  // decide whether to draw graphs in 2D
  //
  if (!fTree) {
    return kFALSE;
  }
  fTree->Lookup();
  Long64_t entries=fTree->GetEntries();
  if (entries==0){
    AliError("No entries found in chain");
    return kFALSE;
  }
  //check whether to draw graphs
  CheckDrawGraph();
  return kTRUE;
}

//______________________________________________________________________________
void AliCalibViewerGUItime::FillRunTypes()
{
  //
  //Loop over the tree entries and fill the run types
  //
  if (!fTree) {
    return;
  }
  Int_t id=0;
  fComboRunType->RemoveAll();
  fComboRunType->AddEntry("ALL",id++);
  fComboRunType->Select(0,kFALSE);
  if (!fTree->GetBranch("runType.")) {
    return;
  }
  TObjString *runType=0x0;
  Int_t nevets=fTree->GetEntries();
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("runType.*",1);
  fTree->SetBranchAddress("runType.",&runType);
  for (Int_t iev=0;iev<nevets;++iev){
    fTree->GetEntry(iev);
    TString type=runType->String();
    if (!type.IsNull()&&!fComboRunType->FindEntry(type)) fComboRunType->AddEntry(type,id++);
  }
  fTree->ResetBranchAddresses();
  fTree->SetBranchStatus("*",1);
}

//______________________________________________________________________________
void AliCalibViewerGUItime::FillCalibTypes()
{
  //
  // loop over configuration and fill calibration types
  //
  Int_t id=0;
  fListCalibType->RemoveAll();
  //  TObject *o=0x0;
  //  fConfigParser->ResetIter();
  TString type="UNSPECIFIED";
  fListCalibType->AddEntry(type.Data(),id);
  fListCalibType->Select(id++);
}

//______________________________________________________________________________
void AliCalibViewerGUItime::CheckDrawGraph()
{
  //
  // Check whether to draw graphs in 2D mode based on the number of entries in the chain
  // GetEstimate() returns the maximum size of the arrays stored in GetV1()...
  //
  if (!fTree) {
    return;
  }
  fNoGraph=kTRUE;
  if (fTree->GetEntries()<fTree->GetEstimate()) fNoGraph=kFALSE;
}

//______________________________________________________________________________
void AliCalibViewerGUItime::Reload(Int_t /*first*/)
{
  //
  // reload the gui contents, this is needed after the input tree has changed
  //
  if ( !fTree ) {
    return;
  }
  
  FillRunTypes();
  
  //activate all branches
  fTree->SetBranchStatus("*",1);
  //reset variables list
  fListVariables->RemoveAll();
  //get selected calibration types
  TList calibTypes;
  fListCalibType->GetSelectedEntries(&calibTypes);
    
  TObjArray *branchList = fTree->GetListOfBranches();
  if ( !branchList ) return;
  TIter nextBranch(branchList);
  Int_t idCount=0,id=0;
  TObject *objBranch=0;
  while ( (objBranch=nextBranch()) ){
    TString branchName(objBranch->GetName());
    TString branchTitle(objBranch->GetName());
    if (branchName == "run" || branchName == "time" || branchName == "runType.") continue;
    Bool_t active=kTRUE;
    TString calibType="UNSPECIFIED";
   
    //check if branch is in selected calibration types
    //if not, don't show it in the list and deactivate the branch.
    Bool_t calibActive=kFALSE;
    TIter nextCalib(&calibTypes);
    TObject *objCalib=0;
    while ( (objCalib=nextCalib()) ) {
      if (calibType==objCalib->GetTitle()) calibActive=kTRUE;
    }
    active&=calibActive;
    if (!active){
      TString s=branchName;
      if (branchName.EndsWith(".")) s+="*";
      fTree->SetBranchStatus(s.Data(),0);
      continue;
    }
    fListVariables->AddEntry(SubstituteUnderscores(branchTitle.Data()),id);
    //fListVariables->Select(id);
    ++idCount;
  }
  //trick to display modifications
  fListVariables->Resize(fListVariables->GetWidth()-1, fListVariables->GetHeight());
  fListCalibType->Resize(fListCalibType->GetWidth()+1, fListCalibType->GetHeight());
}

//______________________________________________________________________________
void AliCalibViewerGUItime::AddReferenceTree(const char* treeFileName, const char* refName)
{
  //
  // map of reference trees that should always be attached to the CalibViewerGUI
  //
  fMapRefTrees->Add(new TObjString(refName), new TObjString(treeFileName));
}

//______________________________________________________________________________
const char* AliCalibViewerGUItime::GetDrawString(){
  //
  // create draw string for ttree by combining the user requestsa
  //
  
  TString selectedVariable="";
  Int_t id=-1;
  if (!fListVariables->GetSelectedEntry()) {
    return "";
  }
  selectedVariable = fListVariables->GetSelectedEntry()->GetTitle();
  id=fListVariables->GetSelectedEntry()->EntryId();
  TString branchName=selectedVariable;
  //  const TObject *key=(*fConfigParser)(id);
  //  if (key) branchName=(*fConfigParser)(id)->GetName();
  //treat case of TVector
  if (branchName.EndsWith(".")){
    Int_t par = (Int_t)(fNmbPar->GetNumber());
    branchName.Append(Form("fElements[%d]",par));
  }

  return branchName.Data();
}
//______________________________________________________________________________
const char* AliCalibViewerGUItime::GetDrawOption() const {
  //
  // get user selected draw options
  //
  TString drawOpt;
  if (fComboAddDrawOpt->GetSelectedEntry()) drawOpt=fComboAddDrawOpt->GetSelectedEntry()->GetTitle();
  if (fChkDrawOptSame->GetState()==kButtonDown && !drawOpt.Contains("same",TString::kIgnoreCase))
    drawOpt+="same";
  return drawOpt.Data();
}
//______________________________________________________________________________
void AliCalibViewerGUItime::GetCutString(TString &cutStr){
  //
  // create cut string
  //
  TCut cuts(fComboCustomCuts->GetTextEntry()->GetText());
  TString runType="";
  if (fComboRunType->GetSelectedEntry()) runType=fComboRunType->GetSelectedEntry()->GetTitle();
  if (runType!="ALL"&&!runType.IsNull()) cuts+=Form("runType.String().Data()==\"%s\"",runType.Data());
  cutStr=cuts.GetTitle();
}
//______________________________________________________________________________
void AliCalibViewerGUItime::UpdateValueArrays(Bool_t withGraph)
{
  //
  // Update arrays of runs
  //
  if (!withGraph){
    fValuesX.ResizeTo(1);
    fValuesY.ResizeTo(1);
    fRunNumbers.ResizeTo(1);
    fTimeStamps.ResizeTo(1);
  } else {
    fValuesX.ResizeTo(fTree->GetSelectedRows());
    fValuesY.ResizeTo(fTree->GetSelectedRows());
    fRunNumbers.ResizeTo(fTree->GetSelectedRows());
    fTimeStamps.ResizeTo(fTree->GetSelectedRows());
    fValuesY.SetElements(fTree->GetV3());
    fRunNumbers.SetElements(fTree->GetV1());
    fTimeStamps.SetElements(fTree->GetV2());
  }
}

//______________________________________________________________________________
void AliCalibViewerGUItime::GetHistogramTitle(TString &title)
{
  //
  // Create string for histogram title
  //

  title=fDrawString;
  Int_t pos=title.First(">>");
  if (pos>0) title=title(0,pos);
  if (!fIsCustomDraw){
    if (fRadioXrun->GetState()==kButtonDown){
      title+=":Run";
    } else if (fRadioXtime->GetState()==kButtonDown){
      title+=":Date";
    }
  }
  TString cuts;
  GetCutString(cuts);
  TObjArray *arr=title.Tokenize(":");
  TObject *o=0x0;
  title+=" {";
  title+=cuts;
  title+="}";
  TIter next(arr,kIterBackward);
  while ( (o=next()) ){
    TString varName=o->GetName();
    title+=";";
    title+=varName;
  }
  delete arr;
}

//______________________________________________________________________________
void AliCalibViewerGUItime::AdjustYRange()
{
  //
  // adjust the range of the Y axis
  //
  TIter nextGraphicObject(fTrashBox);
  TObject *o=0x0;
  Float_t min=0,max=0;
  while ( (o=nextGraphicObject()) ){
    if (o->IsA()==TGraph::Class()){
      TGraph *gr=(TGraph*)o;
      if (min==max) {
        min=TMath::MinElement(gr->GetN(),gr->GetY());
        max=TMath::MaxElement(gr->GetN(),gr->GetY());
      } else {
        Float_t currmin=TMath::MinElement(gr->GetN(),gr->GetY());
        Float_t currmax=TMath::MaxElement(gr->GetN(),gr->GetY());
        if (currmax>max) max=currmax;
        if (currmin<min) min=currmin;
      }
    }
  }
  if (min!=max){
    if (min!=0) min=min-(max-min)/10;
    if (max!=0) max=max+(max-min)/10;
    fCurrentHist->SetMinimum(min);
    fCurrentHist->SetMaximum(max);
  }
}
//______________________________________________________________________________
void AliCalibViewerGUItime::DoDraw() {
  //
  // Draw graphics
  //
  TString drawString=fDrawString;
  TString cutString;
  GetCutString(cutString);
  TString optString  = GetDrawOption();
  Bool_t graphOutput=!fNoGraph;  //ttree buffer for V1, V2... too small
  graphOutput&=(drawString.First(">>")<0); //histogram output in custom draw
  graphOutput&=fRadioXhist->GetState()!=kButtonDown; //histogram drawing selected
  graphOutput&=!(fIsCustomDraw&&!fDrawString.Contains(":")); //custom draw 1D
  Bool_t drawSame=optString.Contains("same",TString::kIgnoreCase);
//   optString+="goff";
  if (graphOutput) {
    drawString.Prepend("run:time:");
    optString="goff";
  }else{
    if (!fIsCustomDraw){
      if (fRadioXrun->GetState()==kButtonDown){
        drawString+=":run";
      } else if (fRadioXtime->GetState()==kButtonDown){
        drawString+=":time";
      }
    }
  }
  TVirtualPad *padsave=gPad;
  fCanvMain->GetCanvas()->cd();
  //delete old histograms and graphs
  if (!drawSame){
    fTrashBox->Delete();
    fCurrentGraph=0x0;
    fCurrentHist=0x0;
  }

  //select data
  fTree->Draw(drawString.Data(),cutString.Data(),optString.Data());
  if (fTree->GetSelectedRows()==-1) {
    return;
  }
  UpdateValueArrays(graphOutput);
  TString title;
  GetHistogramTitle(title);
  Bool_t drawGraph=kFALSE;
  if (graphOutput){
    if (fIsCustomDraw){
      if (fDrawString.Contains(":")){
        fValuesX.SetElements(fTree->GetV4());
        drawGraph=kTRUE;
      } else {
        drawGraph=kFALSE;
      }
    }else{
      drawGraph=kTRUE;
      if (fRadioXrun->GetState()==kButtonDown){
        fValuesX.SetElements(fTree->GetV1());
      } else if (fRadioXtime->GetState()==kButtonDown){
        fValuesX.SetElements(fTree->GetV2());
      } else {
        drawGraph=kFALSE;
      }
    }
  }//create graph according to selection
  if (drawGraph){
    TGraph *graph=new TGraph(fValuesX,fValuesY);
    TString grDraw="p";
    if (!drawSame) grDraw+="a";
    if (!fIsCustomDraw) grDraw+="l";
    graph->Draw(grDraw.Data());
    graph->SetEditable(kFALSE);
    TH1 *hist=graph->GetHistogram();
    hist->SetTitle(title.Data());
    fTrashBox->Add(graph);
    graph->SetLineColor(fTrashBox->GetEntries());
    graph->SetMarkerColor(fTrashBox->GetEntries());
    if (!drawSame) {
      fCurrentGraph=graph;
      fCurrentHist=hist;
    }
  } else {
    TH1 *hist=fTree->GetHistogram();
    hist->SetTitle(title.Data());
    fTrashBox->Add(hist);
    hist->SetLineColor(fTrashBox->GetEntries());
    hist->SetMarkerColor(fTrashBox->GetEntries());
    if (!drawSame) fCurrentHist=hist;
  }
  
  //Set time axis if choosen as x-variables
  if (fRadioXtime->GetState()==kButtonDown&&!fIsCustomDraw&&!drawSame){
    TAxis *xaxis=fCurrentHist->GetXaxis();
    xaxis->SetTimeFormat("#splitline{%d.%m}{%H:%M}");
    xaxis->SetTimeDisplay(1);
    xaxis->SetLabelOffset(xaxis->GetLabelOffset()*3);
    xaxis->SetLabelSize(xaxis->GetLabelSize()/1.3);
  }
  if (!drawSame) {
  //Set title offset
    fCurrentHist->GetYaxis()->SetTitleOffset(1.5);
  } else {
    //adjust y-range
    AdjustYRange();
  }
  gPad->Modified();
  gPad->Update();
  padsave->cd();
}

//______________________________________________________________________________
void AliCalibViewerGUItime::DoDumpRuns()
{
  //
  // Dump the current run numbers to stdout
  //
  Int_t npoints=fRunNumbers.GetNrows();
  Int_t    *sortIndex = new Int_t[npoints];
  TMath::Sort(npoints,fRunNumbers.GetMatrixArray(),sortIndex,kFALSE);
  Int_t run=0, prevRun=-1;
  
  for (Int_t irun=0;irun<npoints;++irun){
    run=(Int_t)fRunNumbers.GetMatrixArray()[sortIndex[irun]];
    if (run!=prevRun) std::cout << Form("%d",run) << std::endl;
    prevRun=run;
  }
  delete [] sortIndex;
}

//______________________________________________________________________________
void AliCalibViewerGUItime::DoParLimitChange()
{
  //
  // DoParLimitChange()
  //
  UpdateParName();
  DoDraw();
}

//______________________________________________________________________________
void AliCalibViewerGUItime::DoNewSelection() {
  //
   // decides whether to redraw if user makes another selection
   //
  UpdateParLimits();
  fDrawString=GetDrawString();
  fIsCustomDraw=kFALSE;
  DoDraw();
}

//______________________________________________________________________________
void AliCalibViewerGUItime::DoCustomDraw()
{
  //
  // Custom draw (TTree::Draw syntax)
  //
  fDrawString=fComboCustomDraw->GetTextEntry()->GetText();
  fNmbPar->SetState(kFALSE);
  fIsCustomDraw=kTRUE;
  DoDraw();
}

//______________________________________________________________________________
void AliCalibViewerGUItime::DoCustomCutsDraw()
{
  //
  // Custom cuts (TTree::Draw syntax)
  //
  if (fIsCustomDraw) DoCustomDraw();
  else {
    fDrawString=GetDrawString();
    fIsCustomDraw=kFALSE;
    DoDraw();
  }
}

//______________________________________________________________________________
void AliCalibViewerGUItime::HandleButtonsDrawSel(Int_t id)
{
  //
  // Draw selection button handling (x-variable)
  //

  if (id == -1) {
    TGButton *btn = (TGButton *) gTQSender;
    id = btn->WidgetId();
  }
  
  Bool_t doDraw=kFALSE;
  switch (id) {
  case (kRadioXhist):
    doDraw=(fRadioXtime->GetState()==kButtonDown||fRadioXrun->GetState()==kButtonDown);
    if (doDraw){
      fRadioXrun->SetState(kButtonUp);
      fRadioXtime->SetState(kButtonUp);
    }
    break;
  case (kRadioXrun):
    doDraw=(fRadioXtime->GetState()==kButtonDown||fRadioXhist->GetState()==kButtonDown);
    if (doDraw){
      fRadioXhist->SetState(kButtonUp);
      fRadioXtime->SetState(kButtonUp);
    }
    break;
  case (kRadioXtime):
    doDraw=(fRadioXhist->GetState()==kButtonDown||fRadioXrun->GetState()==kButtonDown);
    if (doDraw){
      fRadioXrun->SetState(kButtonUp);
      fRadioXhist->SetState(kButtonUp);
    }
    break;
  }
  if (doDraw) DoCustomCutsDraw();
}
//______________________________________________________________________________
void AliCalibViewerGUItime::UpdateParName()
{
  //
  // change parameter name
  //
  
  Int_t par = (Int_t)(fNmbPar->GetNumber());
  TString parName="";
  parName.Form("%d",par);
  fLblPar->SetText(Form("Parameter: %s",parName.Data()));
  fDrawString=GetDrawString();
  fIsCustomDraw=kFALSE;
}

//______________________________________________________________________________
void AliCalibViewerGUItime::UpdateParLimits()
{
  //
  // Adjust limits for TVectorT based variables
  //
  if (!fTree) return;
  TString selectedVariableTitle="";
  Int_t id=-1;
  if (!fListVariables->GetSelectedEntry()) return;
  selectedVariableTitle = fListVariables->GetSelectedEntry()->GetTitle();
  id=fListVariables->GetSelectedEntry()->EntryId();
  TString selectedVariable=selectedVariableTitle;
  //  const TObject *key=(*fConfigParser)(id);
  //  if (key) selectedVariable=(*fConfigParser)(id)->GetName();
  
  if (selectedVariable.IsNull()||!selectedVariable.EndsWith(".")) {
    fNmbPar->SetState(kFALSE);
    fLblPar->SetText("Parameter: none");
    return;
  }
  TVectorD *vD=0x0;
  TVectorF *vF=0x0;
  Int_t maxPar=0;
  fTree->GetEntry(1);
  TBranch *branch=fTree->GetTree()->GetBranch(selectedVariable.Data());
  TString branchClass=branch->GetClassName();
  Int_t event=0;
  if (branchClass=="TVectorT<double>"){
//     branch->SetAddress(&vD);
    fTree->SetBranchAddress(selectedVariable.Data(),&vD);
    while (maxPar<2&&event<fTree->GetEntries()){
      fTree->GetEntry(event++);
      maxPar=vD->GetNrows();
    }
  } else if (branchClass=="TVectorT<float>"){
//     branch->SetAddress(&vF);
    fTree->SetBranchAddress(selectedVariable.Data(),&vF);
    while (maxPar<2&&event<fTree->GetEntries()){
      fTree->GetEntry(event++);
      maxPar=vF->GetNrows();
    }
  } else {
    //class not known
    fNmbPar->SetState(kFALSE);
    return;
  }
  fTree->SetBranchAddress(selectedVariable.Data(),0x0);
  if (fNmbPar->GetNumMax()!=maxPar-1) fNmbPar->SetNumber(0);
  fNmbPar->SetLimitValues(0,maxPar-1);
  fNmbPar->SetState(kTRUE);
  UpdateParName();
}
//______________________________________________________________________________
void AliCalibViewerGUItime::MouseMove(Int_t event, Int_t x, Int_t y, TObject */*selected*/)
{
  //
  // handle mouse events in the draw canvas
  //
  UInt_t dd=0,mm=0,yy=0,HH=0,MM=0,SS=0,run=0;
  Double_t valx=0.,valy=0.;
  if (!fCurrentGraph) {
    fLblRunNumberVal->SetText(Form("%07u",run));
    fLblRunTimeVal->SetText(Form("%02u.%02u.%04u\n%02u:%02u:%02u",dd,mm,yy,HH,MM,SS));
    fLblValueXVal->SetText(Form("%.3e", valx));
    fLblValueYVal->SetText(Form("%.3e", valy));
    return;
  }
  TVirtualPad *padsave=gPad;
  fCanvMain->GetCanvas()->cd();
  Int_t n=fValuesY.GetNrows();
  Double_t *arr=0x0;
  arr=fValuesX.GetMatrixArray();

  Int_t minDist=1000000;
  Int_t minPoint=-1;
  for (Int_t i=0;i<n;++i){
    Int_t pxp = gPad->XtoAbsPixel(gPad->XtoPad(arr[i]));
    Int_t pyp = gPad->YtoAbsPixel(gPad->YtoPad(fValuesY[i]));
    Int_t d   = TMath::Nint(TMath::Sqrt(TMath::Abs(pxp-x) + TMath::Abs(pyp-y)));
    if (d < minDist) {
      minDist  = d;
      minPoint = i;
    }
  }
  if (minDist<2){
    TTimeStamp t((Int_t)fTimeStamps[minPoint]);
    t.GetDate(kTRUE,0,&yy,&mm,&dd);
    t.GetTime(kTRUE,0,&HH,&MM,&SS);
    run=(UInt_t)fRunNumbers[minPoint];
    valx=fValuesX[minPoint];
    valy=fValuesY[minPoint];
  } else {
    dd=0;mm=0;yy=0;HH=0;MM=0;SS=0;run=0;
    valx=0.;
    valy=0.;
  }
  fLblRunNumberVal->SetText(Form("%07u",run));
  fLblRunTimeVal->SetText(Form("%02u.%02u.%04u\n%02u.%02u.%02u",dd,mm,yy,HH,MM,SS));
  if (fIsCustomDraw){
    fLblValueXVal->SetText(Form("%.3e", valx));
  }else{
    if (fRadioXrun->GetState()==kButtonDown){
      fLblValueXVal->SetText("Run");
    } else if (fRadioXtime->GetState()==kButtonDown){
      fLblValueXVal->SetText("Time");
    }
  }
  fLblValueYVal->SetText(Form("%.3e", valy));
  padsave->cd();
  if (run==0) {
    return;
  }
  if (event == kButton1Double ){
    SetGuiTree(run);
  }
  //find closes point of current selection
}

//______________________________________________________________________________
void AliCalibViewerGUItime::SetGuiTree(Int_t run)
{
  //
  // create the AliCalibViewerGUI tree for run
  // cache tree in directory fOutputCacheDir
  // retrieve file from this directory if it already exists
  //

  //try to find file for run in fOutputCacheDir
  TString fileName=fOutputCacheDir;
  if (!fileName.EndsWith("/")) fileName+="/";
  if(fDetector.Contains("TPC")) fileName+=Form("guiTreeRun_%d.root",run);
  else if(fDetector.Contains("TRD")) fileName+=Form("trdOCDBDetails_run%d.root", run);
  else return;
  Bool_t load=kTRUE;
  TFile f(fileName.Data());
  if (!f.IsOpen()){
    if(fDetector.Contains("TRD")) 
      load = fCalibViewerGUI->CreateDetailsTree(run, fileName.Data(), "nothing");
    if(fDetector.Contains("TPC")) 
      load = fCalibViewerGUI->CreateDetailsTree(run, fileName.Data(), "nothing");

    if (!load){
      fCalibViewerGUI->Reset();
      if (fCalibViewerGUItab) fCalibViewerGUItab->SetText(new TGString(Form("Detail - XXXXX")));
      return;
    }
  }
  f.Close();
  if(fDetector.Contains("TPC")) fCalibViewerGUI->Initialize(fileName.Data(), "calPads");
  else fCalibViewerGUI->Initialize(fileName.Data(), "TRDcalibDetails");
  if (fCalibViewerGUItab) fCalibViewerGUItab->SetText(new TGString(Form("Detail - %07d",run)));
  TIter nextRefTree(fMapRefTrees);
  TObject *o=0x0;
  while ( (o=nextRefTree()) ){
    if(fDetector.Contains("TPC"))
      fCalibViewerGUI->GetViewer()->AddReferenceTree(fMapRefTrees->GetValue(o)->GetName(),"calPads",o->GetName());
    else
      fCalibViewerGUI->GetViewer()->AddReferenceTree(fMapRefTrees->GetValue(o)->GetName(),"TRDcalibDetails",o->GetName());
  }
  //if(fDetector.Contains("TPC")) ((AliTPCCalibViewerGUI_new*)fCalibViewerGUI)->Reload();
  //else ((AliTRDCalibViewerGUI*)fCalibViewerGUI)->Reload();
  fCalibViewerGUI->Reload();
}

//______________________________________________________________________________
const char* AliCalibViewerGUItime::SubstituteUnderscores(const char* in)
{
  //
  // Substitute underscores from the branch names
  //
  TString s(in);
  s.ReplaceAll("_{","|{");
  s.ReplaceAll("_"," ");
  s.ReplaceAll("|{","_{");
  return s.Data();
}

//______________________________________________________________________________
void AliCalibViewerGUItime::SavePicture() {
   // 
   // saves the current picture
   // 
   
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
  TString dir(".");
  TGFileInfo fi;
  fi.fFileTypes = kSaveAsTypes;
  // fi.fIniDir    = StrDup(dir);
  fi.fOverwrite = kFALSE;
  new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDSave, &fi);
  if (fi.fFilename && strlen(fi.fFilename)) {
    fCanvMain->GetCanvas()->Print(fi.fFilename);
  }  
}

//______________________________________________________________________________
void AliCalibViewerGUItime::HandleLoadRunButtons() {
  //
  //  Handle the buttons
  //
  Int_t id = ((TGTextButton*)gTQSender)->WidgetId();
  if(id==100) {
    UseFile(fTxtFilename->GetText(), fTxtTreeName->GetText());
    return;
  }
  else
    return;
}

//______________________________________________________________________________
void AliCalibViewerGUItime::HandleLoadRunTextEntry() {
  //
  //  Handle the text entries
  //
  // buttons ID
  // 100 - fTxtFilename
  // 101 - fTxtConfigFile
  
  Int_t id = ((TGTextEntry*)gTQSender)->WidgetId();
  if(id>=100 && id<=103) {
    const char *kTypes[] = {
      "All files",    "*",
      0,              0
    };
    TString dir(".");
    TGFileInfo fi;
    fi.fFileTypes = kTypes;
    fi.fOverwrite = kFALSE;
    new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen, &fi);
    if(fi.fFilename && strlen(fi.fFilename)) {
      if(id==100) {
        fTxtFilename->SetText(fi.fFilename);
      }
      if(id==101) {
	//        fTxtConfigFile->SetText(fi.fFilename);
	//        fConfigFile=fi.fFilename;
      }
      if(id==102) {
	//fTxtRunList->SetText(fi.fFilename);
      }
      if(id==103) {
	//fTxtOutputOCDB->SetText(fi.fFilename);
      }
    }
    return;
  }
  else {
    return;
  }
}

//__________________________________________________________________
TObjArray* AliCalibViewerGUItime::ShowGUI(TChain* chain) {
  //
  //  Launch the time trending GUI. Load the TChain chain into the viewer
  //
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 650);
  frmMain->SetWindowName("Calib Viewer - time trend");
  frmMain->SetCleanup(kDeepCleanup);

  TGTab* tabMain = new TGTab(frmMain, 1000, 600);
  frmMain->AddFrame(tabMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

  TGCompositeFrame* tabCont1 = tabMain->AddTab("Time");

  AliCalibViewerGUItime* calibViewerTime = new AliCalibViewerGUItime(tabCont1, 1000, 650);
  tabCont1->AddFrame(calibViewerTime, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  calibViewerTime->UseChain(chain);

  TObjArray *guiArray = new TObjArray();
  guiArray->Add(calibViewerTime);

  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();
  return guiArray;
}

//__________________________________________________________________
TObjArray* AliCalibViewerGUItime::ShowGUI() {
  //
  //  Launch the time trending GUI. The GUI will be empty but trees can be loaded by using the GUI interface
  //
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 650);
  frmMain->SetWindowName("Calib Viewer - time trend");
  frmMain->SetCleanup(kDeepCleanup);

  TGTab* tabMain = new TGTab(frmMain, 1000, 600);
  frmMain->AddFrame(tabMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

  TGCompositeFrame* tabCont1 = tabMain->AddTab("Time");

  AliCalibViewerGUItime* calibViewerTime = new AliCalibViewerGUItime(tabCont1, 1000, 650);
  tabCont1->AddFrame(calibViewerTime, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  
  TObjArray *guiArray = new TObjArray();
  guiArray->Add(calibViewerTime);

  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();
  return guiArray;
}
