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
//                                                                           //
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
#include "AliTPCCalibViewerGUI.h"
#include "AliTPCCalibViewer.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibDButil.h"
#include "AliTPCConfigParser.h"

#include "AliTPCCalibViewerGUItime.h"

ClassImp(AliTPCCalibViewerGUItime)

AliTPCCalibViewerGUItime::AliTPCCalibViewerGUItime(const TGWindow *p, UInt_t w, UInt_t h) :
TGCompositeFrame(p,w,h),
  fFile(0x0),
  fTree(0x0),
  fCalibViewerGUI(0x0),
  fCalibViewerGUItab(0x0),
  fCurrentHist(0x0),
  fCurrentGraph(0x0),
  fCurrentRunDetails(-1),
  fOutputCacheDir("/tmp"),
  fDrawString(""),
  fConfigFile(""),
  fConfigParser(new AliTPCConfigParser),
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
  fChkDrawOptSparse(0x0),
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
  fContAliases(0x0),
  fListAliases(0x0),
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
  fMapRefTrees->SetOwnerKeyValue();
  fTrashBox->SetOwner();
  DrawGUI(p,w,h);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.5);
  SetInitialValues();
}
//______________________________________________________________________________
AliTPCCalibViewerGUItime::~AliTPCCalibViewerGUItime(){
  //
  // dtor
  //
  delete fConfigParser;
  delete fTrashBox;
  delete fMapRefTrees;
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::DrawGUI(const TGWindow */*p*/, UInt_t w, UInt_t h) {
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
  // draw options
  fContDrawOpt = new TGGroupFrame(fContLeft, "Draw options", kVerticalFrame | kFitWidth | kFitHeight);
  fContLeft->AddFrame(fContDrawOpt, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));

  TGCompositeFrame *cfr = new TGCompositeFrame(fContDrawOpt, 200, 23, kHorizontalFrame | kFitWidth | kFixedHeight);
  fContDrawOpt->AddFrame(cfr, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX , 0, 0, 0, 0));
  
  fChkDrawOptSame = new TGCheckButton(cfr, "Same");
  cfr->AddFrame(fChkDrawOptSame, new TGLayoutHints(kLHintsNormal, 0, 2, 0, 0));
  fChkDrawOptSame->SetToolTipText("Add draw option 'same'");

  fChkDrawOptSparse = new TGCheckButton(cfr, "Sparse");
  cfr->AddFrame(fChkDrawOptSparse, new TGLayoutHints(kLHintsNormal, 0, 2, 0, 0));
  fChkDrawOptSparse->Connect("Clicked()", "AliTPCCalibViewerGUItime", this, "DoNewSelection()");
  fChkDrawOptSparse->SetToolTipText("In case of run trending only plot runs with information");
  
  // additional draw options combo box
  fComboAddDrawOpt = new TGComboBox(fContDrawOpt);
  fComboAddDrawOpt->Resize(0, 22);
  fComboAddDrawOpt->EnableTextInput(kTRUE);
  fContDrawOpt->AddFrame(fComboAddDrawOpt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
//   fComboAddDrawOpt->Connect("ReturnPressed()", "AliTPCCalibViewerGUI", this, "HandleButtonsGeneral(=14)");
//   fComboAddDrawOpt->Connect("Selected(Int_t)", "AliTPCCalibViewerGUI", this, "DoNewSelection()");
  fComboAddDrawOpt->GetTextEntry()->SetText("",kFALSE);
  
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
  fRadioXhist->Connect("Clicked()", "AliTPCCalibViewerGUItime", this, "HandleButtonsDrawSel()");
  fRadioXhist->SetToolTipText("Draw the distribution of the variable");
  
  fRadioXrun = new TGRadioButton(fContDrawSelSubRunTime, ":Run", kRadioXrun);
  fContDrawSelSubRunTime->AddFrame(fRadioXrun, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
  fRadioXrun->Connect("Clicked()", "AliTPCCalibViewerGUItime", this, "HandleButtonsDrawSel()");
  fRadioXrun->SetToolTipText("Use run number as x-value");
  
  // Time radio button
  fRadioXtime = new TGRadioButton(fContDrawSelSubRunTime, ":Time", kRadioXtime);
  fContDrawSelSubRunTime->AddFrame(fRadioXtime, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));
  fRadioXtime->Connect("Clicked()", "AliTPCCalibViewerGUItime", this, "HandleButtonsDrawSel()");
  fRadioXtime->SetToolTipText("Use time stamp number as x-value");
  
  
  // list of variables
  fListVariables = new TGListBox(fContDrawSel);
  fContDrawSel->AddFrame(fListVariables, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  fListVariables->Connect("Selected(Int_t)", "AliTPCCalibViewerGUItime", this, "DoNewSelection()");

  
//-------------------- run type selection ------------------------
  // Parameter label
  fLblRunType = new TGLabel(fContDrawSel, "Run Type:");
  fLblRunType->SetTextJustify(kTextLeft);
  fContDrawSel->AddFrame(fLblRunType, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  
  fComboRunType = new TGComboBox(fContDrawSel);
  fComboRunType->EnableTextInput(kFALSE);
  fContDrawSel->AddFrame(fComboRunType, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  fComboRunType->Connect("Selected(Int_t)", "AliTPCCalibViewerGUItime", this, "DoDraw()");
//   fComboRunType->SetEnabled(kFALSE);
  fComboRunType->Resize(0, 22);
  
  //-------------------- parameter selection ------------------------
  // Parameter label
  fLblPar = new TGLabel(fContDrawSel, "Parameter:");
  fLblPar->SetTextJustify(kTextLeft);
  fContDrawSel->AddFrame(fLblPar, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  
  fNmbPar = new TGNumberEntry(fContDrawSel, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 71);
  fContDrawSel->AddFrame(fNmbPar, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  fNmbPar->Connect("ValueSet(Long_t)", "AliTPCCalibViewerGUItime", this, "DoParLimitChange()");
  fNmbPar->SetState(kFALSE);
  
  //-------------------- calibration type selection ------------------------
  // label
  // draw selection group
  fContCalibType = new TGGroupFrame(fContLeft, "Calib type selection", kVerticalFrame | kFitWidth | kFitHeight);
  fContLeft->AddFrame(fContCalibType, new TGLayoutHints(kLHintsExpandX , 0, 0, 10, 0));
    
    // list of variables
  fListCalibType = new TGListBox(fContCalibType);
  fContCalibType->AddFrame(fListCalibType, new TGLayoutHints(kLHintsNormal | kLHintsExpandX , 0, 0, 0, 0));
  fListCalibType->Connect("Selected(Int_t)", "AliTPCCalibViewerGUItime", this, "DoChangeSelectionList()");
  fListCalibType->Resize(0,88);
  fListCalibType->SetMultipleSelections();
  
  
     // ==========================================================================
   // ************************* content of fContCenter *************************
   // ========================================================================
   // main drawing canvas
  fCanvMain = new TRootEmbeddedCanvas("GUItime_Canvas", fContCenter, 200, 200, kFitWidth | kFitHeight);
  fContCenter->AddFrame(fCanvMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  fCanvMain->GetCanvas()->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)", "AliTPCCalibViewerGUItime", this, "MouseMove(Int_t, Int_t, Int_t, TObject*)");
//   fCanvMain->GetCanvas()->Connect("RangeAxisChanged()", "AliTPCCalibViewerGUItime", this, "GetMinMax()");
  fCanvMain->GetCanvas()->SetToolTipText("The Main_Canvas, here your plots are displayed.");
  fCanvMain->GetCanvas()->SetRightMargin(0.062);
  fCanvMain->GetCanvas()->SetLeftMargin(0.15);
  
   // =========================================================================
   // ************************* content of fContRight *************************
   // ========================================================================
  //group frame for value information
  fContValues = new TGGroupFrame(fContRight, "Data point info", kVerticalFrame | kFitWidth | kFitHeight);
  fContRight->AddFrame(fContValues, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
  //set layout manager
  fContValues->SetLayoutManager(new TGMatrixLayout(fContValues, 0, 2, 5));
  //value information labels

  //run number label
  fLblRunNumber = new TGLabel(fContValues, "Run:");
  fLblRunNumber->SetTextJustify(kTextLeft);
  fContValues->AddFrame(fLblRunNumber, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
  //run number
  fLblRunNumberVal = new TGLabel(fContValues, "000000");
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
  fBtnDumpRuns->Connect("Clicked()", "AliTPCCalibViewerGUItime", this, "DoDumpRuns()");
  fBtnDumpRuns->SetToolTipText("Press to dump the run numbers of the current selection.");
  //group frame for value information
  fContAliases = new TGGroupFrame(fContRight, "Aliases", kVerticalFrame | kFitWidth | kFitHeight);
  fContRight->AddFrame(fContAliases, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 10, 0));
  // list of aliases
  fListAliases = new TGListBox(fContAliases);
  fContAliases->AddFrame(fListAliases, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  fListAliases->Connect("Selected(Int_t)", "AliTPCCalibViewerGUItime", this, "DoNewSelectionAliases()");
//   fListAliases->Resize(0,88);
  //buttons
  TGCompositeFrame *frame1 = new TGCompositeFrame(fContAliases, 200, 23, kHorizontalFrame | kFitWidth | kFitHeight);
  fContAliases->AddFrame(frame1, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX , 0, 0, 0, 0));
  // add button
  TGTextButton *btnAdd = new TGTextButton(frame1, "&Add");
  frame1->AddFrame(btnAdd, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
  btnAdd->Connect("Clicked()", "AliTPCCalibViewerGUItime", this, "DoAddAlias()");
  btnAdd->SetToolTipText("Press to add an alias.");
  btnAdd->Resize(0,22);
  // del button
  TGTextButton *btnDel = new TGTextButton(frame1, "&Del");
  frame1->AddFrame(btnDel, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
  btnDel->Connect("Clicked()", "AliTPCCalibViewerGUItime", this, "DoDelAlias()");
  btnDel->SetToolTipText("Press to delete the selected alias.");
  btnDel->Resize(0,22);
  // update button
  TGTextButton *btnUp = new TGTextButton(frame1, "&Upd");
  frame1->AddFrame(btnUp, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
  btnUp->Connect("Clicked()", "AliTPCCalibViewerGUItime", this, "UpdateAliasList()");
  btnUp->SetToolTipText("Press to update the alias list.");
  btnUp->Resize(0,22);
  

  
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
//   fComboCustomDraw->Resize(0, fLblValueY->GetDefaultHeight());
  fComboCustomDraw->Resize(0, 22);
  fComboCustomDraw->EnableTextInput(kTRUE);
  fContCustom->AddFrame(fComboCustomDraw, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  fComboCustomDraw->Connect("ReturnPressed()", "AliTPCCalibViewerGUItime", this, "DoCustomDraw()");
  fComboCustomDraw->Connect("Selected(Int_t)", "AliTPCCalibViewerGUItime", this, "DoCustomDraw()");
  fComboCustomDraw->GetTextEntry()->SetText("",kFALSE);
  
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
  fComboCustomCuts->Connect("ReturnPressed()", "AliTPCCalibViewerGUItime", this, "DoCustomCutsDraw()");
  fComboCustomCuts->Connect("Selected(Int_t)", "AliTPCCalibViewerGUItime", this, "DoCustomCutsDraw()");
  fComboCustomCuts->GetTextEntry()->SetText("",kFALSE);
  
  SetWindowName("AliTPCCalibViewer GUI - Time");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::SetInitialValues(){
  //
  // Set inital selections of the gui
  //
  fRadioXrun->SetState(kButtonDown);
  fRadioXtime->SetState(kButtonUp);
}

//______________________________________________________________________________
void AliTPCCalibViewerGUItime::UseFile(const char* fileName, const char* treeName) {
  //
  // retrieve tree from file
  //
  TObjArray *arr=0x0;
  TString file(fileName);
  if (file.Contains("://")) {
    if (file.Contains(";")) {
      arr=file.Tokenize(";");
    } else {
      arr=new TObjArray;
      arr->Add(new TObjString(fileName));
    }
  } else {
    TString s=gSystem->GetFromPipe(Form("ls %s",fileName));
    arr=s.Tokenize("\n");
  }

  if (!arr) return;
  TIter next(arr);
  TObject *o=0;
  if (fTree) delete fTree;
  fTree=new TChain(treeName);
  while ( (o=next()) ){
    fTree->AddFile(o->GetName());
  }
  arr->SetOwner();
  delete arr;
  if (!CheckChain()) return;
  UseConfigFile(fConfigFile.Data());
  Reload();
  
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::UseChain(TChain *const chain  = 0)
{
  //
  //
  //
  fTree=chain;
  if (!CheckChain()) return;
  //set configuration file
  UseConfigFile(fConfigFile.Data());
  Reload();
}
//______________________________________________________________________________
Bool_t AliTPCCalibViewerGUItime::CheckChain()
{
  //
  // check whether cahin has entries
  // decide whether to draw graphs in 2D
  //
  if (!fTree) return kFALSE;
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
void AliTPCCalibViewerGUItime::UseConfigFile(const char* file)
{
  //
  // Use 'file' as configuration file
  //
  fConfigFile=file;
  fConfigParser->ParseConfigFileTxt(fConfigFile.Data());
  FillCalibTypes();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::FillRunTypes()
{
  //
  //Loop over the tree entries and fill the run types
  //
  if (!fTree) return;
  Int_t id=0;
  fComboRunType->RemoveAll();
  fComboRunType->AddEntry("ALL",id++);
  fComboRunType->Select(0,kFALSE);
  if (!fTree->GetBranch("runType.")) return;
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
void AliTPCCalibViewerGUItime::FillCalibTypes()
{
  //
  // loop over configuration and fill calibration types
  //
  Int_t id=0;
  fListCalibType->RemoveAll();
  TObject *o=0x0;
  fConfigParser->ResetIter();
  TString type;
  while ( (o=fConfigParser->NextKey()) ){
    type=fConfigParser->GetData(o,kCalibType);
    //remove whitespcaces
    type.Remove(TString::kBoth,' ');
    type.Remove(TString::kBoth,'\t');
    if (type.IsNull()) type="UNSPECIFIED";
//     printf("CalibType: '%s'\n",type.Data());
    if (!fListCalibType->FindEntry(type.Data())) {
      fListCalibType->AddEntry(type.Data(),id);
      fListCalibType->Select(id++);
    }
  }
  //add type for unspecified calibration type
  type="UNSPECIFIED";
  if (!fListCalibType->FindEntry(type.Data())) {
    fListCalibType->AddEntry(type.Data(),id);
    fListCalibType->Select(id++);
  }
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::CheckDrawGraph()
{
  //
  // Check whether to draw graphs in 2D mode based on the number of entries in the chain
  // GetEstimate() returns the maximum size of the arrays stored in GetV1()...
  //
  if (!fTree) return;
  fNoGraph=kTRUE;
  if (fTree->GetEntries()<fTree->GetEstimate()) fNoGraph=kFALSE;
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::Reload(Int_t first)
{
  //
  // reload the gui contents, this is needed after the input tree has changed
  //

  if ( !fTree ) return;
  //in case of the first call create run type and calibration type entries
  if (first){
    FillRunTypes();
    FillCalibTypes();
  }
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
    if (fConfigParser){
      const TObject *key=(*fConfigParser)(branchName.Data());
      if (key){
        //test if branch is active
        active=fConfigParser->GetValue(branchName.Data(),kBranchOnOff);
        id=(*fConfigParser)()->IndexOf(key);
//         branchTitle=fConfigParser->GetData(key,kBranchTitle);
        calibType=fConfigParser->GetData(key,kCalibType);
      }
      else{
        id=1000+idCount;
      }
    } else {
      id=idCount;
    }
    if (calibType.IsNull()) calibType="UNSPECIFIED";
    //check if branch is in selected calibration types
    //if not, don't show it in the list and deactivate the branch.
    Bool_t calibActive=kFALSE;
    TIter nextCalib(&calibTypes);
    TObject *objCalib=0;
    while ( (objCalib=nextCalib()) )
      if (calibType==objCalib->GetTitle()) calibActive=kTRUE;
    active&=calibActive;
    if (!active){
      TString s=branchName;
      if (branchName.EndsWith(".")) s+="*";
//       fTree->SetBranchStatus(s.Data(),0);
      continue;
    }
//     fListVariables->AddEntry(SubstituteUnderscores(branchTitle.Data()),id);
    fListVariables->AddEntry(branchTitle.Data(),id);
    ++idCount;
  }
  //trick to display modifications
  fListVariables->Resize(fListVariables->GetWidth()-1, fListVariables->GetHeight());
  fListVariables->Resize(fListVariables->GetWidth()+1, fListVariables->GetHeight());
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::AddReferenceTree(const char* treeFileName, const char* refName)
{
  //
  // map of reference trees that should always be attached to the CalibViewerGUI
  //
  fMapRefTrees->Add(new TObjString(refName), new TObjString(treeFileName));

}
//______________________________________________________________________________
const TString AliTPCCalibViewerGUItime::GetDrawString(){
  //
  // create draw string for ttree by combining the user requestsa
  //
  
  TString selectedVariable="";
  Int_t id=-1;
  if (!fListVariables->GetSelectedEntry()) return "";
  selectedVariable = fListVariables->GetSelectedEntry()->GetTitle();
  id=fListVariables->GetSelectedEntry()->EntryId();
//   printf("id: %d\n",id);
  TString branchName=selectedVariable;
  if (fConfigParser){
    const TObject *key=(*fConfigParser)(id);
    if (key) branchName=(*fConfigParser)(id)->GetName();
  }
  //treat case of TVector
  if (branchName.EndsWith(".")){
    Int_t par = (Int_t)(fNmbPar->GetNumber());
    branchName.Append(Form("fElements[%d]",par));
  }
//   if (fRadioXrun->GetState()==kButtonDown)
//     selectedVariable.Append(":run");
//   if (fRadioXtime->GetState()==kButtonDown)
//     selectedVariable.Append(":time");
  
  return branchName;
}
//______________________________________________________________________________
const TString AliTPCCalibViewerGUItime::GetDrawOptionString(){
  //
  // get user selected draw options
  //
  TString drawOpt;
  if (fComboAddDrawOpt->GetSelectedEntry()) drawOpt=fComboAddDrawOpt->GetSelectedEntry()->GetTitle();
  if (fChkDrawOptSame->GetState()==kButtonDown && !drawOpt.Contains("same",TString::kIgnoreCase))
    drawOpt+="same";
  return drawOpt;
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::GetCutString(TString &cutStr){
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
void AliTPCCalibViewerGUItime::UpdateValueArrays(Bool_t withGraph, const Double_t *xArr)
{
  //
  //
  //
  if (!withGraph){
    fValuesX.ResizeTo(1);
    fValuesY.ResizeTo(1);
    fRunNumbers.ResizeTo(1);
    fTimeStamps.ResizeTo(1);
  } else {
    const Long64_t nrows=fTree->GetSelectedRows();
    fValuesX.ResizeTo(nrows);
    fValuesY.ResizeTo(nrows);
    fRunNumbers.ResizeTo(nrows);
    fTimeStamps.ResizeTo(nrows);
    long long *index=new long long[nrows];

    //sort data
    Int_t nTime0=0;
    for (Int_t i=0;i<fTree->GetSelectedRows();++i){
      if (fTree->GetV2()[i]<1) ++nTime0;
    }
    
    if (nTime0==fTree->GetSelectedRows()){
        TMath::Sort(nrows,fTree->GetV1(),index,kFALSE);
    } else {
      TMath::Sort(nrows,fTree->GetV2(),index,kFALSE);
    }

    Double_t lastRun=-1.;
    Int_t entries=0;
    const Bool_t drawSparse=(fRadioXrun->GetState()==kButtonDown && fChkDrawOptSparse->GetState()==kButtonDown);
    for (Long64_t i=0; i<nrows; ++i){
      // in case of sparse drawing only use the first entry per run
      Double_t run  = fTree->GetV1()[index[i]];
      Double_t xval = xArr[index[i]];
      
      if (drawSparse){
        if (TMath::Abs(lastRun-run)<.1) {
          lastRun=run;
          continue;
        }
        xval=entries+0.5;
      }
      fValuesX.GetMatrixArray()[entries]=xval;
      fValuesY.GetMatrixArray()[entries]=fTree->GetV3()[index[i]];
      fRunNumbers.GetMatrixArray()[entries]=run;
      fTimeStamps.GetMatrixArray()[entries]=fTree->GetV2()[index[i]];
      lastRun=run;
      ++entries;
    }
    
    if (drawSparse){
      fValuesX.ResizeTo(entries);
      fValuesY.ResizeTo(entries);
      fRunNumbers.ResizeTo(entries);
      fTimeStamps.ResizeTo(entries);
      //      printf("entries: %d\n",entries);
    }
    
    delete [] index;
  }
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::GetHistogramTitle(TString &title)
{
  //
  // Create string for histogram title
  //

  title=fDrawString;
  Int_t pos=title.First(">>");
  if (pos>0) title=title(0,pos);
//   if (!fIsCustomDraw){
    if (fRadioXrun->GetState()==kButtonDown){
      title+=":Run";
    } else if (fRadioXtime->GetState()==kButtonDown){
      title+=":Date";
    }
//   }
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
    //substitue variable names with names in configuration file if available
    if ((*fConfigParser)()->GetEntries()){
      TString branchName=varName;
      Int_t par=0;
      if (branchName.Contains('.')) branchName=branchName(0,branchName.First('.')+1);
      //chek if a configuration for that branch is available
      const TObject *oBranch=(*fConfigParser)(branchName.Data());
      if (oBranch){
        TString branchTitle=fConfigParser->GetData(oBranch,kBranchTitle);
        if (!branchTitle.IsNull()){
          //check for TVectorT type branch
          //add parameter name if available
          if (varName.Contains("fElements")){
            TString parStr=varName(varName.First('[')+1,varName.Length()-varName.First('[')-2);
            par=parStr.Atoi();
            branchTitle+=": ";
            TString yparname=fConfigParser->GetData(oBranch,par+kParamNames);
            if (!yparname.IsNull()){
              branchTitle+=yparname;
            } else {
              branchTitle+="[";
              branchTitle+=par;
              branchTitle+="]";
            }
          }
        }
        varName=branchTitle;
        SubstituteUnderscores(varName);
      }
    }
    title+=varName;
  }
  delete arr;
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::AdjustYRange()
{
  //
  //
  //
  TIter nextGraphicObject(fTrashBox);
  TObject *o=0x0;
  Float_t min=1,max=0;
  while ( (o=nextGraphicObject()) ){
    if (o->IsA()==TGraph::Class()){
      TGraph *gr=(TGraph*)o;
      if (min>max) {
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
void AliTPCCalibViewerGUItime::DoDraw() {
  //
  // Draw graphics
  //
  TString drawString=fDrawString;
  TString cutString;
  GetCutString(cutString);
  TString optString  = GetDrawOptionString();
  Bool_t graphOutput=!fNoGraph;  //ttree buffer for V1, V2... too small
  graphOutput&=(drawString.First(">>")<0); //histogram output in custom draw
  graphOutput&=fRadioXhist->GetState()!=kButtonDown; //histogram drawing selected
//   graphOutput&=!(fIsCustomDraw&&!fDrawString.Contains(":")); //custom draw 1D
//   graphOutput&=fDrawString.CountChar(':')<2; //custom draw 1D
  if (fIsCustomDraw&&fDrawString.Contains(":")) HandleButtonsDrawSel(-kRadioXhist);
  if (fIsCustomDraw&&fDrawString.Contains(":")) HandleButtonsDrawSel(-kRadioXhist);
  Bool_t drawSame=optString.Contains("same",TString::kIgnoreCase);
//   optString+="goff";
  if (graphOutput) {
    drawString.Prepend("run:time:");
    optString="goff";
  }else{
//     if (!fIsCustomDraw){
      if (fRadioXrun->GetState()==kButtonDown){
        drawString+=":run";
      } else if (fRadioXtime->GetState()==kButtonDown){
        drawString+=":time";
      }
//     }
  }
  TVirtualPad *padsave=gPad;
  fCanvMain->GetCanvas()->cd();
  //delete old histograms and graphs
  if (!drawSame){
    fTrashBox->Delete();
    fCurrentGraph=0x0;
    fCurrentHist=0x0;
  }
//   printf("%s (%s) [%s]\n",drawString.Data(), cutString.Data(), optString.Data());
  //select data
  fTree->Draw(drawString.Data(),cutString.Data(),optString.Data());
  if (fTree->GetSelectedRows()==-1) return;
  
  TString title;
  GetHistogramTitle(title);
  Bool_t drawGraph=kFALSE;
  Double_t *xArr=0;
  if (graphOutput){
    drawGraph=kTRUE;
    if (fIsCustomDraw&&fDrawString.Contains(":")){
//       fValuesX.SetElements(fTree->GetV4());
      xArr=fTree->GetV4();
    }else{
      if (fRadioXrun->GetState()==kButtonDown){
//         fValuesX.SetElements(fTree->GetV1());
        xArr=fTree->GetV1();
      } else if (fRadioXtime->GetState()==kButtonDown){
//         fValuesX.SetElements(fTree->GetV2());
        xArr=fTree->GetV2();
      } else {
        drawGraph=kFALSE;
      }
    }
  }
  if (xArr) UpdateValueArrays(graphOutput, xArr);
//   if (graphOutput){
//     if (fIsCustomDraw){
//       if (fDrawString.Contains(":")){
//         fValuesX.SetElements(fTree->GetV4());
//         drawGraph=kTRUE;
//       } else {
//         drawGraph=kFALSE;
//       }
//     }else{
//       drawGraph=kTRUE;
//       if (fRadioXrun->GetState()==kButtonDown){
//         fValuesX.SetElements(fTree->GetV1());
//       } else if (fRadioXtime->GetState()==kButtonDown){
//         fValuesX.SetElements(fTree->GetV2());
//       } else {
//         drawGraph=kFALSE;
//       }
//     }
//   }
//create graph according to selection
  if (drawGraph){
    TGraph *graph=new TGraph(fValuesX,fValuesY);
//     graph->Sort();
    TString grDraw="p";
    if (!drawSame) grDraw+="a";
    if (!fIsCustomDraw) grDraw+="l";
    // sparse drawing, set bin labels
    if (fRadioXrun->GetState()==kButtonDown && fChkDrawOptSparse->GetState()==kButtonDown){
      Int_t nrows=fValuesX.GetNrows();
      Double_t *newBins = new Double_t[nrows+1];
      for (Int_t ibin=0; ibin<nrows+1; ++ibin) newBins[ibin]=ibin;
      graph->GetXaxis()->Set(nrows,newBins);
      graph->GetXaxis()->LabelsOption("v");
      for (Int_t i=0; i<nrows;++i)
        graph->GetXaxis()->SetBinLabel(i+1,Form("%.0f",fRunNumbers.GetMatrixArray()[i]));
      delete [] newBins;
    }
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
//     hist->Draw(optString.Data());
    fTrashBox->Add(hist);
    hist->SetLineColor(fTrashBox->GetEntries());
    hist->SetMarkerColor(fTrashBox->GetEntries());
    if (!drawSame) fCurrentHist=hist;
  }
  
  //Set time axis if choosen as x-variables
//   if (fRadioXtime->GetState()==kButtonDown&&!fIsCustomDraw&&!drawSame){
  if (fRadioXtime->GetState()==kButtonDown&&!drawSame){
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
void AliTPCCalibViewerGUItime::DoDumpRuns()
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
void AliTPCCalibViewerGUItime::DoParLimitChange()
{
  //
  //
  //
  UpdateParName();
  DoDraw();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::DoNewSelection() {
   //
   // decides whether to redraw if user makes another selection
   //
  UpdateParLimits();
  fDrawString=GetDrawString();
  fIsCustomDraw=kFALSE;
  DoDraw();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::DoCustomDraw()
{
  //
  //
  //
  fDrawString=fComboCustomDraw->GetTextEntry()->GetText();
//   if (fDrawString.Contains(">>")){
//     Warning("DoCustomDraw","Currently no user defined histograms allowed!");
//     return;
//   }
  fNmbPar->SetState(kFALSE);
  fIsCustomDraw=kTRUE;
  DoDraw();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::DoCustomCutsDraw()
{
  //
  //
  //
  if (fIsCustomDraw) DoCustomDraw();
  else {
    fDrawString=GetDrawString();
    fIsCustomDraw=kFALSE;
    DoDraw();
  }
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::HandleButtonsDrawSel(Int_t id)
{
  //
  // Draw selection button handling (x-variable)
  //

  if (id == -1) {
    TGButton *btn = (TGButton *) gTQSender;
    id = btn->WidgetId();
  }
  
  Bool_t doDraw=kFALSE;
  Bool_t noDraw=kFALSE;
  if (id<-1){
    noDraw=kTRUE;
    id=TMath::Abs(id);
  }
  switch (id) {
  case (kRadioXhist):
    doDraw=(fRadioXtime->GetState()==kButtonDown||fRadioXrun->GetState()==kButtonDown);
    if (doDraw){
      fRadioXrun->SetState(kButtonUp);
      fRadioXtime->SetState(kButtonUp);
      fRadioXhist->SetState(kButtonDown);
    }
    break;
  case (kRadioXrun):
    doDraw=(fRadioXtime->GetState()==kButtonDown||fRadioXhist->GetState()==kButtonDown);
    if (doDraw){
      fRadioXhist->SetState(kButtonUp);
      fRadioXtime->SetState(kButtonUp);
      fRadioXrun->SetState(kButtonDown);
    }
    break;
  case (kRadioXtime):
    doDraw=(fRadioXhist->GetState()==kButtonDown||fRadioXrun->GetState()==kButtonDown);
    if (doDraw){
      fRadioXrun->SetState(kButtonUp);
      fRadioXhist->SetState(kButtonUp);
      fRadioXtime->SetState(kButtonDown);
    }
    break;
  }
  if (doDraw&&!noDraw) DoCustomCutsDraw();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::UpdateParName()
{
  //
  // change parameter name
  //
  
  Int_t par = (Int_t)(fNmbPar->GetNumber());
  TString parName="";
  Int_t id=fListVariables->GetSelectedEntry()->EntryId();
  if (fConfigParser && (*fConfigParser)(id)) parName=fConfigParser->GetData((*fConfigParser)(id),par+kParamNames);
  if (parName=="") parName.Form("%d",par);
  fLblPar->SetText(Form("Parameter: %s",parName.Data()));
  fDrawString=GetDrawString();
  fIsCustomDraw=kFALSE;
}

//______________________________________________________________________________
void AliTPCCalibViewerGUItime::UpdateParLimits()
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
//   printf("id: %d\n",id);
  TString selectedVariable=selectedVariableTitle;
  const TObject *key=(*fConfigParser)(id);
  if (key) selectedVariable=(*fConfigParser)(id)->GetName();
  
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
//   branch->ResetAddress();
  fTree->SetBranchAddress(selectedVariable.Data(),0x0);
  if (fNmbPar->GetNumMax()!=maxPar-1) fNmbPar->SetNumber(0);
  fNmbPar->SetLimitValues(0,maxPar-1);
  fNmbPar->SetState(kTRUE);
  UpdateParName();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::MouseMove(Int_t event, Int_t x, Int_t y, TObject */*selected*/)
{
  //
  // handle mouse events in the draw canvas
  //
  UInt_t dd=0,mm=0,yy=0,HH=0,MM=0,SS=0,run=0;
  Double_t valx=0.,valy=0.;
  if (!fCurrentGraph) {
    fLblRunNumberVal->SetText(Form("%06u",run));
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
  fLblRunNumberVal->SetText(Form("%06u",run));
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
  if (run==0) return;
  if (event == kButton1Double ){
    SetGuiTree(run);
  }
  //find closes point of current selection
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::SetGuiTree(Int_t run)
{
  //
  // create the AliTPCCalibViewerGUI tree for run
  // cache tree in directory fOutputCacheDir
  // retrieve file from this directory if it already exists
  //

  //
  //Create and set GUI tree  
  //
  //try to find file for run in fOutputCacheDir
  TString fileName=fOutputCacheDir;
  if (!fileName.EndsWith("/")) fileName+="/";
  fileName+=Form("guiTreeRun_%d.root",run);
  Bool_t load=kTRUE;
  if (gSystem->AccessPathName(fileName.Data())){
    load=AliTPCcalibDB::CreateGUITree(run,fileName.Data());
    if (!load){
      fCalibViewerGUI->Reset();
      if (fCalibViewerGUItab) fCalibViewerGUItab->SetText(new TGString(Form("Detail - XXXXX")));
      return;
    }
  }
  fCalibViewerGUI->Initialize(fileName.Data());
  if (fCalibViewerGUItab) fCalibViewerGUItab->SetText(new TGString(Form("Detail - %05d",run)));
  
  //
  //Create and set Reference GUI tree
  //
  AliTPCcalibDButil util;
  util.SetReferenceRun(run);
  fileName=fOutputCacheDir;
  if (!fileName.EndsWith("/")) fileName+="/";
  fileName+=util.GetGUIRefTreeDefaultName();
  //only update if file does not exist
  if (gSystem->AccessPathName(fileName.Data())){
    util.UpdateRefDataFromOCDB();
    util.CreateGUIRefTree(fileName.Data());
  }
  
  fCalibViewerGUI->GetViewer()->AddReferenceTree(fileName.Data(),"calPads","Ref");
  
  //
  // Process additional reference trees
  //
  TIter nextRefTree(fMapRefTrees);
  TObject *o=0x0;
  //Set static reference data
  while ( (o=nextRefTree()) ){
    fCalibViewerGUI->GetViewer()->AddReferenceTree(fMapRefTrees->GetValue(o)->GetName(),"calPads",o->GetName());
  }
  fCalibViewerGUI->Reload();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::SubstituteUnderscores(TString &s)
{
  //
  //
  //
  s.ReplaceAll("_{","|{");
  s.ReplaceAll("_"," ");
  s.ReplaceAll("|{","_{");
}

//______________________________________________________________________________
void AliTPCCalibViewerGUItime::DoNewSelectionAliases()
{
  //
  //
  //
  if (!fTree) return;
  TList *l=fTree->GetListOfAliases();
  if (!l) return;
  TString selectedVariable="";
  if (!fListAliases->GetSelectedEntry()) return;
  selectedVariable = fListAliases->GetSelectedEntry()->GetTitle();
  fDrawString=selectedVariable;
  fIsCustomDraw=kFALSE;
  DoDraw();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::DoAddAlias()
{
  //
  //
  //
  new AliTPCCalibViewerGUItimeAddAliasFrame(gClient->GetRoot(), fContTopBottom, 400, 200, kVerticalFrame, this);
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::DoDelAlias()
{
  //
  //
  //
  if (!fTree) return;
  TList *l=fTree->GetListOfAliases();
  if (!l) return;
  TString selectedVariable="";
  if (!fListAliases->GetSelectedEntry()) return;
  selectedVariable = fListAliases->GetSelectedEntry()->GetTitle();
  l->Remove(l->FindObject(selectedVariable.Data()));
  UpdateAliasList();
}

//______________________________________________________________________________
void AliTPCCalibViewerGUItime::UpdateAliasList()
{
  //
  //
  //
  printf("UpdateAliasList\n");
  if (!fTree) return;
  TList *l=fTree->GetListOfAliases();
  if (!l) return;
  TIter nextAlias(l);
  TObject *o;
  fListAliases->RemoveAll();
  Int_t id=0;
  while( (o=nextAlias()) ){
    fListAliases->AddEntry(o->GetName(),id++);
  }
  fListAliases->Resize(fListAliases->GetWidth()-1, fListAliases->GetHeight());
  fListAliases->Resize(fListAliases->GetWidth()+1, fListAliases->GetHeight());
}
//______________________________________________________________________________
TObjArray* AliTPCCalibViewerGUItime::ShowGUI(const char* fileName, const char* treeName) {
   //
   // Initialize and show GUI for presentation for demonstration purposes
   // or for fast standalone use
   //
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("AliTPCCalibViewer GUItime");
  frmMain->SetCleanup(kDeepCleanup);
  
  TGTab* tabMain = new TGTab(frmMain, 1000, 600);
  frmMain->AddFrame(tabMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  
  TGCompositeFrame* tabCont1 = tabMain->AddTab("Time");
  TGCompositeFrame* tabCont2 = tabMain->AddTab("Detail - XXXXX");
  
  AliTPCCalibViewerGUItime* calibViewerTime = new AliTPCCalibViewerGUItime(tabCont1, 1000, 600);
  tabCont1->AddFrame(calibViewerTime, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  calibViewerTime->SetConfigFileName("$ALICE_ROOT/TPC/CalibMacros/calibVarDescription.txt");
  calibViewerTime->UseFile(fileName, treeName);
  
  AliTPCCalibViewerGUI* calibViewer = new AliTPCCalibViewerGUI(tabCont2, 1000, 600, 0);
  tabCont2->AddFrame(calibViewer, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  calibViewerTime->SetCalibViewerGUI(calibViewer);
  calibViewerTime->SetCalibViewerGUItab(tabMain->GetTabTab(1));
  
  
  TObjArray *guiArray = new TObjArray();
  guiArray->Add(calibViewerTime);
  guiArray->Add(calibViewer);
  
  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();
  
  return guiArray;
}

//______________________________________________________________________________
TObjArray* AliTPCCalibViewerGUItime::ShowGUI(TChain *chain) {
   //
   // Initialize and show GUI for presentation for demonstration purposes
   // or for fast standalone use
   //
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("AliTPCCalibViewer GUItime");
  frmMain->SetCleanup(kDeepCleanup);

  TGTab* tabMain = new TGTab(frmMain, 1000, 600);
  frmMain->AddFrame(tabMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

  TGCompositeFrame* tabCont1 = tabMain->AddTab("Time");

  AliTPCCalibViewerGUItime* calibViewerTime = new AliTPCCalibViewerGUItime(tabCont1, 1000, 600);
  tabCont1->AddFrame(calibViewerTime, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  calibViewerTime->UseChain(chain);

  TObjArray *guiArray = new TObjArray();
  guiArray->Add(calibViewerTime);

  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();

  return guiArray;
}


////////////////////////////////////////////////////////////////////////
//
//   GUI Alias frame
//
////////////////////////////////////////////////////////////////////////


ClassImp(AliTPCCalibViewerGUItimeAddAliasFrame)

AliTPCCalibViewerGUItimeAddAliasFrame::AliTPCCalibViewerGUItimeAddAliasFrame(const TGWindow *p, const TGWindow *main,
                                                                             UInt_t w, UInt_t h, UInt_t options,
                                                                             AliTPCCalibViewerGUItime *gui, TString strAlias) :
  fMain(0x0),
  fTxt1(0x0),
  fTxt2(0x0),
  fGUI(0x0)
{
  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()", "AliTPCCalibViewerGUItimeAddAliasFrame", this, "DoCancel()");
  fMain->DontCallClose(); // to avoid double deletions.
  
   // use hierarchical cleaning
  fMain->SetCleanup(kDeepCleanup);

  //layout
  TGLayoutHints *l1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2);
  TGLayoutHints *l2 = new TGLayoutHints(kLHintsTop | kLHintsRight | kLHintsExpandX, 2, 2, 0, 5);
//   TGLayoutHints *l3 = new TGLayoutHints(kLHintsTop | kLHintsRight, 5, 5, 5, 5);
  
  //input fields
  TGCompositeFrame *f1 = new TGCompositeFrame(fMain, 60, 20, kVerticalFrame);
  fMain->AddFrame(f1, l1);
  TGCompositeFrame *frameName    = new TGCompositeFrame(f1);
  TGCompositeFrame *frameFormula = new TGCompositeFrame(f1);
  f1->AddFrame(frameName,l2);
  f1->AddFrame(frameFormula,l2);
  TGLabel *lblTxt1   = new TGLabel(frameName, "Name:");
  TGLabel *lblTxt2   = new TGLabel(frameFormula, "Formula:");
  fTxt1 = new TGTextEntry(frameName, new TGTextBuffer(1000));
  fTxt2 = new TGTextEntry(frameFormula, new TGTextBuffer(1000));
  
  frameName->AddFrame(lblTxt1, l2);
  frameName->AddFrame(fTxt1, l2);
  frameFormula->AddFrame(lblTxt2, l2);
  frameFormula->AddFrame(fTxt2, l2);
  
  fTxt1->Resize(350, fTxt1->GetDefaultHeight());
  fTxt2->Resize(350, fTxt2->GetDefaultHeight());

  //ok and cancel buttons
  TGHorizontalFrame *frame = new TGHorizontalFrame(fMain, 60, 20, kFixedWidth);
  
  TGTextButton *okButton = new TGTextButton(frame, "&Ok", 1);
  okButton->Connect("Clicked()", "AliTPCCalibViewerGUItimeAddAliasFrame", this, "DoOK()");
  TGTextButton *cancelButton = new TGTextButton(frame, "&Cancel", 2);
  cancelButton->Connect("Clicked()", "AliTPCCalibViewerGUItimeAddAliasFrame", this, "DoCancel()");
  
  
  frame->AddFrame(okButton, l1);
  frame->AddFrame(cancelButton, l1);
  
  frame->Resize(150, okButton->GetDefaultHeight());
  
  fMain->AddFrame(frame, new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1));
  
  fGUI=gui;
  TString aliasName, alias;
  if (!strAlias.IsNull()){
    TChain *c=fGUI->GetChain();
    if (c){
      TList *l=c->GetListOfAliases();
      if (l){
        TNamed *d=(TNamed*)l->FindObject(strAlias);
        if (d){
          aliasName=d->GetName();
          alias=d->GetTitle();
        }
      }
    }
  }else{
    alias=fGUI->GetCustomDrawString();
  }
  fTxt1->SetText(aliasName.Data(),kFALSE);
  fTxt2->SetText(alias.Data(),kFALSE);

  fMain->MapSubwindows();
  fMain->Resize();
  
   // position relative to the parent's window
  fMain->CenterOnParent();
  
  fMain->SetWindowName("Alias Editor");
  
  fMain->MapWindow();
  
}
//______________________________________________________________________________
AliTPCCalibViewerGUItimeAddAliasFrame::~AliTPCCalibViewerGUItimeAddAliasFrame()
{
  //
  //
  //
  fMain->DeleteWindow();  // deletes fMain
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItimeAddAliasFrame::DoOK()
{
  //
  //
  //
  TString aliasName=fTxt1->GetText();
  TString alias=fTxt2->GetText();
  if (!aliasName.IsNull()&&!alias.IsNull()){
    TChain *c=fGUI->GetChain();
    if (c){
      c->SetAlias(aliasName.Data(),alias.Data());
    }
  }
  fGUI->UpdateAliasList();
  delete this;
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItimeAddAliasFrame::DoCancel()
{
  //
  //
  //
  delete this;
}

