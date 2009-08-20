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
#include <TVector.h>
#include <TH1.h>
#include <TCut.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TIterator.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TTimeStamp.h>
#include <TMath.h>
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
#include "AliTPCCalibViewerGUI.h"
#include "AliTPCcalibDB.h"

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
  fIsCustomDraw(kFALSE),
  fRunNumbers(10),
  fTimeStamps(10),
  fValuesX(10),
  fValuesY(10),
  //GUI elements
  //main canvas Top part, bottom part
  fContTopBottom(0x0),
  //top left, centre, right
  fContLCR(0x0),
  //content left
  fContLeft(0x0),
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
  //content bottom
  fContCustom(0x0),
  fContCustomCuts(0x0),
  fLblCustomDraw(0x0),
  fLblCustomCuts(0x0),
  fComboCustomDraw(0x0),
  fComboCustomCuts(0x0)
{
  //
  // ctor
  //
  DrawGUI(p,w,h);
  SetInitialValues();
}
//______________________________________________________________________________
AliTPCCalibViewerGUItime::~AliTPCCalibViewerGUItime(){
  //
  // dtor
  //

}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::DrawGUI(const TGWindow *p, UInt_t w, UInt_t h) {
   //
   // draw the GUI
   //
   // ======================================================================
   // ************************* Display everything *************************
   // ======================================================================
  
  SetCleanup(kDeepCleanup);
  p = p; // to avoid compiler warnings
  
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

  //-------------------- rynType selection ------------------------
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
  
  
     // ==========================================================================
   // ************************* content of fContCenter *************************
   // ========================================================================
   // main drawing canvas
  fCanvMain = new TRootEmbeddedCanvas("GUItime_Canvas", fContCenter, 200, 200, kFitWidth | kFitHeight);
  fContCenter->AddFrame(fCanvMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  fCanvMain->GetCanvas()->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)", "AliTPCCalibViewerGUItime", this, "MouseMove(Int_t, Int_t, Int_t, TObject*)");
//   fCanvMain->GetCanvas()->Connect("RangeAxisChanged()", "AliTPCCalibViewerGUItime", this, "GetMinMax()");
  fCanvMain->GetCanvas()->SetToolTipText("The Main_Canvas, here your plots are displayed.");
  
  
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
  fLblRunNumberVal = new TGLabel(fContValues, "00000");
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
  fLblValueXVal = new TGLabel(fContValues, "00000.000");
  fLblValueXVal->SetTextJustify(kTextRight);
  fContValues->AddFrame(fLblValueXVal, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  //value label y
  fLblValueY = new TGLabel(fContValues, "y-Value:");
  fLblValueY->SetTextJustify(kTextLeft);
  fContValues->AddFrame(fLblValueY, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
  //value y
  fLblValueYVal = new TGLabel(fContValues, "00000.000");
  fLblValueYVal->SetTextJustify(kTextRight);
  fContValues->AddFrame(fLblValueYVal, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
   // draw button
  fBtnDumpRuns = new TGTextButton(fContRight, "&Dump runs");
  fContRight->AddFrame(fBtnDumpRuns, new TGLayoutHints(kLHintsExpandX, 10, 10, 0, 0));
  fBtnDumpRuns->Connect("Clicked()", "AliTPCCalibViewerGUItime", this, "DoDumpRuns()");
  fBtnDumpRuns->SetToolTipText("Press to dump the run numbers of the current selection.");
  
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
  fComboCustomCuts->Connect("ReturnPressed()", "AliTPCCalibViewerGUItime", this, "DoCustomDraw()");
  fComboCustomCuts->Connect("Selected(Int_t)", "AliTPCCalibViewerGUItime", this, "DoCustomDraw()");

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
void AliTPCCalibViewerGUItime::UseFile(const char* fileName) {
  //
  // retrieve tree from file
  //
  TDirectory *save=gDirectory;
  if (fFile) delete fFile;
  fFile = TFile::Open(fileName);
  fTree=(TTree*)fFile->Get("dcs");
  save->cd();
  Reload();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::Reload(){
  //
  // reload the gui contents, this is needed after the input tree has changed
  //

  //reset variables list
  fListVariables->RemoveAll();
  
  if ( !fTree ) return;
  TObjArray *branchList = fTree->GetListOfBranches();
  if ( !branchList ) return;
  TIter nextBranch(branchList);
  Int_t id=0;
  TObject *objBranch=0;
  while ( (objBranch=nextBranch()) ){
    TString branchName(objBranch->GetName());
    if (branchName == "run" || branchName == "time" || branchName == "runType.") continue;
    fListVariables->AddEntry(branchName.Data(),id++);
  }
  //create entris for run types
  id=0;
  fComboRunType->RemoveAll();
  fComboRunType->AddEntry("ALL",id++);
  fComboRunType->Select(0,kFALSE);
  TObjString *runType=0x0;
  Int_t nevets=fTree->GetEntries();
  TBranch *branch=fTree->GetBranch("runType.");
  if (!branch) return;
  branch->SetAddress(&runType);
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("runType.*",1);
  for (Int_t iev=0;iev<nevets;++iev){
    fTree->GetEntry(iev);
    TString type=runType->String();
    if (!type.IsNull()&&!fComboRunType->FindEntry(type)) fComboRunType->AddEntry(type,id++);
  }
  branch->ResetAddress();
  fTree->SetBranchStatus("*",1);
}
//______________________________________________________________________________
const char* AliTPCCalibViewerGUItime::GetDrawString(){
  //
  // create draw string for ttree by combining the user requestsa
  //
  
  TString selectedVariable="";
  if (fListVariables->GetSelectedEntry()) selectedVariable = fListVariables->GetSelectedEntry()->GetTitle();
  //treat case of TVector
  if (selectedVariable.EndsWith(".")){
    Int_t par = (Int_t)(fNmbPar->GetNumber());
    selectedVariable.Append(Form("fElements[%d]",par));
  }
//   if (fRadioXrun->GetState()==kButtonDown)
//     selectedVariable.Append(":run");
//   if (fRadioXtime->GetState()==kButtonDown)
//     selectedVariable.Append(":time");
  
  return selectedVariable.Data();
}
//______________________________________________________________________________
const char* AliTPCCalibViewerGUItime::GetCutString(){
  //
  // create cut string
  //
  TCut cuts(fComboCustomCuts->GetTextEntry()->GetText());
  TString runType="";
  if (fComboRunType->GetSelectedEntry()) runType=fComboRunType->GetSelectedEntry()->GetTitle();
  if (runType!="ALL"&&!runType.IsNull()) cuts+=Form("runType.String().Data()==\"%s\"",runType.Data());
  return cuts.GetTitle();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::DoDraw() {
  TString drawString=fDrawString;
  drawString.Prepend("run:time:");
  TString cutString  = GetCutString();
  TString optString  = "goff";
  TVirtualPad *padsave=gPad;
  fCanvMain->GetCanvas()->cd();
  //delete old histogram and graph
  if (fCurrentGraph) {
    delete fCurrentGraph;
    fCurrentGraph=0x0;
    //fCurrentHist in case of graph is the internal histogram,
    //  which is deletet by the graph itself.
    fCurrentHist=0x0;
  }
  if (fCurrentHist)  delete fCurrentHist;
  //select data
  fTree->Draw(drawString.Data(),cutString.Data(),optString.Data());
  if (fTree->GetSelectedRows()==-1) return;
  fValuesX.ResizeTo(fTree->GetSelectedRows());
  fValuesY.ResizeTo(fTree->GetSelectedRows());
  fRunNumbers.ResizeTo(fTree->GetSelectedRows());
  fTimeStamps.ResizeTo(fTree->GetSelectedRows());
  fValuesY.SetElements(fTree->GetV3());
  fRunNumbers.SetElements(fTree->GetV1());
  fTimeStamps.SetElements(fTree->GetV2());
  TString title="";
  Bool_t drawGraph=kFALSE;
  if (fIsCustomDraw){
    if (fDrawString.Contains(":")){
      fValuesX.SetElements(fTree->GetV4());
      TString yname=fDrawString(0,fDrawString.First(':'));
      TString xname=fDrawString(fDrawString.First(':')+1,fDrawString.Length());
      title=Form("%s;%s;%s",fDrawString.Data(),xname.Data(),yname.Data());
      drawGraph=kTRUE;
    } else {
      drawGraph=kFALSE;
    }
  }else{
    drawGraph=kTRUE;
    if (fRadioXrun->GetState()==kButtonDown){
      fValuesX.SetElements(fTree->GetV1());
      title=Form("%s:Run;Run;%s",fDrawString.Data(),fDrawString.Data());
    } else if (fRadioXtime->GetState()==kButtonDown){
      fValuesX.SetElements(fTree->GetV2());
      title=Form("%s:Time;Time;%s",fDrawString.Data(),fDrawString.Data());
    } else {
      drawGraph=kFALSE;
    }
  }
  //create graph according to selection
  if (drawGraph){
    fCurrentGraph=new TGraph(fValuesX,fValuesY);
    fCurrentGraph->Draw("alp");
    fCurrentHist=fCurrentGraph->GetHistogram();
    fCurrentHist->SetTitle(title.Data());
  } else {
    fCurrentGraph=0x0;
    Float_t add=TMath::Abs(fValuesY.Min()*.05);
    fCurrentHist=new TH1D("hist",Form("%s;%s",fDrawString.Data(),fDrawString.Data()),100,fValuesY.Min()-add,fValuesY.Max()+add);
    fCurrentHist->FillN(fValuesY.GetNrows(),fValuesY.GetMatrixArray(),0);
    fCurrentHist->Draw();
  }
  
//   fCurrentHist=fTree->GetHistogram();
//   fCurrentGraph=(TGraph*)gPad->GetPrimitive("Graph");
  //Set time axis if choosen as x-variables
  if (fRadioXtime->GetState()==kButtonDown&&!fIsCustomDraw){
    TAxis *xaxis=fCurrentHist->GetXaxis();
    xaxis->SetTimeFormat("#splitline{%d.%m}{%H:%M}");
    xaxis->SetTimeDisplay(1);
    xaxis->SetLabelOffset(xaxis->GetLabelOffset()*3);
    xaxis->SetLabelSize(xaxis->GetLabelSize()/1.3);
  }
  if (fCurrentGraph){
    fCurrentGraph->SetEditable(kFALSE);
    fCurrentGraph->SetMarkerStyle(20);
    fCurrentGraph->SetMarkerSize(0.5);
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
  for (Int_t irun=0;irun<fRunNumbers.GetNrows();++irun){
    std::cout << Form("%.0f",fRunNumbers[irun]) << std::endl;
  }
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::DoParLimitChange()
{
  //
  //
  //
  Int_t par = (Int_t)(fNmbPar->GetNumber());
  fLblPar->SetText(Form("Parameter: %02d",par));
  fDrawString=GetDrawString();
  fIsCustomDraw=kFALSE;
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
  if (fDrawString.Contains(">>")){
    Warning("DoCustomDraw","Currently no user defined histograms allowed!");
    return;
  }
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
  if (doDraw) DoNewSelection();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::UpdateParLimits()
{
  //
  // Adjust limits for TVectorT based variables
  //
  if (!fTree) return;
  TString selectedVariable="";
  if (fListVariables->GetSelectedEntry()) selectedVariable = fListVariables->GetSelectedEntry()->GetTitle();
  if (selectedVariable.IsNull()||!selectedVariable.EndsWith(".")) {
    fNmbPar->SetState(kFALSE);
    return;
  }
  TVectorD *vD=0x0;
  TVectorF *vF=0x0;
  Int_t maxPar=0;
  TBranch *branch=fTree->GetBranch(selectedVariable.Data());
  TString branchClass=branch->GetClassName();
  if (branchClass=="TVectorT<double>"){
    branch->SetAddress(&vD);
    fTree->GetEntry(0);
    maxPar=vD->GetNrows();
  } else if (branchClass=="TVectorT<float>"){
    branch->SetAddress(&vF);
    fTree->GetEntry(0);
    maxPar=vF->GetNrows();
  } else {
    //class not known
    fNmbPar->SetState(kFALSE);
    return;
  }
  branch->ResetAddress();
  fNmbPar->SetNumber(0);
  fNmbPar->SetLimitValues(0,maxPar-1);
  fNmbPar->SetState(kTRUE);
}
//______________________________________________________________________________
void AliTPCCalibViewerGUItime::MouseMove(Int_t event, Int_t x, Int_t y, TObject *selected)
{
  //
  // handle mouse events in the draw canvas
  //
  
  //avoid compiler warnings
  selected=selected;
  UInt_t dd=0,mm=0,yy=0,HH=0,MM=0,SS=0,run=0;
  Double_t valx=0.,valy=0.;
  if (!fCurrentGraph) {
    fLblRunNumberVal->SetText(Form("%05u",run));
    fLblRunTimeVal->SetText(Form("%02u.%02u.%04u\n%02u.%02u.%02u",dd,mm,yy,HH,MM,SS));
    fLblValueXVal->SetText(Form("%.3f", valx));
    fLblValueYVal->SetText(Form("%.3f", valy));
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
    Int_t d   = (Int_t)TMath::Sqrt(TMath::Abs(pxp-x) + TMath::Abs(pyp-y));
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
  fLblRunNumberVal->SetText(Form("%05u",run));
  fLblRunTimeVal->SetText(Form("%02u.%02u.%04u\n%02u.%02u.%02u",dd,mm,yy,HH,MM,SS));
  if (fIsCustomDraw){
    fLblValueXVal->SetText(Form("%.3f", valx));
  }else{
    if (fRadioXrun->GetState()==kButtonDown){
      fLblValueXVal->SetText("Run");
    } else if (fRadioXtime->GetState()==kButtonDown){
      fLblValueXVal->SetText("Time");
    }
  }
  fLblValueYVal->SetText(Form("%.3f", valy));
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

  //try to find file for run in fOutputCacheDir
  TString fileName=fOutputCacheDir;
  if (!fileName.EndsWith("/")) fileName+="/";
  fileName+=Form("guiTreeRun_%d.root",run);
  TFile f(fileName.Data());
  if (f.IsOpen()){
    f.Close();
    fCalibViewerGUI->Initialize(fileName.Data());
    if (fCalibViewerGUItab) fCalibViewerGUItab->SetText(new TGString(Form("Detail - %05d",run)));
    return;
  }
  f.Close();
  Bool_t sucess=AliTPCcalibDB::CreateGUITree(run,fileName.Data());
  if (sucess){
    fCalibViewerGUI->Initialize(fileName.Data());
    if (fCalibViewerGUItab) fCalibViewerGUItab->SetText(new TGString(Form("Detail - %05d",run)));
  }else{
    fCalibViewerGUI->Reset();
    if (fCalibViewerGUItab) fCalibViewerGUItab->SetText(new TGString(Form("Detail - XXXXX")));
  }
}
//______________________________________________________________________________
TObjArray* AliTPCCalibViewerGUItime::ShowGUI(const char* fileName) {
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
  calibViewerTime->UseFile(fileName);
  
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

