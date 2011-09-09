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

/* $Id: AliTRDCalibViewerGUI.cxx 40390 2010-04-14 09:43:23Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  GUI for the AliTRDCalibViewerGUI                                         //
//  used for the calibration monitor                                         //
//  All functionalities of the AliTRDCalibViewer are here available          //
//                                                                           //
//  Authors:     Marian Ivanov (Marian.Ivanov@cern.ch)                       //
//               Jens Wiechula (Jens.Wiechula@cern.ch)                       //
//               Ionut Arsene  (iarsene@cern.ch)                             //
//                                                                           //
//  Example usage:                                                           //
/*
  aliroot
  AliTRDCalibViewerGUI::ShowGUI()
  
*/
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTRDCalibViewerGUI.h"
#include <iostream>
#include <string.h>

#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPad.h>

#include <TROOT.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TVector.h>
#include <TH1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGFileDialog.h>
#include <TGInputDialog.h>
#include "AliBaseCalibViewer.h"
#include "AliCalibViewerGUItime.h"
#include "AliTRDCalibViewer.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"

ClassImp(AliTRDCalibViewerGUI)

//________________________________________________________________________________________________
AliTRDCalibViewerGUI::AliTRDCalibViewerGUI(const TGWindow *p, UInt_t w, UInt_t h, char* fileName)
: AliBaseCalibViewerGUI(p, w, h),
    fContLayer(0),
    fContSector(0),
    fContStack(0),
    fLblLayer(0),
    fLblSector(0),
    fLblStack(0),
    fNmbLayer(0),
    fNmbSector(0),
    fNmbStack(0),
    fContLoad(0),
    fContRun(0),
    fLblRun(0),
    fNmbRun(0),
    fContStorage(0),
    fLblStorage(0),
    fTxtStorage(0),
    fContVersion(0),	
    fLblVersion(0),
    fNmbVersion(0),
    fContSubVersion(0),	
    fLblSubVersion(0),
    fNmbSubVersion(0),
    fContChecks(0),
    fChkCalibs(0),
    fChkDCS(0),
    fChkAlign(0),
    fBtnLoad(0),
    fContLoadCalibObjects(0),
    fContCalibInput(0),
    fLblCalibInputFilename(0),
    fTxtCalibInputFilename(0),
    fContCalibOutput(0),
    fLblCalibOutputFilename(0),
    fTxtCalibOutputFilename(0),
    fBtnLoadCalibObjects(0)
{
   //
   // AliTRDCalibViewerGUI constructor; fileName specifies the ROOT tree used for drawing 
   //
   // draw the GUI:
   DrawGUI(p, w, h);
   // initialize
   if (fileName) Initialize(fileName, "TRDcalibDetails");
   // set default button states:
   SetInitialValues();
   // do first drawing: 
   if (fileName) DoDraw();
}

//________________________________________________________________________________________________
AliTRDCalibViewerGUI::AliTRDCalibViewerGUI(const AliTRDCalibViewerGUI &c)
  :AliBaseCalibViewerGUI(c),
   fContLayer(0),
   fContSector(0),
   fContStack(0),
   fLblLayer(0),
   fLblSector(0),
   fLblStack(0),
   fNmbLayer(0),
   fNmbSector(0),
   fNmbStack(0),
   fContLoad(0),
   fContRun(0),
   fLblRun(0),
   fNmbRun(0),
   fContStorage(0),
   fLblStorage(0),
   fTxtStorage(0),
   fContVersion(0),	
   fLblVersion(0),
   fNmbVersion(0),
   fContSubVersion(0),	
   fLblSubVersion(0),
   fNmbSubVersion(0),
   fContChecks(0),
   fChkCalibs(0),
   fChkDCS(0),
   fChkAlign(0),
   fBtnLoad(0),
   fContLoadCalibObjects(0),
   fContCalibInput(0),
   fLblCalibInputFilename(0),
   fTxtCalibInputFilename(0),
   fContCalibOutput(0),
   fLblCalibOutputFilename(0),
   fTxtCalibOutputFilename(0),
   fBtnLoadCalibObjects(0) 
{
  //
  // dummy AliTPCCalibViewerGUI_new copy constructor
  //
}

//________________________________________________________________________________________________
AliTRDCalibViewerGUI & AliTRDCalibViewerGUI::operator =(const AliTRDCalibViewerGUI & /*param*/) {
   //
   // dummy assignment operator
   //
   return (*this);
}

//________________________________________________________________________________________________
AliTRDCalibViewerGUI::~AliTRDCalibViewerGUI() {
   // 
   // Destructor
   // 
  /*
  if (fCanvMain && fCanvMain->GetCanvas()) {
    for (Int_t i = 0; i < fCanvMain->GetCanvas()->GetListOfPrimitives()->GetEntries(); i++) {
      if (strcmp(fCanvMain->GetCanvas()->GetListOfPrimitives()->At(i)->ClassName(), "TFrame") != 0)
	fCanvMain->GetCanvas()->GetListOfPrimitives()->At(i)->Delete();
    }
  } */
  Cleanup();
  if (fViewer) fViewer->Delete();
}

//________________________________________________________________________________________________
Bool_t AliTRDCalibViewerGUI::CreateDetailsTree(Int_t run, const Char_t* outFile, const Char_t* /*ocdbStorage*/) {
  //
  // Get pad level info from OCDB for a given run and dump it into a tree
  // 
  if(!AliCDBManager::Instance()->GetDefaultStorage()){
    std::cout << "AliTRDCalibViewerGUI::CreateDetailsTree(): Default Storage not set. Cannot create Calibration Tree!" << std::endl;
    return kFALSE;
  }
  TString storage = AliCDBManager::Instance()->GetDefaultStorage()->GetURI();
  return ((AliTRDCalibViewer*)fViewer)->DumpOCDBtoTreeDetails("", outFile, run, run, storage.Data());
}

//________________________________________________________________________________________________
void AliTRDCalibViewerGUI::DrawGUI(const TGWindow *p, UInt_t w, UInt_t h) {
  // 
  // draw the GUI
  // 
   
  // draw most of the GUI (all common stuff)
  AliBaseCalibViewerGUI::DrawGUI(p, w, h);
  
  // remove some frames from the virtual class
  fTabRight1->RemoveFrame(fContExport);
  fTabRight1->RemoveFrame(fContTree);
  fTabRight1->RemoveFrame(fContFit);

  // draw and connect slots specific to TRD
  // **************************** content of tabLeft0 *******************************
  // layer options container
  fContLayer = new TGCompositeFrame(fContCuts, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
  fContCuts->AddFrame(fContLayer, new TGLayoutHints(kLHintsExpandX, 5, 0, 0, 0));

    // layer number label
    fLblLayer = new TGLabel(fContLayer, "Layer");
    fLblLayer->SetTextJustify(kTextLeft);
    fContLayer->AddFrame(fLblLayer, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 5, 0, 0, 0));
    // layer number entry
    fNmbLayer = new TGNumberEntry(fContLayer, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 5);
    fContLayer->AddFrame(fNmbLayer, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
    fNmbLayer->SetNumber(0);
    fNmbLayer->Connect("ValueSet(Long_t)", "AliTRDCalibViewerGUI", this, "DoNewSelection()");

  // sector options container
  fContSector = new TGCompositeFrame(fContCuts, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
  fContCuts->AddFrame(fContSector, new TGLayoutHints(kLHintsExpandX, 5, 0, 0, 0));

    // sector number label
    fLblSector = new TGLabel(fContSector, "SM");
    fLblSector->SetTextJustify(kTextLeft);
    fContSector->AddFrame(fLblSector, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 5, 0, 0, 0));
    // sector number entry
    fNmbSector = new TGNumberEntry(fContSector, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, -1, 17);
    fContSector->AddFrame(fNmbSector, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
    fNmbSector->SetNumber(-1);
    fNmbSector->Connect("ValueSet(Long_t)", "AliTRDCalibViewerGUI", this, "DoNewSelection()");

  // stack options container
  fContStack = new TGCompositeFrame(fContCuts, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
  fContCuts->AddFrame(fContStack, new TGLayoutHints(kLHintsExpandX, 5, 0, 0, 0));

    // stack number label
    fLblStack = new TGLabel(fContStack, "Stack");
    fLblStack->SetTextJustify(kTextLeft);
    fContStack->AddFrame(fLblStack, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 5, 0, 0, 0));
    // stack number entry
    fNmbStack = new TGNumberEntry(fContStack, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, -1, 4);
    fContStack->AddFrame(fNmbStack, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
    fNmbStack->SetNumber(-1);
    fNmbStack->Connect("ValueSet(Long_t)", "AliTRDCalibViewerGUI", this, "DoNewSelection()");
    
  // Load run frame
  fContLoad = new TGGroupFrame(fTabRight1, "Load run", kVerticalFrame | kFitWidth | kFitHeight);
  fTabRight1->AddFrame(fContLoad, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
    // Storage
    fContStorage = new TGCompositeFrame(fContLoad, 400, 200, kHorizontalFrame | kFitWidth | kFitHeight);
    fContLoad->AddFrame(fContStorage, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      // label
      fLblStorage = new TGLabel(fContStorage, "OCDB:");
      fLblStorage->SetTextJustify(kTextLeft);
      fContStorage->AddFrame(fLblStorage, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      // text entry
      fTxtStorage = new TGTextEntry(fContStorage, "alien://folder=/alice/data/2010/OCDB/", 111);
      fContStorage->AddFrame(fTxtStorage, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
//      fTxtStorage->Connect("ReturnPressed()", "AliBaseCalibViewerGUI", this, "HandleButtons1D(=111)");
      fTxtStorage->SetToolTipText("Enter the OCDB storage location");
    // Run entry
    fContRun = new TGCompositeFrame(fContLoad, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
    fContLoad->AddFrame(fContRun, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      // label
      fLblRun = new TGLabel(fContRun, "Run:");
      fLblRun->SetTextJustify(kTextLeft);
      fContRun->AddFrame(fLblRun, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      // run number entry
      fNmbRun = new TGNumberEntry(fContRun, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
      fContRun->AddFrame(fNmbRun, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      fNmbRun->SetNumber(-1);
      //fNmbRun->SetToolTipText("Enter the run number");
    // Version entry
    fContVersion = new TGCompositeFrame(fContLoad, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
    fContLoad->AddFrame(fContVersion, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      // label
      fLblVersion = new TGLabel(fContVersion, "Version:");
      fLblVersion->SetTextJustify(kTextLeft);
      fContVersion->AddFrame(fLblVersion, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      // run number entry
      fNmbVersion = new TGNumberEntry(fContVersion, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
      fContVersion->AddFrame(fNmbVersion, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      fNmbVersion->SetNumber(-1);
    // SubVersion entry
    fContSubVersion = new TGCompositeFrame(fContLoad, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
    fContLoad->AddFrame(fContSubVersion, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      // label
      fLblSubVersion = new TGLabel(fContSubVersion, "SubVersion:");
      fLblSubVersion->SetTextJustify(kTextLeft);
      fContSubVersion->AddFrame(fLblSubVersion, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      // run number entry
      fNmbSubVersion = new TGNumberEntry(fContSubVersion, 0, 1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
      fContSubVersion->AddFrame(fNmbSubVersion, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      fNmbSubVersion->SetNumber(-1);

    // Calib & DCS & Align check boxes frame
    fContChecks = new TGCompositeFrame(fContLoad, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
    fContLoad->AddFrame(fContChecks, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
      // Calibs check box
      fChkCalibs = new TGCheckButton(fContChecks, "Calib");
      fContChecks->AddFrame(fChkCalibs, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
      fChkCalibs->SetToolTipText("Get calibration info (gain, pedestal, vdrift, T0 and status)");
      fChkCalibs->SetState(kButtonDown);
      // DCS check box
      fChkDCS = new TGCheckButton(fContChecks, "DCS");
      fContChecks->AddFrame(fChkDCS, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
      fChkDCS->SetToolTipText("Get DCS info");
      fChkDCS->SetState(kButtonDown);
      // Calibs check box
      fChkAlign = new TGCheckButton(fContChecks, "Align");
      fContChecks->AddFrame(fChkAlign, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
      fChkAlign->SetToolTipText("Get alingment info");
      fChkAlign->SetState(kButtonDown);
    // Load button
    fBtnLoad = new TGTextButton(fContLoad, "&Load run");
    fBtnLoad->SetName("loadOCDB");
    fContLoad->AddFrame(fBtnLoad, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
    fBtnLoad->SetToolTipText("Load run from OCDB");
    fBtnLoad->Connect("Clicked()", "AliTRDCalibViewerGUI", this, "SetTree()");

  // Load a file with AliTRDCalPad objects
  // main frame
  fContLoadCalibObjects = new TGGroupFrame(fTabRight1, "Load CalPad from file", kVerticalFrame | kFitWidth | kFitHeight);
  fTabRight1->AddFrame(fContLoadCalibObjects, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
    // container for the input file
    fContCalibInput = new TGCompositeFrame(fContLoadCalibObjects, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
    fContLoadCalibObjects->AddFrame(fContCalibInput, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 0));
      // label
      fLblCalibInputFilename = new TGLabel(fContCalibInput, "Input:");
      fLblCalibInputFilename->SetTextJustify(kTextLeft);
      fContCalibInput->AddFrame(fLblCalibInputFilename, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      // text entry
      fTxtCalibInputFilename = new TGTextEntry(fContCalibInput, "Input file", 200);
      fContCalibInput->AddFrame(fTxtCalibInputFilename, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      fTxtCalibInputFilename->Connect("DoubleClicked()", "AliTRDCalibViewerGUI", this, "HandleFilesystem()");
    // container for the output file
    fContCalibOutput = new TGCompositeFrame(fContLoadCalibObjects, 200, 200, kHorizontalFrame | kFitWidth | kFitHeight);
    fContLoadCalibObjects->AddFrame(fContCalibOutput, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 0));
      // label
      fLblCalibOutputFilename = new TGLabel(fContCalibOutput, "Output:");
      fLblCalibOutputFilename->SetTextJustify(kTextLeft);
      fContCalibOutput->AddFrame(fLblCalibOutputFilename, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      // text entry
      fTxtCalibOutputFilename = new TGTextEntry(fContCalibOutput, "/tmp/output.root", 201);
      fContCalibOutput->AddFrame(fTxtCalibOutputFilename, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
      fTxtCalibOutputFilename->Connect("DoubleClicked()", "AliTRDCalibViewerGUI", this, "HandleFilesystem()");
    // Load button
    fBtnLoadCalibObjects = new TGTextButton(fContLoadCalibObjects, "L&oad calib");
    fBtnLoadCalibObjects->SetName("loadCalPad");
    fContLoadCalibObjects->AddFrame(fBtnLoadCalibObjects, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
    fBtnLoadCalibObjects->SetToolTipText("Extract a tree from an array of AliTRDCalPad objects,\nand load it into this application");
    fBtnLoadCalibObjects->Connect("Clicked()", "AliTRDCalibViewerGUI", this, "SetTree()");
      
  SetWindowName("AliTRDCalibViewer GUI");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
}

//___________________________________________________________________________________
void AliTRDCalibViewerGUI::SetTree() {
  //
  //  handles the loading of a new tree and updates the GUI
  //
  TString fileName;
  TGTextButton *button = ((TGTextButton*)gTQSender);
  TString name=button->GetName();
  if(name.Contains("OCDB")) {
    Int_t run = (Int_t)fNmbRun->GetNumber();
    TString storage(fTxtStorage->GetText());
    Bool_t getCalibs = (fChkCalibs->GetState()==kButtonDown ? kTRUE : kFALSE);
    Bool_t getDCS = (fChkDCS->GetState()==kButtonDown ? kTRUE : kFALSE);
    Bool_t getAlign = (fChkAlign->GetState()==kButtonDown ? kTRUE : kFALSE);
    Int_t version = (Int_t)fNmbVersion->GetNumber();
    Int_t subVersion = (Int_t)fNmbSubVersion->GetNumber();
    
    fileName = Form("trdOCDBDetails_run%d.root", run);
    Bool_t load=kTRUE;
    load=((AliTRDCalibViewer*)fViewer)->DumpOCDBtoTreeDetails("", fileName.Data(), run, run, storage.Data(),
                                                              version, subVersion, getCalibs, getDCS, getAlign);
    if(!load) {
      Reset();
      return;
    }
  }
  if(name.Contains("CalPad")) {
    ((AliTRDCalibViewer*)fViewer)->DumpCalibToTree(fTxtCalibInputFilename->GetText(),
                                                   fTxtCalibOutputFilename->GetText());
    fileName = fTxtCalibOutputFilename->GetText();
  }
  Initialize(fileName.Data(), "TRDcalibDetails");
  Reload();
  DoDraw();
}

//________________________________________________________________________________________________
void AliTRDCalibViewerGUI::HandleFilesystem() {
  //
  //  Slot used by the text buttons to trigger the file system dialog
  //
  Int_t id = ((TGTextEntry*)gTQSender)->WidgetId();
  const char *kTypes[] = {
      "All files",    "*",
       0,              0};
  TString dir(".");
  TGFileInfo fi;
  fi.fFileTypes = kTypes;
  fi.fOverwrite = kFALSE;
  new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen, &fi);
  if(fi.fFilename && strlen(fi.fFilename)) {
    if(id==200)
      fTxtCalibInputFilename->SetText(fi.fFilename);
    if(id==201)
      fTxtCalibOutputFilename->SetText(fi.fFilename);
  }
  return;
}

//________________________________________________________________________________________________
void AliTRDCalibViewerGUI::Initialize(const char* fileName, const char* treeName) {
  // 
  // initialize the GUI with a calibrationTree from fileName
  // 
  // create AliTRDCalibViewer object, which will be used for generating all drawings
  if (fViewer) delete fViewer;
  fViewer = new AliTRDCalibViewer(fileName, treeName);
  Initialize(fViewer);   
}

//________________________________________________________________________________________________
void AliTRDCalibViewerGUI::Initialize(AliBaseCalibViewer *viewer) {
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
  TObjArray* arr = ((AliTRDCalibViewer*)fViewer)->GetListOfVariables(1);
  if (!arr) {
    return;
  }
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
  TObjArray *arrNorm = ((AliTRDCalibViewer*)fViewer)->GetListOfNormalizationVariables();
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
  
  if (fInitialized) Info("Initialize", "AliTRDCalibViewerGUI new initialized.");
  fInitialized = kTRUE;
}

//________________________________________________________________________________________________
void AliTRDCalibViewerGUI::Reset(){
  //
  // reset variables, delete calib viewer
  //
  if (fViewer) delete fViewer;
  fListVariables->RemoveAll();
  fListNormalization->RemoveAll();
  fInitialized = kFALSE;
}

//________________________________________________________________________________________________
TString* AliTRDCalibViewerGUI::GetDrawString() {
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
    desiredData = TString(((AliTRDCalibViewer*)fViewer)->AddAbbreviations((char*)desiredData.Data()));
  return new TString(desiredData.Data());
}

//________________________________________________________________________________________________
TString* AliTRDCalibViewerGUI::GetCutString() {
   // 
   // create the cut string out of selection
   // 
  
   TString cutsStr("");
      
   // try to add forgotten '~'
   if(fChkAutoAppend->GetState() == kButtonDown) 
      cutsStr = TString(((AliTRDCalibViewer*)fViewer)->AddAbbreviations((char*)cutsStr.Data()));
   return new TString(cutsStr.Data());
}

//________________________________________________________________________________________________
TString* AliTRDCalibViewerGUI::GetSectorString() {
  // 
  // create the sector string out of selection
  // 

  Int_t layerNo = (Int_t)(fNmbLayer->GetNumber());
  Int_t sectorNo = (Int_t)(fNmbSector->GetNumber());
  Int_t stackNo = (Int_t)(fNmbStack->GetNumber());

  TString sectorStr("");
  sectorStr = Form("Layer%dSector%dStack%d", layerNo, sectorNo, stackNo);

  return new TString(sectorStr.Data());
}   

//________________________________________________________________________________________________
void AliTRDCalibViewerGUI::DoDraw() {
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
      entries = ((AliTRDCalibViewer*)fViewer)->EasyDraw1D(desiredData.Data(), sectorStr.Data(), cutsStr.Data(), addDrawOpt.Data());
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
    entries = ((AliTRDCalibViewer*)fViewer)->EasyDraw(desiredData.Data(), sectorStr.Data(), cutsStr.Data(), addDrawOpt.Data());
  }
  if (entries == -1) return; // nothing was drawn, there is no histogram to get min and max
   
  SetMinMaxLabel();
  fCanvMain->GetCanvas()->Update();
}

//________________________________________________________________________________________________
void AliTRDCalibViewerGUI::MouseMove(Int_t event, Int_t x, Int_t y, TObject *selectedObject) { 
  //
  // mouse move
  // zoom to chamber works ONLY in 2D mode
  // 
  if(event != kButton1Double )
    return;
  if(!selectedObject->InheritsFrom("TH2")) return;

  //Int_t layerNo = (Int_t)(fNmbLayer->GetNumber());
  Int_t sectorNo = (Int_t)(fNmbSector->GetNumber());
  Int_t stackNo = (Int_t)(fNmbStack->GetNumber());

  // zoom out to the current layer if a chamber is viewed now
  if(sectorNo!=-1 && stackNo!=-1) {
    fNmbSector->SetNumber(-1);
    fNmbStack->SetNumber(-1);
    DoNewSelection();
    return;
  }

  // check what kind of parameter we visualize
  TString drawStr(GetDrawString()->Data());
  Int_t viewedParamClass = -1;     // -1 nothing, 0 calibration, 1 FEE params
  if(drawStr.Contains("Status") || drawStr.Contains("Gain") || drawStr.Contains("Noise") ||
     drawStr.Contains("Vdrift") || drawStr.Contains("T0") ||
     drawStr.Contains("gain") || drawStr.Contains("chiSquare"))
    viewedParamClass = 0;
  if(drawStr.Contains("SORandEOR") || 
     drawStr.Contains("gsmSOR") || drawStr.Contains("gsmDelta") ||
     drawStr.Contains("nimSOR") || drawStr.Contains("nimDelta") ||
     drawStr.Contains("nevSOR") || drawStr.Contains("nevDelta") ||
     drawStr.Contains("nptSOR") || drawStr.Contains("nptDelta")) {
    viewedParamClass = 1;
  }
  if(viewedParamClass==-1) return;

  // some constants refering to the TRD geometry
  const Int_t gkNRows[ 5] = {16, 16, 12, 16, 16};  // number of pad rows in the chambers from each of the 5 stacks
  const Int_t gkNCols = 144;    // number of pad cols per chamber

  // get the coordinate of the clicked point in physical coordinates
  Float_t upy = gPad->AbsPixeltoY(y);
  Float_t upx = gPad->AbsPixeltoX(x);
  Float_t gy  = gPad->PadtoY(upy);
  Float_t gx  = gPad->PadtoX(upx);
  Int_t selectedStack = -1;
  Int_t selectedSector = -1;

  // retrieve the double-clicked chamber 
  if(sectorNo==-1 && stackNo==-1) {
    // get the selected stack
    Float_t rowLowBound = -0.5;
    Float_t rowHighBound = -0.5;
    for(Int_t i=0; i<5; i++) {
      if(i>0) rowLowBound += gkNRows[i-1];
      rowHighBound += gkNRows[i];
      if(gx>=rowLowBound && gx<=rowHighBound)
	selectedStack = i;
    }
    // get the selected sector
    if(viewedParamClass==0) {   // calibration params
      selectedSector = (Int_t)TMath::Floor((gy+0.5)/Float_t(gkNCols));
    }
    if(viewedParamClass==1) {   // FEE params
      selectedSector = (Int_t)TMath::Floor((gy+0.5)/8.0);   // 8 MCMs per chamber in pad cols direction
    }
  }
  if(sectorNo!=-1 && stackNo==-1) {
    // get the selected stack
    Float_t rowLowBound = -0.5;
    Float_t rowHighBound = -0.5;
    for(Int_t i=0; i<5; i++) {
      if(i>0) rowLowBound += gkNRows[i-1];
      rowHighBound += gkNRows[i];
      if(gx>=rowLowBound && gx<=rowHighBound)
	selectedStack = i;
    }
    // get the selected sector
    selectedSector = sectorNo;
  }
  if(sectorNo==-1 && stackNo!=-1) {
    // get the selected stack
    selectedStack = stackNo;
    // get the selected sector
    if(viewedParamClass==0) {   // calibration params
      selectedSector = (Int_t)TMath::Floor((gy+0.5)/144.0);
    }
    if(viewedParamClass==1) {   // FEE params
      selectedSector = (Int_t)TMath::Floor((gy+0.5)/8.0);
    }
  }

  fNmbSector->SetNumber(selectedSector);
  fNmbStack->SetNumber(selectedStack);
  DoNewSelection();
  return;
}

//___________________________________________________________________________
void AliTRDCalibViewerGUI::ShowGUI() {
  //
  //   Draw the graphical user interface
  //
  TGMainFrame* mainWindow = new TGMainFrame(gClient->GetRoot(), 1000, 700);
  mainWindow->SetWindowName("Run OCDB details");
  mainWindow->SetCleanup(kDeepCleanup);
  AliTRDCalibViewerGUI *calibViewer = new AliTRDCalibViewerGUI(mainWindow, 1000, 650, 0);
  mainWindow->AddFrame(calibViewer, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  mainWindow->MapSubwindows();
  mainWindow->Resize();
  mainWindow->MapWindow();
}

//___________________________________________________________________________
void AliTRDCalibViewerGUI::ShowGUI(const Char_t* treeFile, const Char_t* treeName) {
  //
  //   Draw the graphical user interface
  //
  TGMainFrame* mainWindow = new TGMainFrame(gClient->GetRoot(), 1000, 700);
  mainWindow->SetWindowName("Run OCDB details");
  mainWindow->SetCleanup(kDeepCleanup);
  AliBaseCalibViewerGUI *calibViewer = new AliTRDCalibViewerGUI(mainWindow, 1000, 650, 0);
  mainWindow->AddFrame(calibViewer, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

  calibViewer->Initialize(treeFile, treeName);
  calibViewer->Reload();
  calibViewer->DoDraw();

  mainWindow->MapSubwindows();
  mainWindow->Resize();
  mainWindow->MapWindow();
}

//___________________________________________________________________________
void AliTRDCalibViewerGUI::ShowGUIwithTrending() {
  //
  // Draw a GUI application containing 2 tabs:
  //    -- tab for time/run trending
  //    -- tab for run details
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 700);
  frmMain->SetCleanup(kDeepCleanup);
  
  TGTab* tabMain = new TGTab(frmMain, 1000, 700);
  frmMain->AddFrame(tabMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  TGCompositeFrame* tabCont1 = tabMain->AddTab("Time");
  TGCompositeFrame* tabCont2 = tabMain->AddTab("Detail - XXXXX");
  
  AliCalibViewerGUItime* calibViewerTime = new AliCalibViewerGUItime(tabCont1, 1000, 650, "TRD");
  tabCont1->AddFrame(calibViewerTime, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  
  AliTRDCalibViewerGUI *calibViewer = new AliTRDCalibViewerGUI(tabCont2, 1000, 700, 0);
  tabCont2->AddFrame(calibViewer, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  calibViewerTime->SetCalibViewerGUI(calibViewer);
  calibViewerTime->SetCalibViewerGUItab(tabMain->GetTabTab(1));
  
  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();
}
