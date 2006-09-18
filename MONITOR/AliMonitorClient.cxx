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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class receives monitor histograms from a monitor process and        //
//  provides a graphical user interface for browsing and analysing the       //
//  monitor histograms.                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliMonitorClient.h"
#include "AliMonitorProcess.h"
#include "AliLog.h"
#include <TGMsgBox.h>
#include <TGFileDialog.h>
#include <TMessage.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TGMenu.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TGToolBar.h>
#include <TG3DLine.h>
#include <TGNumberEntry.h>
#include <TGCanvas.h>
#include <TGSplitter.h>
#include <TGListTree.h>
#include <TRootEmbeddedCanvas.h>
#include <TGTextView.h>
#include <TFolder.h>
#include <TSocket.h>
#include <TTimer.h>
#include <TFile.h>
#include "AliMonitorHisto.h"


ClassImp(AliMonitorClient) 


const char* AliMonitorClient::fgSettingsFileName = ".AliMonitorClient";


//_____________________________________________________________________________
AliMonitorClient::AliMonitorStringDlg::AliMonitorStringDlg(TString& string, 
							   TGFrame* main, 
							   const char* title, 
							   const char* label) :
  AliMonitorDialog(main, 300, 80),
  fStringLayout(new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 2, 2)),
  fStringLabel(new TGLabel(fFrame, label)),
  fStringEntry(new TGTextEntry(fFrame, string.Data())),
  fString(string)
{
// create a dialog for connecting to a monitor process

  fFrame->AddFrame(fStringLabel, fStringLayout);
  fStringEntry->Resize(100, fStringEntry->GetDefaultHeight());
  fFrame->AddFrame(fStringEntry, fStringLayout);

  fString = "";

  fMain->SetWindowName(title);
  fMain->MapSubwindows();
  fMain->Layout();
  gClient->WaitFor(fMain);
}

//_____________________________________________________________________________
AliMonitorClient::AliMonitorStringDlg::~AliMonitorStringDlg()
{
// clean up

  delete fStringLabel;
  delete fStringLayout;
  delete fStringEntry;
}

//_____________________________________________________________________________
void AliMonitorClient::AliMonitorStringDlg::OnOkClicked()
{
  fString = fStringEntry->GetText();
}


//_____________________________________________________________________________
AliMonitorClient::AliMonitorNumberDlg::AliMonitorNumberDlg(Float_t& value, 
							   TGFrame* main, 
							   const char* title, 
							   const char* label, 
							   Float_t min) :
  AliMonitorDialog(main, 250, 80),
  fNumberLayout(new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 2, 2)),
  fNumberLabel(new TGLabel(fFrame, label)),
  fNumberEntry(new TGNumberEntry(fFrame, value, 4, -1, TGNumberFormat::kNESRealOne)),
  fNumber(value)
{
// create a dialog for getting a number

  fFrame->AddFrame(fNumberLabel, fNumberLayout);
  fNumberEntry->SetLimits(TGNumberFormat::kNELLimitMin, min, 0);
  fFrame->AddFrame(fNumberEntry, fNumberLayout);

  fNumber = -1;

  fMain->SetWindowName(title);
  fMain->MapSubwindows();
  fMain->Layout();
  gClient->WaitFor(fMain);
}

//_____________________________________________________________________________
AliMonitorClient::AliMonitorNumberDlg::~AliMonitorNumberDlg()
{
// clean up

  delete fNumberLabel;
  delete fNumberLayout;
  delete fNumberEntry;
}

//_____________________________________________________________________________
void AliMonitorClient::AliMonitorNumberDlg::OnOkClicked()
{
  fNumber = fNumberEntry->GetNumber();
}



// constants for menu entries
enum {kMenuFileConnect, kMenuFileDisconnect, kMenuFileOpen, kMenuFileExit,
      kMenuViewToolBar, kMenuViewTree, kMenuViewDescription, 
      kMenuViewReference, kMenuViewStatistics,
      kMenuFavAdd, kMenuFavDelete, 
      kMenuFavLoad, kMenuFavSave,  kMenuFavSaveAs, kMenuFavSaveOnExit,
      kMenuRefLoad, kMenuRefThreshold, 
      kMenuRefTakeCurrent, kMenuRefSave, kMenuRefSaveAs,
      kMenuOptLoop, kMenuOptPrint, 
      kMenuOptSave, kMenuOptSaveOnExit,
      kMenuHelpDoc, kMenuHelpAbout
};

//_____________________________________________________________________________
AliMonitorClient::AliMonitorClient():
  TGMainFrame(gClient->GetRoot(), 500, 300),
  fQObject(),
  fMenuBarLayout(new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1)),
  fMenuBarItemLayout(new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0)),
  fMenuBarHelpLayout(new TGLayoutHints(kLHintsTop | kLHintsRight)),
  fMenuFile(new TGPopupMenu(gClient->GetRoot())),
  fMenuView(new TGPopupMenu(gClient->GetRoot())),
  fMenuFavorites(new TGPopupMenu(gClient->GetRoot())),
  fMenuReference(new TGPopupMenu(gClient->GetRoot())),
  fMenuOptions(new TGPopupMenu(gClient->GetRoot())),
  fMenuHelp(new TGPopupMenu(gClient->GetRoot())),
  fMenuBar(new TGMenuBar(this, 1, 1, kHorizontalFrame)),
  fToolBarLayout(new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 2, 2)),
  fToolBarSep(new TGHorizontal3DLine(this)),
  fToolBar(new TGToolBar(this, 60, 20, kHorizontalFrame)),
  fEventNumberLayout(new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 10, 2, 10, 5)),
  fEventNumber(new TGNumberEntry(fToolBar, 1, 4, 10, TGNumberFormat::kNESInteger)),
  fEventButton(NULL),
  fSumButton(NULL),
  fRunButton(NULL),
  fLoopButton(NULL),
  fLoopOnPicture(NULL),
  fLoopOffPicture(NULL),
  fPreviousButton(NULL),
  fNextButton(NULL),
  fCopyButton(NULL),
  fSaveButton(NULL),
  fPrintButton(NULL),
  fBottomLayout(new TGLayoutHints(kLHintsExpandX | kLHintsBottom)),
  fLeftLayout(new TGLayoutHints(kLHintsLeft | kLHintsExpandY)),
  fExpandLayout(new TGLayoutHints(kLHintsExpandX | kLHintsExpandY)),
  fVerticalFrame(new TGVerticalFrame(this, 10, 10)),
  fHorizontalFrame(new TGHorizontalFrame(fVerticalFrame, 10, 10)),
  fTreeFrame(new TGCompositeFrame(fHorizontalFrame, 10, 10, kSunkenFrame | kFixedWidth)),
  fTreeCanvas(new TGCanvas(fTreeFrame, 10, 10)),
  fTree(new TGListTree(fTreeCanvas, kHorizontalFrame)),
  fHistoPicture(fClient->GetPicture("h1_t.xpm")),
  fAllItem(fTree->AddItem(NULL, "All")),
  fFavoritesItem(fTree->AddItem(NULL, "Favorites")),
  fComparisonItem(fTree->AddItem(NULL, "Comparison")),
  fTreeSplitter(new TGVSplitter(fHorizontalFrame, 4)),
  fDrawFrame(new TGCompositeFrame(fHorizontalFrame, 10, 10, kSunkenFrame)),
  fDrawCanvas(new TRootEmbeddedCanvas("current monitor histogram", fDrawFrame, 10, 10)),
  fDescriptionSplitter(new TGHSplitter(fVerticalFrame, 4, 4)),
  fDescriptionFrame(new TGCompositeFrame(fVerticalFrame, 10, 60, kSunkenFrame | kFixedHeight)),
  fDescription(new TGTextView(fDescriptionFrame, 10, 60, "")),
  fServerName("localhost"),
  fSocket(NULL),
  fSocketHandler(NULL),
  fFolder(CreateTopFolder()),
  fCurrentItem(NULL),
  fBaseItem(NULL),
  fLoopTimer(NULL),
  fLoopInterval(1000),
  fFavoritesFileName(""),
  fReferenceFileName(""),
  fReference(CreateTopFolder()),
  fPrintCommand("gv")
{
// initialize the monitoring client window

  // File menu
  fMenuFile->AddEntry("&Connect...", kMenuFileConnect);
  fMenuFile->AddEntry("&Disconnect...", kMenuFileDisconnect);
  fMenuFile->HideEntry(kMenuFileDisconnect);
  fMenuFile->AddEntry("&Open...", kMenuFileOpen);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("E&xit", kMenuFileExit);
  fMenuFile->Connect("Activated(Int_t)", "AliMonitorClient", this,
		     "OnMenuActivated(Int_t)");

  // View menu
  fMenuView->AddEntry("Tool&bar", kMenuViewToolBar);
  fMenuView->CheckEntry(kMenuViewToolBar);
  fMenuView->AddEntry("&Tree", kMenuViewTree);
  fMenuView->CheckEntry(kMenuViewTree);
  fMenuView->AddEntry("&Description", kMenuViewDescription);
  fMenuView->CheckEntry(kMenuViewDescription);
  fMenuView->AddSeparator();
  fMenuView->AddEntry("&Reference", kMenuViewReference);
  fMenuView->CheckEntry(kMenuViewReference);
  AliMonitorHisto::SetDrawRef(kTRUE);
  fMenuView->AddEntry("&Statistics", kMenuViewStatistics);
  fMenuView->CheckEntry(kMenuViewStatistics);
  gStyle->SetOptStat(1110);
  fMenuView->Connect("Activated(Int_t)", "AliMonitorClient", this,
		     "OnMenuActivated(Int_t)");

  // Favorites menu
  fMenuFavorites->AddEntry("&Add", kMenuFavAdd);
  fMenuFavorites->DisableEntry(kMenuFavAdd);
  fMenuFavorites->AddEntry("&Delete", kMenuFavDelete);
  fMenuFavorites->DisableEntry(kMenuFavDelete);
  fMenuFavorites->AddSeparator();
  fMenuFavorites->AddEntry("&Load...", kMenuFavLoad);
  fMenuFavorites->AddEntry("&Save", kMenuFavSave);
  fMenuFavorites->DisableEntry(kMenuFavSave);
  fMenuFavorites->AddEntry("&Save As...", kMenuFavSaveAs);
  fMenuFavorites->DisableEntry(kMenuFavSaveAs);
  fMenuFavorites->AddEntry("Save On &Exit", kMenuFavSaveOnExit);
  fMenuFavorites->CheckEntry(kMenuFavSaveOnExit);
  fMenuFavorites->Connect("Activated(Int_t)", "AliMonitorClient", this,
			  "OnMenuActivated(Int_t)");

  // Reference menu
  fMenuReference->AddEntry("&Load...", kMenuRefLoad);
  fMenuReference->AddEntry("&Threshold...", kMenuRefThreshold);
  fMenuReference->AddSeparator();
  fMenuReference->AddEntry("Take &Current", kMenuRefTakeCurrent);
  fMenuReference->DisableEntry(kMenuRefTakeCurrent);
  fMenuReference->AddEntry("&Save", kMenuRefSave);
  fMenuReference->DisableEntry(kMenuRefSave);
  fMenuReference->AddEntry("Save &As...", kMenuRefSaveAs);
  fMenuReference->DisableEntry(kMenuRefSaveAs);
  fMenuReference->Connect("Activated(Int_t)", "AliMonitorClient", this,
			  "OnMenuActivated(Int_t)");

  // Options menu
  fMenuOptions->AddEntry("&Loop Interval...", kMenuOptLoop);
  fMenuOptions->AddEntry("&Print Command...", kMenuOptPrint);
  fMenuOptions->AddSeparator();
  fMenuOptions->AddEntry("&Save Settings", kMenuOptSave);
  fMenuOptions->AddEntry("Save Settings on &Exit", kMenuOptSaveOnExit);
  fMenuOptions->CheckEntry(kMenuOptSaveOnExit);
  fMenuOptions->Connect("Activated(Int_t)", "AliMonitorClient", this,
			"OnMenuActivated(Int_t)");

  // Help menu
  fMenuHelp->AddEntry("&Documentation...", kMenuHelpDoc);
  fMenuHelp->AddEntry("A&bout...", kMenuHelpAbout);
  fMenuHelp->Connect("Activated(Int_t)", "AliMonitorClient", this,
		     "OnMenuActivated(Int_t)");

  // menu bar
  fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
  fMenuBar->AddPopup("&View", fMenuView, fMenuBarItemLayout);
  fMenuBar->AddPopup("F&avorites", fMenuFavorites, fMenuBarItemLayout);
  fMenuBar->AddPopup("&Reference", fMenuReference, fMenuBarItemLayout);
  fMenuBar->AddPopup("&Options", fMenuOptions, fMenuBarItemLayout);
  fMenuBar->AddPopup("&Help", fMenuHelp, fMenuBarHelpLayout);

  AddFrame(fMenuBar, fMenuBarLayout);


  // *** tool bar ***
  AddFrame(fToolBarSep, fToolBarLayout);
  AddFrame(fToolBar, fToolBarLayout);

  // event number field
  fEventNumber->SetLimits(TGNumberFormat::kNELLimitMin, 1, 0);
  fToolBar->AddFrame(fEventNumber, fEventNumberLayout);
  fEventNumber->GetNumberEntry()->Connect("ReturnPressed()", 
					  "AliMonitorClient", this,
					  "OnEventNumberChanged()");
  fEventNumber->GetButtonUp()->Connect("Clicked()", 
				       "AliMonitorClient", this,
				       "OnEventNumberChanged()");
  fEventNumber->GetButtonDown()->Connect("Clicked()", 
					 "AliMonitorClient", this,
					 "OnEventNumberChanged()");

  // draw event button
  ToolBarData_t eventButton = {"h1_s.xpm", "draw event histogram", 
			       kTRUE, 11, NULL};
  fToolBar->AddButton(this, &eventButton, 5);
  fEventButton = eventButton.fButton;
  fEventButton->Connect("Pressed()", "AliMonitorClient", this,
			"OnEventButtonPressed()");

  // draw sum button
  ToolBarData_t sumButton = {"h2_s.xpm", "draw sum histogram", 
			     kTRUE, 12, NULL};
  fToolBar->AddButton(this, &sumButton, 5);
  fSumButton = sumButton.fButton;
  fSumButton->Connect("Pressed()", "AliMonitorClient", this,
		      "OnSumButtonPressed()");

  // draw run button
  ToolBarData_t runButton = {"h3_s.xpm", "draw run histogram", 
			     kTRUE, 13, NULL};
  fToolBar->AddButton(this, &runButton, 5);
  fRunButton = runButton.fButton;
  fRunButton->SetDown(kTRUE);
  fRunButton->Connect("Pressed()", "AliMonitorClient", this,
		      "OnRunButtonPressed()");

  // loop button
  char fileName[256];
  sprintf(fileName, "%s/MONITOR/loop_on.xpm", gSystem->Getenv("ALICE_ROOT"));
  ToolBarData_t loopButton = {fileName, "loop over histograms", 
			      kFALSE, 21, NULL};
  fToolBar->AddButton(this, &loopButton, 25);
  fLoopButton = loopButton.fButton;
  fLoopButton->Connect("Clicked()", "AliMonitorClient", this,
		       "OnLoopButtonClicked()");
  fLoopOnPicture = gClient->GetPicture(fileName);
  sprintf(fileName, "%s/MONITOR/loop_off2.xpm", gSystem->Getenv("ALICE_ROOT"));
  fLoopOffPicture = gClient->GetPicture(fileName);

  // previous button
  sprintf(fileName, "%s/MONITOR/previous_s.xpm", 
	  gSystem->Getenv("ALICE_ROOT"));
  ToolBarData_t previousButton = {fileName, "go to previous histogram", 
				  kFALSE, 22, NULL};
  fToolBar->AddButton(this, &previousButton, 5);
  fPreviousButton = previousButton.fButton;
  fPreviousButton->Connect("Clicked()", "AliMonitorClient", this,
			   "OnPreviousButtonClicked()");

  // next button
  sprintf(fileName, "%s/MONITOR/next_s.xpm", gSystem->Getenv("ALICE_ROOT"));
  ToolBarData_t nextButton = {fileName, "go to next histogram", 
			      kFALSE, 23, NULL};
  fToolBar->AddButton(this, &nextButton, 5);
  fNextButton = nextButton.fButton;
  fNextButton->Connect("Clicked()", "AliMonitorClient", this,
		       "OnNextButtonClicked()");

  // copy button
  sprintf(fileName, "%s/MONITOR/copy_s.xpm", gSystem->Getenv("ALICE_ROOT"));
  ToolBarData_t copyButton = {fileName, 
			      "copy the current histogram to a new canvas", 
			      kFALSE, 31, NULL};
  fToolBar->AddButton(this, &copyButton, 25);
  fCopyButton = copyButton.fButton;
  fCopyButton->Connect("Clicked()", "AliMonitorClient", this,
		       "OnCopyButtonClicked()");

  // save button
  sprintf(fileName, "%s/MONITOR/save_s.xpm", gSystem->Getenv("ALICE_ROOT"));
  ToolBarData_t saveButton = {fileName, "save the current histogram", 
			      kFALSE, 32, NULL};
  fToolBar->AddButton(this, &saveButton, 5);
  fSaveButton = saveButton.fButton;
  fSaveButton->Connect("Clicked()", "AliMonitorClient", this,
		       "OnSaveButtonClicked()");

  // print button
  ToolBarData_t printButton = {"printer_s.xpm", "print the current histogram", 
			     kFALSE, 33, NULL};
  fToolBar->AddButton(this, &printButton, 5);
  fPrintButton = printButton.fButton;
  fPrintButton->Connect("Clicked()", "AliMonitorClient", this,
			"OnPrintButtonClicked()");

  // *** frames ***
  AddFrame(fVerticalFrame, fExpandLayout);

  fVerticalFrame->AddFrame(fHorizontalFrame, fExpandLayout);

  // tree frame
  fHorizontalFrame->AddFrame(fTreeFrame, fLeftLayout);
  fTreeFrame->AddFrame(fTreeCanvas, fExpandLayout);
  fTreeCanvas->AddFrame(fTree, fExpandLayout);
  fTree->Connect("Clicked(TGListTreeItem*,Int_t)", "AliMonitorClient",
		 this, "OnTreeClicked(TGListTreeItem*,Int_t)");
  fTree->Connect("ReturnPressed(TGListTreeItem*)", "AliMonitorClient",
		 this, "OnTreeReturnPressed(TGListTreeItem*)");

  // tree items
  fTreeFrame->Resize(100, fTreeFrame->GetDefaultHeight());

  // tree / histogram splitter
  fTreeSplitter->SetFrame(fTreeFrame, kTRUE);
  fHorizontalFrame->AddFrame(fTreeSplitter, fLeftLayout);

  // histogram frame
  fHorizontalFrame->AddFrame(fDrawFrame, fExpandLayout);
  fDrawFrame->AddFrame(fDrawCanvas, fExpandLayout);

  // description frame
  fVerticalFrame->AddFrame(fDescriptionFrame, fBottomLayout);
  fDescriptionFrame->AddFrame(fDescription, fExpandLayout);

  // histogram / description splitter
  fVerticalFrame->AddFrame(fDescriptionSplitter, fBottomLayout);
  fDescriptionSplitter->SetFrame(fDescriptionFrame, kFALSE);

  // main window
  Connect("CloseWindow()", "AliMonitorClient", this, "CloseWindow()");
  SetWindowName("Monitor Client");
  SetWMSize(GetWidth(), GetHeight());
  Move(100, 100);
  SetWMPosition(100, 100);
  MapSubwindows();
  Layout();
  MapWindow();


  // default data members

  fFolder = CreateTopFolder();

  fReference = CreateTopFolder();
  AliMonitorHisto::SetThreshold(5.0);

  // load saved settings
  LoadSettings();
}

//_____________________________________________________________________________
AliMonitorClient::~AliMonitorClient()
{
// clean up

  delete fMenuBarLayout;
  delete fMenuBarItemLayout;
  delete fMenuBarHelpLayout;
  delete fMenuFile;
  delete fMenuView;
  delete fMenuFavorites;
  delete fMenuReference;
  delete fMenuOptions;
  delete fMenuHelp;
  delete fMenuBar;

  delete fEventNumberLayout;
  delete fEventNumber;
  delete fToolBarLayout;
  delete fToolBarSep;
  delete fToolBar;

  delete fBottomLayout;
  delete fLeftLayout;
  delete fExpandLayout;

  delete fTree;
  delete fTreeCanvas;
  delete fTreeFrame;
  delete fTreeSplitter;

  delete fDrawCanvas;
  delete fDrawFrame;

  delete fDescriptionSplitter;
  delete fDescription;
  delete fDescriptionFrame;

  delete fVerticalFrame;
  delete fHorizontalFrame;

  if (fSocket) {
    fSocketHandler->Remove();
    fSocket->Send("disconnect");
    fSocket->Close();
    delete fSocket;
    delete fSocketHandler;
  }

  if (fFolder) delete fFolder;

  if (fLoopTimer) {
    fLoopTimer->TurnOff();
    delete fLoopTimer;
  }

  if (fReference) delete fReference;
}


//_____________________________________________________________________________
void AliMonitorClient::CloseWindow()
{
// terminate the application when the window is closed

  if (fMenuFavorites->IsEntryChecked(kMenuFavSaveOnExit)) {
    SaveFavorites();
  }
  if (fMenuOptions->IsEntryChecked(kMenuOptSaveOnExit)) {
    SaveSettings();
  }
  if (fSocket) {
    fSocketHandler->Remove();
    fSocket->Send("disconnect");
    fSocket->Close();
    delete fSocket;
    fSocket = NULL;
    delete fSocketHandler;
    fSocketHandler = NULL;
  }

  gApplication->Terminate(0);
}


//_____________________________________________________________________________
void AliMonitorClient::OnNewData()
{
// called when data has arrived from the monitor server

  if (CheckForNewData()) UpdateAll();
}


//_____________________________________________________________________________
void AliMonitorClient::OnMenuActivated(Int_t id)
{
// called when a menu item was selected

  switch (id) {

  case kMenuFileConnect: {
    if (ConnectToServer()) UpdateAll();
    break;
  }

  case kMenuFileDisconnect: {
    DisconnectFromServer();
    break;
  }

  case kMenuFileOpen: {
    if (OpenFile()) UpdateAll();
    break;
  }

  case kMenuFileExit: {
    CloseWindow();
    break;
  }

  case kMenuViewToolBar: {
    ViewToolBar(!fMenuView->IsEntryChecked(kMenuViewToolBar));
    break;
  }

  case kMenuViewTree: {
    ViewTree(!fMenuView->IsEntryChecked(kMenuViewTree));
    break;
  }

  case kMenuViewDescription: {
    ViewDescription(!fMenuView->IsEntryChecked(kMenuViewDescription));
    break;
  }

  case kMenuViewReference: {
    ViewReference(!fMenuView->IsEntryChecked(kMenuViewReference));
    UpdateHisto();
    break;
  }

  case kMenuViewStatistics: {
    ViewStatistics(!fMenuView->IsEntryChecked(kMenuViewStatistics));
    UpdateHisto();
    break;
  }

  case kMenuFavAdd: {
    if (AddFavorite()) {
      fMenuFavorites->EnableEntry(kMenuFavSave);
      fMenuFavorites->EnableEntry(kMenuFavSaveAs);
    }
    break;
  }

  case kMenuFavDelete: {
    if (DeleteFavorite()) {
      UpdateHisto();
      UpdateDescription();
      fMenuFavorites->DisableEntry(kMenuFavDelete);
      fMenuFavorites->EnableEntry(kMenuFavSave);
      fMenuFavorites->EnableEntry(kMenuFavSaveAs);
    }
    break;
  }

  case kMenuFavLoad: {
    if (LoadFavorites()) {
      UpdateHisto();
      UpdateDescription();
      fMenuFavorites->DisableEntry(kMenuFavSave);
      fMenuFavorites->EnableEntry(kMenuFavSaveAs);
    }
    break;
  }

  case kMenuFavSave: {
    if (SaveFavorites()) {
      fMenuFavorites->DisableEntry(kMenuFavSave);
    }
    break;
  }

  case kMenuFavSaveAs: {
    if (SaveFavoritesAs()) {
      fMenuFavorites->DisableEntry(kMenuFavSave);
    }
    break;
  }

  case kMenuFavSaveOnExit: {
    if (fMenuFavorites->IsEntryChecked(kMenuFavSaveOnExit)) {
      fMenuFavorites->UnCheckEntry(kMenuFavSaveOnExit);
    } else {
      fMenuFavorites->CheckEntry(kMenuFavSaveOnExit);
    }
    break;
  }

  case kMenuRefLoad: {
    if (LoadReference()) {
      SetReference();
      UpdateHisto();
      UpdateComparisonTree();
      fMenuReference->EnableEntry(kMenuRefSaveAs);
    }
    break;
  }

  case kMenuRefThreshold: {
    Float_t threshold = AliMonitorHisto::GetThreshold();
    new AliMonitorNumberDlg(threshold, this, 
			    "Comparison with Reference Histograms",
			    "threshold for comparison:", 0.);
    if (threshold < 0) break;

    AliMonitorHisto::SetThreshold(threshold);
    UpdateHisto();
    UpdateComparisonTree();
    break;
  }

  case kMenuRefTakeCurrent: {
    if (TakeCurrentReference()) {
      UpdateHisto();
      UpdateComparisonTree();
      fMenuReference->EnableEntry(kMenuRefSave);
      fMenuReference->EnableEntry(kMenuRefSaveAs);
    }
    break;
  }

  case kMenuRefSave: {
    if (SaveReference()) {
      fMenuReference->DisableEntry(kMenuRefSave);
    }
    break;
  }

  case kMenuRefSaveAs: {
    if (SaveReferenceAs()) {
      fMenuReference->DisableEntry(kMenuRefSave);
    }
    break;
  }

  case kMenuOptLoop: {
    Float_t interval = fLoopInterval * 0.001;
    new AliMonitorNumberDlg(interval, this, "Loop Interval",
			    "loop time in seconds:", 0.1);
    if (interval < 0) break;

    fLoopInterval = Int_t(1000 * interval);
    if (fLoopTimer) {
      fLoopTimer->Stop();
      fLoopTimer->Start(fLoopInterval);
    }
    break;
  }

  case kMenuOptPrint: {
    TString printCommand(fPrintCommand);
    new AliMonitorStringDlg(printCommand, this, "Print Command",
			    "shell command for printing:");
    if (printCommand.IsNull()) break;

    fPrintCommand = printCommand;
    break;
  }

  case kMenuOptSaveOnExit: {
    if (fMenuOptions->IsEntryChecked(kMenuOptSaveOnExit)) {
      fMenuOptions->UnCheckEntry(kMenuOptSaveOnExit);
    } else {
      fMenuOptions->CheckEntry(kMenuOptSaveOnExit);
    }
    break;
  }

  case kMenuOptSave: {
    SaveSettings();
    break;
  }

  case kMenuHelpAbout: {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), this, 
		 "About", "AliMonitorClient $Revision$", 
		 kMBIconAsterisk, kMBOk, &result);
    break;
  }

  default: {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), this, 
		 "AliMonitorClient", "not yet implemented", 
		 kMBIconExclamation, kMBOk, &result);
  }
  }
}


//_____________________________________________________________________________
void AliMonitorClient::OnEventNumberChanged()
{
// called when the event number in the tool button was changed

  if (!fRunButton->IsDown()) {
    UpdateHisto();
    UpdateComparisonTree();
  }
}

//_____________________________________________________________________________
void AliMonitorClient::OnEventButtonPressed()
{
// called when the event tool button was clicked

  fSumButton->SetDown(kFALSE);
  fRunButton->SetDown(kFALSE);
  UpdateHisto();
  UpdateComparisonTree();
}

//_____________________________________________________________________________
void AliMonitorClient::OnSumButtonPressed()
{
// called when the sum tool button was clicked

  fEventButton->SetDown(kFALSE);
  fRunButton->SetDown(kFALSE);
  UpdateHisto();
  UpdateComparisonTree();
}

//_____________________________________________________________________________
void AliMonitorClient::OnRunButtonPressed()
{
// called when the run tool button was clicked

  fEventButton->SetDown(kFALSE);
  fSumButton->SetDown(kFALSE);
  UpdateHisto();
  UpdateComparisonTree();
}

//_____________________________________________________________________________
void AliMonitorClient::OnLoopButtonClicked()
{
// called when the loop tool button was clicked

  // if loop is not running, start the loop timer and 
  // set the stop loop button picture
  if (!fLoopTimer) {  
    if (!fBaseItem) {
      if (!GetBaseItem()) return;
    }
    fLoopTimer = new TTimer(fLoopInterval);
    fLoopTimer->Connect("Timeout()", "AliMonitorClient", this, 
			"OnLoopTimer()");
    ((TGPictureButton*) fLoopButton)->SetPicture(fLoopOffPicture);
    fLoopTimer->TurnOn();

  // if the loop is running, stop it
  } else {
    StopLoop();
  }
}

//_____________________________________________________________________________
void AliMonitorClient::OnPreviousButtonClicked()
{
// called when the previous tool button was clicked

  if (!fBaseItem) {
    if (!GetBaseItem()) return;
  }
  StopLoop();
  GoToPreviousItem();
}

//_____________________________________________________________________________
void AliMonitorClient::OnNextButtonClicked()
{
// called when the next tool button was clicked

  if (!fBaseItem) {
    if (!GetBaseItem()) return;
  }
  StopLoop();
  GoToNextItem();
}

//_____________________________________________________________________________
void AliMonitorClient::OnCopyButtonClicked()
{
// called when the copy tool button was clicked

  fDrawCanvas->GetCanvas()->DrawClone();
}

//_____________________________________________________________________________
void AliMonitorClient::OnSaveButtonClicked()
{
// called when the save tool button was clicked

  // display a file save dialog
  static TGFileInfo fileInfo;
  static const char* fileTypes[] = {"PostScript",   "*.ps",
				    "Encapsulated PostScript", "*.eps",
				    "SVG",          "*.svg",
				    "Gif files",    "*.gif",
				    "Macro files",  "*.C",
				    "ROOT files", "*.root", 
				    "All files",  "*", 
				    NULL,         NULL};
  fileInfo.fFileTypes = fileTypes;
  new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fileInfo);
  if (!fileInfo.fFilename) return;

  fDrawCanvas->GetCanvas()->SaveAs(fileInfo.fFilename);
}

//_____________________________________________________________________________
void AliMonitorClient::OnPrintButtonClicked()
{
// called when the print tool button was clicked

  // save the canvas to a temporary postscript file
  char fileName[L_tmpnam+4];
  sprintf(fileName, "%s.ps", tmpnam(NULL));
  fDrawCanvas->GetCanvas()->SaveAs(fileName);

  // call the print command and delete the temporary file
  char command[256];
  sprintf(command, "%s %s", fPrintCommand.Data(), fileName);
  gSystem->Exec(command);
  gSystem->Unlink(fileName);
}


//_____________________________________________________________________________
void AliMonitorClient::OnTreeClicked(TGListTreeItem* item, Int_t)
{
// called when an item in the histogram tree is clicked

  OnTreeReturnPressed(item);
}

//_____________________________________________________________________________
void AliMonitorClient::OnTreeReturnPressed(TGListTreeItem* item)
{
// called when return is pressed at an item in the histogram tree

  fCurrentItem = item;
  fBaseItem = NULL;
  StopLoop();
  UpdateItem(kFALSE);
}


//_____________________________________________________________________________
void AliMonitorClient::OnLoopTimer()
{
// called by the loop timer when a new histogram should be displayed

  if (!fBaseItem) {
    if (!GetBaseItem()) return;
  }
  GoToNextItem();
}


//_____________________________________________________________________________
TFolder* AliMonitorClient::CreateTopFolder() const
{
// create a top folder for monitor histograms

  return (new TFolder("Monitor", "monitor histograms"));
}

//_____________________________________________________________________________
AliMonitorHisto* AliMonitorClient::GetHisto(const char* folderName, 
					    const char* histoName)
{
// find the monitor histogram with the given name in the given folder

  TGListTreeItem* folderItem = fTree->FindChildByName(fAllItem, folderName);
  if (folderItem) {
    TGListTreeItem* histoItem = fTree->FindChildByName(folderItem, histoName);
    if (histoItem) return (AliMonitorHisto*) histoItem->GetUserData();
  }
  return NULL;
}

//_____________________________________________________________________________
TGListTreeItem* AliMonitorClient::GetItem(TGListTreeItem* base, 
					  const char* folderName, 
					  const char* histoName, 
					  Bool_t create)
{
// find the tree item with given name in the given folder
// if create is kTRUE it is created if it is not there

  // get or create the folder
  TGListTreeItem* folderItem = fTree->FindChildByName(base, folderName);
  if (!folderItem) {
    if (!create) return NULL;
    folderItem = fTree->AddItem(base, folderName);
  }

  // get or create the histo
  TGListTreeItem* histoItem = fTree->FindChildByName(folderItem, histoName);
  if (!histoItem) {
    if (!create) return NULL;
    histoItem = fTree->AddItem(folderItem, histoName,
			       fClient->GetPicture("h1_t.xpm"),
			       fClient->GetPicture("h1_t.xpm"));
  }
  return histoItem;
}


//_____________________________________________________________________________
Bool_t AliMonitorClient::ConnectToServer()
{
// display the dialog for the server name or ip and try to connect to it

  TString serverName(fServerName);

  do {
    // ask for the server name or ip
    new AliMonitorStringDlg(serverName, this, "Connection to monitor process",
			    "monitor server name or ip:");
    if (serverName.IsNull()) return kFALSE;

    // connect to the server
    fSocket = new TSocket(serverName, AliMonitorProcess::GetPort());
    if (!fSocket || !fSocket->IsValid() || (fSocket->Send("client") <= 0)) {
      if (fSocket) delete fSocket;
      fSocket = NULL;
      Int_t result;
      new TGMsgBox(gClient->GetRoot(), this, "Connect",
		   "connection to monitor server failed", 
		   kMBIconExclamation, kMBOk, &result);

    } else {  // set up a handler for notifying when new data is there
      fServerName = serverName;
      fSocketHandler = new TFileHandler(fSocket->GetDescriptor(), 
					TFileHandler::kRead);
      fSocketHandler->Connect("Notified()", "AliMonitorClient", this,
			      "OnNewData()");
      fSocketHandler->Add();
      TInetAddress adr = fSocket->GetInetAddress();
      AliInfo(Form("connected to server: %s (%s), port %d",
		   adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
      fMenuFile->HideEntry(kMenuFileConnect);
      fMenuFile->HideEntry(kMenuFileOpen);
      fMenuFile->EnableEntry(kMenuFileDisconnect);
    };

  } while (!fSocket);

  // delete the old monitor histos
  if (fFolder) delete fFolder;
  fFolder = CreateTopFolder();
  return kTRUE;
}

//_____________________________________________________________________________
void AliMonitorClient::DisconnectFromServer()
{
// disconnect from the monitor server

  // are-you-really-sure-dialog
  Int_t result;
  new TGMsgBox(gClient->GetRoot(), this, "Disconnect",
	       "Do you really want to disconnect from the monitor server?", 
	       kMBIconQuestion, kMBYes | kMBNo, &result);
  if (result == kMBNo) return;

  // disconnect from the server
  fSocketHandler->Remove();
  fSocket->Send("disconnect");
  fSocket->Close();
  delete fSocket;
  fSocket = NULL;
  fMenuFile->HideEntry(kMenuFileDisconnect);
  fMenuFile->EnableEntry(kMenuFileConnect);
  fMenuFile->EnableEntry(kMenuFileOpen);
}

//_____________________________________________________________________________
Bool_t AliMonitorClient::OpenFile()
{
// open a file with monitor histograms

  // display a file open dialog
  static TGFileInfo fileInfo;
  static const char* fileTypes[] = {"ROOT files", "*.root", 
				    "All files",  "*", 
				    NULL,         NULL};
  fileInfo.fFileTypes = fileTypes;
  new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fileInfo);
  if (!fileInfo.fFilename) return kFALSE;

  // open the root file
  TFile* file = TFile::Open(fileInfo.fFilename);
  if (!file || !file->IsOpen()) {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), this, "Open",
		 "The file could not be opened.", 
		 kMBIconExclamation, kMBOk, &result);
    return kFALSE;
  }

  // get the folder with the monitor histograms
  TFolder* folder = (TFolder*) file->Get("Monitor");
  if (!folder || !folder->InheritsFrom(TFolder::Class())) {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), this, "Open",
		 "The file does not contain monitor histograms.", 
		 kMBIconExclamation, kMBOk, &result);
    file->Close();
    delete file;
    return kFALSE;
  }

  // replace the old folder with the new one
  file->Clear();
  file->Close();
  delete file;
  if (fFolder) delete fFolder;
  fFolder = folder;
  return kTRUE;
}


//_____________________________________________________________________________
void AliMonitorClient::ViewToolBar(Bool_t visible)
{
// en- or disable the view of the tool bar

  if (visible) {
    ShowFrame(fToolBarSep);
    ShowFrame(fToolBar);
    fMenuView->CheckEntry(kMenuViewToolBar);
  } else {
    HideFrame(fToolBarSep);
    HideFrame(fToolBar);
    fMenuView->UnCheckEntry(kMenuViewToolBar);
  }
}

//_____________________________________________________________________________
void AliMonitorClient::ViewTree(Bool_t visible)
{
// en- or disable the view of the tree

  if (visible) {
    fHorizontalFrame->ShowFrame(fTreeFrame);
    fHorizontalFrame->ShowFrame(fTreeSplitter);
    fMenuView->CheckEntry(kMenuViewTree);
  } else {
    fHorizontalFrame->HideFrame(fTreeFrame);
    fHorizontalFrame->HideFrame(fTreeSplitter);
    fMenuView->UnCheckEntry(kMenuViewTree);
  }
}

//_____________________________________________________________________________
void AliMonitorClient::ViewDescription(Bool_t visible)
{
// en- or disable the view of the histogram description

  if (visible) {
    fVerticalFrame->ShowFrame(fDescriptionFrame);
    fVerticalFrame->ShowFrame(fDescriptionSplitter);
    fMenuView->CheckEntry(kMenuViewDescription);
  } else {
    fVerticalFrame->HideFrame(fDescriptionFrame);
    fVerticalFrame->HideFrame(fDescriptionSplitter);
    fMenuView->UnCheckEntry(kMenuViewDescription);
  }
}

//_____________________________________________________________________________
void AliMonitorClient::ViewReference(Bool_t visible)
{
// en- or disable the view of the reference histos

  if (visible) {
    AliMonitorHisto::SetDrawRef(kTRUE);
    fMenuView->CheckEntry(kMenuViewReference);
  } else {
    AliMonitorHisto::SetDrawRef(kFALSE);
    fMenuView->UnCheckEntry(kMenuViewReference);
  }
}

//_____________________________________________________________________________
void AliMonitorClient::ViewStatistics(Bool_t visible)
{
// en- or disable the view of the statistics box

  if (visible) {
    gStyle->SetOptStat(1110);
    fMenuView->CheckEntry(kMenuViewStatistics);
  } else {
    gStyle->SetOptStat(0);
    fMenuView->UnCheckEntry(kMenuViewStatistics);
  }
}


//_____________________________________________________________________________
Bool_t AliMonitorClient::AddFavorite()
{
// add the current histogram or folder to the list of favorites

  if (!fCurrentItem || !fCurrentItem->GetParent()) return kFALSE;

  // get the folder item
  TGListTreeItem* folderItem = fCurrentItem->GetParent();
  if (fCurrentItem->GetFirstChild()) folderItem = fCurrentItem;

  Bool_t result = kFALSE;

  // add a folder
  if (fCurrentItem->GetFirstChild()) {
    TGListTreeItem* histoItem = fCurrentItem->GetFirstChild();
    while (histoItem) {
      if (!GetItem(fFavoritesItem, folderItem->GetText(), 
		   histoItem->GetText(), kFALSE)) result = kTRUE;
      TGListTreeItem* item = GetItem(fFavoritesItem, folderItem->GetText(), 
				     histoItem->GetText(), kTRUE);
      item->SetUserData(histoItem->GetUserData());
      histoItem = histoItem->GetNextSibling();
    }

  // add a histo
  } else {
    if (!GetItem(fFavoritesItem, folderItem->GetText(), 
		 fCurrentItem->GetText(), kFALSE)) result = kTRUE;
    TGListTreeItem* item = GetItem(fFavoritesItem, folderItem->GetText(), 
				   fCurrentItem->GetText(), kTRUE);
    item->SetUserData(fCurrentItem->GetUserData());
  }

  if (result) gClient->NeedRedraw(fTree);
  return result;
}

//_____________________________________________________________________________
Bool_t AliMonitorClient::DeleteFavorite()
{
// delete the current histogram or folder from the list of favorites

  // check that the current item is a favorites histo or folder
  if (!fCurrentItem || !fCurrentItem->GetParent()) return kFALSE;
  TGListTreeItem* base = fCurrentItem;
  while (base->GetParent()) base = base->GetParent();
  if (base != fFavoritesItem) return kFALSE;

  // delete it
  TGListTreeItem* parent = fCurrentItem->GetParent();
  fTree->DeleteItem(fCurrentItem);
  fCurrentItem = NULL;

  // delete the parent folder if it is empty now
  if (parent->GetParent() != NULL) {
    if (!parent->GetFirstChild()) fTree->DeleteItem(parent);
  }

  gClient->NeedRedraw(fTree);
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorClient::LoadFavorites(Bool_t dialog)
{
// load the list of favorites from a file

  // display a file open dialog
  TGFileInfo fileInfo;
  static const char* fileTypes[] = {"Favorites files", "*.fav", 
				    "All files",  "*", 
				    NULL,         NULL};
  fileInfo.fFileTypes = fileTypes;
  fileInfo.fIniDir = StrDup(".");
  fileInfo.fFilename = StrDup(fFavoritesFileName.Data());
  if (dialog) {
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fileInfo);
    if (!fileInfo.fFilename) return kFALSE;
  }

  // open the text file
  FILE* file = fopen(fileInfo.fFilename, "rt");
  if (!file) {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), this, "Load Favorites",
		 "The file could not be opened.", 
		 kMBIconExclamation, kMBOk, &result);
    return kFALSE;
  }

  // delete the old favorites
  TGListTreeItem* favFolderItem = fFavoritesItem->GetFirstChild();
  while (favFolderItem) {
    TGListTreeItem* deleteItem = favFolderItem;
    favFolderItem = favFolderItem->GetNextSibling();
    fTree->DeleteItem(deleteItem);
  }

  // scan the text file and add the favorites histos
  char buffer[256];
  while (!feof(file)) {
    if (fgets(buffer, 255, file) == NULL) break;
    char* folder = strtok(buffer, "/");
    char* item = strtok(NULL, "\n");
    if (item[strlen(item)-1] == '\n') item[strlen(item)-1] = 0;
    if (!folder || !item) continue;

    AliMonitorHisto* histo = GetHisto(folder, item);
    TGListTreeItem* histoItem = GetItem(fFavoritesItem, folder, item, kTRUE);
    histoItem->SetUserData(histo);
  }
  fclose(file);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorClient::SaveFavorites()
{
// save the list of favorites in a file with the name fFavoritesFileName

  // if no file name is given use a default one
  if (fFavoritesFileName.IsNull()) {
    char* fileName = gSystem->ConcatFileName(gSystem->HomeDirectory(), 
					     "AliMonitorClient.fav");
    fFavoritesFileName = fileName;
    free(fileName);
  }

  // open the text file
  FILE* file = fopen(fFavoritesFileName.Data(), "wt");
  if (!file) {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), this, "Save Favorites",
		 "The file could not be opened.", 
		 kMBIconExclamation, kMBOk, &result);
    return kFALSE;
  }

  // loop over folders
  TGListTreeItem* folderItem = fFavoritesItem->GetFirstChild();
  while (folderItem) {

    // loop over histos
    TGListTreeItem* histoItem = folderItem->GetFirstChild();
    while (histoItem) {

      // write the favorites
      if (fprintf(file, "%s/%s\n", folderItem->GetText(), 
		  histoItem->GetText()) <= 0) {
	Int_t result;
	new TGMsgBox(gClient->GetRoot(), this, "Save Favorites",
		     "An error occured while sving the favorites.", 
		     kMBIconExclamation, kMBOk, &result);
	fclose(file);
	return kFALSE;
      }
      histoItem = histoItem->GetNextSibling();
    }

    folderItem = folderItem->GetNextSibling();
  }

  fclose(file);
  AliInfo(Form("favorites saved to file %s", fFavoritesFileName.Data()));
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorClient::SaveFavoritesAs()
{
// ask for a file name and save the favorites to this file

  // display a save file dialog
  TGFileInfo fileInfo;
  static const char* fileTypes[] = {"Favorites files", "*.fav", 
				    "All files",  "*", 
				    NULL,         NULL};
  fileInfo.fFileTypes = fileTypes;
  fileInfo.fIniDir = StrDup(".");
  fileInfo.fFilename = StrDup(fFavoritesFileName.Data());
  new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fileInfo);
  if (!fileInfo.fFilename) return kFALSE;

  // save the favorites
  fFavoritesFileName = fileInfo.fFilename;
  return SaveFavorites();
}


//_____________________________________________________________________________
Bool_t AliMonitorClient::LoadReference(Bool_t dialog)
{
// load reference histograms from a file

  // display a file open dialog
  TGFileInfo fileInfo;
  static const char* fileTypes[] = {"ROOT files", "*.root", 
				    "All files",  "*", 
				    NULL,         NULL};
  fileInfo.fFileTypes = fileTypes;
  fileInfo.fIniDir = StrDup(".");
  fileInfo.fFilename = StrDup(fReferenceFileName.Data());
  if (dialog) {
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fileInfo);
    if (!fileInfo.fFilename) return kFALSE;
  }

  // open the root file
  TFile* file = TFile::Open(fileInfo.fFilename);
  if (!file || !file->IsOpen()) {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), this, "Load Reference",
		 "The file could not be opened.", 
		 kMBIconExclamation, kMBOk, &result);
    return kFALSE;
  }

  // get the folder with the monitor histograms
  TFolder* folder = (TFolder*) file->Get("Monitor");
  if (!folder || !folder->InheritsFrom(TFolder::Class())) {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), this, "Load Reference",
		 "The file does not contain monitor histograms.", 
		 kMBIconExclamation, kMBOk, &result);
    file->Close();
    delete file;
    return kFALSE;
  }

  // replace the old reference folder with the new one
  fReferenceFileName = fileInfo.fFilename;
  file->Clear();
  file->Close();
  delete file;
  if (fReference) delete fReference;
  fReference = folder;
  return kTRUE;
}

//_____________________________________________________________________________
void AliMonitorClient::SetReference()
{
// set the reference histograms for all monitor histograms

  // loop over folder
  TIterator* iFolder = fFolder->GetListOfFolders()->MakeIterator();
  while (TFolder* folder = (TFolder*) iFolder->Next()) {
    TFolder* refFolder = (TFolder*) fReference->FindObject(folder->GetName());
    if (!refFolder) continue;

    // loop over histos
    TIterator* iHisto = folder->GetListOfFolders()->MakeIterator();
    while (AliMonitorHisto* histo = (AliMonitorHisto*) iHisto->Next()) {
      AliMonitorHisto* refHisto = 
	(AliMonitorHisto*) refFolder->FindObject(histo->GetName());
      if (!refHisto) continue;
      histo->SetReference(refHisto);
    }
    delete iHisto;

  }
  delete iFolder;
}

//_____________________________________________________________________________
Bool_t AliMonitorClient::TakeCurrentReference()
{
// take the current monitor histogram or folder as reference

  // check that a histo or folder is selected
  if (!fCurrentItem) return kFALSE;
  AliMonitorHisto* histo = (AliMonitorHisto*) fCurrentItem->GetUserData();
  if (!histo && !fCurrentItem->GetFirstChild()) return kFALSE;

  // confirm-dialog
  char message[256];
  if (histo) {
    sprintf(message, "Do you want to take the current histogram %s/%s "
	    "as reference histogram?", fCurrentItem->GetParent()->GetText(),
	    histo->GetName());
  } else {
    sprintf(message, "Do you want to take all current histogram of the "
	    "folder %s as reference histogram?", fCurrentItem->GetText());
  }
  Int_t result;
  new TGMsgBox(gClient->GetRoot(), this, "Take Current as Reference",
	       message, kMBIconQuestion, kMBYes | kMBNo, &result);
  if (result != kMBYes) return kFALSE;

  // take ...
  if (histo) {   // ... a histo
    TakeReferenceHisto(fCurrentItem->GetParent()->GetText(), histo);
  } else if (fCurrentItem->GetParent()) {  // ... a folder
    TakeReferenceFolder(fCurrentItem);
  } else {  // ... a top folder
    TGListTreeItem* folderItem = fCurrentItem->GetFirstChild();
    while (folderItem) {
      TakeReferenceFolder(folderItem);
      folderItem = folderItem->GetNextSibling();
    }
  }
  return kTRUE;
}

//_____________________________________________________________________________
void AliMonitorClient::TakeReferenceHisto(const char* folderName,
					  AliMonitorHisto* histo)
{
// take the given monitor histogram as reference histogram

  // get or create the reference folder
  TFolder* refFolder = (TFolder*) fReference->FindObject(folderName);
  if (!refFolder) refFolder = fReference->AddFolder(folderName, folderName);

  // delete the old reference histo
  AliMonitorHisto* refHisto = 
    (AliMonitorHisto*) refFolder->FindObject(histo->GetName());
  if (refHisto) {
    refFolder->Remove(refHisto);
    delete refHisto;
  }

  // add the new one and use it as reference
  refFolder->Add(new AliMonitorHisto(*histo));
  histo->SetReference(histo);
}

//_____________________________________________________________________________
void AliMonitorClient::TakeReferenceFolder(TGListTreeItem* item)
{
// take all monitor histogram in the given folder as reference histograms

  // loop over histos
  TGListTreeItem* histoItem = item->GetFirstChild();
  while (histoItem) {
    AliMonitorHisto* histo = (AliMonitorHisto*) histoItem->GetUserData();
    if (histo) TakeReferenceHisto(item->GetText(), histo);
    histoItem = histoItem->GetNextSibling();
  }
}

//_____________________________________________________________________________
Bool_t AliMonitorClient::SaveReference()
{
// save the reference histograms to the file with name fReferenceFileName

  // if no file name is given use a default one
  if (fFavoritesFileName.IsNull()) {
    char* fileName = gSystem->ConcatFileName(gSystem->HomeDirectory(), 
					     "AliMonitorClientRef.root");
    fFavoritesFileName = fileName;
    free(fileName);
  }

  // open the root file
  TFile* file = TFile::Open(fReferenceFileName, "RECREATE");
  if (!file || !file->IsOpen()) {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), this, "Save Reference",
		 "The file could not be opened.", 
		 kMBIconExclamation, kMBOk, &result);
    if (file) delete file;
    return kFALSE;
  }

  // write the reference folder
  fReference->Write();
  file->Close();
  delete file;
  AliInfo(Form("reference histograms saved to file %s", 
	       fReferenceFileName.Data()));
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorClient::SaveReferenceAs()
{
// ask for a file name and save the reference histograms to this file

  // display a save file dialog
  TGFileInfo fileInfo;
  static const char* fileTypes[] = {"ROOT files", "*.root", 
				    "All files",  "*", 
				    NULL,         NULL};
  fileInfo.fFileTypes = fileTypes;
  fileInfo.fIniDir = StrDup(".");
  fileInfo.fFilename = StrDup(fReferenceFileName.Data());
  new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fileInfo);
  if (!fileInfo.fFilename) return kFALSE;

  // save the references
  fReferenceFileName = fileInfo.fFilename;
  return SaveReference();
}


//_____________________________________________________________________________
void AliMonitorClient::StopLoop()
{
// stop looping over histograms

  // stop the timer and reset the loop button picture
  if (fLoopTimer) {
    fLoopTimer->TurnOff();
    delete fLoopTimer;
    fLoopTimer = NULL;
    ((TGPictureButton*) fLoopButton)->SetPicture(fLoopOnPicture);
  }
}

//_____________________________________________________________________________
void AliMonitorClient::LoadSettings()
{
// load setting from the file with the name fgSettingFileName and apply them

  // open the text file
  char* fileName = gSystem->ConcatFileName(gSystem->HomeDirectory(), 
					   fgSettingsFileName);
  FILE* file = fopen(fileName, "rt");
  if (!file) {
    free(fileName);
    return;
  }

  // scan the text file and apply the settings
  char buffer[256];
  while (!feof(file)) {
    if (fgets(buffer, 255, file) == NULL) break;
    char* token = strtok(buffer, "=");
    char* value = strtok(NULL, "\n");
    if (!token || !value) continue;
    if (value[strlen(value)-1] == '\n') value[strlen(value)-1] = 0;

    if (strcmp(token, "ServerName") == 0) {
      fServerName = value;

    } else if (strcmp(token, "ViewToolBar") == 0) {
      ViewToolBar(strcmp(value, "0") != 0);
    } else if (strcmp(token, "ViewTree") == 0) {
      ViewTree(strcmp(value, "0") != 0);
    } else if (strcmp(token, "ViewDescription") == 0) {
      ViewDescription(strcmp(value, "0") != 0);
    } else if (strcmp(token, "ViewReference") == 0) {
      ViewReference(strcmp(value, "0") != 0);
    } else if (strcmp(token, "ViewStatistics") == 0) {
      ViewStatistics(strcmp(value, "0") != 0);

    } else if (strcmp(token, "FavoritesFileName") == 0) {
      fFavoritesFileName = value;
      LoadFavorites(kFALSE);
    } else if (strcmp(token, "FavoritesSaveOnExit") == 0) {
      if (strcmp(value, "0") != 0) {
	fMenuFavorites->CheckEntry(kMenuFavSaveOnExit);
      } else {
	fMenuFavorites->UnCheckEntry(kMenuFavSaveOnExit);
      }

    } else if (strcmp(token, "ReferenceFileName") == 0) {
      fReferenceFileName = value;
      LoadReference(kFALSE);
    } else if (strcmp(token, "ReferenceThreshold") == 0) {
      AliMonitorHisto::SetThreshold(atof(value));

    } else if (strcmp(token, "LoopInterval") == 0) {
      fLoopInterval = atoi(value);
    } else if (strcmp(token, "PrintCommand") == 0) {
      fPrintCommand = value;
    } else if (strcmp(token, "SettingsSaveOnExit") == 0) {
      if (strcmp(value, "0") != 0) {
	fMenuOptions->CheckEntry(kMenuOptSaveOnExit);
      } else {
	fMenuOptions->UnCheckEntry(kMenuOptSaveOnExit);
      }
    }
  }

  fclose(file);
  AliDebug(1, Form("settings from file %s applied", fileName));
  free(fileName);
}

//_____________________________________________________________________________
void AliMonitorClient::SaveSettings()
{
// save setting to the file with the name fgSettingFileName

  // open the text file
  char* fileName = gSystem->ConcatFileName(gSystem->HomeDirectory(), 
					   fgSettingsFileName);
  FILE* file = fopen(fileName, "wt");
  if (!file) {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), this, "Save Settings",
		 "The file could not be opened.", 
		 kMBIconExclamation, kMBOk, &result);
    free(fileName);
    return;
  }

  // write the settings
  fprintf(file, "ServerName=%s\n", fServerName.Data());

  fprintf(file, "ViewToolBar=%d\n", 
	  fMenuView->IsEntryChecked(kMenuViewToolBar));
  fprintf(file, "ViewTree=%d\n", 
	  fMenuView->IsEntryChecked(kMenuViewTree));
  fprintf(file, "ViewDescription=%d\n", 
	  fMenuView->IsEntryChecked(kMenuViewDescription));
  fprintf(file, "ViewReference=%d\n", 
	  fMenuView->IsEntryChecked(kMenuViewReference));
  fprintf(file, "ViewStatistics=%d\n", 
	  fMenuView->IsEntryChecked(kMenuViewStatistics));

  if (!fFavoritesFileName.IsNull()) {
    fprintf(file, "FavoritesFileName=%s\n", fFavoritesFileName.Data());
  }
  fprintf(file, "FavoritesSaveOnExit=%d\n", 
	  fMenuFavorites->IsEntryChecked(kMenuFavSaveOnExit));

  if (!fReferenceFileName.IsNull()) {
    fprintf(file, "ReferenceFileName=%s\n", fReferenceFileName.Data());
  }
  fprintf(file, "ReferenceThreshold=%.1f\n", AliMonitorHisto::GetThreshold());

  fprintf(file, "LoopInterval=%d\n", fLoopInterval);
  fprintf(file, "PrintCommand=%s\n", fPrintCommand.Data());
  fprintf(file, "SettingsSaveOnExit=%d\n", 
	  fMenuOptions->IsEntryChecked(kMenuOptSaveOnExit));

  fclose(file);
  AliInfo(Form("settings saved to file %s", fileName));
  free(fileName);
}


//_____________________________________________________________________________
Bool_t AliMonitorClient::GetBaseItem()
{
// get the base item for looping over histograms

  if (fCurrentItem) {
    // the base item is a folder
    fBaseItem = fCurrentItem;
    // if the current item is a histo, its parent is the base folder
    if (fBaseItem->GetParent() && fBaseItem->GetParent()->GetParent()) {
      fBaseItem = fBaseItem->GetParent();
    }

  } else {  // if no item is selected the All item is the base item
    fBaseItem = fAllItem;
    fCurrentItem = fBaseItem->GetFirstChild();
    if (!fCurrentItem) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorClient::GoToNextItem()
{
// go to the next histogram in the loop

  if (!fCurrentItem) return kFALSE;
  if (!fBaseItem || !fBaseItem->GetFirstChild()) return kFALSE;

  // remember the start item to avoid an endless loop
  TGListTreeItem* startItem = fCurrentItem;

  do {
    // folder -> first child
    if (fCurrentItem->GetFirstChild()) {
      fCurrentItem = fCurrentItem->GetFirstChild();

    // histo -> next histo
    } else if ((fCurrentItem != fBaseItem) &&
	       (fCurrentItem->GetNextSibling())) {
      fCurrentItem = fCurrentItem->GetNextSibling();

    // last histo in folder -> next folder
    } else if ((fCurrentItem != fBaseItem) &&
	       (fCurrentItem->GetParent() != fBaseItem) &&
	       fCurrentItem->GetParent()->GetNextSibling()) {
      fCurrentItem = fCurrentItem->GetParent()->GetNextSibling();

    // last histo in last folder -> first folder
    } else {
      fCurrentItem = fBaseItem->GetFirstChild();
    }

    // abort if no next item found
    if (fCurrentItem == startItem) return kFALSE;

  // end loop if an item with a monitor histo was found
  } while (!fCurrentItem->GetUserData());

  UpdateItem(kTRUE);
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorClient::GoToPreviousItem()
{
// go to the previous histogram in the loop

  if (!fCurrentItem) return kFALSE;
  if (!fBaseItem || !fBaseItem->GetFirstChild()) return kFALSE;

  // remember the start item to avoid an endless loop
  TGListTreeItem* startItem = fCurrentItem;

  do {
    // folder -> last child
    if (fCurrentItem->GetFirstChild()) {
      fCurrentItem = fCurrentItem->GetFirstChild();
      while (fCurrentItem->GetNextSibling()) {
	fCurrentItem = fCurrentItem->GetNextSibling();
      }

    // histo -> previous histo
    } else if ((fCurrentItem != fBaseItem) && 
	       (fCurrentItem->GetPrevSibling())) {
      fCurrentItem = fCurrentItem->GetPrevSibling();

    // first histo in folder -> previous folder
    } else if ((fCurrentItem != fBaseItem) && 
	       (fCurrentItem->GetParent() != fBaseItem) &&
	       fCurrentItem->GetParent()->GetPrevSibling()) {
      fCurrentItem = fCurrentItem->GetParent()->GetPrevSibling();

    // first histo in first folder -> last folder
    } else {
      fCurrentItem = fBaseItem->GetFirstChild();
      while (fCurrentItem->GetNextSibling()) {
	fCurrentItem = fCurrentItem->GetNextSibling();
      }
    }

    // abort if no next item found
    if (fCurrentItem == startItem) return kFALSE;

  // end loop if an item with a monitor histo was found
  } while (!fCurrentItem->GetUserData());

  UpdateItem(kTRUE);
  return kTRUE;
}

//_____________________________________________________________________________
void AliMonitorClient::UpdateItem(Bool_t highlight)
{
// update the displayed tree for currently selected item

  if (highlight) {  // highlight the selected item and open its parent folders
    TGListTreeItem* parent = fCurrentItem->GetParent();
    while (parent) {
      if (!parent->IsOpen()) fTree->OpenItem(parent);
      parent = parent->GetParent();
    }
    fTree->HighlightItem(fCurrentItem);
    gClient->NeedRedraw(fTree);
  }

  // update displayed histo
  UpdateDescription();
  UpdateHisto();

  // en- or disable the "Reference/Take Current" menu item
  TGListTreeItem* base = fCurrentItem;
  while (base->GetParent()) base = base->GetParent();
  if (base != fComparisonItem) {
    fMenuReference->EnableEntry(kMenuRefTakeCurrent);
  } else {
    fMenuReference->DisableEntry(kMenuRefTakeCurrent);
  }

  // en- or disable the "Favorites/Add" and "Favorites/Delete" menu items
  if (fCurrentItem->GetParent()) {
    if (base == fFavoritesItem) {
      fMenuFavorites->DisableEntry(kMenuFavAdd);
      fMenuFavorites->EnableEntry(kMenuFavDelete);
    } else {
      fMenuFavorites->EnableEntry(kMenuFavAdd);
      fMenuFavorites->DisableEntry(kMenuFavDelete);
    }
  } else {
    fMenuFavorites->DisableEntry(kMenuFavAdd);
    fMenuFavorites->DisableEntry(kMenuFavDelete);
  }
}


//_____________________________________________________________________________
Bool_t AliMonitorClient::CheckForNewData()
{
// check whether the monitor process server sent new data

  // disable the socket handler in this method
  if (!fSocket || !fSocket->IsValid()) return kFALSE;
  fSocketHandler->Remove();

  // receive a control message from the server
  char controlMessage[256];
  if (fSocket->Recv(controlMessage, 255) <= 0) {
    fSocketHandler->Add();
    return kFALSE;
  }

  // if it is new histogram data, send ok
  if ((strcmp(controlMessage, "histograms") != 0) ||
      (fSocket->Send("ok") <= 0)) {
    fSocketHandler->Add();
    return kFALSE;
  }

  // get the histogram data
  TMessage* message = NULL;
  if (fSocket->Recv(message) <= 0) {
    fSocketHandler->Add();
    return kFALSE;
  }

  // replace the old folder of monitor histos with the new one
  if (message->GetClass()->InheritsFrom(TFolder::Class())) {
    if (fFolder) delete fFolder;
    fFolder = (TFolder*) message->ReadObject(message->GetClass());
    delete message;
    fSocketHandler->Add();
    return kTRUE;
  }

  delete message;
  fSocketHandler->Add();
  return kFALSE;
}

//_____________________________________________________________________________
void AliMonitorClient::ClearItems(TGListTreeItem* base) const
{
// remove the references to the histograms from all subitems of the 
// given tree item

  // loop over folders
  TGListTreeItem* folderItem = base->GetFirstChild();
  while (folderItem) {

    // loop over histos
    TGListTreeItem* histoItem = folderItem->GetFirstChild();
    while (histoItem) {
      histoItem->SetUserData(NULL);
      histoItem = histoItem->GetNextSibling();
    }

    folderItem = folderItem->GetNextSibling();
  }
}

//_____________________________________________________________________________
void AliMonitorClient::CleanUpTree(TGListTreeItem* base)
{
// remove items without monitor histograms and 
// folders without monitor histograms

  // loop over folders
  TGListTreeItem* folderItem = base->GetFirstChild();
  while (folderItem) {

    // loop over histos
    TGListTreeItem* histoItem = folderItem->GetFirstChild();
    while (histoItem) {
      TGListTreeItem* deleteItem = NULL;
      if (!histoItem->GetUserData()) deleteItem = histoItem;
      histoItem = histoItem->GetNextSibling();
      if (fCurrentItem == deleteItem) fCurrentItem = NULL;
      if (deleteItem) fTree->DeleteItem(deleteItem);
    }

    folderItem = folderItem->GetNextSibling();
  }

  // loop over folders and remove empty folders
  folderItem = base->GetFirstChild();
  while (folderItem) {
    TGListTreeItem* deleteItem = NULL;
    if (!folderItem->GetFirstChild()) deleteItem = folderItem;
    folderItem = folderItem->GetNextSibling();
    if (fCurrentItem == deleteItem) fCurrentItem = NULL;
    if (deleteItem) fTree->DeleteItem(deleteItem);
  }
}

//_____________________________________________________________________________
void AliMonitorClient::UpdateTree()
{
// update the tree of monitor histograms

  // delete references to old monitor histograms
  ClearItems(fAllItem);

  // loop over folder
  TIterator* iFolder = fFolder->GetListOfFolders()->MakeIterator();
  while (TFolder* folder = (TFolder*) iFolder->Next()) {

    // loop over histos
    TIterator* iHisto = folder->GetListOfFolders()->MakeIterator();
    while (AliMonitorHisto* histo = (AliMonitorHisto*) iHisto->Next()) {

      // add new monitor histograms
      TGListTreeItem* histoItem = GetItem(fAllItem, folder->GetName(),
					  histo->GetName(), kTRUE);
      histoItem->SetUserData(histo);
    }
    delete iHisto;

  }
  delete iFolder;

  // remove items and folders without monitor histograms
  CleanUpTree(fAllItem);

  gClient->NeedRedraw(fTree);
}

//_____________________________________________________________________________
void AliMonitorClient::UpdateFavoritesTree()
{
// update the tree of favorite monitor histograms

  // loop over folders
  TGListTreeItem* folderItem = fFavoritesItem->GetFirstChild();
  while (folderItem) {

    // loop over histos
    TGListTreeItem* histoItem = folderItem->GetFirstChild();
    while (histoItem) {

      // set monitor histo
      histoItem->SetUserData(GetHisto(folderItem->GetText(), 
				      histoItem->GetText()));
      histoItem = histoItem->GetNextSibling();
    }

    folderItem = folderItem->GetNextSibling();
  }
}

//_____________________________________________________________________________
void AliMonitorClient::UpdateComparisonTree()
{
// update the tree of monitor histograms with significant deviation
// from the reference histograms

  if (!fFolder) return;

  // delete references to old monitor histograms
  ClearItems(fComparisonItem);

  // add monitor histograms where the comparison returns a deviation
  // loop over folders
  TIterator* iFolder = fFolder->GetListOfFolders()->MakeIterator();
  while (TFolder* folder = (TFolder*) iFolder->Next()) {

    // loop over histos
    TIterator* iHisto = folder->GetListOfFolders()->MakeIterator();
    while (AliMonitorHisto* histo = (AliMonitorHisto*) iHisto->Next()) {

      // compare the histo to its reference
      Bool_t comparison = kTRUE;
      if (fEventButton->IsDown()) {
	comparison = histo->CompareEvent(fEventNumber->GetIntNumber());
      } else if (fSumButton->IsDown()) {
	comparison = histo->CompareSum(fEventNumber->GetIntNumber());
      } else {
	comparison = histo->CompareRun();
      }

      // add it to the comparison tree in case of a bad comparison result
      if (!comparison) {
	TGListTreeItem* histoItem = GetItem(fComparisonItem, folder->GetName(),
					    histo->GetName(), kTRUE);
	histoItem->SetUserData(histo);
      }
    }
    delete iHisto;

  }
  delete iFolder;

  // remove items and folders without monitor histograms
  CleanUpTree(fComparisonItem);

  gClient->NeedRedraw(fTree);
}

//_____________________________________________________________________________
void AliMonitorClient::UpdateDescription()
{
// update the description of the current monitor histogram

  fDescription->Clear();
  AliMonitorHisto* histo = NULL;
  if (fCurrentItem) histo = (AliMonitorHisto*) fCurrentItem->GetUserData();
  if (histo) fDescription->LoadBuffer(histo->GetDescription().Data());

  gClient->NeedRedraw(fDescription);
}

//_____________________________________________________________________________
void AliMonitorClient::UpdateHisto()
{
// update the current monitor histogram

  // clear the canvas if no histo is selected
  AliMonitorHisto* histo = NULL;
  if (fCurrentItem) histo = (AliMonitorHisto*) fCurrentItem->GetUserData();
  if (!histo) {
    fDrawCanvas->GetCanvas()->Clear();
    fDrawCanvas->GetCanvas()->Update();
    return;
  }

  // draw the histo for a single event or summed over several events or a run
  fDrawCanvas->GetCanvas()->cd();
  if (fEventButton->IsDown()) {
    histo->DrawEvent(fEventNumber->GetIntNumber());
  } else if (fSumButton->IsDown()) {
    histo->DrawSum(fEventNumber->GetIntNumber());
  } else {
    histo->DrawRun();
  }
}

//_____________________________________________________________________________
void AliMonitorClient::UpdateAll()
{
// update the trees, the histogram and the description

  UpdateTree();
  UpdateFavoritesTree();
  UpdateComparisonTree();
  UpdateDescription();
  UpdateHisto();
}

