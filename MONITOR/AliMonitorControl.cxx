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
//  This class provides a graphical user interface for the monitor process   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliMonitorControl.h"
#include "AliMonitorHisto.h"
#include "AliLog.h"
#include <TGNumberEntry.h>
#include <TGTextView.h>
#include <TGMsgBox.h>
#include <TSystem.h>
#include <TSocket.h>
#include <TGFrame.h>
#include <TGMenu.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TTimer.h>
#include <TApplication.h>
#include "AliMonitorProcess.h"


ClassImp(AliMonitorControl) 



//_____________________________________________________________________________
AliMonitorControl::AliMonitorBufferDlg::AliMonitorBufferDlg(Int_t& size, 
							    TGFrame* main) :
  AliMonitorDialog(main, 250, 80), fSize(size)
{
// create a dialog for setting the size of the buffer for monitor histos

  fBufferLayout = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 2, 2);
  fBufferLabel = new TGLabel(fFrame, "size of histogram buffer:");
  fFrame->AddFrame(fBufferLabel, fBufferLayout);
  fBufferEntry = new TGNumberEntry(fFrame, size, 2, -1,
				   TGNumberFormat::kNESInteger);
  fBufferEntry->SetLimits(TGNumberFormat::kNELLimitMinMax, 1, 99);
  fFrame->AddFrame(fBufferEntry, fBufferLayout);

  fSize = -1;

  fMain->SetWindowName("Buffer Size");
  fMain->MapSubwindows();
  fMain->Layout();
  gClient->WaitFor(fMain);
}

//_____________________________________________________________________________
AliMonitorControl::AliMonitorBufferDlg::~AliMonitorBufferDlg()
{
// clean up

  delete fBufferLabel;
  delete fBufferLayout;
  delete fBufferEntry;
}

//_____________________________________________________________________________
void AliMonitorControl::AliMonitorBufferDlg::OnOkClicked()
{
  fSize = fBufferEntry->GetIntNumber();
}


//_____________________________________________________________________________
AliMonitorControl::AliMonitorClientsDlg::AliMonitorClientsDlg(TObjArray* clients, 
							      TGFrame* main) :
  AliMonitorDialog(main, 450, 300, kFALSE)
{
// create a dialog to display the list of clients

  delete fFrameLayout;
  fFrameLayout = new TGLayoutHints(kLHintsCenterX | kLHintsTop, 
				   10, 10, 15, 15);
  ((TGFrameElement*)(fMain->GetList()->First()))->fLayout = fFrameLayout;
  fClientsLayout = new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 
				     0, 0, 0, 0);
  fClients = new TGTextView(fFrame, 420, 230);
  fFrame->AddFrame(fClients, fClientsLayout);

  char line[256];
  const char* format = "%-30s%-20s%-8d";
  sprintf(line, format, "name", "address", 0);
  strncpy(&line[50], "port", 4);
  fClients->AddLine(line);
  for (Int_t i = 0; i < (Int_t) strlen(line); i++) line[i] = '-';
  fClients->AddLine(line);

  for (Int_t iClient = 0; iClient < clients->GetEntriesFast(); iClient++) {
    TSocket* socket = (TSocket*) clients->At(iClient);
    if (!socket) continue;
    TInetAddress adr = socket->GetInetAddress();
    sprintf(line, format,
	    adr.GetHostName(), adr.GetHostAddress(), adr.GetPort());
    fClients->AddLine(line);
  }

  fMain->SetWindowName("List of Clients");
  fMain->MapSubwindows();
  fMain->Layout();
  gClient->WaitFor(fMain);
}

//_____________________________________________________________________________
AliMonitorControl::AliMonitorClientsDlg::~AliMonitorClientsDlg()
{
// clean up

  delete fClientsLayout;
  delete fClients;
}


// constants for menu entries
enum {kMenuFileExit, kMenuFileAbort,
      kMenuOptBuffer, kMenuOptClients,
      kMenuHelpDoc, kMenuHelpAbout
};


//_____________________________________________________________________________
AliMonitorControl::AliMonitorControl(AliMonitorProcess* process)
{
// initialize the monitoring control window

  fMonitorProcess = process;


  // colors
  gClient->GetColorByName("lightblue", fColorStatus);
  gClient->GetColorByName("green", fColorStart);
  gClient->GetColorByName("red", fColorStop);

  // main window
  fMain = new TGMainFrame(gClient->GetRoot(), 380, 200);

  // menu bar
  fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
                                      0, 0, 1, 1);
  fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
  fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

  fMenuFile = new TGPopupMenu(gClient->GetRoot());
  fMenuFile->AddEntry("E&xit", kMenuFileExit);
  fMenuFile->AddEntry("&Abort", kMenuFileAbort);
  fMenuFile->DisableEntry(kMenuFileAbort);
  fMenuFile->Connect("Activated(Int_t)", "AliMonitorControl", this,
		     "HandleMenu(Int_t)");

  fMenuOptions = new TGPopupMenu(gClient->GetRoot());
  fMenuOptions->AddEntry("&Histogram buffer...", kMenuOptBuffer);
  fMenuOptions->AddEntry("List of &Clients...", kMenuOptClients);
  fMenuOptions->Connect("Activated(Int_t)", "AliMonitorControl", this,
			"HandleMenu(Int_t)");

  fMenuHelp = new TGPopupMenu(gClient->GetRoot());
  fMenuHelp->AddEntry("&Documentation...", kMenuHelpDoc);
  fMenuHelp->AddEntry("A&bout...", kMenuHelpAbout);
  fMenuHelp->Connect("Activated(Int_t)", "AliMonitorControl", this,
		     "HandleMenu(Int_t)");

  fMenuBar = new TGMenuBar(fMain, 1, 1, kHorizontalFrame);
  fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
  fMenuBar->AddPopup("&Options", fMenuOptions, fMenuBarItemLayout);
  fMenuBar->AddPopup("&Help", fMenuHelp, fMenuBarHelpLayout);

  fMain->AddFrame(fMenuBar, fMenuBarLayout);


  // status frame
  fFrameLayout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 
				   2, 2, 2, 2);
  fFrame = new TGVerticalFrame(fMain, 0, 0, kChildFrame | kSunkenFrame);
  fMain->AddFrame(fFrame, fFrameLayout);

  fStatusLayout = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2, 2, 2);
  fStatusFrameLayout = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 2, 2);


  // run and event number
  fStatus1Frame = new TGHorizontalFrame(fFrame, 0, 0);
  fFrame->AddFrame(fStatus1Frame, fStatusFrameLayout);

  fRunNumberLabel = new TGLabel(fStatus1Frame, "current run:");
  fStatus1Frame->AddFrame(fRunNumberLabel, fStatusLayout);
  fRunNumber = new TGTextEntry(fStatus1Frame, "-");
  fRunNumber->Resize(60, fRunNumber->GetDefaultHeight());
  fRunNumber->SetAlignment(kTextRight);
  fRunNumber->SetEnabled(kFALSE);
  fRunNumber->SetBackgroundColor(fColorStatus);
  fStatus1Frame->AddFrame(fRunNumber, fStatusLayout);

  fEventNumberLabel = new TGLabel(fStatus1Frame, 
				  "  event:");
  fStatus1Frame->AddFrame(fEventNumberLabel, fStatusLayout);
  fEventNumber = new TGTextEntry(fStatus1Frame, "-/-/-");
  fEventNumber->Resize(100, fEventNumber->GetDefaultHeight());
  fEventNumber->SetAlignment(kTextRight);
  fEventNumber->SetEnabled(kFALSE);
  fEventNumber->SetBackgroundColor(fColorStatus);
  fStatus1Frame->AddFrame(fEventNumber, fStatusLayout);


  // process status
  fStatus2Frame = new TGHorizontalFrame(fFrame, 0, 0);
  fFrame->AddFrame(fStatus2Frame, fStatusFrameLayout);

  fStatusLabel = new TGLabel(fStatus2Frame, "current status:");
  fStatus2Frame->AddFrame(fStatusLabel, fStatusLayout);
  fStatus = new TGTextEntry(fStatus2Frame, "stopped");
  fStatus->Resize(250, fStatus->GetDefaultHeight());
  fStatus->SetAlignment(kTextLeft);
  fStatus->SetEnabled(kFALSE);
  fStatus->SetBackgroundColor(fColorStatus);
  fStatus2Frame->AddFrame(fStatus, fStatusLayout);


  // process status
  fStatus3Frame = new TGHorizontalFrame(fFrame, 0, 0);
  fFrame->AddFrame(fStatus3Frame, fStatusFrameLayout);

  fEventsLabel = new TGLabel(fStatus3Frame, "monitored events:");
  fStatus3Frame->AddFrame(fEventsLabel, fStatusLayout);
  fEvents = new TGTextEntry(fStatus3Frame, "-");
  fEvents->Resize(60, fEvents->GetDefaultHeight());
  fEvents->SetAlignment(kTextRight);
  fEvents->SetEnabled(kFALSE);
  fEvents->SetBackgroundColor(fColorStatus);
  fStatus3Frame->AddFrame(fEvents, fStatusLayout);

  fClientsLabel = new TGLabel(fStatus3Frame, " number of clients:");
  fStatus3Frame->AddFrame(fClientsLabel, fStatusLayout);
  fClients = new TGTextEntry(fStatus3Frame, "-");
  fClients->Resize(40, fClients->GetDefaultHeight());
  fClients->SetAlignment(kTextRight);
  fClients->SetEnabled(kFALSE);
  fClients->SetBackgroundColor(fColorStatus);
  fStatus3Frame->AddFrame(fClients, fStatusLayout);


  // buttons
  fButtonFrameLayout = new TGLayoutHints(kLHintsExpandX | kLHintsBottom, 
					 50, 50, 10, 10);
  fButtonFrame = new TGHorizontalFrame(fMain, 0, 0);
  fMain->AddFrame(fButtonFrame, fButtonFrameLayout);

  fButtonLayout = new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 
				    2, 2, 2, 2);
  fResetButton = new TGTextButton(fButtonFrame, " &Reset ", 1);
  fResetButton->Connect("Clicked()", "AliMonitorControl", this, "DoReset()");
  fButtonFrame->AddFrame(fResetButton, fButtonLayout);

  fStartStopButton = new TGTextButton(fButtonFrame, " &Start ", 2);
  fStartStopButton->SetBackgroundColor(0x00FF00);
  fStartStopButton->Connect("Clicked()", "AliMonitorControl", this, 
			    "DoStartStop()");
  fButtonFrame->AddFrame(fStartStopButton, fButtonLayout);
  fStartButtonStatus = kTRUE;


  // main window
  fMain->SetWindowName("Monitor Process Control");
  fMain->DontCallClose();
  fMain->SetWMSize(fMain->GetWidth(), fMain->GetHeight());
  fMain->SetWMSizeHints(fMain->GetWidth(), fMain->GetHeight(), 
			fMain->GetWidth(), fMain->GetHeight(), 0, 0);
  fMain->SetWMPosition(100, 100);
  fMain->MapSubwindows();
  fMain->Layout();
  fMain->MapWindow();


  fTerminating = kFALSE;

  fTimer = new TTimer(this, 10, kTRUE);
  fTimer->TurnOn();
}

//_____________________________________________________________________________
AliMonitorControl::AliMonitorControl(const AliMonitorControl& control) :
  TObject(control)
{
  AliFatal("copy constructor not implemented");
}

//_____________________________________________________________________________
AliMonitorControl& AliMonitorControl::operator = (const AliMonitorControl& 
						  /*control*/)
{
  AliFatal("assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliMonitorControl::~AliMonitorControl()
{
// clean up

  delete fMenuBarLayout;
  delete fMenuBarItemLayout;
  delete fMenuBarHelpLayout;
  delete fMenuFile;
  delete fMenuOptions;
  delete fMenuHelp;
  delete fMenuBar;

  delete fFrameLayout;
  delete fFrame;
  delete fStatusLayout;
  delete fStatusFrameLayout;

  delete fStatus1Frame;
  delete fRunNumberLabel;
  delete fRunNumber;
  delete fEventNumberLabel;
  delete fEventNumber;

  delete fStatus2Frame;
  delete fStatusLabel;
  delete fStatus;

  delete fStatus3Frame;
  delete fEventsLabel;
  delete fEvents;
  delete fClientsLabel;
  delete fClients;

  delete fButtonFrameLayout;
  delete fButtonFrame;
  delete fButtonLayout;
  delete fResetButton;
  delete fStartStopButton;

  delete fTimer;
}


//_____________________________________________________________________________
void AliMonitorControl::HandleMenu(Int_t id)
{
// called when a menu item was selected

  switch (id) {

  case kMenuFileExit: {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), fMain, "Exit", 
		 "Do you really want to terminate the monitor process?", 
		 kMBIconQuestion, kMBYes | kMBNo, &result);
    if (result == kMBYes) {
      fMenuFile->EnableEntry(kMenuFileAbort);
      if (fMonitorProcess->IsStopped()) {
	delete fMonitorProcess;
	gApplication->Terminate(0);
      }
      fMonitorProcess->Stop();
      fTerminating = kTRUE;
    }
    break;
  }

  case kMenuFileAbort:
    exit(1);

  case kMenuOptBuffer: {
    Int_t size = AliMonitorHisto::GetNHistosMax();
    new AliMonitorBufferDlg(size, fMain);
    if (size < 0) break;
    AliMonitorHisto::SetNHistosMax(size);
    break;
  }

  case kMenuOptClients: {
    new AliMonitorClientsDlg(fMonitorProcess->GetListOfClients(), fMain);
    break;
  }

  case kMenuHelpAbout: {
    Int_t result;
    char text[256];
    sprintf(text, "AliMonitorControl $Revision$\nrunning\n"
	    "AliMonitorProcess %s", AliMonitorProcess::GetRevision());
    new TGMsgBox(gClient->GetRoot(), fMain, 
		 "About", text, kMBIconAsterisk, kMBOk, &result);
    break;
  }

  default: {
    Int_t result;
    new TGMsgBox(gClient->GetRoot(), fMain, 
		 "AliMonitorControl", "not yet implemented", 
		 kMBIconExclamation, kMBOk, &result);
  }
  }
}

//_____________________________________________________________________________
void AliMonitorControl::DoReset()
{
// called when the reset button was clicked

  Int_t result;
  new TGMsgBox(gClient->GetRoot(), fMain, "Reset", 
	       "Do you really want to reset the monitor histograms?", 
	       kMBIconQuestion, kMBYes | kMBNo, &result);
  if (result == kMBYes) {
    fMonitorProcess->Reset();
  }
}

//_____________________________________________________________________________
void AliMonitorControl::DoStartStop()
{
// called when the start/stop button was clicked

  if (fStartButtonStatus) {
    if (fMonitorProcess->GetStatus() == AliMonitorProcess::kStopped) {
      fStartStopButton->SetBackgroundColor(fColorStop);
      fStartStopButton->SetText(" &Stop ");
      fStartButtonStatus = kFALSE;
      fMonitorProcess->Run();
      fStartStopButton->SetBackgroundColor(fColorStart);
      fStartStopButton->SetText(" &Start ");
      fStartButtonStatus = kTRUE;
    }

  } else {
    if (fMonitorProcess->GetStatus() != AliMonitorProcess::kStopped) {
      fMonitorProcess->Stop();
    }
  }
}



//_____________________________________________________________________________
Bool_t AliMonitorControl::HandleTimer(TTimer* timer)
{
// update the displayed information

  timer->TurnOff();
  if (fTerminating && fMonitorProcess->IsStopped()) {
    delete fMonitorProcess;
    gApplication->Terminate(0);
  }
  UpdateStatus();
  gSystem->ProcessEvents();
  timer->TurnOn();

  return kFALSE;
}


//_____________________________________________________________________________
void AliMonitorControl::UpdateStatus()
{
// update the displayed status information

  char text[256];

  if (fMonitorProcess->IsStopped()) {
    fRunNumber->SetText("-");
    fEventNumber->SetText("-/-/-");
    fEvents->SetText("-");
    fClients->SetText("-");

  } else {
    sprintf(text, "%d", fMonitorProcess->GetRunNumber());
    fRunNumber->SetText(text);
    sprintf(text, "%d/%d/%d", fMonitorProcess->GetEventPeriodNumber(),
	    fMonitorProcess->GetEventOrbitNumber(),
	    fMonitorProcess->GetEventBunchNumber());
    fEventNumber->SetText(text);
    sprintf(text, "%d", fMonitorProcess->GetNumberOfEvents());
    fEvents->SetText(text);
    sprintf(text, "%d", fMonitorProcess->GetNumberOfClients());
    fClients->SetText(text);
  }

  const char* status = NULL;
  switch (fMonitorProcess->GetStatus()) {
  case AliMonitorProcess::kStopped: 
    status = "stopped"; break;
  case AliMonitorProcess::kWaiting: 
    status = "waiting for new data"; break;
  case AliMonitorProcess::kReading: 
    status = "reading raw data"; break;
  case AliMonitorProcess::kRecTPC: 
    status = "running TPC reconstruction"; break;
  case AliMonitorProcess::kRecITS: 
    status = "running ITS reconstruction"; break;
  case AliMonitorProcess::kRecV0s: 
    status = "running V0 reconstruction"; break;
  case AliMonitorProcess::kRecHLT: 
    status = "running HLT reconstruction"; break;
  case AliMonitorProcess::kFilling: 
    status = "filling monitor histograms"; break;
  case AliMonitorProcess::kUpdating: 
    status = "updating monitor histograms"; break;
  case AliMonitorProcess::kWriting: 
    status = "writing monitor histograms"; break;
  case AliMonitorProcess::kResetting: 
    status = "resetting monitor histograms"; break;
  case AliMonitorProcess::kConnecting: 
    status = "checking for new clients"; break;
  case AliMonitorProcess::kBroadcasting: 
    status = "broadcasting monitor histograms"; break;
  }

  if (fMonitorProcess->WillStop()) {
    sprintf(text, "stopping... (%s)", status);
  } else {
    sprintf(text, "%s", status);
  }
  if (strcmp(text, fStatus->GetText()) != 0) {
    fStatus->SetText(text);
  }

  /*
  if (fStartButtonStatus != fMonitorProcess->IsStopped()) {
    fStartButtonStatus = fMonitorProcess->IsStopped();
    if (fStartButtonStatus) {
      fStartStopButton->SetBackgroundColor(fColorStart);
      fStartStopButton->SetText(" &Start ");
    } else {
      fStartStopButton->SetBackgroundColor(fColorStop);
      fStartStopButton->SetText(" &Stop ");
    }
  }
  */
}

