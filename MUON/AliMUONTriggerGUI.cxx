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

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONTriggerGUI
/// Graphical User Interface utility class for the MUON trigger detector
/// It creates, after initialisation with a data file, a sensitive map
/// of the trigger boards
/// \author Bogdan Vulpescu, LPC Clermont-Ferrand
//-----------------------------------------------------------------------------

#include "AliMUONTriggerGUI.h"
#include "AliMUONTriggerGUIboard.h"
#include "AliMUONTriggerGUIdimap.h"
#include "AliMUONTriggerGUIbdmap.h"

#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpPad.h"
#include "AliMpIntPair.h"
#include "AliMUON.h"
#include "AliMpDEIterator.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONCalibrationData.h"

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"

#include <TSystem.h>
#include <TGLabel.h>
#include <TGFrame.h>
#include <TApplication.h>
#include <TGDimension.h>
#include <TString.h>
#include <TGMenu.h>
#include <TGTextEntry.h>
#include <TGButton.h>
#include <TFile.h>
#include <TGImageMap.h>
#include <TGTextBuffer.h>

/// \cond CLASSIMP
ClassImp(AliMUONTriggerGUI)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerGUI::AliMUONTriggerGUI()
  : TObject(),
    fMain(0),
    fImageMap(0),
    fTxtBuffer1(0),
    fTxtBuffer2(0),
    fTxtCircuit(0),
    fRunInput(0),
    fError(0),
    fControl(0),
    fCircuit(0),
    fSkipToEventTxt(0),
    fFileName(0),
    fPath(0),
    fEvString(0),
    fChamber(0),
    fEvent(0),
    fEventsPerRun(0),
    fLoader(0),
    fRunLoader(0),
    fCDBManager(0),
    fCalibrationData(0),
    fBoardsInit(0),
    fDiMap(0),
    fTriggerProcessor(0),
    fBoards(0)
{
  /// main GUI frame of the trigger monitor

  fChamber = 0;
  fEvent   = 0;
  fEventsPerRun = 0;
  fRunLoader = 0;
  fDiMap = 0;

  fBoards = 0;
  fFileName = new TString("");
  fEvString = new TString("");
  fPath = new TString("");

  fTriggerProcessor = 0;

  fCDBManager = AliCDBManager::Instance();
  fCDBManager->SetDefaultStorage("local://$ALICE_ROOT");
  fCalibrationData = new AliMUONCalibrationData(0);  // runnumber zero
  
  // Main frame

  fMain = new TGMainFrame(gClient->GetRoot(), 750, 420);
  fMain->Connect("CloseWindow()", "AliMUONTriggerGUI", this, "CloseWindow()");

  // Menu bar
  
  TGMenuBar *menuBar = new TGMenuBar(fMain);
  
  // File menu

  TGPopupMenu *menuFile = new TGPopupMenu(gClient->GetRoot());
  //menuFile->AddLabel("");

  menuFile->AddEntry("Run",     kMFILERUN);
  menuFile->AddEntry("Control", kMFILECNTRL);
  menuFile->AddEntry("Exit",    kMFILEEXIT);

  menuFile->Connect("Activated(Int_t)", "AliMUONTriggerGUI", this, "HandleMenu(Int_t)");

  // Circuit menu

  TGPopupMenu *menuCircuit = new TGPopupMenu(gClient->GetRoot());
  //menuCircuit->AddLabel("");

  menuCircuit->AddEntry("Open",     kMCIRCUITOPEN);

  menuCircuit->Connect("Activated(Int_t)", "AliMUONTriggerGUI", this, "HandleMenu(Int_t)");

  // Digits map menu

  TGPopupMenu *menuMap = new TGPopupMenu(gClient->GetRoot());
  //menuMap->AddLabel("");

  menuMap->AddEntry("Digits map",   kMMAPDIGITS);
  menuMap->AddEntry("Reset digits", kMRESETDIGITS);

  menuMap->Connect("Activated(Int_t)", "AliMUONTriggerGUI", this,
		    "HandleMenu(Int_t)");

  // Trigger menu

  TGPopupMenu *menuTrigger = new TGPopupMenu(gClient->GetRoot());
  //menuTrigger->AddLabel("");

  menuTrigger->AddEntry("Trigger DSET",     kMTRIGGERDSET);

  menuTrigger->Connect("Activated(Int_t)", "AliMUONTriggerGUI", this, "HandleMenu(Int_t)");

  // Add menus to the menu bar

  menuBar->AddPopup("File", menuFile, 
		    new TGLayoutHints(kLHintsTop | kLHintsLeft, 5,5,2,2)
		    );

  menuBar->AddPopup("Maps", menuMap, 
		    new TGLayoutHints(kLHintsTop | kLHintsLeft, 5,5,2,2)
		    );

  menuBar->AddPopup("Circuit", menuCircuit, 
		    new TGLayoutHints(kLHintsTop | kLHintsLeft, 5,5,2,2)
		    );

  menuBar->AddPopup("Trigger", menuTrigger, 
		    new TGLayoutHints(kLHintsTop | kLHintsLeft, 5,5,2,2)
		    );

  // Add menu bar to the main frame

  fMain->AddFrame(menuBar,
		  new TGLayoutHints(kLHintsTop | 
				    kLHintsLeft | 
				    kLHintsExpandX, 
				    0, 0, 1, 1)
		  );
 
  // The image map

  fImageMap = new TGImageMap(fMain,"$ALICE_ROOT/MUON/data/guimap.gif");
  
  fImageMap->Connect("RegionClicked(Int_t)", "AliMUONTriggerGUI", this, "OpenBoard(Int_t)");

  fImageMap->SetToolTipText("Map of the local boards as seen from the I.P.");

  // Add image map to the main frame

  fMain->AddFrame(fImageMap);
  fMain->SetWindowName("Map of the local boards as seen from the I.P.");

  // Resize the main frame
  
  TGDimension size = fMain->GetDefaultSize();
  fMain->Resize(size);

  fMain->MapSubwindows();

  fMain->MapWindow();

  fBoardsInit = kFALSE;

  //HandleMenu(kMFILERUN);  // temporary

  //InitBoards();
  //Init();
  
}

//__________________________________________________________________________
void AliMUONTriggerGUI::HandleMenu(Int_t id)
{
  /// handles entry numbers in the available menus (EMenuIdentifiers)

  TGCompositeFrame *runInput1, *runInput2, *runInput3;
  TGCompositeFrame *control1, *control2, *circuit1, *circuit2;
  TGLabel *runL1, *runL2, *circuitL1;
  TGTextEntry *runText1, *runText2, *circuitText1;
  TGTextButton *runApply, *runCancel; 
  TGTextButton *controlClose, *nextEvent, *previousEvent, *skipToEvent;
  TGTextButton *circuitCancel, *circuitOpen;

  //Int_t trigInfo[kNBoards*6] = {-1};
  //Int_t nLocalTrigger = 0;

  TString error = TString("");
  if (id != kMFILEEXIT && id != kMFILERUN && fRunLoader == 0) {
    error.Append("Run not initialized (Menu: File/Run)");
    DoErrorGUI(error.Data());
    return;
  }

  switch (id) {

  case kMFILEEXIT:

    printf("\nBye. \n");
    CloseWindow();
    break;
    
  case kMFILERUN:
    
    // input main frame

    fRunInput = new TGTransientFrame(gClient->GetRoot(), fMain, 400, 200);
    fRunInput->Connect("CloseWindow()", "AliMUONTriggerGUI", this, "CloseRunInput()");
    fRunInput->DontCallClose(); // to avoid double deletions.

    // use hierarchical cleaning
    fRunInput->SetCleanup(kDeepCleanup);

    fRunInput->SetWindowName("Input file and event number");

    // input galice and event number frames

    runInput1 = new TGCompositeFrame(fRunInput, 400, 200, kHorizontalFrame);
    runInput2 = new TGCompositeFrame(fRunInput, 400, 200, kHorizontalFrame);

    // .. with labels

    runL1 = new TGLabel(runInput1, new TGString("Path to gAlice:"));
    runL2 = new TGLabel(runInput2, new TGString("Event number:"));

    // galice text entry

    runText1 = new TGTextEntry(runInput1, fTxtBuffer1 = new TGTextBuffer(100));

    //fPath->Append("dset");                 // temporary
    //fTxtBuffer1->AddText(0,fPath->Data());  // temporary

    runText1->SetToolTipText("Enter the path to galice.root");
    runText1->Resize(300, runText1->GetDefaultHeight());

    // event number text entry

    runText2 = new TGTextEntry(runInput2, fTxtBuffer2 = new TGTextBuffer(5));

    fEvString->Form("%d",0);               
    fTxtBuffer2->AddText(0,fEvString->Data());

    runText2->SetToolTipText("Enter the first event number to start with");
    runText2->Resize(300, runText2->GetDefaultHeight());

    // add to galice frame

    runInput1->AddFrame(runL1,
			 new TGLayoutHints(kLHintsLeft | 
					   kLHintsCenterY, 
					   3, 5, 0, 0)
			);

    runInput1->AddFrame(runText1, 
			new TGLayoutHints(kLHintsRight | 
					  kLHintsCenterY, 
					  0, 2, 0, 0)
			);

    // add to event number frame

    runInput2->AddFrame(runL2,
			 new TGLayoutHints(kLHintsLeft | 
					   kLHintsCenterY, 
					   3, 5, 0, 0)
			);

    runInput2->AddFrame(runText2, 
			new TGLayoutHints(kLHintsRight | 
					  kLHintsCenterY, 
					  0, 2, 0, 0)
			);

    // add input frames to main input frame

    fRunInput->AddFrame(runInput1, 
			new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 3, 0));
    
    fRunInput->AddFrame(runInput2, 
			new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 3, 0));
    
    // frame with buttons

    runInput3 = new TGCompositeFrame(fRunInput, 400, 200, kHorizontalFrame);

    // buttons

    runApply = new TGTextButton(runInput3, "Apply", 1);
    runApply->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoRunApply()");
    runApply->SetToolTipText("Apply changes");

    runCancel = new TGTextButton(runInput3, "Cancel", 2);
    runCancel->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoRunCancel()");    runCancel->SetToolTipText("Cancel changes");

    // add buttons

    runInput3->AddFrame(runApply, 
			 new TGLayoutHints(kLHintsTop | 
					   kLHintsLeft, 
					   3, 3, 2, 2)
			);
    
    runInput3->AddFrame(runCancel, 
			 new TGLayoutHints(kLHintsTop | 
					   kLHintsLeft, 
					   3, 3, 2, 2)
			);
    
    // add to the main input frame

    fRunInput->AddFrame(runInput3, 
			new TGLayoutHints(kLHintsTop | 
					  kLHintsExpandX, 
					  2, 2, 3, 0)
			);
    
    fRunInput->MapSubwindows();
    fRunInput->Resize();
    fRunInput->MapWindow();

    //DoRunApply();   // temporary

    break;
    
  case kMFILECNTRL:

    // control main frame

    fControl = new TGTransientFrame(gClient->GetRoot(), fMain, 50, 50);
    fControl->Connect("CloseWindow()", "AliMUONTriggerGUI", this, "CloseControl()");
    fControl->DontCallClose(); // to avoid double deletions.
  
    // use hierarchical cleaning
    fControl->SetCleanup(kDeepCleanup);
  
    fControl->SetWindowName("Run controls");

    // frame to hold buttons

    control1 = new TGCompositeFrame(fControl, 50, 50, kVerticalFrame);
  
    // buttons

    controlClose = new TGTextButton(control1, "Close", 1);
    controlClose->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoControlClose()");
    //controlClose->Resize(300, controlClose->GetDefaultHeight());  

    nextEvent = new TGTextButton(control1, "Next event", 2);
    nextEvent->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoNextEvent()");

    previousEvent = new TGTextButton(control1, "Previous event", 3);
    previousEvent->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoPreviousEvent()");
  
    // frame to hold event skip

    control2 = new TGCompositeFrame(fControl, 50, 50, kHorizontalFrame);
  
    // skip to event text field

    fSkipToEventTxt = new TGTextEntry(control2, fTxtBuffer2 = new TGTextBuffer(5));

    fTxtBuffer2->AddText(0,fEvString->Data());

    // skip to event button

    skipToEvent = new TGTextButton(control2, "Skip to event", 1);
    skipToEvent->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoSkipToEvent()");

    // add to event skip frame

    control2->AddFrame(fSkipToEventTxt, 
		       new TGLayoutHints(kLHintsTop | 
					 kLHintsRight, 
					 10, 10, 5, 5)
		       );
  
    control2->AddFrame(skipToEvent, 
		       new TGLayoutHints(kLHintsTop | 
					 kLHintsRight, 
					 10, 10, 5, 5)
		       );
  
    // add buttons
  
    control1->AddFrame(controlClose, 
		 new TGLayoutHints(kLHintsBottom | 
				   kLHintsExpandX, 
				   10, 10, 5, 5)
		 );
  
    control1->AddFrame(nextEvent, 
		 new TGLayoutHints(kLHintsBottom | 
				   kLHintsExpandX, 
				   10, 10, 5, 5)
		 );
  
    control1->AddFrame(previousEvent, 
		 new TGLayoutHints(kLHintsBottom | 
				   kLHintsExpandX, 
				   10, 10, 5, 5)
		 );
  
    // add to the main frame

    fControl->AddFrame(control1,
		       new TGLayoutHints(kLHintsBottom |
					 kLHintsLeft | 
					 kLHintsCenterY, 
					 3, 5, 0, 0)
		       );
    
    fControl->AddFrame(control2,
		       new TGLayoutHints(kLHintsTop | 
					 kLHintsCenterX, 
					 3, 5, 0, 0)
		       );
    
    fControl->MapSubwindows();
    fControl->Resize();
    fControl->MapWindow();
   
    break;

  case kMMAPDIGITS:
    
    if (fDiMap == 0) {
      fDiMap = new AliMUONTriggerGUIdimap(fLoader,fBoards,gClient->GetRoot(), fMain, 400, 200);
    } else if (!fDiMap->IsOn()) {
      fDiMap = new AliMUONTriggerGUIdimap(fLoader,fBoards,gClient->GetRoot(), fMain, 400, 200);
    }

    break;
    
  case kMRESETDIGITS:
    
    Int_t number, over, pos;
    for (Int_t ib = 0; ib < kNBoards; ib++) {
      AliMUONTriggerGUIboard *board = GetBoard(ib);
      board->ClearXDigits();
      board->ClearYDigits();
      // extended y-strips
      number = board->GetNumber();
      pos    = board->GetPosition();
      over   = board->GetYOver();
      for (Int_t io = 1; io <= over; io++) {
	if (io == pos) continue;
	board = GetBoard(number+io-pos);
	board->ClearYDigits();
      }
    }

    break;
    
  case kMCIRCUITOPEN:

    fCircuit = new TGTransientFrame(gClient->GetRoot(), fMain, 50, 50);
    fCircuit->Connect("CloseWindow()", "AliMUONTriggerGUI", this, "CloseCircuit()");
    fCircuit->DontCallClose(); // to avoid double deletions.
  
    // use hierarchical cleaning
    fCircuit->SetCleanup(kDeepCleanup);
  
    fCircuit->SetWindowName("Board circuit");

    // sub-frames

    circuit1 = new TGCompositeFrame(fCircuit, 400, 200, kHorizontalFrame);
    circuit2 = new TGCompositeFrame(fCircuit, 400, 200, kHorizontalFrame);

    // labels

    circuitL1 = new TGLabel(circuit1, new TGString("Circuit number:"));
    
    // text entry
    
    circuitText1 = new TGTextEntry(circuit1, fTxtCircuit = new TGTextBuffer(10));
    // buttons

    circuitCancel = new TGTextButton(circuit2, "Cancel", 1);
    circuitCancel->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoCircuitCancel()");
    //circuitCancel->Resize(100, circuitCancel->GetDefaultHeight());  

    circuitOpen = new TGTextButton(circuit2, "Open", 2);
    circuitOpen->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoCircuitOpen()");
    //circuitOpen->Resize(100, circuitOpen->GetDefaultHeight());  

    // adding

    circuit1->AddFrame(circuitL1,
		       new TGLayoutHints(kLHintsLeft | 
					 kLHintsCenterY, 
					 5, 5, 2, 2)
		       );
    
    circuit1->AddFrame(circuitText1, 
		       new TGLayoutHints(kLHintsRight | 
					 kLHintsCenterY, 
					 0, 2, 2, 2)
		       );
    
    circuit2->AddFrame(circuitCancel, 
		 new TGLayoutHints(kLHintsBottom | 
				   kLHintsCenterY, 
				   5, 5, 2, 2)
		 );
  
    circuit2->AddFrame(circuitOpen, 
		 new TGLayoutHints(kLHintsBottom | 
				   kLHintsCenterY, 
				   5, 5, 2, 2)
		 );
  
    fCircuit->AddFrame(circuit1, 
		       new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 3, 0));
    
    fCircuit->AddFrame(circuit2, 
		       new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 3, 0));
    
    fCircuit->MapSubwindows();
    fCircuit->Resize();
    fCircuit->MapWindow();

    break;

  case kMTRIGGERDSET:
    /*
    cout << "Trigger with boards digits....." << endl;
    fTriggerProcessor->FeedBoardsGUI(Boards());

    nLocalTrigger = fTriggerProcessor->TriggerGUI(trigInfo,kTRUE);

    cout << "Trigger done with " << nLocalTrigger << " local decisions !" << endl;

    for (Int_t ilo = 0; ilo < nLocalTrigger; ilo++) {

      cout << "Local decision " << ilo << endl;
      cout << "Circuit = "  << trigInfo[6*ilo+0] << endl;
      cout << "LoStripX = " << trigInfo[6*ilo+1] << endl;
      cout << "LoStripY = " << trigInfo[6*ilo+2] << endl;
      cout << "LoDev = "    << trigInfo[6*ilo+3] << endl;
      cout << "LoLpt = "    << trigInfo[6*ilo+4] << endl;
      cout << "LoHpt = "    << trigInfo[6*ilo+5] << endl;
      cout                                       << endl;
    }
    */
    break;

  default:
    printf("Menu item %d selected\n", id);
    break;
    
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::CloseRunInput() const
{
  /// close the run input frame

  delete fRunInput;

}

//__________________________________________________________________________
void AliMUONTriggerGUI::CloseError() const
{
  /// close the error frame

  delete fError;

}

//__________________________________________________________________________
void AliMUONTriggerGUI::CloseControl() const
{
  /// close the event control frame

  delete fControl;

  //gSystem->Exec("cd $PWD/mtrun; aliroot -b -q mtsim.C");

}

//__________________________________________________________________________
void AliMUONTriggerGUI::CloseCircuit() const
{
  /// close the circuit frame

  delete fCircuit;

}

//__________________________________________________________________________
void AliMUONTriggerGUI::CloseWindow() 
{
  /// close the main frame and exit aplication

  delete fRunLoader;
  printf("\n"); 
  gApplication->Terminate(); 

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoRunApply()
{
  /// apply changes in the run control frame

  printf("Input 1 = %s \n",fTxtBuffer1->GetString());

  TString es = TString(fTxtBuffer2->GetString());
  fEvent = es.Atoi();

  printf("Input 2 = %s event = %d \n",fTxtBuffer2->GetString(),fEvent);

  TString error = TString("");;

  fPath->Form("%s",fTxtBuffer1->GetString());
  fFileName->Form("%s",fTxtBuffer1->GetString());
  printf("File location: %s \n",fFileName->Data());

  if (gSystem->AccessPathName(fFileName->Data()) || !fFileName->EndsWith(".root")) {

    error.Append("No galice file: ");
    error.Append(fFileName->Data());
    DoErrorGUI(error.Data());

  } else {

    TFile *ftest = new TFile(fFileName->Data(),"READ");
    AliRun *galice = (AliRun*)ftest->Get("gAlice");
    if (galice == 0) {

      ftest->Close();
      delete ftest;
      
      error.Append("No gAlice in file: ");
      error.Append(fFileName->Data());
      DoErrorGUI(error.Data());

      return;

    } else {
      ftest->Close();
      delete ftest;
    }

    if (fRunLoader) delete fRunLoader;

    fRunLoader = AliRunLoader::Open(fFileName->Data(),"MUONFolder","READ");

    if (fRunLoader == 0x0) {

      DoErrorGUI("No run loader !");

    } else {
      
      fRunLoader->LoadgAlice();
      gAlice = fRunLoader->GetAliRun();
      fEventsPerRun = gAlice->GetEventsPerRun();
      
      fLoader = fRunLoader->GetLoader("MUONLoader");
      fLoader->LoadDigits("READ");
      fRunLoader->GetEvent(fEvent);
           
      fRunInput->SendCloseMessage();
      
      if (!fBoardsInit) {
	InitBoards();
      }
      
      if (fDiMap) {
	if (fDiMap->IsOn()) {
	  fDiMap->SetLoader(fLoader);
	}
      }
      
      //fTriggerProcessor = new AliMUONTriggerElectronics(fCalibrationData);
    }

  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoRunCancel()
{
  /// cancel the changes in the run control frame

  printf("Input 1 = %s \n",fTxtBuffer1->GetString());

  TString es = TString(fTxtBuffer2->GetString());
  fEvent = es.Atoi();

  printf("Input 2 = %s event = %d \n",fTxtBuffer2->GetString(),fEvent);

  fRunInput->SendCloseMessage();

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoErrorGUI(const Char_t *wt)
{
  /// show an error message in a new frame

  fError = new TGTransientFrame(gClient->GetRoot(), fMain, 50, 50);
  fError->Connect("CloseWindow()", "AliMUONTriggerGUI", this, "CloseError()");
  fError->DontCallClose(); // to avoid double deletions.
  
  // use hierarchical cleaning
  fError->SetCleanup(kDeepCleanup);
  
  fError->SetWindowName("Error !");

  TGCompositeFrame *fW = new TGCompositeFrame(fError, 50, 50, kVerticalFrame);
  
  TGTextButton *fErrorOK = new TGTextButton(fW, "&Ok", 1);
  fErrorOK->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoErrorOK()");
  
  fW->AddFrame(fErrorOK, 
	       new TGLayoutHints(kLHintsBottom | 
				 kLHintsCenterX, 
				 5, 5, 5, 5)
	       );
  
  TGLabel *fWL = new TGLabel(fW, new TGString(wt));
  
  fW->AddFrame(fWL,
	       new TGLayoutHints(kLHintsTop | 
				 kLHintsCenterX, 
				 5, 5, 5, 5)
	       );
  
  fError->AddFrame(fW,
		   new TGLayoutHints(kLHintsLeft | 
				     kLHintsCenterY, 
				     5, 5, 5, 5)
		   );
  
  fError->MapSubwindows();
  fError->Resize();
  fError->MapWindow();
   
}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoNextEvent()
{
  /// load next event from the file

  TString error = TString("");

  if (fEvent < (fEventsPerRun-1)) {
    fEvent++;
    fRunLoader->GetEvent(fEvent);

    fEvString->Form("%d",fEvent);               
    fTxtBuffer2->RemoveText(0,5);
    fTxtBuffer2->AddText(0,fEvString->Data());
    fSkipToEventTxt->SetFocus();

  } else {
    error.Form("Only %d event(s) in the run !",fEventsPerRun);
    DoErrorGUI(error.Data());
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoPreviousEvent()
{
  /// load previous event from the input file

  TString error = TString("");

  if (fEvent > 0) {
    fEvent--;
    fRunLoader->GetEvent(fEvent);

    fEvString->Form("%d",fEvent);               
    fTxtBuffer2->RemoveText(0,5);
    fTxtBuffer2->AddText(0,fEvString->Data());
    fSkipToEventTxt->SetFocus();

  } else {
    error.Form("Already at event 0 !");
    DoErrorGUI(error.Data());
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoSkipToEvent()
{
  /// skip to event -input- from the input file

  TString error = TString("");

  TString es = TString(fTxtBuffer2->GetString());
  fEvent = es.Atoi();

  if (fEvent < 0 || fEvent > (fEventsPerRun-1)) {
    error.Form("Event number out of range !");
    DoErrorGUI(error.Data());
  } else {
    fRunLoader->GetEvent(fEvent);
    /*
    fRunLoader->LoadKinematics();

    AliStack* stack = gAlice->Stack();
    Int_t nParticles = stack->GetNtrack();
    Int_t nPrimaries = stack->GetNprimary();
    
    TParticle *part;
    Int_t nMuons = 0;
    Int_t pdgCode;
    Double_t px, py, pz, theta, phi;
    for (Int_t i = 0; i < nPrimaries; i++) {
      part = stack->Particle(i);
      if (!part) continue;
      if (TMath::Abs(part->GetPdgCode()) == 13) {
	nMuons++;
	pdgCode = part->GetPdgCode();
	px = part->Px();
	py = part->Py();
	pz = part->Pz();
	theta = part->Theta();
	phi = part->Phi();
	printf("Kine %d px %f py %f pz %f th %f ph %f \n",pdgCode,px,py,pz,theta,phi); 
      }
    }
    */    
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoErrorOK()
{
  /// close the error frame

  fError->SendCloseMessage();

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoControlClose()
{
  /// close the event control frame

  fControl->SendCloseMessage();

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoCircuitCancel()
{
  /// close the circuit frame

  fCircuit->SendCloseMessage();

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoCircuitOpen()
{
  /// open a circuit

  TString cs = TString(fTxtCircuit->GetString());
  Int_t icirc = cs.Atoi();

  AliMUONTriggerGUIboard *board;

  for (Int_t ib = 0; ib < kNBoards; ib++) {

    board = GetBoard(ib);

    if (board->GetIdCircuit() == icirc) {

      OpenBoard(ib);

      if (fDiMap) {
	if (fDiMap->IsOn()) {
	  fDiMap->SelectBoard(ib);
	}
      }
      
      fCircuit->SendCloseMessage();
      return;
    }

  }

}

//__________________________________________________________________________
AliMUONTriggerGUIboard* AliMUONTriggerGUI::GetBoard(Int_t id) const
{
  /// return board with "id"

  if (fBoards == 0) return 0;
  void * b = fBoards->UncheckedAt(id);
  if (b == 0) return 0;

  return (AliMUONTriggerGUIboard*)b;

}

//__________________________________________________________________________
void AliMUONTriggerGUI::OpenBoard(Int_t id)
{
  /// open board with "id" in a new frame
   
  AliMUONTriggerGUIboard *board = GetBoard(id);
  UShort_t status = board->GetStatus();
  board->SetOpen(kTRUE);

  AliMUONTriggerGUIbdmap *bf;
  Char_t text[200];

  bf = new AliMUONTriggerGUIbdmap(gClient->GetRoot(), fMain, 400, 200);

  if (status & kGood) {
    sprintf(text,"%s (Circuit %4d) status : working",
	    board->GetBoardName(),board->GetIdCircuit());
  }

  if (status & kWithProblems) {
    sprintf(text,"%s (Circuit %4d)  status : has problems...",
	    board->GetBoardName(),board->GetIdCircuit());
  }

  if (status & kNotWorking) {
    sprintf(text,"%s (Circuit %4d)  status : not working",
	    board->GetBoardName(),board->GetIdCircuit());
  }

  if (status & kUnknown) {
    sprintf(text,"%s (Circuit %4d)  status : unknown",
	    board->GetBoardName(),board->GetIdCircuit());
  }

  bf->SetName(text);
  bf->SetBoard(Boards(),id);
  bf->SetLoader(fLoader);
  bf->Init();

  bf->Show();

  if (fDiMap) {
    if (fDiMap->IsOn()) {
      fDiMap->SelectBoard(id);
    }
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::Init()
{
  /// initialize the main GUI frame
  
  if (!fBoardsInit) {
    InitBoards();
  }

  /*
  AliMUONTriggerGUIboard *board;
  for (Int_t ib = 0; ib < kNBoards; ib++) {
    board = GetBoard(ib);
    board->Dump();
  }
  */
}

//__________________________________________________________________________
void AliMUONTriggerGUI::InitBoards()
{
  /// create board objects and define the sensitive regions in the image map
  
  fBoardsInit = kTRUE;

  AliMUONTriggerCrateStore* crateManager = new AliMUONTriggerCrateStore();
  crateManager->ReadFromFile(0x0);

  Int_t nPixelX = 700;
  Int_t nPixelY = 676;

  Int_t nPixelBorderX = 40;  // to guess...
  Int_t nPixelBorderY = 40;  // to guess...

  Float_t boardsX = 2*257.00;  // cm
  Float_t boardsY = 2*306.61;  // cm

  FILE *fmap;

  Int_t side, col, line, nbx, detElemId = 0;
  UShort_t status = 1;
  Float_t xCenter, yCenter, zCenter, xWidth, yWidth;
  Float_t xc, yc;
  Char_t name[8], text[200];
  Int_t x, y;
  UInt_t w, h;
  Int_t xp[5];
  Int_t yp[5];
  TString mapspath = gSystem->Getenv("ALICE_ROOT");
  mapspath.Append("/MUON/data");

  TGRegion *reg;
  AliMUONTriggerGUIboard *board;

  // regions for the image map

  sprintf(text,"%s/guimapp11.txt",mapspath.Data());
  fmap = fopen(text,"r");

  for (Int_t ib = 0; ib < kNBoards; ib++) {

    fscanf(fmap,"%d   %d   %d   %d   %f   %f   %f   %f   %f   %s   \n",&side,&col,&line,&nbx,&xCenter,&yCenter,&xWidth,&yWidth,&zCenter,&name[0]);

    //printf("%d   %d   %d   %d   %f   %f   %f   %f   %f   %s   \n",side,col,line,nbx,xCenter,yCenter,xWidth,yWidth,zCenter,name);

    board = new AliMUONTriggerGUIboard(ib,name);

    status = 1;
    board->SetStatus(status);

    // calculate detElemId%100
    // side=0 left
    // side=1 right
    // ALICE SC
    if (side == 0)              detElemId = 14 - line;
    if (side == 1 && line <  5) detElemId = 13 + line;
    if (side == 1 && line >= 5) detElemId = line - 5;
    
    board->SetDetElemId(detElemId);

    Boards()->Add(board);

    xc = xCenter;
    yc = yCenter;

    x = (Int_t)(nPixelX/2 + xc * (nPixelX - 2*nPixelBorderX)/boardsX);
    y = (Int_t)(nPixelY/2 - yc * (nPixelY - 2*nPixelBorderY)/boardsY);

    if (x < 0) x = 0;
    if (y < 0) y = 0;

    w = (UInt_t)(xWidth*(nPixelX-2*nPixelBorderX)/boardsX);
    h = (UInt_t)(yWidth*(nPixelY-2*nPixelBorderY)/boardsY);
    
    xp[0] = x-w/2;
    xp[1] = x+w/2;
    xp[2] = x+w/2;
    xp[3] = x-w/2;
    xp[4] = x-w/2;
    
    yp[0] = y-h/2;
    yp[1] = y-h/2;
    yp[2] = y+h/2;
    yp[3] = y+h/2;
    yp[4] = y-h/2;

    reg = new TGRegion(5,xp,yp);
    fImageMap->AddRegion(*reg, ib);

    if (status & kGood) {
      sprintf(text,"%s working",name);
    }
    if (status & kWithProblems) {
      sprintf(text,"%s has problems...",name);
    }
    if (status & kNotWorking) {
      sprintf(text,"%s not working",name);
    }
    if (status & kUnknown) {
      sprintf(text,"%s status unknown",name);
    }
    
    //fImageMap->SetToolTipText(ib, text);
    
  }

  fclose(fmap);
  
  // MT position and dimension in board
  
  for (Int_t imt = 0; imt < kNMT; imt++) {

    sprintf(text,"%s/guimapp%2d.txt",mapspath.Data(),11+imt);

    fmap = fopen(text,"r");

    for (Int_t ib = 0; ib < kNBoards; ib++) {

      fscanf(fmap,"%d   %d   %d   %d   %f   %f   %f   %f   %f   %s   \n",&side,&col,&line,&nbx,&xCenter,&yCenter,&xWidth,&yWidth,&zCenter,&name[0]);

      board = GetBoard(ib);
      board->SetDimensions(imt,xCenter,yCenter,zCenter,xWidth,yWidth);

    }

    fclose(fmap);

  }

  // MT x-strips indices and circuit number
  Int_t ix, iy1, iy2, sIx, sIy1, cathode, icirc;
  sprintf(text,"%s/guimapix11.txt",mapspath.Data());

  fmap = fopen(text,"r");
  
  for (Int_t ib = 0; ib < kNBoards; ib++) {
    
    fscanf(fmap,"%d   %d   %d   %d   %d   %d   %d   %s   \n",&side,&col,&line,&nbx,&ix,&iy1,&iy2,&name[0]);
    
    board = GetBoard(ib);
    board->SetXSindex(ix,iy1,iy2);

    // set the circuit number
    detElemId = board->GetDetElemId();
    detElemId += 100 * 11;
    cathode = 0;
    sIx  = board->GetXSix();
    sIy1 = board->GetXSiy1();
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));
    
    AliMpPad pad = seg->PadByIndices(AliMpIntPair(sIx,sIy1),kTRUE);
    AliMpIntPair location = pad.GetLocation(0);

    Int_t nboard = location.GetFirst();
    AliMUONLocalTriggerBoard* b = crateManager->LocalBoard(nboard);
    icirc = b->GetNumber();
    board->SetBoardName((Char_t*)b->GetName());
    board->SetIdCircuit(icirc);

    TString crateName = b->GetCrate();
    
    sprintf(text,"%s (crate %s circuit %3d, number %3d)",board->GetBoardName(),crateName.Data(),icirc,ib);
    fImageMap->SetToolTipText(ib, text);
        
  }
  
  fclose(fmap);

  // MT y-strips indices
  Int_t ix1, ix2, iy;

  sprintf(text,"%s/guimapiy11.txt",mapspath.Data());
  
  fmap = fopen(text,"r");
  
  for (Int_t ib = 0; ib < kNBoards; ib++) {
    
    fscanf(fmap,"%d   %d   %d   %d   %d   %d   %d   %s   \n",&side,&col,&line,&nbx,&ix1,&ix2,&iy,&name[0]);
    
    board = GetBoard(ib);
    board->SetYSindex(ix1,ix2,iy);
    
  }
  
  fclose(fmap);

  // Extended y-strips over neighbouring boards

  sprintf(text,"%s/guimapp11.txt",mapspath.Data());
  
  fmap = fopen(text,"r");
  
  for (Int_t ib = 0; ib < kNBoards; ib++) {
    
    fscanf(fmap,"%d   %d   %d   %d   %f   %f   %f   %f   %f   %s   \n",&side,&col,&line,&nbx,&xCenter,&yCenter,&xWidth,&yWidth,&zCenter,&name[0]);
    
    board = GetBoard(ib);

    board->SetPosition(nbx);
    
    if ((col == 2 || col == 3) && (line >= 4 && line <= 6)) {
      board->SetYOver(4);
    } else
      
    if (col == 1 && (line == 4 || line == 6)) {
      board->SetYOver(3);
    } else
	
    if (col == 7 || line == 1 || line == 9) {
      board->SetYOver(1);
    } else
	  
    {
      board->SetYOver(2);
    }
    
  }
  
  fclose(fmap);

  // Set coordinates of strips boxes

  for (Int_t ib = 0; ib < kNBoards; ib++) {
    
    board = GetBoard(ib);
    SetStripBoxes(board);

  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::SetStripBoxes(AliMUONTriggerGUIboard *board) 
{
  /// set coordinates of strip boxes

  gAlice = fRunLoader->GetAliRun();
  AliMUON *pMUON = (AliMUON*)gAlice->GetModule("MUON");
  const AliMUONGeometryTransformer* kGeomTransformer = pMUON->GetGeometryTransformer();

  const AliMpVSegmentation* seg;
  AliMpPad pad;
  AliMpDEIterator it;

  Int_t chamber, detElemId, maxX, maxY, ic;
  Float_t xpmin, xpmax, ypmin, ypmax;
  Float_t xg1, xg2, yg1, yg2, zg1;
  Float_t xlocal1, xlocal2, ylocal1, ylocal2;
  Float_t xCenter, yCenter, xWidth, yWidth;

  for (Int_t i = 0; i < kNMT; i++) {

    chamber = 11+i;

    xCenter = board->GetXCenter(i);
    yCenter = board->GetYCenter(i);
    xWidth  = board->GetXWidth(i);
    yWidth  = board->GetYWidth(i);

    for ( it.First(chamber-1); ! it.IsDone(); it.Next() ) {
    
      detElemId = it.CurrentDEId();
    
      if (detElemId%100 != board->GetDetElemId()) continue;

      /*---------- y-pads ic = 2 ----------*/
      ic = 2;
      
      seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(ic-1));

      maxX = seg->MaxPadIndexX();
      maxY = seg->MaxPadIndexY();
      
      for (Int_t ix = 0; ix <= maxX; ix++) {
	for (Int_t iy = 0; iy <= maxY; iy++) {
	  
	  pad = seg->PadByIndices(AliMpIntPair(ix,iy),kFALSE);
	  
	  if (!pad.IsValid()) continue;
	  
	  // get the pad position and dimensions
	  xlocal1 = pad.Position().X();
	  ylocal1 = pad.Position().Y();
	  xlocal2 = pad.Dimensions().X();
	  ylocal2 = pad.Dimensions().Y();
	  
	  kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
	  // (no transformation for pad dimensions)
	  xg2 = xlocal2;
	  yg2 = ylocal2;
	  
	  // transform in the monitor coordinate system
	  //xpmin = -(xg1+xg2);
	  //xpmax = -(xg1-xg2);
	  //ypmin = -(yg2-yg1);
	  //ypmax = +(yg2+yg1);
	  // ALICE SC
	  xpmin = +(xg1-xg2);
	  xpmax = +(xg1+xg2);
	  ypmin = -(yg2-yg1);
	  ypmax = +(yg2+yg1);
	  
	  xpmin -= xCenter;
	  xpmax -= xCenter;
	  ypmin -= yCenter;
	  ypmax -= yCenter;
	  
	  Int_t iX1, iX2, iY, ixDig;
	  iX1 = board->GetYSix1();
	  iX2 = board->GetYSix2();
	  iY  = board->GetYSiy();
	  if (ix >= iX1 && ix <= iX2 && iy == iY) {
	    
	    ypmin = -0.5*yWidth;
	    ypmax = +0.5*yWidth;

	    ixDig = ix - iX1;

	    board->SetYDigBox(i,ixDig,(Double_t)xpmin,(Double_t)ypmin,
			              (Double_t)xpmax,(Double_t)ypmax);

	  }

	}  // end maxY
      }  // end maxX

      /*---------- x-pads ic = 1 ----------*/
      ic = 1;
      
      seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(ic-1));

      maxX = seg->MaxPadIndexX();
      maxY = seg->MaxPadIndexY();
            
      for (Int_t ix = 0; ix <= maxX; ix++) {
	for (Int_t iy = 0; iy <= maxY; iy++) {
	  
	  pad = seg->PadByIndices(AliMpIntPair(ix,iy),kFALSE);
	  
	  if (!pad.IsValid()) continue;
	  
	  // get the pad position and dimensions
	  xlocal1 = pad.Position().X();
	  ylocal1 = pad.Position().Y();
	  xlocal2 = pad.Dimensions().X();
	  ylocal2 = pad.Dimensions().Y();
	  
	  kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
	  // (no transformation for pad dimensions)
	  xg2 = xlocal2;
	  yg2 = ylocal2;
	  
	  // transform in the monitor coordinate system
	  //xpmin = -(xg1+xg2);
	  //xpmax = -(xg1-xg2);
	  //ypmin = -(yg2-yg1);
	  //ypmax = +(yg2+yg1);
	  // ALICE SC
	  xpmin = +(xg1+xg2);
	  xpmax = +(xg1-xg2);
	  ypmin = -(yg2-yg1);
	  ypmax = +(yg2+yg1);
	  
	  xpmin -= xCenter;
	  xpmax -= xCenter;
	  ypmin -= yCenter;
	  ypmax -= yCenter;

	  Int_t iX, iY1, iY2, iyDig;
	  iX  = board->GetXSix();
	  iY1 = board->GetXSiy1();
	  iY2 = board->GetXSiy2();
	  if (ix == iX && iy >= iY1 && iy <= iY2) {

	    iyDig = iy - iY1;

	    board->SetXDigBox(i,iyDig,(Double_t)xpmin,(Double_t)ypmin,
			              (Double_t)xpmax,(Double_t)ypmax);

	  }

	}  // end maxY
      }  // end maxX

    }  // end detElemId

  }  // end kNMT loop
	
}

