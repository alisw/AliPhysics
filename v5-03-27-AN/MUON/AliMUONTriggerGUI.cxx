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

#include "AliMpDDLStore.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpTriggerSegmentation.h"
#include "AliMpEncodePair.h"
#include "AliMpCDB.h"
#include "AliMpDEManager.h"
#include "AliMpDEIterator.h"

#include "AliMUONGeometryTransformer.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONDigitStoreV1.h"
#include "AliMUONDigitStoreV2R.h"
#include "AliMUONTriggerStoreV1.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRawWriter.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONGlobalTrigger.h"

#include "AliRun.h"
#include "AliDAQ.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"
#include "AliRawDataHeaderSim.h"
#include "AliRawReader.h"

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

#include <cstdio>

/// \cond CLASSIMP
ClassImp(AliMUONTriggerGUI)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerGUI::AliMUONTriggerGUI(Int_t runNumber)
  : TObject(),
    fMain(0),
    fImageMap(0),
    fTxtBuffer1(0),
    fTxtBuffer2(0),
    fTxtCircuit(0),
    fTxtFETRegOn(0),
    fTxtFETRegOff(0),
    fRunInput(0),
    fError(0),
    fControl(0),
    fCircuit(0),
    fFETRegOn(0),
    fFETRegOff(0),
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
    fCrateManager(0),
    fMCDataInterface(0),
    fBoardsInit(0),
    fControlOn(kFALSE),
    fDiMap(0),
    fTriggerProcessor(0),
    fBoards(0),
    fDigitStore(0),
    fTriggerStore(0),
    fTStoreOn(kFALSE),
    fRUNRAW(kFALSE),
    fRawReader(0),
    fCurrentRawEvent(-1),
    fRawDigitStore(0),
    fRawTriggerStore(0)
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

  fDigitStore = new AliMUONDigitStoreV2R;
  fDigitStore->Create();

  fTriggerStore = new AliMUONTriggerStoreV1;
  fTriggerStore->Create();

  fRawDigitStore = new AliMUONDigitStoreV1;
  fRawDigitStore->Create();

  fRawTriggerStore = new AliMUONTriggerStoreV1;
  fRawTriggerStore->Create();

  fCDBManager = AliCDBManager::Instance();
  fCDBManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  fCDBManager->SetRun(runNumber);
  AliMpCDB::LoadDDLStore();
  fCalibrationData = new AliMUONCalibrationData(runNumber);
  
  fTriggerProcessor = new AliMUONTriggerElectronics(fCalibrationData);

  // Main frame

  fMain = new TGMainFrame(gClient->GetRoot(), 750, 420);
  fMain->Connect("CloseWindow()", "AliMUONTriggerGUI", this, "CloseWindow()");

  // Menu bar
  
  TGMenuBar *menuBar = new TGMenuBar(fMain);
  
  // File menu

  TGPopupMenu *menuFile = new TGPopupMenu(gClient->GetRoot());
  //menuFile->AddLabel("");

  menuFile->AddEntry("Run input", kMFILERUN);
  menuFile->AddEntry("Control",   kMFILECNTRL);
  menuFile->AddEntry("Exit",      kMFILEEXIT);

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

  TGPopupMenu *menuTriggerD = new TGPopupMenu(gClient->GetRoot());
  TGPopupMenu *menuTriggerT = new TGPopupMenu(gClient->GetRoot());
  TGPopupMenu *menuTriggerF = new TGPopupMenu(gClient->GetRoot());

  menuTrigger->AddPopup("Digit store",       menuTriggerD);
  menuTrigger->AddSeparator();
  menuTrigger->AddPopup("Trigger store",menuTriggerT);
  menuTrigger->AddSeparator();
  menuTrigger->AddPopup("Front End Test",    menuTriggerF);
  menuTrigger->AddSeparator();

  menuTriggerD->AddEntry("Create", kMDSTORE);
  menuTriggerD->AddEntry("Print",  kMDSTOREP);
  menuTriggerD->AddEntry("Clear",  kMDSTORECL);
  menuTriggerT->AddEntry("Create", kMTSTORE);
  menuTriggerT->AddEntry("Print",  kMTSTOREP);
  menuTriggerT->AddEntry("Clear",  kMTSTORECL);
  menuTriggerF->AddEntry("On",     kMFETON);
  menuTriggerF->AddEntry("Off",    kMFETOFF);
  menuTriggerF->AddEntry("Reg On", kMFETREGON);
  menuTriggerF->AddEntry("Reg Off",kMFETREGOFF);
  menuTrigger->AddEntry("Write raw data", kMTRAWDATA);

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

  menuBar->AddPopup("TriggerDSET", menuTrigger, 
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

  InitBoards();
  
}

//__________________________________________________________________________
void AliMUONTriggerGUI::HandleMenu(Int_t id)
{
  /// handles entry numbers in the available menus (EMenuIdentifiers)

  TGCompositeFrame *runInput1, *runInput2, *runInput3;
  TGCompositeFrame *control1, *control2, *circuit1, *circuit2;
  TGCompositeFrame *fetregon1, *fetregon2;
  TGCompositeFrame *fetregoff1, *fetregoff2;
  TGLabel *runL1, *runL2, *circuitL1, *fetregonL1, *fetregoffL1;
  TGTextEntry *runText1, *runText2, *circuitText1;
  TGTextEntry *fetregonText1, *fetregoffText1;
  TGTextButton *runApply1, *runApply2, *runCancel; 
  TGTextButton *controlClose, *nextEvent, *previousEvent, *skipToEvent;
  TGTextButton *circuitCancel, *circuitOpen;
  TGTextButton *fetregonCancel, *fetregoffCancel;
  TGTextButton *fetregonRun, *fetregoffRun;

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

    runL1 = new TGLabel(runInput1, new TGString("Full file path:"));
    runL2 = new TGLabel(runInput2, new TGString("Event number:"));

    // galice text entry

    runText1 = new TGTextEntry(runInput1, fTxtBuffer1 = new TGTextBuffer(100));

    runText1->SetToolTipText("Enter the path to galice.root or the raw data file (root)");
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

    runApply1 = new TGTextButton(runInput3, "Apply (galice)", 1);
    runApply1->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoRunGalApply()");
    runApply1->SetToolTipText("Apply changes (galice input)");

    runApply2 = new TGTextButton(runInput3, "Apply (raw)", 1);
    runApply2->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoRunRawApply()");
    runApply2->SetToolTipText("Apply changes (raw data input)");

    runCancel = new TGTextButton(runInput3, "Cancel", 2);
    runCancel->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoRunCancel()");    
    runCancel->SetToolTipText("Cancel changes");

    // add buttons

    runInput3->AddFrame(runApply1, 
			 new TGLayoutHints(kLHintsTop | 
					   kLHintsLeft, 
					   3, 3, 2, 2)
			);
    
    runInput3->AddFrame(runApply2, 
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

    fControlOn = kTRUE;
   
    break;

  case kMMAPDIGITS:
    
    if (fDiMap == 0) {
      fDiMap = new AliMUONTriggerGUIdimap(fBoards,gClient->GetRoot(), fMain, 400, 200);
      fDiMap->SetLoader(fLoader);
      if (fRUNRAW) {
	fDiMap->SetRawDigitStore(fRawDigitStore);
      } else {
	fDiMap->SetMCDataInterface(fMCDataInterface);
      }
      fDiMap->DrawAllMaps();
    } else if (!fDiMap->IsOn()) {
      fDiMap = new AliMUONTriggerGUIdimap(fBoards,gClient->GetRoot(), fMain, 400, 200);
      fDiMap->SetLoader(fLoader);
      if (fRUNRAW) {
	fDiMap->SetRawDigitStore(fRawDigitStore);
      } else {
	fDiMap->SetMCDataInterface(fMCDataInterface);
      }
      fDiMap->DrawAllMaps();
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

  case kMDSTORE:
    CreateDigitStore();
    break;

  case kMDSTOREP:
    PrintDigitStore();
    break;

  case kMDSTORECL:
    ClearDigitStore();
    break;

  case kMTSTORE:
    CreateTriggerStore();
    break;

  case kMTSTOREP:
    PrintTriggerStore();
    break;

  case kMTSTORECL:
    ClearTriggerStore();
    break;

  case kMTRAWDATA:
    WriteTriggerRawData();
    break;

  case kMFETON:
    FET(1);
    break;

  case kMFETOFF:
    FET(0);
    break;

  case kMFETREGON:

    fFETRegOn = new TGTransientFrame(gClient->GetRoot(), fMain, 50, 50);
    fFETRegOn->Connect("CloseWindow()", "AliMUONTriggerGUI", this, "CloseFETRegOn()");
    fFETRegOn->DontCallClose(); // to avoid double deletions.
  
    // use hierarchical cleaning
    fFETRegOn->SetCleanup(kDeepCleanup);
  
    fFETRegOn->SetWindowName("FET ON regional crate");

    // sub-frames

    fetregon1 = new TGCompositeFrame(fFETRegOn, 400, 200, kHorizontalFrame);
    fetregon2 = new TGCompositeFrame(fFETRegOn, 400, 200, kHorizontalFrame);

    // labels

    fetregonL1 = new TGLabel(fetregon1, new TGString("Regional crate name:"));
    
    // text entry
    
    fetregonText1 = new TGTextEntry(fetregon1, fTxtFETRegOn = new TGTextBuffer(10));
    // buttons

    fetregonCancel = new TGTextButton(fetregon2, "Cancel", 1);
    fetregonCancel->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoFETRegOnCancel()");
    //fetregonCancel->Resize(100, fetregonCancel->GetDefaultHeight());  

    fetregonRun = new TGTextButton(fetregon2, "Run FET", 2);
    fetregonRun->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoFETRegOnRun()");
    //fetregonRun->Resize(100, fetregonRun->GetDefaultHeight());  

    // adding

    fetregon1->AddFrame(fetregonL1,
		       new TGLayoutHints(kLHintsLeft | 
					 kLHintsCenterY, 
					 5, 5, 2, 2)
		       );
    
    fetregon1->AddFrame(fetregonText1, 
		       new TGLayoutHints(kLHintsRight | 
					 kLHintsCenterY, 
					 0, 2, 2, 2)
		       );
    
    fetregon2->AddFrame(fetregonCancel, 
		 new TGLayoutHints(kLHintsBottom | 
				   kLHintsCenterY, 
				   5, 5, 2, 2)
		 );
  
    fetregon2->AddFrame(fetregonRun, 
		 new TGLayoutHints(kLHintsBottom | 
				   kLHintsCenterY, 
				   5, 5, 2, 2)
		 );
  
    fFETRegOn->AddFrame(fetregon1, 
		       new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 3, 0));
    
    fFETRegOn->AddFrame(fetregon2, 
		       new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 3, 0));
    
    fFETRegOn->MapSubwindows();
    fFETRegOn->Resize();
    fFETRegOn->MapWindow();

    break;

  case kMFETREGOFF:

    fFETRegOff = new TGTransientFrame(gClient->GetRoot(), fMain, 50, 50);
    fFETRegOff->Connect("CloseWindow()", "AliMUONTriggerGUI", this, "CloseFETRegOff()");
    fFETRegOff->DontCallClose(); // to avoid double deletions.
  
    // use hierarchical cleaning
    fFETRegOff->SetCleanup(kDeepCleanup);
  
    fFETRegOff->SetWindowName("FET OFF regional crate");

    // sub-frames

    fetregoff1 = new TGCompositeFrame(fFETRegOff, 400, 200, kHorizontalFrame);
    fetregoff2 = new TGCompositeFrame(fFETRegOff, 400, 200, kHorizontalFrame);

    // labels

    fetregoffL1 = new TGLabel(fetregoff1, new TGString("Regional crate name:"));
    
    // text entry
    
    fetregoffText1 = new TGTextEntry(fetregoff1, fTxtFETRegOff = new TGTextBuffer(10));
    // buttons

    fetregoffCancel = new TGTextButton(fetregoff2, "Cancel", 1);
    fetregoffCancel->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoFETRegOffCancel()");
    //fetregoffCancel->Resize(100, fetregoffCancel->GetDefaultHeight());  

    fetregoffRun = new TGTextButton(fetregoff2, "Run FET", 2);
    fetregoffRun->Connect("Clicked()", "AliMUONTriggerGUI", this, "DoFETRegOffRun()");
    //fetregoffRun->Resize(100, fetregoffRun->GetDefaultHeight());  

    // adding

    fetregoff1->AddFrame(fetregoffL1,
		       new TGLayoutHints(kLHintsLeft | 
					 kLHintsCenterY, 
					 5, 5, 2, 2)
		       );
    
    fetregoff1->AddFrame(fetregoffText1, 
		       new TGLayoutHints(kLHintsRight | 
					 kLHintsCenterY, 
					 0, 2, 2, 2)
		       );
    
    fetregoff2->AddFrame(fetregoffCancel, 
		 new TGLayoutHints(kLHintsBottom | 
				   kLHintsCenterY, 
				   5, 5, 2, 2)
		 );
  
    fetregoff2->AddFrame(fetregoffRun, 
		 new TGLayoutHints(kLHintsBottom | 
				   kLHintsCenterY, 
				   5, 5, 2, 2)
		 );
  
    fFETRegOff->AddFrame(fetregoff1, 
		       new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 3, 0));
    
    fFETRegOff->AddFrame(fetregoff2, 
		       new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 3, 0));
    
    fFETRegOff->MapSubwindows();
    fFETRegOff->Resize();
    fFETRegOff->MapWindow();

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
void AliMUONTriggerGUI::CloseFETRegOn() const
{
  /// close the FET regional on frame

  delete fFETRegOn;

}

//__________________________________________________________________________
void AliMUONTriggerGUI::CloseFETRegOff() const
{
  /// close the FET regional off frame

  delete fFETRegOff;

}

//__________________________________________________________________________
void AliMUONTriggerGUI::CloseWindow() 
{
  /// close the main frame and exit aplication

  delete fRunLoader;
  delete fMCDataInterface;
  printf("\n"); 
  gApplication->Terminate(); 

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoRunGalApply()
{
  /// apply changes in the run control frame for galice input

  if (fRUNRAW) {
    printf("This is a run with raw data input.\n");
    return;
  }

  TString es = TString(fTxtBuffer2->GetString());
  fEvent = es.Atoi();

  printf("Galice input = %s event = %d \n",fTxtBuffer1->GetString(),fEvent);

  TString error = TString("");;

  fPath->Form("%s",fTxtBuffer1->GetString());
  fFileName->Form("%s",fTxtBuffer1->GetString());

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

    if (fRunLoader) {
      delete fRunLoader;
      delete fMCDataInterface;
    }

    fRunLoader = AliRunLoader::Open(fFileName->Data(),"MUONFolder","READ");

    if (fRunLoader == 0x0) {

      DoErrorGUI("No run loader !");

    } else {
      
      fRunLoader->LoadgAlice();
      gAlice = fRunLoader->GetAliRun();
      fEventsPerRun = AliRunLoader::Instance()->GetNumberOfEvents();
      
      fLoader = fRunLoader->GetLoader("MUONLoader");
      fRunLoader->GetEvent(fEvent);
           
      fMCDataInterface = new AliMUONMCDataInterface(fFileName->Data());

      fRunInput->SendCloseMessage();
      
      if (!fBoardsInit) {
	InitBoards();
      }
      
      if (fDiMap) {
	if (fDiMap->IsOn()) {
	  fDiMap->SetLoader(fLoader);
	  fDiMap->SetMCDataInterface(fMCDataInterface);
	}
      }
      
    }

  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoRunRawApply()
{
  /// apply changes in the run control frame for raw date input

  if (fControlOn) DoControlClose();

  if (fRunLoader != 0) {
    printf("This is a run with galice input.\n");
    return;
  }

  TString es = TString(fTxtBuffer2->GetString());
  fEvent = es.Atoi();

  printf("Raw data input = %s event = %d \n",fTxtBuffer1->GetString(),fEvent);

  TString error = TString("");;

  fPath->Form("%s",fTxtBuffer1->GetString());
  fFileName->Form("%s",fTxtBuffer1->GetString());

  if (gSystem->AccessPathName(fFileName->Data()) || !fFileName->EndsWith(".root")) {

    error.Append("No raw data file: ");
    error.Append(fFileName->Data());
    DoErrorGUI(error.Data());

  } else {

    if (fRawReader == 0) {
      fRawReader = AliRawReader::Create(fFileName->Data());
    } else {
      delete fRawReader;
      fRawReader = AliRawReader::Create(fFileName->Data());
    }
    if (fRawReader == 0) {
      error.Append("Not a valid raw data file: ");
      error.Append(fFileName->Data());
      DoErrorGUI(error.Data());
      return;
    }

    fRawReader->GotoEvent(fEvent);
    fCurrentRawEvent = fEvent;

    AliMUONDigitMaker digitMaker;
    digitMaker.SetMakeTriggerDigits(kTRUE);
    digitMaker.Raw2Digits(fRawReader,fRawDigitStore,fRawTriggerStore);

    fRunInput->SendCloseMessage();
      
    if (!fBoardsInit) {
      InitBoards();
    }
    
    if (fDiMap) {
      if (fDiMap->IsOn()) {
	fDiMap->SetRawDigitStore(fRawDigitStore);
      }
    }

    fRUNRAW = kTRUE;

  }
    
}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoRunCancel()
{
  /// cancel the changes in the run control frame

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

  if (!fRUNRAW) {
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
  } else {
    if (fRawReader->NextEvent()) {
      fCurrentRawEvent++;

      AliMUONDigitMaker digitMaker;
      digitMaker.SetMakeTriggerDigits(kTRUE);
      digitMaker.Raw2Digits(fRawReader,fRawDigitStore,fRawTriggerStore);
      
      fEvString->Form("%d",fCurrentRawEvent);               
      fTxtBuffer2->RemoveText(0,5);
      fTxtBuffer2->AddText(0,fEvString->Data());
      fSkipToEventTxt->SetFocus();
      
    }

  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoPreviousEvent()
{
  /// load previous event from the input file

  TString error = TString("");

  if (!fRUNRAW) {
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
  } else {
    if (fCurrentRawEvent > 0) {
      fCurrentRawEvent--;
      fRawReader->GotoEvent(fCurrentRawEvent);

      AliMUONDigitMaker digitMaker;
      digitMaker.SetMakeTriggerDigits(kTRUE);
      digitMaker.Raw2Digits(fRawReader,fRawDigitStore,fRawTriggerStore);
      
      fEvString->Form("%d",fCurrentRawEvent);               
      fTxtBuffer2->RemoveText(0,5);
      fTxtBuffer2->AddText(0,fEvString->Data());
      fSkipToEventTxt->SetFocus();
    
    }

  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoSkipToEvent()
{
  /// skip to event -input- from the input file

  TString error = TString("");

  TString es = TString(fTxtBuffer2->GetString());
  fEvent = es.Atoi();

  if (!fRUNRAW) {
    if (fEvent < 0 || fEvent > (fEventsPerRun-1)) {
      error.Form("Event number out of range !");
      DoErrorGUI(error.Data());
    } else {
      fRunLoader->GetEvent(fEvent);
    }
  } else {
    if (fRawReader->GotoEvent(fEvent)) {
      fCurrentRawEvent = fEvent;

      AliMUONDigitMaker digitMaker;
      digitMaker.SetMakeTriggerDigits(kTRUE);
      digitMaker.Raw2Digits(fRawReader,fRawDigitStore,fRawTriggerStore);
      
      fEvString->Form("%d",fCurrentRawEvent);               
      fTxtBuffer2->RemoveText(0,5);
      fTxtBuffer2->AddText(0,fEvString->Data());
      fSkipToEventTxt->SetFocus();
    }

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
  fControlOn = kFALSE;

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoCircuitCancel()
{
  /// close the circuit frame

  fCircuit->SendCloseMessage();

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoFETRegOnCancel()
{
  /// close the FET regional on window

  fFETRegOn->SendCloseMessage();

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoFETRegOffCancel()
{
  /// close the FET regional off window

  fFETRegOff->SendCloseMessage();

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoFETRegRun(Int_t onoff)
{
  /// FET ON/OFF for the regional crate

  TString crateName;
  if (onoff == 1) 
    crateName = TString(fTxtFETRegOn->GetString());
  if (onoff == 0) 
    crateName = TString(fTxtFETRegOff->GetString());

  AliMUONTriggerGUIboard *board;
  Int_t amp = onoff;

  for (Int_t ib = 0; ib < kNBoards; ib++) {
    board = (AliMUONTriggerGUIboard*)fBoards->At(ib);
    if (strcmp(board->GetCrateName(),crateName.Data()) != 0) continue;
    FETboard(ib,amp);
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoFETRegOnRun()
{
  /// FET ON for the regional crate

  DoFETRegRun(1);
  fFETRegOn->SendCloseMessage();

}

//__________________________________________________________________________
void AliMUONTriggerGUI::DoFETRegOffRun()
{
  /// FET ON for the regional crate

  DoFETRegRun(0);
  fFETRegOff->SendCloseMessage();

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
    snprintf(text,200,"%s (Circuit %4d) status : working",
	    board->GetBoardName(),board->GetIdCircuit());
  }

  if (status & kWithProblems) {
    snprintf(text,200,"%s (Circuit %4d)  status : has problems...",
	    board->GetBoardName(),board->GetIdCircuit());
  }

  if (status & kNotWorking) {
    snprintf(text,200,"%s (Circuit %4d)  status : not working",
	    board->GetBoardName(),board->GetIdCircuit());
  }

  if (status & kUnknown) {
    snprintf(text,200,"%s (Circuit %4d)  status : unknown",
	    board->GetBoardName(),board->GetIdCircuit());
  }

  bf->SetName(text);
  bf->SetBoard(Boards(),id);
  bf->SetLoader(fLoader);
  bf->SetCrateManager(fCrateManager);
  if (fRUNRAW) {
    bf->SetRawDigitStore(fRawDigitStore);
    bf->SetRawTriggerStore(fRawTriggerStore);
  } else {
    bf->SetMCDataInterface(fMCDataInterface);
  }
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

  // used by bdmap
  if (fCrateManager == 0x0) {
    fCrateManager = new AliMUONTriggerCrateStore();
    fCrateManager->ReadFromFile(fCalibrationData);
  }

  // create boards geometry from the mapping
  AliMUONTriggerGUIboard *board;
  for (Int_t ib = 0; ib < kNBoards; ib++) {
    Boards()->Add(new AliMUONTriggerGUIboard());
  }
  
  // circuit number to board number in array
  Int_t cIdtobId[235];
  for (Int_t i = 0; i < 235; i++) cIdtobId[i] = -1;

  AliMpDEIterator it;
  Int_t boardId = -1;
  Int_t manuIdPrev, ich, idet, boardIdTmp = -1;
  for (Int_t chamber = 0; chamber < kNMT; chamber++) {
    for ( it.First(); ! it.IsDone(); it.Next() ) {
      
      if (it.CurrentDEId()/100 < 11) continue;
      
      ich = it.CurrentDEId()/100 - 11;
      if (ich != chamber) continue;
      idet = it.CurrentDEId()%100;
      
      const AliMpVSegmentation* seg0 = AliMpSegmentation::Instance()->GetMpSegmentation(it.CurrentDEId(), AliMp::kCath0);
      const AliMpVSegmentation* seg1 = AliMpSegmentation::Instance()->GetMpSegmentation(it.CurrentDEId(), AliMp::kCath1);

      // x-strips
      manuIdPrev = 0;
      for (Int_t ix = 0; ix <= seg0->MaxPadIndexX(); ix++) {
	for (Int_t iy = 0; iy <= seg0->MaxPadIndexY(); iy++) {
	  AliMpPad pad = seg0->PadByIndices(ix,iy,kFALSE);
	  if (pad.IsValid()) {
	    Int_t manuId = pad.GetLocalBoardId(0);
	    if (manuId != manuIdPrev) {
	      AliMpLocalBoard *mpboard = AliMpDDLStore::Instance()->GetLocalBoard(manuId);
	      manuIdPrev = manuId;
	      if (ich == 0) {
		boardId++;
	      } else {
		boardId = cIdtobId[manuId];
	      }
	      board = GetBoard(boardId);
	      if (board->GetNumber() == -1)
		board->SetNumber(boardId);
	      if (board->GetDetElemId() == -1) 
		board->SetDetElemId(idet);
	      if (!strcmp(board->GetBoardName(),""))
		board->SetBoardName(mpboard->GetName());
	      cIdtobId[manuId] = boardId;
	    }
	    GetBoard(boardId)->AddPadX(pad,ich);
	  }
	}
      }  // end plane 0 (x-strips)
         
      // y-strips
      manuIdPrev = 0;
      for (Int_t ix = 0; ix <= seg1->MaxPadIndexX(); ix++) {
	for (Int_t iy = 0; iy <= seg1->MaxPadIndexY(); iy++) {
	  AliMpPad pad = seg1->PadByIndices(ix,iy,kFALSE);
	  if (pad.IsValid()) {
	    Int_t nloc = pad.GetNofLocations();
	    for (Int_t iloc = 0; iloc < nloc; iloc++) {
	      Int_t manuId = pad.GetLocalBoardId(iloc);
	      if (manuId != manuIdPrev) {
		manuIdPrev = manuId;
		boardIdTmp = cIdtobId[manuId];
	      }
	      GetBoard(boardIdTmp)->AddPadY(pad,ich);
	    }
	  }
	}
      }  // end plane 1 (y-strips)
      
    } // end det elem loop
  } // end chamber loop
  
  for (Int_t ib = 0; ib < kNBoards; ib++) {
    board = GetBoard(ib);
    board->MakeGeometry();
    AliMUONLocalTriggerBoard* b = fCrateManager->LocalBoard(board->GetIdCircuit());
    board->SetCrateName((b->GetCrate()).Data());
  }
  
  // create the sensitive map

  Int_t nPixelX = 700;
  Int_t nPixelY = 676;

  Int_t nPixelBorderX = 40;  // to guess...
  Int_t nPixelBorderY = 40;  // to guess...

  // boards limits
  Float_t boardsX = 2*257.00;  // cm
  Float_t boardsY = 2*306.61;  // cm

  UShort_t status = 1;
  Float_t xc, yc, xw, yw;
  Char_t text[256];
  Int_t x, y;
  UInt_t w, h;
  Int_t xp[5];
  Int_t yp[5];

  TGRegion *reg;

  // regions for the image map (from MT11)

  for (Int_t ib = 0; ib < kNBoards; ib++) {

    board = GetBoard(ib);

    status = 1;
    board->SetStatus(status);

    xc = board->GetXCenter(0);
    yc = board->GetYCenter(0);
    xw = board->GetXWidth(0);
    yw = board->GetYWidth(0);

    x = (Int_t)(nPixelX/2 + xc * (nPixelX - 2*nPixelBorderX)/boardsX);
    y = (Int_t)(nPixelY/2 - yc * (nPixelY - 2*nPixelBorderY)/boardsY);

    if (x < 0) x = 0;
    if (y < 0) y = 0;

    w = (UInt_t)(xw*(nPixelX-2*nPixelBorderX)/boardsX);
    h = (UInt_t)(yw*(nPixelY-2*nPixelBorderY)/boardsY);
    
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

    snprintf(text,256,"%s (crate %s circuit %3d, number %3d)",board->GetBoardName(),board->GetCrateName(),board->GetIdCircuit(),board->GetNumber());
    fImageMap->SetToolTipText(ib, text);

    // Set coordinates of strips boxes
    
    SetStripBoxes(board);
    
    //board->PrintBoard();
    
  }
  
}

//__________________________________________________________________________
void AliMUONTriggerGUI::SetStripBoxes(AliMUONTriggerGUIboard *board) 
{
  /// set coordinates of strip boxes

  AliMUONGeometryTransformer transformer;
  transformer.LoadGeometryData("transform.dat");

  const AliMpVSegmentation* seg;
  AliMpPad pad;
  AliMpDEIterator it;

  Int_t chamber, detElemId, maxX, maxY;
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

      /*---------- y-pads cath = 1 ----------*/
      
      seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::kCath1);

      maxX = seg->MaxPadIndexX();
      maxY = seg->MaxPadIndexY();
      
      for (Int_t ix = 0; ix <= maxX; ix++) {
	for (Int_t iy = 0; iy <= maxY; iy++) {
	  
	  pad = seg->PadByIndices(ix,iy,kFALSE);
	  
	  if (!pad.IsValid()) continue;
	  
	  // get the pad position and dimensions
	  xlocal1 = pad.GetPositionX();
	  ylocal1 = pad.GetPositionY();
	  xlocal2 = pad.GetDimensionX();
	  ylocal2 = pad.GetDimensionY();
	  
	  transformer.Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
	  // (no transformation for pad dimensions)
	  xg2 = xlocal2;
	  yg2 = ylocal2;
	  
	  // transform in the monitor coordinate system
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

      /*---------- x-pads cath = 0 ----------*/
      
      seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::kCath0);

      maxX = seg->MaxPadIndexX();
      maxY = seg->MaxPadIndexY();
            
      for (Int_t ix = 0; ix <= maxX; ix++) {
	for (Int_t iy = 0; iy <= maxY; iy++) {
	  
	  pad = seg->PadByIndices(ix,iy,kFALSE);
	  
	  if (!pad.IsValid()) continue;
	  
	  // get the pad position and dimensions
	  xlocal1 = pad.GetPositionX();
	  ylocal1 = pad.GetPositionY();
	  xlocal2 = pad.GetDimensionX();
	  ylocal2 = pad.GetDimensionY();
	  
	  transformer.Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
	  // (no transformation for pad dimensions)
	  xg2 = xlocal2;
	  yg2 = ylocal2;
	  
	  // transform in the monitor coordinate system
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

//__________________________________________________________________________
void AliMUONTriggerGUI::CreateDigitStore() 
{
  /// create memory resident digits store with set strips

  Int_t nstripX, nstripY, detElemId, charge, ix, iy, iX1, iY1;
  Int_t cathode, maxX, maxY;
  Int_t manuId, manuChannel;
  const AliMpVSegmentation* seg;
  AliMpPad pad;
  AliMUONTriggerGUIboard* board;
  AliMUONVDigit *dig;

  for (Int_t ib = 0; ib < kNBoards; ib++) {

    board = (AliMUONTriggerGUIboard*)fBoards->At(ib);
    if (board == 0) continue;

    nstripX = board->GetNStripX();
    nstripY = board->GetNStripY();
    for (Int_t ichamber = 11; ichamber <= 14; ichamber++) {
      
      detElemId = ichamber * 100 + board->GetDetElemId();
      // x strips
      cathode = 0;
      for (Int_t isx = 0; isx < nstripX; isx++) {
	
	charge = (Int_t)board->GetXDig(ichamber-11,isx);
	if (charge == 0) continue;
        ix  = board->GetXSix();
        iY1 = board->GetXSiy1();
        iy  = isx + iY1;
	seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(cathode));
	maxX = seg->MaxPadIndexX();
	maxY = seg->MaxPadIndexY();
	if (ix > maxX) printf("Index x > maximum!\n");
	if (iy > maxY) printf("Index y > maximum!\n");
        pad = seg->PadByIndices(ix,iy,kTRUE);
	manuId = pad.GetLocalBoardId(0);
	manuChannel = pad.GetLocalBoardChannel(0);

	dig = fDigitStore->Add(detElemId,manuId,manuChannel,cathode,AliMUONVDigitStore::kAllow);
	dig->SetCharge(charge);
	dig->SetPadXY(ix,iy);
	//printf("Cathode 0: ix %3d iy %3d manuId %3d manuChannel %3d \n",ix,iy,manuId,manuChannel);

      }  // end x strips

      // y strips
      cathode = 1;
      for (Int_t isy = 0; isy < nstripY; isy++) {
	
	charge = board->GetYDig(ichamber-11,isy);
	if (charge == 0) continue;
        iX1 = board->GetYSix1();
        ix  = isy + iX1;
        iy  = board->GetYSiy();
	seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(cathode));
	maxX = seg->MaxPadIndexX();
	maxY = seg->MaxPadIndexY();
	if (ix > maxX) printf("Index x > maximum!\n");
	if (iy > maxY) printf("Index y > maximum!\n");
        pad = seg->PadByIndices(ix,iy,kTRUE);
	manuId = pad.GetLocalBoardId(0);
	manuChannel = pad.GetLocalBoardChannel(0);

	dig = fDigitStore->Add(detElemId,manuId,manuChannel,cathode,AliMUONVDigitStore::kAllow);
	dig->SetCharge(charge);
	dig->SetPadXY(ix,iy);
	//printf("Cathode 1: ix %3d iy %3d manuId %3d manuChannel %3d \n",ix,iy,manuId,manuChannel);

      }  // end y strips

    }  // end chambers

  }  // end boards

}

//__________________________________________________________________________
void AliMUONTriggerGUI::PrintDigitStore() const
{
  /// Print the digits created in the GUI

  const AliMpVSegmentation* seg;
  AliMpPad pad;
  Int_t ix, iy, charge, detElemId, cathode;

  Int_t chamber;
  for (Int_t i = 0; i < kNMT; i++) {

    chamber = 11+i;

    MpPair_t deRange = AliMpDEManager::GetDetElemIdRange(chamber-1);
    TIter next(fDigitStore->CreateIterator(AliMp::PairFirst(deRange),AliMp::PairSecond(deRange)));
    AliMUONVDigit *mdig;

    while ( ( mdig = static_cast<AliMUONVDigit*>(next())) )
    {
      cathode = mdig->Cathode();
      
      ix = mdig->PadX();
      iy = mdig->PadY();
      detElemId = mdig->DetElemId(); 
      charge = (Int_t)mdig->Charge();

      seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(cathode));  
      pad = seg->PadByIndices(ix,iy,kTRUE);
 
      printf("Digit: detElemId %4d cath %1d ix %2d iy %3d charge %1d \n",detElemId,cathode,ix,iy,charge);
      
    }

  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::CreateTriggerStore() 
{
  /// Process the DSET digit store and fill the trigger store

  if (fDigitStore->GetSize() == 0) {
    printf("The digit store is empty...\n");
    return;
  }
  fTStoreOn = kTRUE;

  AliMUONVDigitStore *digitStore = static_cast<AliMUONVDigitStore*>(fDigitStore);

  fTriggerProcessor->Digits2Trigger(*digitStore,*fTriggerStore);

}

//__________________________________________________________________________
void AliMUONTriggerGUI::PrintTriggerStore() const
{
  /// Print the trigger output for DSET digits store

  if (!fTStoreOn) return;

  TIter next(fTriggerStore->CreateLocalIterator());
  
  UShort_t x2m, x2u, x2d;
  Int_t loStripX, loStripY, loDev, loCircuit, iStripX, iStripY, loLpt, loHpt;
  AliMUONLocalTrigger *mlt;
  while ( ( mlt = static_cast<AliMUONLocalTrigger*>(next()) ) )
  {    
    loCircuit = mlt->LoCircuit();

    AliMUONLocalTriggerBoard* ltb = fCrateManager->LocalBoard(loCircuit);
    x2d = ltb->GetSwitch(0);
    x2m = ltb->GetSwitch(1);
    x2u = ltb->GetSwitch(2);
    
    loStripX = mlt->LoStripX();
    loStripY = mlt->LoStripY();
    loDev    = mlt->LoDev();
    loLpt    = mlt->LoLpt();
    loHpt    = mlt->LoHpt();

    iStripX = loStripX/2;
    if ((x2u == 1 || x2m == 1 || x2d == 1) && x2m == 1) {
      iStripY = loStripY/2;
    } else {
      iStripY = loStripY;
    }
    
    printf("Circ %3d Xs %2d Ys %2d Dev %2d Lpt %1d Hpt %1d \n",loCircuit,loStripX,loStripY,loDev,loLpt,loHpt);

    AliMUONGlobalTrigger *globalTrigger = fTriggerStore->Global();
    globalTrigger->Print();

  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::WriteTriggerRawData() 
{
  /// Write raw data (DATE and ROOT) for the trigger store from 
  /// the DSET digit store

  if (!fTStoreOn) {
    printf("The trigger store is empty... \n");
    return;
  }

  // raw data (ddl)

  AliMUONRawWriter *rawWriter = new AliMUONRawWriter();
  AliRawDataHeaderSim header;
  rawWriter->SetHeader(header);
  rawWriter->Digits2Raw(0,fTriggerStore);

  delete rawWriter;

  // raw data (ddl) to date
  // AliSimulation::ConvertRawFilesToDate

  char command[256];
  char dateFileName[256];
  snprintf(dateFileName,256,"TriggerGUI.date");
  snprintf(command, 256,"dateStream -c -s -D -o %s -# %d -C -run %d",
          dateFileName, 1, 0);
  FILE* pipe = gSystem->OpenPipe(command, "w");

  UInt_t detectorPattern = 0;
  fprintf(pipe, "GDC DetectorPattern %u\n", detectorPattern);
  Float_t ldc = 0;
  Int_t prevLDC = -1;
  
  // loop over detectors and DDLs
  for (Int_t iDet = 0; iDet < AliDAQ::kNDetectors; iDet++) {
    for (Int_t iDDL = 0; iDDL < AliDAQ::NumberOfDdls(iDet); iDDL++) {
      
      Int_t ddlID = AliDAQ::DdlID(iDet,iDDL);
      Int_t ldcID = Int_t(ldc + 0.0001);
      ldc += AliDAQ::NumberOfLdcs(iDet) / AliDAQ::NumberOfDdls(iDet);
      
      char rawFileName[256];
      snprintf(rawFileName, 256,"%s",AliDAQ::DdlFileName(iDet,iDDL));
      // check existence and size of raw data file
      FILE* file = fopen(rawFileName, "rb");
      if (!file) continue;
      fseek(file, 0, SEEK_END);
      unsigned long size = ftell(file);
      fclose(file);
      if (!size) continue;

      if (ldcID != prevLDC) {
	fprintf(pipe, " LDC Id %d\n", ldcID);
	prevLDC = ldcID;
      }
      fprintf(pipe, "  Equipment Id %d Payload %s\n", ddlID, rawFileName);
    }
  }
  Int_t result = gSystem->ClosePipe(pipe);

  // raw data (date) to root
  // AliSimulation::ConvertDateToRoot

  char rootFileName[256];
  snprintf(rootFileName,256,"TriggerGUI.root");

  // ALIMDC setup
  const Int_t kDBSize = 2000000000;
  const Int_t kTagDBSize = 1000000000;
  const Bool_t kFilter = kFALSE;
  const Int_t kCompression = 1;

  char* path = gSystem->Which(gSystem->Getenv("PATH"), "alimdc");
  if (!path) {
    printf("the program alimdc was not found\n");
    return;
  } else {
    delete[] path;
  }

  printf("converting DATE file %s to root file %s \n",
               dateFileName, rootFileName);

  const char* rawDBFS[2] = { "/tmp/mdc1", "/tmp/mdc2" };
  const char* tagDBFS    = "/tmp/mdc1/tags";

  // User defined file system locations
  if (gSystem->Getenv("ALIMDC_RAWDB1"))
    rawDBFS[0] = gSystem->Getenv("ALIMDC_RAWDB1");
  if (gSystem->Getenv("ALIMDC_RAWDB2"))
    rawDBFS[1] = gSystem->Getenv("ALIMDC_RAWDB2");
  if (gSystem->Getenv("ALIMDC_TAGDB"))
    tagDBFS = gSystem->Getenv("ALIMDC_TAGDB");

  gSystem->Exec(Form("rm -rf %s",rawDBFS[0]));
  gSystem->Exec(Form("rm -rf %s",rawDBFS[1]));
  gSystem->Exec(Form("rm -rf %s",tagDBFS));

  gSystem->Exec(Form("mkdir %s",rawDBFS[0]));
  gSystem->Exec(Form("mkdir %s",rawDBFS[1]));
  gSystem->Exec(Form("mkdir %s",tagDBFS));

  result = gSystem->Exec(Form("alimdc %d %d %d %d %s",
			      kDBSize, kTagDBSize, kFilter, kCompression,
			      dateFileName));
  gSystem->Exec(Form("mv %s/*.root %s", rawDBFS[0], rootFileName));

  gSystem->Exec(Form("rm -rf %s",rawDBFS[0]));
  gSystem->Exec(Form("rm -rf %s",rawDBFS[1]));
  gSystem->Exec(Form("rm -rf %s",tagDBFS));

}

//__________________________________________________________________________
void AliMUONTriggerGUI::FETboard(Int_t ib, Int_t amp) 
{
  /// Front End test set all strips for board with index "ib"
  /// AliMUONTriggerGUIbdmap::DoDigits()

  AliMUONTriggerGUIboard *board;
  AliMUONTriggerGUIboard *btmp;
  Int_t pos, over, number, nStripX, nStripY;

  board = (AliMUONTriggerGUIboard*)fBoards->At(ib);
  if (board == 0) return;
  
  nStripX = board->GetXSiy2() - board->GetXSiy1() + 1;
  nStripY = board->GetYSix2() - board->GetYSix1() + 1;
  
  number = board->GetNumber();
  pos    = board->GetPosition();
  over   = board->GetYOver();
  
  for (Int_t imt = 0; imt < kNMT; imt++) {
    
    for (Int_t ix = 0; ix < nStripX; ix++) {
      board->SetDigitX(imt,ix,amp);
    }
    
    for (Int_t iy = 0; iy < nStripY; iy++) {
      board->SetDigitY(imt,iy,amp);
      
      // extended y-strips
      for (Int_t io = 1; io <= over; io++) {
	if (io == pos) continue;
	btmp = (AliMUONTriggerGUIboard*)fBoards->UncheckedAt(number+io-pos);
	btmp->SetDigitY(imt,iy,amp);
      }
      
    }
    
  }
  
}

//__________________________________________________________________________
void AliMUONTriggerGUI::FET(Int_t onoff) 
{
  /// Front End test set all strips for all boards
  /// AliMUONTriggerGUIbdmap::DoDigits()

  Int_t amp = onoff;

  for (Int_t ib = 0; ib < kNBoards; ib++) {
    FETboard(ib,amp);
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUI::ClearDigitStore() 
{
  /// Clear the DSET digit store

  fDigitStore->Clear();

}

//__________________________________________________________________________
void AliMUONTriggerGUI::ClearTriggerStore() 
{
  /// Clear the trigger store from the DSET digit store

  fTriggerStore->Clear();
  fTStoreOn = kFALSE;

}
