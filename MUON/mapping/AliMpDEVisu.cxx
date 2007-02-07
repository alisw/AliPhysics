#include <TObject.h>
#include <TList.h>
#include <TString.h>
#include <TVector2.h>
#include <TCanvas.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>
#include <TGLabel.h>
#include <TGComboBox.h>
#include <TGNumberEntry.h>
#include <TGTextView.h>
#include <TGTextEntry.h>

#include "AliMpSlatMotifMap.h"
#include "AliMpSt345Reader.h"
#include "AliMpSectorReader.h"
#include "AliMpSlat.h"
#include "AliMpPCB.h"
#include "AliMpPCBPainter.h"
#include "AliMpSectorReader.h"
#include "AliMpSector.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpRowPainter.h"
#include "AliMpVPainter.h"
#include "AliMpMotifPainter.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifMap.h"

#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpStationType.h"
#include "AliMpSegmentation.h"
#include "AliMpPad.h"
#include "AliMpDDLStore.h"

#include "AliMpDEVisu.h"

// Category: graphics
//
// Class AliMpDEVisu
// -----------------------
// GUI for drawing segmentation
// motif manu and associated channels
// date: 2007/01/26
// Author: Ch. Finck

/// \cond CLASSIMP
ClassImp(AliMpDEVisu)
/// \endcond

//_________________________________________________________________
AliMpDEVisu::AliMpDEVisu(UInt_t w, UInt_t h) 
    : TGFrame(gClient->GetRoot(), w, h),
      fkMainWindow(gClient->GetRoot()),
      fMain(new TGMainFrame(gClient->GetRoot(), w, h)),
      fEcanvas(0x0),
      fChamberCombo(0x0),
      fDECombo(0x0),
      fNumberEntry(0x0),
      fPlaneButton(0x0),
      fNameDEView(0x0),
      fLogMessage(0x0),
      fLogFile(0x0),
      fTrashList(0x0),
      fDEComboIdx(),
      fCurrentPlane(AliMp::kBendingPlane),
      fCurrentDetElem(100),
      fCurrentDEName(),
      fSegmentation(),
      fDDLStore(AliMpDDLStore::Instance())

{

// Create canvas widget

    Int_t width  = Int_t(w*0.99);
    Int_t height = Int_t(h*0.99);

    fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain, width, height);
    fEcanvas->GetCanvas()->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)",
				    "AliMpDEVisu",
				    this,
				    "HandleMovement(Int_t, Int_t, Int_t, TObject*)");

// Create a horizontal frame widget with buttons
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,width,height/6);

    TGTextButton *draw = new TGTextButton(hframe,"&Draw");
    draw->Connect("Clicked()","AliMpDEVisu",this,"DrawDE()");
    hframe->AddFrame(draw, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,15,5,3,4));

    TGTextButton *exit = new TGTextButton(hframe,"&Exit","gApplication->Terminate(0)");
    hframe->AddFrame(exit, new TGLayoutHints(kLHintsLeft| kLHintsCenterY,5,5,3,4));


    Int_t i = 0;
    Char_t text[20];

// chamber label
    TGLabel* chamberLabel = new TGLabel(hframe, "Chamber :");
    hframe->AddFrame(chamberLabel, new TGLayoutHints(kLHintsCenterX | kLHintsLeft | kLHintsCenterY, 10, 0, 20, 0));
    fChamberCombo = new TGComboBox(hframe, kChamberCombo);

    fDEComboIdx.Set(26);
    for(i = 0; i < 10; i++)
    {
      sprintf(text,"%d",i+1);
      fChamberCombo->AddEntry(text,i);
    }
    fChamberCombo->Resize(40,20);
    fChamberCombo->Select(0);
    fChamberCombo->Associate(this);
    hframe->AddFrame(fChamberCombo, new TGLayoutHints(kLHintsLeft, 10, 0, 9, 0));

// DE label
    TGLabel*  detElemLabel = new TGLabel(hframe, "DE :");
    hframe->AddFrame(detElemLabel, new TGLayoutHints(kLHintsCenterX | kLHintsLeft | kLHintsCenterY, 10, 0, 20, 0));
    fDECombo = new TGComboBox(hframe, kDECombo);
    UpdateComboDE();

    fDECombo->Resize(80,20);
    fDECombo->Select(0);
    fDECombo->Associate(this);
    hframe->AddFrame(fDECombo, new TGLayoutHints(kLHintsLeft, 10, 0, 9, 0));


// DE name
    TGLabel*  detElemName = new TGLabel(hframe, "Name :");
    hframe->AddFrame(detElemName, new TGLayoutHints(kLHintsCenterX | kLHintsLeft | kLHintsCenterY, 10, 0, 20, 0));

    AliMpDetElement* detElem = AliMpDEManager::GetDetElement(fCurrentDetElem);
    fCurrentDEName = detElem->GetDEName();
    fNameDEView = new TGTextView(hframe, 180, 25, fCurrentDEName.Data(), kDEName);
    hframe->AddFrame(fNameDEView, new TGLayoutHints(kLHintsLeft, 10, 0, 9, 0));

// plane type
    fPlaneButton = new TGCheckButton(hframe, "NB Plane", kPlaneType);
    fPlaneButton->SetState(kButtonUp);
    fPlaneButton->Associate(this);
    hframe->AddFrame(fPlaneButton, new TGLayoutHints(kLHintsLeft, 10, 0, 9, 0));


// button motif
    TGTextButton* drawManu = new TGTextButton(hframe,"&Search");
    drawManu->Connect("Clicked()","AliMpDEVisu",this,"DrawManuMotif(Bool_t)");
    drawManu->SetToolTipText("Search for a given manu number");
    hframe->AddFrame(drawManu, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,15,5,3,4));

// entry manu
    fNumberEntry  = new TGNumberEntry(hframe, 0, 4, -1, 
					TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
					TGNumberFormat::kNELLimitMinMax, 1, 1500);
    fNumberEntry->Resize(60,20);
    hframe->AddFrame(fNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 0, 9, 0));

// reset button
    TGTextButton *resetManu = new TGTextButton(hframe,"&Reset");
    resetManu->Connect("Clicked()","AliMpDEVisu",this,"ResetManu()");
    hframe->AddFrame(resetManu, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,10,5,3,4));

// delete window button
    TGTextButton* deletePopup = new TGTextButton(hframe,"&Delete_Popup");
    deletePopup->Connect("Clicked()","AliMpDEVisu",this,"DeletePopUp()");
    deletePopup->SetToolTipText("Delete motif popup window");
    hframe->AddFrame(deletePopup, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,15,5,3,4));


// log mesg
    TGHorizontalFrame *logFrame = new TGHorizontalFrame(fMain, width/2, height/6);

    TGLabel*  logName = new TGLabel(logFrame, "Log MSG:");
    logFrame->AddFrame(logName, new TGLayoutHints(kLHintsCenterX | kLHintsLeft | kLHintsCenterY, 10, 0, 2, 0));

    fLogMessage = new TGTextView(logFrame, width/2, 60);
    fLogMessage->ShowBottom();

    logFrame->AddFrame(fLogMessage, new TGLayoutHints(kLHintsLeft, 10, 0, 2, 0));

// clear log mesg
  TGTextButton* clearLog = new TGTextButton(logFrame,"&Clear");
    clearLog->Connect("Clicked()","AliMpDEVisu",this,"ClearLogMessage()");
    clearLog->SetToolTipText("Clear log message");
    logFrame->AddFrame(clearLog, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,15,5,3,4));

// save log mesg
    TGTextButton* saveLog = new TGTextButton(logFrame,"&Save");
    saveLog->Connect("Clicked()","AliMpDEVisu",this,"SaveLogMessage()");
    saveLog->SetToolTipText("Save log message into file");
    logFrame->AddFrame(saveLog, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,5,5,3,4));

// log file name
    fLogFile = new TGTextEntry(logFrame,"AliMpDEVisu.log");
    fLogFile->SetToolTipText("Default log file name");
    logFrame->AddFrame(fLogFile, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 0, 9, 0));

// frame
    fMain->AddFrame(hframe, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,2,2,10,10));
    fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsCenterX,	10,10,10,10));

    fMain->AddFrame(logFrame, new TGLayoutHints(kLHintsLeft,2,2,2,10));

// Set a name to the main frame
    fMain->SetWindowName("MUON Detection Element Visualization");

// Map all subwindows of main frame
    fMain->MapSubwindows();

// Initialize the layout algorithm
    fMain->Resize(fMain->GetDefaultSize());

// Map main frame
    fMain->MapWindow();

// instance segmentation
    fSegmentation = AliMpSegmentation::Instance()->GetMpSegmentation(fCurrentDetElem,
								     AliMp::GetCathodType(fCurrentPlane));
    fLogMessage->AddLine("Segmentation loaded");
    fLogMessage->ShowBottom();
}

//____________________________________________________________
AliMpDEVisu::~AliMpDEVisu() 
{
// Clean up used widgets: frames, buttons, layouthints
    
    fChamberCombo->Delete();
    fDECombo->Delete();
    fNumberEntry->Delete();
    fPlaneButton->Delete();
    fNameDEView->Delete();
    fLogMessage->Delete();
    fLogFile->Delete();
    fMain->Cleanup();
    delete fMain;
}

//__________________________________________________________
void AliMpDEVisu::HandleMovement(Int_t eventType, Int_t eventX, Int_t eventY, TObject* /*select*/)
{
/// handle cursor mouvement

    if (eventType == 11) {// 61) {// double click
	
      TCanvas *canvas = fEcanvas->GetCanvas();
      canvas->cd(1);

// estimate graphic pad sizes
      TVector2 ul(gPad->XtoPixel(0.02), gPad->YtoPixel(0.98));
      TVector2 br(gPad->XtoPixel(1.00), gPad->YtoPixel(0.02));

      TVector2 padDim = br - ul;
     
      fSegmentation = AliMpSegmentation::Instance()
	  ->GetMpSegmentation(fCurrentDetElem, AliMp::GetCathodType(fCurrentPlane));

// get DE dimension (half size)
      TVector2 deDim =  fSegmentation-> Dimensions();
      deDim *= 2.;

// calculated in pixel 
      Double_t x = 0.;
      Double_t y = 0.;

     if (fCurrentDetElem >= 500) {
       // origin at center
       x = (eventX - ul.X() - padDim.X()/2.)/(padDim.X())*deDim.X(); 
       y = (br.Y() - eventY - padDim.Y()/2.)/(padDim.Y())*deDim.Y(); 

     } else {
       // origin at bl
       x = (eventX - ul.X())/(padDim.X())*deDim.X(); 
       y = (br.Y() - eventY)/(padDim.Y())*deDim.Y(); 

     }

      TVector2 padReal(x,y);
  

// get manu
      AliMpPad pad = fSegmentation->PadByPosition(padReal, false);

      Char_t log[255];
      if (!pad.IsValid()) {
	sprintf(log, "PopupManuMotif: No manu for DE: %d at position (%5.2f, %5.2f)",
			fCurrentDetElem, x, y);
	fLogMessage->AddLine(log);
	fLogMessage->ShowBottom();
	DrawDE();
	return;
      }

      Int_t manu = pad.GetLocation().GetFirst();
	
      fNumberEntry->SetNumber(manu);

      sprintf(log, "PopupManuMotif: DE: %d, manu: %d at position: %5.2f, %5.2f", fCurrentDetElem, manu, x, y);
      fLogMessage->AddLine(log);
      fLogMessage->ShowBottom();

      DrawManuMotif(true);

    }
}

//__________________________________________________________
void AliMpDEVisu::DrawDE() 
{
/// Draws function graphics in randomly choosen interval

    InfoDE();

    if (fCurrentDetElem >= 500) {

      DrawSlat("PMCI");

    } else {

      DrawQuadrant("RSMCI");
    
    }
    DeletePopUp();
}

//__________________________________________________________
void AliMpDEVisu::DrawManuMotif(Bool_t popup) 
{
//  Draw manu motif in yellow and popup channel motif window
 

    if (!fNumberEntry->GetIntNumber()) return;

    Char_t command[255];

    if (fCurrentDetElem >= 500) {
 
      sprintf(command, "%s%d", "PMCI:", (Int_t)fNumberEntry->GetIntNumber());

      DrawSlat(command, popup);

    } else {   

      sprintf(command, "%s%d", "RSMCI:", (Int_t)fNumberEntry->GetIntNumber());
  
      DrawQuadrant(command, popup);

    }
}

//__________________________________________________________
void AliMpDEVisu::DrawSlat(Option_t* option, Bool_t popup) 
{
/// draw slat segmentation

    TCanvas *canvas = fEcanvas->GetCanvas();
    canvas->Clear();

    AliMpDetElement* detElem = AliMpDEManager::GetDetElement(fCurrentDetElem);
    TString nameType =  detElem->GetSegType();

    AliMpSlatMotifMap mm;
    AliMpSt345Reader reader(mm);  
    AliMpSlat* slatCurrent = reader.ReadSlat(nameType.Data(), fCurrentPlane);
    AliMpSlat* slatOther   = reader.ReadSlat(nameType.Data(), AliMp::OtherPlaneType(fCurrentPlane));

    canvas->Divide(1);
    canvas->cd(1);

    AliMpVPainter* painter = AliMpVPainter::CreatePainter(slatCurrent);
    painter->Draw(option);

    canvas->Update();

    delete painter;

    Int_t numberOfManu =  slatCurrent->GetNofElectronicCards() + 
					     slatOther->GetNofElectronicCards();

    Char_t log[255];
  
    sprintf(log, "DrawSlat: number of manus: %d", numberOfManu);

    fLogMessage->AddLine(log);
    fLogMessage->ShowBottom();   

    if (popup)
	PopUpManuMotif(slatCurrent);

}

//__________________________________________________________
void AliMpDEVisu:: DrawQuadrant(Option_t* option, Bool_t popup)
{
/// draw quadrant segmentation

    TCanvas *canvas = fEcanvas->GetCanvas();
    canvas->Clear();

    AliMpDetElement* detElem = AliMpDEManager::GetDetElement(fCurrentDetElem);
    AliMp::StationType  station = detElem->GetStationType();

    AliMpSectorReader readerCurrent(station, fCurrentPlane);
    AliMpSector* sectorCurrent = readerCurrent.BuildSector();

    AliMpSectorReader readerOther(station, AliMp::OtherPlaneType(fCurrentPlane));
    AliMpSector* sectorOther = readerOther.BuildSector();
     
    canvas->Divide(1);
    canvas->cd(1);

    AliMpVPainter *painter = AliMpVPainter::CreatePainter(sectorCurrent);
    painter->Draw(option);

    canvas->Update();

    delete painter;

    Int_t numberOfManu =  
      sectorCurrent->GetMotifMap()->GetNofMotifPositions() + 
      sectorOther->GetMotifMap()->GetNofMotifPositions();

    Char_t log[255];
  
    sprintf(log, "DrawQuadrant: number of manus: %d", numberOfManu);

    fLogMessage->AddLine(log);
    fLogMessage->ShowBottom();   

    if (popup)
	PopUpManuMotif(sectorCurrent);
}

//__________________________________________________________
void AliMpDEVisu::ResetManu() 
{
/// reset manu search entry 
    fLogMessage->AddLine("Reset Motif Search Entry:");
    fLogMessage->ShowBottom();
    fNumberEntry->SetNumber(0);

}

//__________________________________________________________
void AliMpDEVisu::DeletePopUp() 
{
/// delete motif popup windows 
  
    if (fTrashList.GetEntries() > 0) {
      fLogMessage->AddLine("Delete Motif PopUp Windows:");
      fLogMessage->ShowBottom();

//       for (Int_t i = 0; i < fTrashList.GetEntries(); ++i) {

// 	TGTransientFrame* trans = (TGTransientFrame*)fTrashList.At(i);
// 	if (trans)
// 	    delete trans;
//       }

      fTrashList.Delete();
    }
}

//__________________________________________________________
void AliMpDEVisu::SaveLogMessage() 
{
/// save log message into log file
    
    TString logFile = fLogFile->GetDisplayText();
    fLogMessage->GetText()->Save(logFile.Data());

    Char_t log[255];
    sprintf(log, "SaveLogMessage: saving log message into logfile: %s", logFile.Data());
    fLogMessage->AddLine(log);
    fLogMessage->ShowBottom();  
}

//__________________________________________________________
void AliMpDEVisu::ClearLogMessage() 
{
/// clear log message 
    fLogMessage->GetText()->Clear();
    Char_t log[255];
    sprintf(log, "ClearLogMessage: clear log messages:");
    fLogMessage->AddLine(log);
    fLogMessage->ShowBottom();  
}

//__________________________________________________________
void AliMpDEVisu::InfoDE() 
{
/// info concerning the whole DE

    AliMpDetElement* detElem = fDDLStore->GetDetElement(fCurrentDetElem);
    Int_t ddlId        = detElem->GetDdlId();
    Int_t numberOfBus  = detElem->GetNofBusPatches();
    Int_t firstBus     = detElem->GetBusPatchId(0);
    Int_t lastBus      = detElem->GetBusPatchId(numberOfBus - 1); // expect a continuous numbering 

    detElem = AliMpDEManager::GetDetElement(fCurrentDetElem);
    Int_t numberOfSerialManu = detElem->GetNofManus(); // number of manu with an identified serial number

    Char_t log[255];
    sprintf(log, "DrawDE: detection element: %d, name: %s", fCurrentDetElem, fCurrentDEName.Data());
    fLogMessage->AddLine(log);
    fLogMessage->ShowBottom();   

    
    sprintf(log, "DrawDE: DDL: %d, number of buspatches %d from %d to %d", 
	    ddlId, numberOfBus, firstBus, lastBus);
    fLogMessage->AddLine(log);
    fLogMessage->ShowBottom();

    if (numberOfSerialManu != 0) { // not available yet for all DE 
	sprintf(log, "DrawDE: number of manus with serial number: %d", numberOfSerialManu);
	fLogMessage->AddLine(log);
	fLogMessage->ShowBottom();
    }

}

//__________________________________________________________
Bool_t AliMpDEVisu::ProcessMessage(Long_t msg, Long_t parm1, Long_t /*parm2*/)
{
/// process message from widgets actions/entries

    switch(GET_MSG(msg)) 
    {
    case kC_COMMAND: 
      switch (GET_SUBMSG(msg)) 
      {
      case kCM_COMBOBOX: 

	switch (parm1) 
	{
	case kChamberCombo: 
	    UpdateComboDE();
	    UpdateNameView();
	    break;
   
	case kDECombo:
	    UpdateNameView();
	    break; 
	}
	break;
  
      case kCM_CHECKBUTTON:
 	  if (fPlaneButton->GetState() == kButtonDown) {
	    fCurrentPlane = AliMp::kNonBendingPlane;
	    if (fNumberEntry->GetIntNumber() && fNumberEntry->GetIntNumber() <= 1024)
		fNumberEntry->SetNumber(fNumberEntry->GetIntNumber() + 1024);
	  } else {
	    fCurrentPlane = AliMp::kBendingPlane;
	    if (fNumberEntry->GetIntNumber() && fNumberEntry->GetIntNumber() > 1024)
		fNumberEntry->SetNumber(fNumberEntry->GetIntNumber() - 1024);
	  }
	break;

      }
      break;
    }
    return true;
}

//__________________________________________________________
void AliMpDEVisu::UpdateComboDE()
{
/// update DE in respect to selected chamber
 
    fDECombo->RemoveAll();

    AliMpDEIterator it;
    Int_t i = 0;
    Char_t text[20];

    for ( it.First(fChamberCombo->GetSelected()); ! it.IsDone(); it.Next() ) {
      sprintf(text,"%d",it.CurrentDE()->GetId());
      fDECombo->AddEntry(text,i);
      fDEComboIdx[i++] = it.CurrentDE()->GetId();
    }
    fDECombo->Select(0);
}

//__________________________________________________________
void AliMpDEVisu::UpdateNameView()
{
/// update DE name in respect to selected DE id.

    fNameDEView->Clear();

    fCurrentDetElem = fDEComboIdx[fDECombo->GetSelected()];
    AliMpDetElement* detElem = AliMpDEManager::GetDetElement(fCurrentDetElem);
    fCurrentDEName = detElem->GetDEName();

    fNameDEView->AddLine(fCurrentDEName.Data());
    fNameDEView->ShowBottom();
}

//__________________________________________________________
void AliMpDEVisu::PopUpManuMotif(AliMpSlat* slat)
{
/// pop up manu window motif painter for slat

// Create transient frame
    TGTransientFrame* trans = new TGTransientFrame(fkMainWindow, fMain, 400, 400);
    trans->DontCallClose();
    trans->CenterOnParent();

// fill trash
    fTrashList.Add(trans);

    Char_t title[255];
    sprintf(title,"Manu Motif: %d", (Int_t)fNumberEntry->GetIntNumber()); 
    trans->SetWindowName(title);

// Create canvas widget
    TRootEmbeddedCanvas* eTransCanvas = 
	new TRootEmbeddedCanvas("ETransCanvas",trans,400,400);

    trans->AddFrame(eTransCanvas, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10,10,1,10));

 
// motif painter
    AliMpMotifPosition* motifPosFound = 0x0;

    TCanvas *canvas = eTransCanvas->GetCanvas();
    canvas->Clear();

    for ( AliMpSlat::Size_t i = 0; i < slat->GetSize(); ++i ) {

    AliMpPCB* pcb = slat->GetPCB(i);
    AliMpPCBPainter* pcbPainter = new AliMpPCBPainter(pcb);

    for ( AliMpPCB::Size_t j = 0; j < slat->GetPCB(i)->GetSize(); ++j ) {
      
      AliMpMotifPosition*motifPos = pcb->GetMotifPosition(j);
      
      Int_t manuId = motifPos->GetID();
      if (manuId == (Int_t)fNumberEntry->GetIntNumber()) {
	  motifPosFound = motifPos;
	  break;
      }
    }
    if (motifPosFound)
	break;
    delete pcbPainter;
  }  

  if(motifPosFound) {
// maps
    trans->MapSubwindows();
    trans->MapWindow();
// painter
    AliMpVPainter* painter = AliMpVPainter::CreatePainter(motifPosFound);
    painter->Draw("ZT");
// canvas
    canvas->Update();

  }
 

}
//__________________________________________________________
void AliMpDEVisu::PopUpManuMotif(AliMpSector* sector)
{

/// pop up manu window motif painter for sector

// Create transient frame
    TGTransientFrame* trans = new TGTransientFrame(fkMainWindow, fMain, 400, 400);
    trans->DontCallClose();
    trans->CenterOnParent();

// fill trash
    fTrashList.Add(trans);

    Char_t title[255];
    sprintf(title,"Manu Motif: %d", (Int_t)fNumberEntry->GetIntNumber()); 
    trans->SetWindowName(title);

// Create canvas widget
    TRootEmbeddedCanvas* eTransCanvas = 
	new TRootEmbeddedCanvas("ETransCanvas",trans,400,400);

    trans->AddFrame(eTransCanvas, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,
						10,10,1,10));

// motif painter
    AliMpMotifPosition* motifPosFound = 0x0;

    TCanvas *canvas = eTransCanvas->GetCanvas();
    canvas->Clear();

    for (Int_t iRow = 0; iRow < sector->GetNofRows(); ++iRow) {

      AliMpRow* row = sector->GetRow(iRow);
      AliMpRowPainter* rowPainter = new  AliMpRowPainter(row);

      for (Int_t iRowSeg = 0; iRowSeg < sector->GetRow(iRow)->GetNofRowSegments(); ++iRowSeg){

	for (Int_t iRowSeg = 0; iRowSeg < row->GetNofRowSegments(); ++iRowSeg) {
	  AliMpVRowSegment *rowSegment = row->GetRowSegment(iRowSeg);

	  for (Int_t iMotif = 0; iMotif < rowSegment->GetNofMotifs(); ++iMotif){

	    Int_t motifPositionId = rowSegment->GetMotifPositionId(iMotif);
	    AliMpMotifPosition *motifPos = rowSegment->GetRow()
		->GetMotifMap()->FindMotifPosition(motifPositionId);
		
	    Int_t manuId = motifPos->GetID();
	    if (manuId == (Int_t)fNumberEntry->GetIntNumber()) {
	      motifPosFound = motifPos;
	      break;
	    }
	  }
	  if (motifPosFound)
	      break;
	}
      }
      if (motifPosFound)
	  break;
      delete rowPainter;
    }

    if(motifPosFound) {
// map
      trans->MapSubwindows();
      trans->MapWindow();
// painter
      AliMpVPainter* painter = AliMpVPainter::CreatePainter(motifPosFound);
      painter->Draw("ZT");
// canvas
      canvas->Update();
    }
}
