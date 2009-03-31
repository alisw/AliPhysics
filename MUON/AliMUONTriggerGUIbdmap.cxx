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
/// \class AliMUONTriggerGUIbdmap
///
/// The gui map of a single board object, strips/digits
///
/// \author Bogdan Vulpescu, LPC Clermont-Ferrand
//-----------------------------------------------------------------------------

#include <TPolyLine.h>
#include <TStyle.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TGTextEdit.h>
#include <TGText.h>
#include <TTimer.h>
#include <TObjArray.h>
#include <TCanvas.h>

#include "AliRunLoader.h"
#include "AliMUONVDigit.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONTriggerStoreV1.h"
#include "AliMUONDigitStoreV1.h"
#include "AliMUONTriggerGUIboard.h"
#include "AliMUONTriggerGUIbdmap.h"

#include "AliMpDEManager.h"
#include "AliMpDEIterator.h"
#include "AliMpEncodePair.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"


/// \cond CLASSIMP
ClassImp(AliMUONTriggerGUIbdmap)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerGUIbdmap::AliMUONTriggerGUIbdmap(const TGWindow *p, const TGWindow *mainWindow, UInt_t w, UInt_t h)
  : TGFrame(0),
    fMain(0),
    fLocTrigE(0),
    fBoard(0),
    fLoader(0),
    fMCDataInterface(0),
    fRawDigitStore(0),
    fRawTriggerStore(0),
    fXStrips(0),
    fYStrips(0),
    fEditStrips(0),
    fXOn(0),
    fYOn(0),
    fLabelX(0),
    fLabelY(0),
    fIsEditable(0),
    fCanvasSize(0),
    fNStripX(0),
    fNStripY(0),
    fBoards(0),
    fCalibrationData(0x0),
    fCrateManager(0)
{
  /// frame constructor

  gStyle->SetPadLeftMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.05);

  fMain = new TGTransientFrame(p, mainWindow, w, h, kVerticalFrame);
  fMain->Connect("CloseWindow()", "AliMUONTriggerGUIbdmap", this, "CloseWindow()");
  fMain->DontCallClose(); // to avoid double deletions.
  
  // use hierarchical cleaning
  fMain->SetCleanup(kDeepCleanup);

  fMain->ChangeOptions((fMain->GetOptions() & ~kVerticalFrame) | kHorizontalFrame);

  // strips/digits canvases
  
  fCanvasSize = 200;

  TRootEmbeddedCanvas *recanvas[4];
  Char_t ecname[1];
  for (Int_t i = 0; i < kNMT; i++) {

    sprintf(ecname,"%1d",i+1);
    recanvas[i] = new TRootEmbeddedCanvas(ecname,fMain,fCanvasSize,fCanvasSize);
    
    fCanvas[i] = recanvas[i]->GetCanvas();
    
    fCanvas[i]->SetBorderMode(0);
    fCanvas[i]->SetBit(kNoContextMenu);
    
    fCanvas[i]->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",
			"AliMUONTriggerGUIbdmap",this,
			"EditStrips(Int_t,Int_t,Int_t,TObject*)");

    fMain->AddFrame(recanvas[i],
		    /*
		    new TGTableLayoutHints(2,5,2,6,
					   kLHintsExpandX|kLHintsExpandY |
					   kLHintsShrinkX|kLHintsShrinkY |
					   kLHintsFillX|kLHintsFillY)
		    */
		    new TGLayoutHints(kLHintsTop | 
				      kLHintsLeft, 
				      2, 2, 2, 2)
		    );

  }

  //TGDimension size = fMain->GetDefaultSize();
  //fMain->Resize(size);
    
  // the other frames

  TGCompositeFrame *cf0 = new TGCompositeFrame(fMain, 60, 20, kVerticalFrame);

  TGCompositeFrame *cf1 = new TGCompositeFrame(cf0, 60, 20, kVerticalFrame | kFixedWidth);

  cf1->AddFrame(fXStrips = new TGCheckButton(cf1, "Draw X strips and digits", 1),
		new TGLayoutHints(kLHintsTop | 
				  kLHintsLeft, 
				  2, 2, 2, 2)
		);

  cf1->Resize(fXStrips->GetDefaultWidth()+10, fMain->GetDefaultHeight());

  cf1->AddFrame(fYStrips = new TGCheckButton(cf1, "Draw Y strips and digits", 2),
		new TGLayoutHints(kLHintsTop | 
				  kLHintsLeft, 
				  2, 2, 2, 2)	       
		);

  cf1->Resize(fYStrips->GetDefaultWidth()+10, fMain->GetDefaultHeight());

  cf1->AddFrame(fEditStrips = new TGCheckButton(cf1, "Set/unset strips", 3),
		new TGLayoutHints(kLHintsTop | 
				  kLHintsLeft, 
				  2, 2, 2, 2)	       
		);

  fXStrips->Connect("Clicked()", "AliMUONTriggerGUIbdmap", this, "HandleButtons()");
  fYStrips->Connect("Clicked()", "AliMUONTriggerGUIbdmap", this, "HandleButtons()");
  fEditStrips->Connect("Clicked()", "AliMUONTriggerGUIbdmap", this, "HandleEditButton()");

  cf0->AddFrame(cf1, 
		new TGLayoutHints(kLHintsTop | 
				  kLHintsExpandX, 
				  2, 2, 2, 2)	       
		);
  
  TGCompositeFrame *cf2 = new TGCompositeFrame(cf0, 60, 20, kHorizontalFrame | kFixedWidth);

  TGTextButton *digitsButton = new TGTextButton(cf2, "Digits", 4);
  digitsButton->Connect("Clicked()", "AliMUONTriggerGUIbdmap", this, "DoDigits()");

  //cf2->Resize(digitsButton->GetDefaultWidth()+40, fMain->GetDefaultHeight());

  cf2->AddFrame(digitsButton, 
		new TGLayoutHints(kLHintsTop | 
				  kLHintsLeft |
				  kLHintsExpandX, 
				  2, 2, 5, 5)	       
		);

  TGTextButton *dresetButton = new TGTextButton(cf2, "Reset", 5);
  dresetButton->Connect("Clicked()", "AliMUONTriggerGUIbdmap", this, "ResetDigits()");

  //cf2->Resize(dresetButton->GetDefaultWidth()+40, fMain->GetDefaultHeight());

  cf2->AddFrame(dresetButton, 
		new TGLayoutHints(kLHintsTop | 
				  kLHintsLeft |
				  kLHintsExpandX, 
				  2, 2, 5, 5)	       
		);

  TGTextButton *closeButton = new TGTextButton(cf2, "Close", 6);
  closeButton->Connect("Clicked()", "AliMUONTriggerGUIbdmap", this, "DoClose()");

  //cf2->Resize(closeButton->GetDefaultWidth()+40, fMain->GetDefaultHeight());

  cf2->AddFrame(closeButton, 
		new TGLayoutHints(kLHintsTop | 
				  kLHintsLeft |
				  kLHintsExpandX, 
				  2, 2, 5, 5)	       
		);

  cf0->AddFrame(cf2, 
		  new TGLayoutHints(kLHintsTop | 
				    kLHintsExpandX, 
				    2, 2, 2, 2)	       
		  );
  
  // editor window for the local trigger

  TGCompositeFrame *cf3 = new TGCompositeFrame(cf0, 60, 20, kVerticalFrame | kFixedWidth);

  fLocTrigE = new TGTextEdit(cf3, 100, 100, kSunkenFrame | kDoubleBorder);
  cf3->AddFrame(fLocTrigE, 
		new TGLayoutHints(kLHintsExpandX | 
				  kLHintsExpandY, 
				  3, 3, 3, 3)
		);

  cf0->AddFrame(cf3, 
		  new TGLayoutHints(kLHintsTop | 
				    kLHintsExpandX, 
				    2, 2, 2, 2)	       
		  );
  
  fMain->AddFrame(cf0, 
		  new TGLayoutHints(kLHintsTop | 
				    kLHintsExpandX, 
				    2, 2, 2, 2)	       
		  );
  
  fIsEditable = kFALSE;
  fLabelX = kFALSE;
  fLabelY = kFALSE;

  fMain->MapSubwindows();
  fMain->Resize();
  
  fMain->CenterOnParent();

  //fMain->MapWindow();
  //gClient->WaitFor(fMain);
  
}

//__________________________________________________________________________
AliMUONTriggerGUIbdmap::~AliMUONTriggerGUIbdmap()
{
  /// frame destructor

  fMain->DeleteWindow();

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::Show()
{
  /// map the main frame

  fMain->MapWindow();
  LocalTriggerInfo();

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::LocalTriggerInfo()
{
  /// print the local trigger

  TGText txt;
  Char_t buffer[20];

  sprintf(buffer,"Local trigger info\n");
  fLocTrigE->LoadBuffer(buffer);

  AliMUONVTriggerStore *triggerStore;

  if (fLoader != 0x0) {
    AliRunLoader *runLoader = fLoader->GetRunLoader();
    triggerStore = fMCDataInterface->TriggerStore(runLoader->GetEventNumber());
  } else if (fRawTriggerStore != 0x0) {
    triggerStore = static_cast<AliMUONVTriggerStore*>(fRawTriggerStore);
  } else {
    sprintf(buffer,"No data loaded yet...\n");
    txt.LoadBuffer(buffer);
    fLocTrigE->AddText(&txt);
    return;
  }

  Int_t circuitNumber = fBoard->GetIdCircuit();

  UShort_t x2m, x2u, x2d;

  Int_t loStripX, loStripY, loDev, loCircuit, iStripX, iStripY, loLpt, loHpt;
  AliMUONLocalTrigger *mlt;

  TIter next(triggerStore->CreateLocalIterator());
  
  while ( ( mlt = static_cast<AliMUONLocalTrigger*>(next()) ) )
  {    
    loCircuit = mlt->LoCircuit();

    if (loCircuit == circuitNumber) {

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

      sprintf(buffer,"Circuit = %03d",loCircuit);
      txt.LoadBuffer(buffer);
      fLocTrigE->AddText(&txt);

      sprintf(buffer,"LoStripX = %2d",loStripX);
      txt.LoadBuffer(buffer);
      fLocTrigE->AddText(&txt);

      sprintf(buffer,"LoStripY = %2d",loStripY);
      txt.LoadBuffer(buffer);
      fLocTrigE->AddText(&txt);

      sprintf(buffer,"LoDev = %2d",loDev);
      txt.LoadBuffer(buffer);
      fLocTrigE->AddText(&txt);

      sprintf(buffer,"--------------------");
      txt.LoadBuffer(buffer);
      fLocTrigE->AddText(&txt);

      sprintf(buffer,"X-strip = %2d ( %2d )",iStripX,(loStripX+loDev+1)/2);
      txt.LoadBuffer(buffer);
      fLocTrigE->AddText(&txt);

      sprintf(buffer,"Y-strip = %2d",iStripY);
      txt.LoadBuffer(buffer);
      fLocTrigE->AddText(&txt);

      sprintf(buffer,"--------------------");
      txt.LoadBuffer(buffer);
      fLocTrigE->AddText(&txt);

      sprintf(buffer,"LoLpt = %2d",loLpt);
      txt.LoadBuffer(buffer);
      fLocTrigE->AddText(&txt);

      sprintf(buffer,"LoHpt = %2d",loHpt);
      txt.LoadBuffer(buffer);
      fLocTrigE->AddText(&txt);

      break;

    }

  }
  
}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::Init()
{
  /// initialize the boards in the four MT

  for (Int_t i = 0; i < kNMT; i++) {

    fXWidth[i]  = fBoard->GetXWidth(i);
    fYWidth[i]  = fBoard->GetYWidth(i); 
    fXCenter[i] = fBoard->GetXCenter(i);
    fYCenter[i] = fBoard->GetYCenter(i); 
    
  }

  Float_t xw, yw;
  for (Int_t i = 0; i < kNMT; i++) {
    xw = 1.20*fXWidth[i];
    yw = 1.20*fYWidth[i];
    fCanvas[i]->Range(-xw/2,-yw/2,+xw/2,+yw/2);
    for (Int_t is = 0; is < kNS; is++) {

      fXDigPL[i][is] = new TPolyLine(5);
      fYDigPL[i][is] = new TPolyLine(5);
      fXDigPL[i][is]->SetBit(kCannotPick);
      fYDigPL[i][is]->SetBit(kCannotPick);
      fXDigPL[i][is]->SetLineColor(4);
      fYDigPL[i][is]->SetLineColor(4);

      fXDigBox[i][is] = new TBox(0,0,0,0);
      fYDigBox[i][is] = new TBox(0,0,0,0);
      fXDigBox[i][is]->SetBit(kCannotPick);
      fYDigBox[i][is]->SetBit(kCannotPick);
      fXDigBox[i][is]->SetFillStyle(0);
      fYDigBox[i][is]->SetFillStyle(0);
      fXDigBox[i][is]->SetLineColor(4);
      fYDigBox[i][is]->SetLineColor(4);

      fXLabelL[i][is] = 0;
      fXLabelR[i][is] = 0;
      fYLabelL[i][is] = 0;
      fYLabelR[i][is] = 0;

    }
  }


  fXOn = kFALSE;
  fYOn = kFALSE;

  fNStripX = fBoard->GetXSiy2() - fBoard->GetXSiy1() + 1;
  fNStripY = fBoard->GetYSix2() - fBoard->GetYSix1() + 1;

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::HandleEditButton()
{
  /// handle the editable check button

  if (((fXOn && !fYOn) || (!fXOn && fYOn)) && 
      (fEditStrips->GetState() == kButtonDown)) {
    fIsEditable = kTRUE;
  } else {
    fIsEditable = kFALSE;
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::EditStrips(Int_t event, Int_t x, Int_t y, TObject *sel)
{
  /// set/unset strips

  TString cs;
  Int_t iMT;
  Double_t *px, *py;
  Double_t xf1, yf1, xf2, yf2;
  Int_t np = 5;
  Double_t xMin, xMax, yMin, yMax;
  Float_t xd, yd, fxDim, fyDim, cDim;
  Char_t cln[2];

  cDim = (Float_t)fCanvasSize;
  
  // the (x,y) event does not really go up to the canvas size...
  cDim = 196.0;

  if (fIsEditable) {

    if (event == kButton1Down) {
      
      cs = TString(sel->GetName());
      iMT = cs.Atoi()-1;

      fCanvas[iMT]->cd();

      fCanvas[iMT]->GetRange(xf1,yf1,xf2,yf2);      
      fxDim = (Float_t)xf2;
      fyDim = (Float_t)yf2;

      //xd = fxDim*(+2.0*(Float_t)(x)/cDim - 1.0);
      //yd = fyDim*(-2.0*(Float_t)(y)/cDim + 1.0);

      xd = +(2.0*fxDim*(Float_t)(x))/cDim - fxDim;
      yd = -(2.0*fyDim*(Float_t)(y))/cDim + fyDim;

      if (fXOn) {

	for (Int_t ix = 0; ix < fNStripX; ix++) {
	  
	  px = fXDigPL[iMT][ix]->GetX();
	  py = fXDigPL[iMT][ix]->GetY();

	  xMin = +9999.;
	  xMax = -9999.;
	  yMin = +9999.;
	  yMax = -9999.;
	  for (Int_t ip = 0; ip < np; ip++) {
	    xMin = TMath::Min(xMin,px[ip]);
	    xMax = TMath::Max(xMax,px[ip]);
	    yMin = TMath::Min(yMin,py[ip]);
	    yMax = TMath::Max(yMax,py[ip]);
	  }
	    
	  if (yd > (Float_t)yMin && yd < (Float_t)yMax) {

	    if (fXDigBox[iMT][ix]->GetFillStyle() == 0) {

	      fXDigBox[iMT][ix]->SetFillStyle(1001);
	      fXDigBox[iMT][ix]->SetFillColor(2);

	      sprintf(cln,"%2d",ix);

	      fXLabelL[iMT][ix]->Clear();
	      fXLabelL[iMT][ix]->AddText(cln);
	      fXLabelL[iMT][ix]->Draw();

	      fXLabelR[iMT][ix]->Clear();
	      fXLabelR[iMT][ix]->AddText(cln);
	      fXLabelR[iMT][ix]->Draw();

	      fXDigBox[iMT][ix]->SetX1(xMin);
	      fXDigBox[iMT][ix]->SetY1(yMin);
	      fXDigBox[iMT][ix]->SetX2(xMax);
	      fXDigBox[iMT][ix]->SetY2(yMax);

	      fXDigBox[iMT][ix]->Draw();

	    } else if (fXDigBox[iMT][ix]->GetFillStyle() == 1001) {

	      fXDigBox[iMT][ix]->SetFillStyle(0);

	      fXLabelL[iMT][ix]->Clear();
	      fXLabelL[iMT][ix]->Draw();

	      fXLabelR[iMT][ix]->Clear();
	      fXLabelR[iMT][ix]->Draw();

	      fXDigBox[iMT][ix]->SetX1(-fBoard->GetXCenter(iMT));
	      fXDigBox[iMT][ix]->SetY1(-fBoard->GetYCenter(iMT));
	      fXDigBox[iMT][ix]->SetX2(-fBoard->GetXCenter(iMT));
	      fXDigBox[iMT][ix]->SetY2(-fBoard->GetYCenter(iMT));

	      fXDigBox[iMT][ix]->Draw();

	    }

	    if (!fXDigBox[iMT][ix]->TestBit(kObjInCanvas))
	      fXDigBox[iMT][ix]->Draw();

	    fCanvas[iMT]->Modified();
	    fCanvas[iMT]->Update();

	    break;
	    
	  } 
	}

      }

      if (fYOn) {

	for (Int_t iy = 0; iy < fNStripY; iy++) {

	  px = fYDigPL[iMT][iy]->GetX();
	  py = fYDigPL[iMT][iy]->GetY();

	  xMin = +9999.;
	  xMax = -9999.;
	  yMin = +9999.;
	  yMax = -9999.;
	  for (Int_t ip = 0; ip < np; ip++) {
	    xMin = TMath::Min(xMin,px[ip]);
	    xMax = TMath::Max(xMax,px[ip]);
	    yMin = TMath::Min(yMin,py[ip]);
	    yMax = TMath::Max(yMax,py[ip]);
	  }

	  if (xd > (Float_t)xMin && xd < (Float_t)xMax) {

	    if (fYDigBox[iMT][iy]->GetFillStyle() == 0) {

	      fYDigBox[iMT][iy]->SetFillStyle(1001);
	      fYDigBox[iMT][iy]->SetFillColor(2);

	      sprintf(cln,"%2d",iy);

	      fYLabelL[iMT][iy]->Clear();
	      fYLabelL[iMT][iy]->AddText(cln);
	      fYLabelL[iMT][iy]->Draw();

	      fYLabelR[iMT][iy]->Clear();
	      fYLabelR[iMT][iy]->AddText(cln);
	      fYLabelR[iMT][iy]->Draw();

	      fYDigBox[iMT][iy]->SetX1(xMin);
	      fYDigBox[iMT][iy]->SetY1(yMin);
	      fYDigBox[iMT][iy]->SetX2(xMax);
	      fYDigBox[iMT][iy]->SetY2(yMax);
	      
	      fYDigBox[iMT][iy]->Draw();

	    } else if (fYDigBox[iMT][iy]->GetFillStyle() == 1001) {

	      fYDigBox[iMT][iy]->SetFillStyle(0);

	      fYLabelL[iMT][iy]->Clear();
	      fYLabelL[iMT][iy]->Draw();

	      fYLabelR[iMT][iy]->Clear();
	      fYLabelR[iMT][iy]->Draw();

	      fYDigBox[iMT][iy]->SetX1(-fBoard->GetXCenter(iMT));
	      fYDigBox[iMT][iy]->SetY1(-fBoard->GetYCenter(iMT));
	      fYDigBox[iMT][iy]->SetX2(-fBoard->GetXCenter(iMT));
	      fYDigBox[iMT][iy]->SetY2(-fBoard->GetYCenter(iMT));

	      fYDigBox[iMT][iy]->Draw();

	    }

	    if (!fYDigBox[iMT][iy]->TestBit(kObjInCanvas))
	      fYDigBox[iMT][iy]->Draw();

	    fCanvas[iMT]->Modified();
	    fCanvas[iMT]->Update();

	    break;
	    
	  } 
	}

      }

    }  // end button event

  }  // end IsEditable

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::DoDigits()
{
  /// set the board digits from GUI input (set/unset)

  Int_t amp = 0;
  Int_t number = fBoard->GetNumber();
  Int_t pos, over;
  pos  = fBoard->GetPosition();
  over = fBoard->GetYOver();
  AliMUONTriggerGUIboard *board;

  for (Int_t imt = 0; imt < kNMT; imt++) {

    for (Int_t ix = 0; ix < fNStripX; ix++) {
      if (fXDigBox[imt][ix]->GetFillStyle() ==    0) {
	amp = 0;
      }
      if (fXDigBox[imt][ix]->GetFillStyle() == 1001) {
	amp = 1;
      }
      fBoard->SetDigitX(imt,ix,amp);
    }

    for (Int_t iy = 0; iy < fNStripY; iy++) {
      if (fYDigBox[imt][iy]->GetFillStyle() ==    0) {
	amp = 0;
      }
      if (fYDigBox[imt][iy]->GetFillStyle() == 1001) {
	amp = 1;
      }
      fBoard->SetDigitY(imt,iy,amp);

      // extended y-strips
      for (Int_t io = 1; io <= over; io++) {
	if (io == pos) continue;
	board = (AliMUONTriggerGUIboard*)fBoards->UncheckedAt(number+io-pos);
	board->SetDigitY(imt,iy,amp);
      }

    }

  }

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::ResetDigits()
{
  /// set the board digits from GUI input (set/unset)

  Int_t amp = 0;
  Int_t number = fBoard->GetNumber();
  AliMUONTriggerGUIboard *board;
  Int_t pos, over;
  pos  = fBoard->GetPosition();
  over = fBoard->GetYOver();

  for (Int_t imt = 0; imt < kNMT; imt++) {
    for (Int_t ix = 0; ix < fNStripX; ix++) {
      fXDigBox[imt][ix]->SetFillStyle(0);
      fXDigBox[imt][ix]->SetX1(-fBoard->GetXCenter(imt));
      fXDigBox[imt][ix]->SetY1(-fBoard->GetYCenter(imt));
      fXDigBox[imt][ix]->SetX2(-fBoard->GetXCenter(imt));
      fXDigBox[imt][ix]->SetY2(-fBoard->GetYCenter(imt));

      fBoard->SetDigitX(imt,ix,amp);
    }
    for (Int_t iy = 0; iy < fNStripY; iy++) {
      fYDigBox[imt][iy]->SetFillStyle(0);
      fYDigBox[imt][iy]->SetX1(-fBoard->GetXCenter(imt));
      fYDigBox[imt][iy]->SetY1(-fBoard->GetYCenter(imt));
      fYDigBox[imt][iy]->SetX2(-fBoard->GetXCenter(imt));
      fYDigBox[imt][iy]->SetY2(-fBoard->GetYCenter(imt));

      fBoard->SetDigitY(imt,iy,amp);

      // extended y-strips
      for (Int_t io = 1; io <= over; io++) {
	if (io == pos) continue;
	board = (AliMUONTriggerGUIboard*)fBoards->UncheckedAt(number+io-pos);

	board->SetDigitY(imt,iy,amp);

      }

    }
  }

  fBoard->ClearXDigits();
  fBoard->ClearYDigits();

  // extended y-strips
  for (Int_t io = 1; io <= over; io++) {
    if (io == pos) continue;
    board = (AliMUONTriggerGUIboard*)fBoards->UncheckedAt(number+io-pos);
    board->ClearYDigits();
  }
  
  fXStrips->SetState(kButtonUp);
  fYStrips->SetState(kButtonUp);
  fEditStrips->SetState(kButtonUp);
  fIsEditable = kFALSE;
		      
  DrawClear();
  fXOn = kFALSE;
  fYOn = kFALSE;

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::HandleButtons(Int_t id)
{
  /// handle the check buttons

  if (id == -1) {
    TGButton *btn = (TGButton *) gTQSender;
    id = btn->WidgetId();
  }

  // draw x and y 
  if(fXStrips->GetState() == kButtonDown && 
     fYStrips->GetState() == kButtonDown    ) {
    
    DrawClear();
    DrawStrips(kTRUE,kTRUE);
    DrawDigits(kTRUE,kTRUE);
    fXOn = kTRUE;
    fYOn = kTRUE;

  }

  // draw only x
  if(fXStrips->GetState() == kButtonDown && 
     fYStrips->GetState() == kButtonUp      ) {
    
    DrawClear();
    DrawStrips(kTRUE,kFALSE);
    DrawDigits(kTRUE,kFALSE);
    fXOn = kTRUE;
    fYOn = kFALSE;

  }

  // draw only y
  if(fXStrips->GetState() == kButtonUp   && 
     fYStrips->GetState() == kButtonDown    ) {
    
    DrawClear();
    DrawStrips(kFALSE,kTRUE);
    DrawDigits(kFALSE,kTRUE);
    fXOn = kFALSE;
    fYOn = kTRUE;

  }

  // clear
  if(fXStrips->GetState() == kButtonUp   && 
     fYStrips->GetState() == kButtonUp      ) {
    
    DrawClear();
    fXOn = kFALSE;
    fYOn = kFALSE;

  }

  HandleEditButton();

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::DrawClear()
{
  /// draw the frame of the board image in the canvas
  
  for (Int_t i = 0; i < kNMT; i++) {

    fCanvas[i]->cd();
    fCanvas[i]->Clear();

    fCanvas[i]->Modified();
    fCanvas[i]->Update();

  }

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::DrawDigits(Bool_t bx, Bool_t by)
{
  /// draw the digits in "x" or/and "y"

  Bool_t drawDigits = kTRUE;
  Bool_t drawDigitsRaw = kTRUE;
  if (fLoader == 0x0) {
    drawDigits = kFALSE;
  }
  if (fRawDigitStore == 0x0) {
    drawDigitsRaw = kFALSE;
  }
  
  AliMUONTriggerGUIboard *board;
  Int_t over, pos, number;
  const AliMpVSegmentation* seg;
  AliMpPad pad;
  Int_t cathode, detElemId, ix, iy, charge;
  Int_t chamber, np = 5;
  AliMpDEIterator it;
  Float_t xpmin, xpmax, ypmin, ypmax;
  Float_t xg1, xg2, yg1, yg2, zg1;
  Float_t xdw, ydw, xcw, ycw;
  Double_t xc1, xc2, yc1, yc2;
  Char_t cln[2];
  TBox *boxd;
  Double_t xMin, xMax, yMin, yMax;
  Double_t *px, *py;
    
  number = fBoard->GetNumber();
  pos    = fBoard->GetPosition();
  over   = fBoard->GetYOver();
  
  if (drawDigits || drawDigitsRaw) {

  for (Int_t i = 0; i < kNMT; i++) {
    
    fCanvas[i]->cd();
    
    fCanvas[i]->GetRange(xc1,yc1,xc2,yc2);
    xcw = (Float_t)(0.5*(xc2-xc1));
    ycw = (Float_t)(0.5*(yc2-yc1));

    AliMUONVDigitStore *digitStore = 0x0;
    AliMUONGeometryTransformer transformer;
    transformer.LoadGeometryData("transform.dat");

    if (drawDigits) {
      AliRunLoader *runLoader = fLoader->GetRunLoader();
      digitStore = fMCDataInterface->DigitStore(runLoader->GetEventNumber());
    }
    if (drawDigitsRaw) {
      digitStore = static_cast<AliMUONVDigitStore*>(fRawDigitStore);
    }

    chamber = 11+i;
    
    MpPair_t deRange = AliMpDEManager::GetDetElemIdRange(chamber-1);
    TIter next(digitStore->CreateIterator(AliMp::PairFirst(deRange),AliMp::PairSecond(deRange)));
    AliMUONVDigit *mdig;
    
    while ( ( mdig = static_cast<AliMUONVDigit*>(next())) )
      {
	cathode = mdig->Cathode();
	
	ix = mdig->PadX();
	iy = mdig->PadY();
	detElemId = mdig->DetElemId(); 
	charge = (Int_t)mdig->Charge();
	
	Bool_t triggerBgn = kFALSE;
	Int_t schg = (Int_t)(charge + 0.5);
	// APPLY CONDITION ON SOFT BACKGROUND	
	Int_t tchg = schg - (Int_t(schg/10))*10;	
	if (schg<=10 || tchg>0) {
	  triggerBgn = kFALSE;
	} else {
	  triggerBgn = kTRUE;
	}
	
	if (detElemId%100 != fBoard->GetDetElemId()) continue;
	
	seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(cathode));  
	
	pad = seg->PadByIndices(ix,iy,kTRUE);
	
	//if (cathode == 0) printf("GUI x:  ix %d iy %d \n",ix,iy);
	//if (cathode == 1) printf("GUI y:  ix %d iy %d \n",ix,iy);
	
	// get the pad position and dimensions
	Float_t xlocal1 = pad.Position().X();
	Float_t ylocal1 = pad.Position().Y();
	Float_t xlocal2 = pad.Dimensions().X();
	Float_t ylocal2 = pad.Dimensions().Y();
	
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
	
	xpmin -= fXCenter[i];
	xpmax -= fXCenter[i];
	ypmin -= fYCenter[i];
	ypmax -= fYCenter[i];
	
	xdw = xpmax-xpmin;
	ydw = ypmax-ypmin;
	
	// x-strips
	if ((xdw > ydw) && bx) {
	  
	  //printf("X strips mdig->Cathode() = %1d \n",mdig->Cathode());
	  
	  Int_t iX, iY1, iY2;
	  iX  = fBoard->GetXSix();
	  iY1 = fBoard->GetXSiy1();
	  iY2 = fBoard->GetXSiy2();
	  if (ix == iX && iy >= iY1 && iy <= iY2) {
	    //printf("Digit indices (x-strip) ix = %3d iy = %3d board = %s %d chamber = %2d \n",ix,iy,fBoard->GetBoardName(),fBoard->GetNumber(),chamber); 
	    /*
	      xpmin += 0.01*fXWidth[i];
	      xpmax -= 0.01*fXWidth[i];
	      ypmin += 0.1*ydw;
	      ypmax -= 0.1*ydw;
	    */
	    boxd = new TBox(xpmin,ypmin,xpmax,ypmax);
	    boxd->SetFillStyle(1001);
	    if (triggerBgn) boxd->SetFillColor(6);
	    else boxd->SetFillColor(5);
	    boxd->SetBit(kCannotPick);
	    boxd->Draw();
	    
	    fXDigBox[i][iy-iY1]->SetFillStyle(1001);
	    fXDigBox[i][iy-iY1]->SetFillColor(2);
	    fXDigBox[i][iy-iY1]->SetX1(xpmin);
	    fXDigBox[i][iy-iY1]->SetY1(ypmin);
	    fXDigBox[i][iy-iY1]->SetX2(xpmax);
	    fXDigBox[i][iy-iY1]->SetY2(ypmax);
	    fXDigBox[i][iy-iY1]->Draw();
	    
	    sprintf(cln,"%2d",(iy-iY1));
	    fXLabelL[i][iy-iY1]->Clear();
	    fXLabelL[i][iy-iY1]->AddText(cln);
	    fXLabelL[i][iy-iY1]->Draw();
	    fXLabelR[i][iy-iY1]->Clear();
	    fXLabelR[i][iy-iY1]->AddText(cln);
	    fXLabelR[i][iy-iY1]->Draw();
	    
	    fBoard->SetDigitX(i,iy-iY1,charge);
	    
	  }
	  
	}
	
	// y-strips
	if ((xdw < ydw) && by) {
	  
	  //printf("Y strips mdig->Cathode() = %1d \n",mdig->Cathode());
	  
	  Int_t iX1, iX2, iY;
	  iX1 = fBoard->GetYSix1();
	  iX2 = fBoard->GetYSix2();
	  iY  = fBoard->GetYSiy();
	  if (ix >= iX1 && ix <= iX2 && iy == iY) {
	    //printf("Digit indices (y-strip) ix = %3d iy = %3d board = %s chamber = %2d \n",ix,iy,fBoard->GetBoardName(),chamber); 
	    ypmin = -0.5*fYWidth[i];
	    ypmax = +0.5*fYWidth[i];
	    /*
	      ypmin = -0.5*fYWidth[i];
	      ypmax = +0.5*fYWidth[i];
	      ypmin += 0.01*fYWidth[i];
	      ypmax -= 0.01*fYWidth[i];	  
	      xpmin += 0.1*xdw;
	      xpmax -= 0.1*xdw;
	    */
	    boxd = new TBox(xpmin,ypmin,xpmax,ypmax);
	    boxd->SetFillStyle(1001);
	    if (triggerBgn) boxd->SetFillColor(6);
	    else boxd->SetFillColor(5);
	    boxd->SetBit(kCannotPick);
	    boxd->Draw();
	    
	    fYDigBox[i][ix-iX1]->SetFillStyle(1001);
	    fYDigBox[i][ix-iX1]->SetFillColor(2);
	    fYDigBox[i][ix-iX1]->SetX1(xpmin);
	    fYDigBox[i][ix-iX1]->SetY1(ypmin);
	    fYDigBox[i][ix-iX1]->SetX2(xpmax);
	    fYDigBox[i][ix-iX1]->SetY2(ypmax);
	    fYDigBox[i][ix-iX1]->Draw();
	    
	    sprintf(cln,"%2d",(ix-iX1));
	    fYLabelL[i][ix-iX1]->Clear();
	    fYLabelL[i][ix-iX1]->AddText(cln);
	    fYLabelL[i][ix-iX1]->Draw();
	    fYLabelR[i][ix-iX1]->Clear();
	    fYLabelR[i][ix-iX1]->AddText(cln);
	    fYLabelR[i][ix-iX1]->Draw();
	    
	    fBoard->SetDigitY(i,ix-iX1,charge);
	    
	    // extended y-strips
	    for (Int_t io = 1; io <= over; io++) {
	      if (io == pos) continue;
	      board = (AliMUONTriggerGUIboard*)fBoards->UncheckedAt(number+io-pos);
	      board->SetDigitY(i,ix-iX1,charge);
	    }
	    
	  }
	}
	
      }  // end digits loop

  }  // end chamber loop

  }  // end drawDigits

  // DSET digits
  for (Int_t i = 0; i < kNMT; i++) {
    
    fCanvas[i]->cd();
    
    fCanvas[i]->GetRange(xc1,yc1,xc2,yc2);
    xcw = (Float_t)(0.5*(xc2-xc1));
    ycw = (Float_t)(0.5*(yc2-yc1));
    
    // x-strips DSET
    if (bx) {
      
      for (ix = 0; ix < fNStripX; ix++) {
	//if (fBoard->GetXDig(i,ix) > 0) {
	if (fXDigBox[i][ix]->GetFillStyle() == 1001 || 
	    fBoard->GetXDig(i,ix) > 0) {
	  
	  px = fXDigPL[i][ix]->GetX();
	  py = fXDigPL[i][ix]->GetY();
	  
	  xMin = +9999.;
	  xMax = -9999.;
	  yMin = +9999.;
	  yMax = -9999.;
	  for (Int_t ip = 0; ip < np; ip++) {
	    xMin = TMath::Min(xMin,px[ip]);
	    xMax = TMath::Max(xMax,px[ip]);
	    yMin = TMath::Min(yMin,py[ip]);
	    yMax = TMath::Max(yMax,py[ip]);
	  }
	  
	  if (fXDigBox[i][ix]->GetFillStyle() == 0) {
	    fXDigBox[i][ix]->SetFillStyle(1001);
	    fXDigBox[i][ix]->SetFillColor(2);
	    fXDigBox[i][ix]->SetX1(xMin);
	    fXDigBox[i][ix]->SetY1(yMin);
	    fXDigBox[i][ix]->SetX2(xMax);
	    fXDigBox[i][ix]->SetY2(yMax);
	  }
	  
	  sprintf(cln,"%2d",ix);
	  
	  fXLabelL[i][ix]->Clear();
	  fXLabelL[i][ix]->AddText(cln);
	  fXLabelL[i][ix]->Draw();
	  
	  fXLabelR[i][ix]->Clear();
	  fXLabelR[i][ix]->AddText(cln);
	  fXLabelR[i][ix]->Draw();
	  
	  fXDigBox[i][ix]->Draw();
	  
	}
	
      }

    }
      
    // y-strips DSET
    if (by) {
      
      for (iy = 0; iy < fNStripY; iy++) {
	//if (fBoard->GetYDig(i,iy) > 0) {
	if (fYDigBox[i][iy]->GetFillStyle() == 1001 || 
	    fBoard->GetYDig(i,iy) > 0) {
	  
	  px = fYDigPL[i][iy]->GetX();
	  py = fYDigPL[i][iy]->GetY();
	  
	  xMin = +9999.;
	  xMax = -9999.;
	  yMin = +9999.;
	  yMax = -9999.;
	  for (Int_t ip = 0; ip < np; ip++) {
	    xMin = TMath::Min(xMin,px[ip]);
	    xMax = TMath::Max(xMax,px[ip]);
	    yMin = TMath::Min(yMin,py[ip]);
	    yMax = TMath::Max(yMax,py[ip]);
	  }
	  
	  if (fYDigBox[i][iy]->GetFillStyle() == 0) {
	    fYDigBox[i][iy]->SetFillStyle(1001);
	    fYDigBox[i][iy]->SetFillColor(2);
	    fYDigBox[i][iy]->SetX1(xMin);
	    fYDigBox[i][iy]->SetY1(yMin);
	    fYDigBox[i][iy]->SetX2(xMax);
	    fYDigBox[i][iy]->SetY2(yMax);
	  }
	  
	  sprintf(cln,"%2d",iy);
	  
	  fYLabelL[i][iy]->Clear();
	  fYLabelL[i][iy]->AddText(cln);
	  fYLabelL[i][iy]->Draw();
	  
	  fYLabelR[i][iy]->Clear();
	  fYLabelR[i][iy]->AddText(cln);
	  fYLabelR[i][iy]->Draw();
	  
	  fYDigBox[i][iy]->Draw();
	  
	}
	
      }

    }
      
    fCanvas[i]->Modified();
    fCanvas[i]->Update();
      
  }  // end canvas (chamber) loop

  fMain->MapWindow();

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::DrawStrips(Bool_t bx, Bool_t by)
{
  /// draw the "x" or/and "y" strips

  AliMUONGeometryTransformer transformer;
  transformer.LoadGeometryData("transform.dat");

  const AliMpVSegmentation* seg;
  AliMpPad pad;
  AliMpDEIterator it;
  Int_t chamber;
  Float_t xg1, xg2, yg1, yg2, zg1;
  Float_t xlocal1, xlocal2, ylocal1, ylocal2;
  Int_t detElemId, maxX, maxY;
  Float_t xdw, ydw, xpmin, xpmax, ypmin, ypmax;
  Float_t ptx1, ptx2, pty1, pty2;
  Char_t cln[2];
  Double_t xc1, xc2, yc1, yc2;
  Float_t xcw, ycw;

  Bool_t makeLabelsX = kFALSE;
  Bool_t makeLabelsY = kFALSE;

  if (bx && !fLabelX) {
    makeLabelsX = kTRUE;
    fLabelX = kTRUE;
  }

  if (by && !fLabelY) {
    makeLabelsY = kTRUE;
    fLabelY = kTRUE;
  }

  for (Int_t i = 0; i < kNMT; i++) {

    fCanvas[i]->cd();
    
    fCanvas[i]->GetRange(xc1,yc1,xc2,yc2);
    xcw = (Float_t)(0.5*(xc2-xc1));
    ycw = (Float_t)(0.5*(yc2-yc1));

    chamber = 11+i;

    for ( it.First(chamber-1); ! it.IsDone(); it.Next() ) {
    
      detElemId = it.CurrentDEId();
    
      if (detElemId%100 != fBoard->GetDetElemId()) continue;

      /*---------- y-pads ic = 1 ----------*/
      
      seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::kCath1);  

      maxX = seg->MaxPadIndexX();
      maxY = seg->MaxPadIndexY();
      
      for (Int_t ix = 0; ix <= maxX; ix++) {
	for (Int_t iy = 0; iy <= maxY; iy++) {
	  
	  pad = seg->PadByIndices(ix,iy,kFALSE);
	  
	  if (!pad.IsValid()) continue;
	  
	  // get the pad position and dimensions
	  xlocal1 = pad.Position().X();
	  ylocal1 = pad.Position().Y();
	  xlocal2 = pad.Dimensions().X();
	  ylocal2 = pad.Dimensions().Y();
	  
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
	  
	  xpmin -= fXCenter[i];
	  xpmax -= fXCenter[i];
	  ypmin -= fYCenter[i];
	  ypmax -= fYCenter[i];
	  
	  xdw = xpmax-xpmin;
	  ydw = ypmax-ypmin;
	  
	  // y-strips
	  if (by) {
	    
	    Int_t iX1, iX2, iY, ixDig;
	    iX1 = fBoard->GetYSix1();
	    iX2 = fBoard->GetYSix2();
	    iY  = fBoard->GetYSiy();
	    if (ix >= iX1 && ix <= iX2 && iy == iY) {
	      
	      ypmin = -0.5*fYWidth[i];
	      ypmax = +0.5*fYWidth[i];
	      /*
	      ypmin += 0.01*fYWidth[i];
	      ypmax -= 0.01*fYWidth[i];
	      xpmin += 0.1*xdw;
	      xpmax -= 0.1*xdw;
	      */
	      ixDig = ix - iX1;
	      fYDigPL[i][ixDig]->SetPoint(0,xpmin,ypmin);
	      fYDigPL[i][ixDig]->SetPoint(1,xpmax,ypmin);
	      fYDigPL[i][ixDig]->SetPoint(2,xpmax,ypmax);
	      fYDigPL[i][ixDig]->SetPoint(3,xpmin,ypmax);
	      fYDigPL[i][ixDig]->SetPoint(4,xpmin,ypmin);
	      fYDigPL[i][ixDig]->Draw();
	      /*
	      fYDigBox[i][ixDig]->SetFillStyle(1001);
	      fYDigBox[i][ixDig]->SetFillColor(5);
	      fYDigBox[i][ixDig]->DrawBox(xpmin,ypmin,xpmax,ypmax);
	      */
	      if (makeLabelsY) {
		sprintf(cln,"%2d",(ix-iX1));
		ptx1 = xpmin;
		ptx2 = xpmax;
		
		pty1 = 1.065*ypmin - 0.04*ycw;
		pty2 = 1.065*ypmin + 0.04*ycw;
		fYLabelL[i][ix-iX1] = new TPaveText(ptx1,pty1,ptx2,pty2);
		fYLabelL[i][ix-iX1]->SetBorderSize(0);
		fYLabelL[i][ix-iX1]->SetBit(kCannotPick);
		
		pty1 = 1.065*ypmax - 0.04*ycw;
		pty2 = 1.065*ypmax + 0.04*ycw;
		fYLabelR[i][ix-iX1] = new TPaveText(ptx1,pty1,ptx2,pty2);
		fYLabelR[i][ix-iX1]->SetBorderSize(0);
		fYLabelR[i][ix-iX1]->SetBit(kCannotPick);
	      }

	    }
	    
	  }
	  
	}  // end maxY
	
      }  // end maxX
      
      /*---------- x-pads ic = 0 ----------*/
      
      seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::kCath0); 

      maxX = seg->MaxPadIndexX();
      maxY = seg->MaxPadIndexY();
      
      for (Int_t ix = 0; ix <= maxX; ix++) {
	for (Int_t iy = 0; iy <= maxY; iy++) {
	  
	  pad = seg->PadByIndices(ix,iy,kFALSE);
	  
	  if (!pad.IsValid()) continue;
	  
	  // get the pad position and dimensions
	  xlocal1 = pad.Position().X();
	  ylocal1 = pad.Position().Y();
	  xlocal2 = pad.Dimensions().X();
	  ylocal2 = pad.Dimensions().Y();
	  
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
	  
	  xpmin -= fXCenter[i];
	  xpmax -= fXCenter[i];
	  ypmin -= fYCenter[i];
	  ypmax -= fYCenter[i];
	  xdw = xpmax-xpmin;
	  ydw = ypmax-ypmin;
	  
	  // x-strips
	  if (bx) {
	    
	    Int_t iX, iY1, iY2, iyDig;
	    iX  = fBoard->GetXSix();
	    iY1 = fBoard->GetXSiy1();
	    iY2 = fBoard->GetXSiy2();
	    if (ix == iX && iy >= iY1 && iy <= iY2) {
	      /*		  
	      xpmin += 0.01*fXWidth[i];
	      xpmax -= 0.01*fXWidth[i];
	      ypmin += 0.1*ydw;
	      ypmax -= 0.1*ydw;
	      */
	      iyDig = iy - iY1;
	      fXDigPL[i][iyDig]->SetPoint(0,xpmin,ypmin);
	      fXDigPL[i][iyDig]->SetPoint(1,xpmax,ypmin);
	      fXDigPL[i][iyDig]->SetPoint(2,xpmax,ypmax);
	      fXDigPL[i][iyDig]->SetPoint(3,xpmin,ypmax);
	      fXDigPL[i][iyDig]->SetPoint(4,xpmin,ypmin);
	      fXDigPL[i][iyDig]->Draw();
	      /*	      
	      fXDigBox[i][iyDig]->SetFillStyle(1001);
	      fXDigBox[i][iyDig]->SetFillColor(5);
	      fXDigBox[i][iyDig]->DrawBox(xpmin,ypmin,xpmax,ypmax);
	      */
	      if (makeLabelsX) {
		sprintf(cln,"%2d",(iy-iY1));
		pty1 = ypmin;
		pty2 = ypmax;
		
		ptx1 = 1.065*xpmin - 0.04*xcw;
		ptx2 = 1.065*xpmin + 0.04*xcw;
		fXLabelL[i][iy-iY1] = new TPaveText(ptx1,pty1,ptx2,pty2);
		fXLabelL[i][iy-iY1]->SetBorderSize(0);
		fXLabelL[i][iy-iY1]->SetBit(kCannotPick);
		
		ptx1 = 1.065*xpmax - 0.04*xcw;
		ptx2 = 1.065*xpmax + 0.04*xcw;
		fXLabelR[i][iy-iY1] = new TPaveText(ptx1,pty1,ptx2,pty2);
		fXLabelR[i][iy-iY1]->SetBorderSize(0);
		fXLabelR[i][iy-iY1]->SetBit(kCannotPick);
	      }

	    }
	    
	  }
	  
	}  // end maxY
	
      }  // end maxX
      
    }  // end DEIterator
    
    fCanvas[i]->Modified();
    fCanvas[i]->Update();

  }

  fMain->MapWindow();

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::CloseWindow() const
{
  /// close dialog in response to window manager close.

  delete this;

}

//__________________________________________________________________________
void AliMUONTriggerGUIbdmap::DoClose()
{
  /// handle Close button.

  fBoard->SetOpen(kFALSE);
  TTimer::SingleShot(150,"AliMUONTriggerGUIbdmap",this,"CloseWindow()");

}


