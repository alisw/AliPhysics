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

// --- MUON header files ---
#include "AliMUONTriggerDisplay.h"

#include "AliMpDDLStore.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpConstants.h"
#include "AliMpPad.h"
#include "AliMpLocalBoard.h"
#include "AliMpCDB.h"

// --- AliRoot header files ---
#include "AliLog.h"

// --- ROOT system ---
#include <TH1.h> 
#include <TH2.h>
#include <TGraph.h>

//-----------------------------------------------------------------------------
/// \class AliMUONTriggerDisplay
///
/// MUON base class for converting histos as a function of strip/board/slat number
/// into display histos showing the detection element position in the 
/// trigger chamber.
///
/// Input histos can be given as:
///  - TH2 with x -> board # [1-234]  and   y -> strip # in board [0-15].  Option: kDisplayStrips
///  - TH1 with x -> board # [1-234]                                       Option: kDisplayBoards
///  - TH1 with x -> slat #  [0-17]                                        Option: kDisplaySlats
///
/// \author D. Stocco

/// \cond CLASSIMP
ClassImp(AliMUONTriggerDisplay)
/// \endcond
           
//____________________________________________________________________________ 
AliMUONTriggerDisplay::AliMUONTriggerDisplay() : 
TObject()
{
  /// ctor
}

//__________________________________________________________________
AliMUONTriggerDisplay::~AliMUONTriggerDisplay()
{
  /// dtor
}


//____________________________________________________________________________ 
TH2* AliMUONTriggerDisplay::GetEmptyDisplayHisto(TString displayHistoName, EDisplayType displayType,
						 Int_t cathode, Int_t chamber,
						 TString displayHistoTitle)
{
  //
  /// Return the display histogram with optimized binning
  /// but do not fill it.
  //
  TH2F* displayHisto = new TH2F();
  
  InitOrDisplayTriggerInfo(0x0, displayHisto, displayType,
			   cathode, chamber, displayHistoName, displayHistoTitle);
  
  return displayHisto;
}


//____________________________________________________________________________ 
TH2* AliMUONTriggerDisplay::GetBoardNumberHisto(TString displayHistoName,
						Int_t chamber,
						TString displayHistoTitle)
{
  //
  /// Return the display histogram with optimized binning
  /// and fill it with the board number
  //

  const Int_t kNboards = AliMpConstants::NofLocalBoards();
  TH1F* inputHisto = new TH1F("boardNumbers","Board Numbers",kNboards,0.5,(Float_t)kNboards + 0.5);
  for(Int_t ibin=1; ibin<=kNboards; ibin++){
    inputHisto->Fill(ibin,ibin);
  }

  TH2F* displayHisto = (TH2F*)GetEmptyDisplayHisto(displayHistoName, kDisplayBoards, 0, chamber, displayHistoTitle);
  FillDisplayHistogram(inputHisto,displayHisto,kDisplayBoards,0,chamber,kNumbered);
  
  delete inputHisto;

  displayHisto->SetStats(kFALSE);
  return displayHisto;
}


//____________________________________________________________________________ 
TH2* AliMUONTriggerDisplay::GetDisplayHistogram(TH1* inputHisto, TString displayHistoName,
						EDisplayType displayType, Int_t cathode,
						Int_t chamber, TString displayHistoTitle,
						EDisplayOption displayOpt)
{
  //
  /// Get histogram displaying the information contained in the input histogram
  //
  TH2* displayHisto = GetEmptyDisplayHisto(displayHistoName, displayType,
					   cathode, chamber, displayHistoTitle);

  FillDisplayHistogram(inputHisto, displayHisto, displayType, cathode, chamber, displayOpt);

  return displayHisto;
}

//____________________________________________________________________________ 
TH2* AliMUONTriggerDisplay::GetDisplayHistogram(TGraph* inputGraph, TString displayHistoName,
						EDisplayType displayType, Int_t cathode,
						Int_t chamber, TString displayHistoTitle,
						EDisplayOption displayOpt)
{
  //
  /// Get histogram displaying the information contained in the input graph
  //
  TH2* displayHisto = GetEmptyDisplayHisto(displayHistoName, displayType,
					   cathode, chamber, displayHistoTitle);

  FillDisplayHistogram(inputGraph, displayHisto, displayType, cathode, chamber, displayOpt);

  return displayHisto;
}

//____________________________________________________________________________ 
Bool_t AliMUONTriggerDisplay::FillDisplayHistogram(TH1* inputHisto, TH2* displayHisto,
						   EDisplayType displayType, Int_t cathode,
						   Int_t chamber, EDisplayOption displayOpt)
{
  //
  /// Fill a previously initialized display histogram 
  /// with the information contained in inputHisto.
  /// To get initialized display, please use GetEmptyDisplayHisto method
  //
  return InitOrDisplayTriggerInfo(inputHisto, displayHisto, displayType,
				  cathode, chamber, "", "",displayOpt);
}

//____________________________________________________________________________ 
Bool_t AliMUONTriggerDisplay::FillDisplayHistogram(TGraph* inputGraph, TH2* displayHisto,
						   EDisplayType displayType, Int_t cathode,
						   Int_t chamber, EDisplayOption displayOpt)
{
  //
  /// Fill a previously initialized display histogram 
  /// with the information contained in inputGraph.
  /// To get initialized display, please use GetEmptyDisplayHisto method
  //
  return InitOrDisplayTriggerInfo(inputGraph, displayHisto, displayType,
				  cathode, chamber, "", "",displayOpt);
}

//____________________________________________________________________________ 
Bool_t AliMUONTriggerDisplay::InitOrDisplayTriggerInfo(TObject* inputObject, TH2* displayHisto,
						       EDisplayType displayType,
						       Int_t cathode, Int_t chamber,
						       TString displayHistoName, TString displayHistoTitle,
						       EDisplayOption displayOpt)
{
  //
  /// Initialize trigger information display histograms using mapping
  /// Trigger information is displayed in a user-friendly way:
  /// from local board and strip numbers to their position on chambers.
  //

  // Load mapping
  if ( ! AliMpSegmentation::Instance(kFALSE) ) {
    /// Load mapping
    if ( ! AliMpCDB::LoadDDLStore() ) {
      AliError("Could not access mapping from OCDB !");
      AliError("Histograms are not initialized !");
      return kFALSE;
    }
  }

  Int_t iCh = (chamber > AliMpConstants::NofTrackingChambers()) ? chamber - 11 : chamber;
  Int_t iCath = cathode;

  TArrayD xAxisStrip;
  TArrayD yAxisStrip;
  TArrayD xAxisBoard;
  TArrayD yAxisBoard;

  Float_t yOffsetLine, xOffsetLine = 0.;

  const Float_t kResetValue=1234567.;

  if(!inputObject){
    xAxisBoard.Set(55);
    xAxisBoard.Reset(kResetValue);
    yAxisBoard.Set(50);
    yAxisBoard.Reset(kResetValue);

    xAxisStrip.Set(420);
    xAxisStrip.Reset(kResetValue);
    yAxisStrip.Set(710);
    yAxisStrip.Reset(kResetValue);
  }
  else if(!displayHisto){
      AliWarning("Display histogram not initialized. Please initialize it first!");
      return kFALSE;
  }
  else {
    TH1* inputHisto = dynamic_cast<TH1*>(inputObject);    
    if ( inputHisto ) {
      if ( inputHisto->GetEntries() == 0 ) {
	return kTRUE;
      }
    }
    else {
      TGraph* inputGraph = dynamic_cast<TGraph*>(inputObject);
      if ( inputGraph ) {
        if ( inputGraph->GetN() == 0 ){
	  return kTRUE;
        }
      }  
      else {
        AliWarning("The object should inherit from TH1 or TGraph!");
        return kFALSE;
      }
    }
  }

  Float_t xWidth, yWidth, yWidthSlat=0., xWidthCol=0.;
  Float_t x1,x2,y1,y2;
  Float_t x1b=0., x2b=0., y1b=0., y2b=0.;
  Float_t xcPad, ycPad;
  Int_t line=0, slat;
  Float_t sign = 1.;

  const Float_t kShiftB = 0.5;
  const Float_t kShiftS = 0.1;

  const Float_t kShiftX = (iCath==0) ? kShiftB : kShiftS;
  const Float_t kShiftY = (iCath==0) ? kShiftS : kShiftB;
  const Float_t kShiftEl = (displayType==kDisplaySlats) ? 0.01 : kShiftB;

  Int_t iChamber = iCh + AliMpConstants::NofTrackingChambers();
  for(Int_t iLoc = 0; iLoc < AliMpConstants::NofLocalBoards(); iLoc++) {  
    Int_t iBoard = iLoc+1;
    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromLocalBoard(iBoard, iChamber);

    if (!detElemId) continue;

    AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(iBoard, kFALSE);

    // skip copy cards
    if( !localBoard->IsNotified()) 
      continue;

    // get segmentation
    const AliMpVSegmentation* seg[2] = {
      AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(0)),
      AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(1))};

    if(iLoc==0){
      AliMpPad pad1 = seg[1]->PadByLocation(iBoard,0,kFALSE);
      yWidthSlat = pad1.GetDimensionY();
      AliMpPad pad0 = seg[0]->PadByLocation(iBoard,0,kFALSE);
      xOffsetLine = TMath::Abs(pad0.GetPositionX()) + pad0.GetDimensionX();
      xWidthCol = 2.* pad0.GetDimensionX();
    }

    // Get ideal global position of DetElemId center
    slat = detElemId%100;
    line = (4 + slat)%18;
    sign = 1.;
    if(line>8) {
      line = 17 - line;
      sign = -1.;
    }
    yOffsetLine = (Float_t)(line - 4) * 2. * yWidthSlat;
	  
    Int_t nLocations = 1;

    for(Int_t cath=0; cath<AliMpConstants::NofCathodes(); cath++){
      // Loop on cathodes: 
      // necessary because strip info is read for each cathode
      // board info is read only from cathode 0

      // loop over strips
      for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy) {
	// get pad from electronics

	Int_t offset = 0;
	if (cath && localBoard->GetSwitch(AliMpLocalBoard::kZeroAllYLSB)) offset = -8;

	AliMpPad pad = seg[cath]->PadByLocation(iBoard,ibitxy+offset,kFALSE);

	if (!pad.IsValid()) continue;
        
  // For non-bending plane fill only the first board covered by the strip
  // i.e. avoide filling many times the same information
  if ( cath == 1 && pad.GetLocalBoardId(0) != iBoard ) continue;

	xWidth = pad.GetDimensionX();
	yWidth = pad.GetDimensionY();
	xcPad = sign * (pad.GetPositionX() + xOffsetLine);
	if(line==4) xcPad += 0.75 * sign * xWidthCol;
	ycPad = pad.GetPositionY() + yOffsetLine;
	nLocations = pad.GetNofLocations();
	Int_t iStrip = pad.GetLocalBoardChannel(0);

	if(cath==iCath){
	  x1 = xcPad - xWidth + kShiftX;
	  y1 = ycPad - yWidth + kShiftY;
	  x2 = xcPad + xWidth - kShiftX;
	  y2 = ycPad + yWidth - kShiftY;

	  if(!inputObject){
	    AddSortedPoint(x1, xAxisStrip, kResetValue);
	    AddSortedPoint(x2, xAxisStrip, kResetValue);

	    AddSortedPoint(y1, yAxisStrip, kResetValue);
	    AddSortedPoint(y2, yAxisStrip, kResetValue);
	  }
	  else if(displayType==kDisplayStrips) 
	    FillBins(inputObject, displayHisto, iBoard, iStrip, x1, x2, y1, y2, kShiftX, kShiftY, displayOpt);
	}

	if(cath==0){
	  if(iStrip==0) {
	    x1b = xcPad - xWidth + kShiftEl;
	    y1b = ycPad - yWidth + kShiftEl;
	  }
	  x2b = xcPad + xWidth - kShiftEl;
	  y2b = ycPad + yWidth - kShiftEl;
	}

      } // loop on strips

      // if iCath==0 strip info and board info are both filled -> break!
      // if iCath==1 board info is filled at cath==0. Strip info to be filled at cath==1
      if(iCath==0) break;
    } // loop on cathodes

    if(!inputObject){
      // Per board
      AddSortedPoint(x1b, xAxisBoard, kResetValue);
      AddSortedPoint(x2b, xAxisBoard, kResetValue);

      AddSortedPoint(y1b, yAxisBoard, kResetValue);
      AddSortedPoint(y2b, yAxisBoard, kResetValue);
    }
    else if(displayType==kDisplayBoards) 
      FillBins(inputObject, displayHisto, iBoard, -nLocations, x1b, x2b, y1b, y2b, kShiftEl, kShiftEl, displayOpt);
    else if(displayType==kDisplaySlats) 
      FillBins(inputObject, displayHisto, slat, -1, x1b, x2b, y1b, y2b, kShiftEl, kShiftEl, displayOpt);
  } // loop on local boards

  if ( inputObject ) return kTRUE;

  displayHisto->Reset();

  // Book histos
  const Float_t kMinDiff = 0.1;

  TArrayD* currArray[4] = {&xAxisStrip, &yAxisStrip, 
			   &xAxisBoard, &yAxisBoard};
  for(Int_t iaxis=0; iaxis<4; iaxis++){
    Int_t ipoint=0;
    while(TMath::Abs((*currArray[iaxis])[ipoint]-kResetValue)>kMinDiff) { 
      ipoint++;
    }
    if(ipoint>currArray[iaxis]->GetSize()-2) 
      AliWarning(Form("Array size (%i) lower than the number of points!", currArray[iaxis]->GetSize()));
    currArray[iaxis]->Set(ipoint);
  }

  switch(displayType){
  case kDisplayStrips:
    displayHisto->SetBins(xAxisStrip.GetSize()-1, xAxisStrip.GetArray(),
			  yAxisStrip.GetSize()-1, yAxisStrip.GetArray());
    break;
  case kDisplayBoards:
  case kDisplaySlats:
    displayHisto->SetBins(xAxisBoard.GetSize()-1, xAxisBoard.GetArray(),
			  yAxisBoard.GetSize()-1, yAxisBoard.GetArray());
    break;
  }

  displayHisto->SetName(displayHistoName.Data());
  displayHisto->SetTitle(displayHistoTitle.Data());
  displayHisto->SetXTitle("X (cm)");
  displayHisto->SetYTitle("Y (cm)");
  //displayHisto->SetStats(kFALSE);

  return kTRUE;
}


//____________________________________________________________________________ 
Bool_t AliMUONTriggerDisplay::AddSortedPoint(Float_t currVal, TArrayD& position, const Float_t kResetValue)
{
  //
  /// Add sorted point in array according to an increasing order.
  /// Used to build display histograms axis.
  //
  Int_t nEntries = position.GetSize()-1;
  Float_t tmp1, tmp2;
  const Float_t kMinDiff = 0.1;
  for(Int_t i=0; i<nEntries; i++){
    if(TMath::Abs(position[i]-currVal)<kMinDiff) return kFALSE;
    if(TMath::Abs(position[i]-kResetValue)<kMinDiff) {
      position[i] = currVal;
      return kTRUE;
    }
    if(currVal>position[i]) continue;
    tmp1 = position[i];
    position[i] = currVal;
    for(Int_t j=i+1; j<nEntries; j++){
      tmp2 = position[j];
      position[j] = tmp1;
      tmp1 = tmp2;
      if(tmp1==kResetValue) break;
    }
    return kTRUE;
  }
  return kFALSE;
}


//____________________________________________________________________________ 
void AliMUONTriggerDisplay::FillBins(TObject* inputObject, TH2* displayHisto,
				     Int_t iElement1, Int_t iElement2,
				     Float_t x1, Float_t x2, Float_t y1, Float_t y2,
				     const Float_t kShiftX, const Float_t kShiftY,
				     EDisplayOption displayOpt)
{
  //
  /// Given the bin in inputHisto, search the corresponding bins
  /// in display histo and fill it.
  //
  Int_t binY=0;
  Float_t binContent=0;
  TH1* inputHisto = dynamic_cast<TH1*>(inputObject);
  if ( inputHisto ) {
    Int_t binX = inputHisto->GetXaxis()->FindBin(iElement1);
    if ( inputObject->IsA()->InheritsFrom("TH2") ) {
      binY = inputHisto->GetYaxis()->FindBin(iElement2);
      binContent = inputHisto->GetBinContent(binX, binY);
    }
    else binContent = inputHisto->GetBinContent(binX);
  }
  else {
    TGraph* inputGraph = dynamic_cast<TGraph*>(inputObject);
    if ( inputGraph ) {
      Double_t xpt, ypt;
      for ( Int_t ipt=0; ipt<inputGraph->GetN(); ipt++){
        inputGraph->GetPoint(ipt, xpt, ypt);
        if ( TMath::Abs(xpt - iElement1) < 0.1 ) {
	  binContent = ypt;
	  break;
        }
      }
    }
    else return;
  }  

  if(binContent==0) {
    if(displayOpt==kShowZeroes) binContent = 1e-5;
    else return;
  }

  Int_t binX1 = displayHisto->GetXaxis()->FindBin(x1 + 0.01*kShiftX);
  Int_t binX2 = displayHisto->GetXaxis()->FindBin(x2 - 0.01*kShiftX);
  Int_t binY1 = displayHisto->GetYaxis()->FindBin(y1 + 0.01*kShiftY);
  Int_t binY2 = displayHisto->GetYaxis()->FindBin(y2 - 0.01*kShiftY);

  if(displayOpt==kNumbered) {
    Int_t meanBin = (binX1+binX2)/2;
    binX1 = meanBin;
    binX2 = meanBin;
    
    meanBin = (binY1+binY2)/2;
    binY1 = meanBin;
    binY2 = meanBin;
  }

  Float_t elementArea = 1.;
  if(displayOpt == kNormalizeToArea) {
    elementArea = (x2 - x1 + 2*kShiftX) * 
      (y2 - y1 + 2*kShiftY); // In InitOrDisplayTriggerInfo: 
                             // x2 = x_c + xHalfWidth - kShiftX
                             // x1 = x_c - xHalfWidth + kShiftX
                             // so x2 - x1 + 2*kShiftX returns the element width.
	  
	  // If iElement2 is less than 0, then its meaning is the 
	  // number of boards covered by one strip
	  // the area has therefore to be multiplied accordingly
	  // This fixes the problem when filling the trigger rate per boards in non-bending plane
	  // which is overestimated since the segmentation is always given by the bending-plane
	  if ( iElement2 < 0 ) elementArea *= -(Double_t)iElement2;

  }

  for(Int_t ibinx=binX1; ibinx<=binX2; ibinx++){
    for(Int_t ibiny=binY1; ibiny<=binY2; ibiny++){
      displayHisto->SetBinContent(ibinx,ibiny,binContent/elementArea);
    }
  }
}

