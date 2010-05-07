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

#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMpConstants.h"
#include "AliMUONConstants.h"

// Classes for display
#include "AliMUONTriggerDisplay.h"
#include "AliCDBManager.h"
#include "AliMpDDLStore.h"

#include "AliLog.h"

#include "TRandom.h"
#include "Riostream.h"
#include "TH1F.h"
#include "TObjArray.h"
#include "TGraphAsymmErrors.h"

#include "TH2F.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "AliMUONTriggerChamberEfficiency.h"

//-----------------------------------------------------------------------------
/// \class AliMUONTriggerChamberEfficiency
/// A class to store and give access to the trigger chamber efficiency.
///
/// Efficiency is stored per cathode on local boards
///
/// The main method of this class is IsTriggered().
///
/// \author Diego Stocco; INFN Torino
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONTriggerChamberEfficiency)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerChamberEfficiency::AliMUONTriggerChamberEfficiency(AliMUONTriggerEfficiencyCells* effCells)
:
TObject(),
fIsOwner(kFALSE),
fEfficiencyMap(effCells),
fEfficiencyObjects(0x0),
fDisplayList(0x0)
{
///  Default constructor.
  FillFromList();
}

//__________________________________________________________________________
AliMUONTriggerChamberEfficiency::AliMUONTriggerChamberEfficiency(const Char_t* filename, const Char_t* listname)
:
TObject(),
fIsOwner(kTRUE),
fEfficiencyMap(0x0),
fEfficiencyObjects(0x0),
fDisplayList(0x0)
{
///  Constructor using an ASCII file.
  fEfficiencyMap = new AliMUONTriggerEfficiencyCells(filename, listname);
  FillFromList();
}


//_____________________________________________________________________________
AliMUONTriggerChamberEfficiency::AliMUONTriggerChamberEfficiency(const AliMUONTriggerChamberEfficiency& other)
:
TObject(other),
fIsOwner(other.fIsOwner),
fEfficiencyMap(other.fEfficiencyMap),
fEfficiencyObjects(other.fEfficiencyObjects),
fDisplayList(other.fDisplayList)
{
/// Copy constructor
}

//_____________________________________________________________________________
AliMUONTriggerChamberEfficiency& AliMUONTriggerChamberEfficiency::operator=(const AliMUONTriggerChamberEfficiency& other)
{
  /// Asignment operator
  // check assignement to self
  if (this == &other)
    return *this;

  fIsOwner = other.fIsOwner;
  fEfficiencyMap = other.fEfficiencyMap;
  fEfficiencyObjects = other.fEfficiencyObjects;
  fDisplayList = other.fDisplayList;
    
  return *this;
}

//__________________________________________________________________________
AliMUONTriggerChamberEfficiency::~AliMUONTriggerChamberEfficiency()
{
///  Destructor.
  if ( fIsOwner )
    delete fEfficiencyMap;
  delete fEfficiencyObjects;
  delete fDisplayList;
}


//__________________________________________________________________________
Float_t AliMUONTriggerChamberEfficiency::GetCellEfficiency(Int_t detElemId, Int_t localBoard, Int_t hType) const
{
///  Get the efficiencies of the 2 cathodes at a given local board

  Int_t chamber = FindChamberIndex(detElemId);
  Int_t index = GetIndex(kHboardEff, hType, chamber);
  TGraphAsymmErrors* effGraph = ((TGraphAsymmErrors*)fEfficiencyObjects->At(index));

  // Some graphs are not available in the old implementation
  if ( ! effGraph ) return -1.;

  Double_t xpt, ypt;
  effGraph->GetPoint(localBoard-1, xpt, ypt);
  return ypt;
}


//__________________________________________________________________________
void 
AliMUONTriggerChamberEfficiency::IsTriggered(Int_t detElemId, Int_t localBoard, Bool_t &trigBend, Bool_t &trigNonBend) const
{
///  Whether or not a given local board has a chance to trig, on each cathode.

  // P(B) : probability to fire bending plane
  Float_t effBend = GetCellEfficiency(detElemId, localBoard, AliMUONTriggerEfficiencyCells::kBendingEff);

  // P(BN) : probability to fire bending and non-bending plane
  Float_t effBoth = GetCellEfficiency(detElemId, localBoard, AliMUONTriggerEfficiencyCells::kBothPlanesEff);

  trigBend =  ( gRandom->Rndm() > effBend ) ? kFALSE : kTRUE;

  // P(N) : probability to fire non-bending plane
  Float_t effNonBend = GetCellEfficiency(detElemId, localBoard, AliMUONTriggerEfficiencyCells::kNonBendingEff);

  if ( effBoth > 0 ) {
    effNonBend = ( trigBend ) ? 
      effBoth / effBend :   // P(N|B) = P(BN) / P(B)
      ( effNonBend - effBoth ) / ( 1. - effBend );  // P(N|!B) = ( P(N) - P(BN) ) / ( 1 - P(B) )
  }

  trigNonBend =  ( gRandom->Rndm() > effNonBend ) ? kFALSE : kTRUE;

  AliDebug(2,Form("Ch %i  board %i  resp (%i, %i)  prob (%.2f, %.2f)  effNB %.2f  effBoth %.2f\n", detElemId/100, localBoard, trigBend, trigNonBend, effBend, effNonBend, GetCellEfficiency(detElemId, localBoard, AliMUONTriggerEfficiencyCells::kNonBendingEff), effBoth));


}


//__________________________________________________________________________
Int_t AliMUONTriggerChamberEfficiency::FindChamberIndex(Int_t detElemId) const
{
///  From detElemId to chamber number

  // Int_t iChamber = AliMpDEManager::GetChamberId(detElemId);
  Int_t iChamber = detElemId/100 - 1;
  return iChamber-AliMpConstants::NofTrackingChambers();
}


//__________________________________________________________________________
void
AliMUONTriggerChamberEfficiency::FillFromList(Bool_t useMeanValues)
{
///  Fills internal histos from list.

  if ( fEfficiencyObjects )
    delete fEfficiencyObjects;

  const Int_t kNeffHistos = 
    2 * ( AliMUONTriggerEfficiencyCells::kNcounts - 1 ) * AliMUONConstants::NTriggerCh();

  fEfficiencyObjects = new TObjArray(kNeffHistos);
  fEfficiencyObjects->SetOwner();

  TH1F *histoNum = 0x0, *histoDen=0x0;
  TString histoName = "";
  Int_t deType[2] = {AliMUONTriggerEfficiencyCells::kHboardCount,
		     AliMUONTriggerEfficiencyCells::kHslatCount};
  Int_t deTypeEff[2] = {kHboardEff, kHslatEff};
  Int_t index = -1;

  Bool_t rebuildEfficiency = kTRUE;

  for ( Int_t ide=0; ide<2; ide++){
    Int_t currDe = deType[ide];

    if ( useMeanValues && currDe == AliMUONTriggerEfficiencyCells::kHboardCount ) 
      continue;

    for(Int_t ich=0; ich<AliMUONConstants::NTriggerCh(); ich++){
      histoName = fEfficiencyMap->GetHistoName(currDe, AliMUONTriggerEfficiencyCells::kAllTracks, ich);
      if ( fEfficiencyMap->GetHistoList() ) {
	histoDen = (TH1F*)fEfficiencyMap->GetHistoList()->FindObject(histoName.Data());
	if ( !histoDen ) {
	  AliWarning(Form("Histogram %s not found. Efficiency won't be re-build", histoName.Data()));
	  rebuildEfficiency = kFALSE;
	}
      }
      else {
	AliWarning("Histogram list not present: efficiency won't be re-build");
	rebuildEfficiency = kFALSE;
      }

      Int_t nTypes = ( rebuildEfficiency ) ? AliMUONTriggerEfficiencyCells::kNcounts-1 : 2; 
      for(Int_t hType=0; hType<nTypes; hType++){
	histoName = fEfficiencyMap->GetHistoName(currDe, hType, ich);

	histoNum = ( rebuildEfficiency ) ? 
	  (TH1F*)fEfficiencyMap->GetHistoList()->FindObject(histoName.Data()) :
	  fEfficiencyMap->GetOldEffHisto(currDe, ich, hType);
      
	if ( !histoNum ) {
	  AliWarning(Form("Histogram %s not found. Skip to next", histoName.Data()));
	  continue;
	}

	index = GetIndex(deTypeEff[ide], hType, ich);
	TGraphAsymmErrors* effGraph = GetEfficiencyGraph(histoNum,histoDen);
	histoName.ReplaceAll("Count","Eff");
	effGraph->SetName(histoName.Data());
	fEfficiencyObjects->AddAt(effGraph, index);
	AliDebug(5,Form("Adding object %s (%s/%s) at index %i",effGraph->GetName(),histoNum->GetName(),histoDen->GetName(),index));

	if ( useMeanValues ){
	  Int_t currChamber = ich + AliMpConstants::NofTrackingChambers();
	  histoName = fEfficiencyMap->GetHistoName(AliMUONTriggerEfficiencyCells::kHboardCount, hType, ich);
	  TH1F* auxHistoNum = (TH1F*)fEfficiencyMap->GetHistoList()->FindObject(histoName.Data())->Clone("tempHistoNum");
	  TH1F* auxHistoDen = (TH1F*)fEfficiencyMap->GetHistoList()->FindObject(histoName.Data())->Clone("tempHistoDen");
	  for ( Int_t iBinBoard = 1; iBinBoard<=AliMpConstants::NofLocalBoards(); iBinBoard++){
	    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromLocalBoard(iBinBoard, currChamber);
	    Int_t iBin = histoNum->FindBin(detElemId%100);

	    auxHistoNum->SetBinContent(iBinBoard, histoNum->GetBinContent(iBin));
	    auxHistoDen->SetBinContent(iBinBoard, histoDen->GetBinContent(iBin));
	  }
	  index = GetIndex(kHboardEff, hType, ich);
	  effGraph = GetEfficiencyGraph(auxHistoNum,auxHistoDen);
	  histoName.ReplaceAll("Count","Eff");
	  effGraph->SetName(histoName.Data());
	  fEfficiencyObjects->AddAt(effGraph, index);
	  AliDebug(5,Form("Adding object %s (%s/%s) at index %i",effGraph->GetName(),histoNum->GetName(),histoDen->GetName(),index));
	  delete auxHistoNum;
	  delete auxHistoDen;
	} // if (useMeanValues)
      } // loop on count type
    } // loop on chamber
  } // loop on detection element histogram
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEfficiency::DisplayEfficiency(Bool_t perSlat, Bool_t show2Dhisto)
{
  //
  /// Display calculated efficiency.
  //

  if ( !AliCDBManager::Instance()->GetDefaultStorage() ){
    AliWarning("Please set default CDB storage (needed for mapping).");
    return;
  }
  if ( AliCDBManager::Instance()->GetRun() < 0 ){
    AliWarning("Please set CDB run number (needed for mapping).");
    return;
  }

  TString baseCanName = "MTRtrigChEffCan";
  TString histoName;

  // Remove previously created canvases
  TCanvas* can = 0x0;
  TIter next(gROOT->GetListOfCanvases());
  while ((can = (TCanvas *)next())) {
    histoName = can->GetName();
    if ( histoName.Contains(baseCanName.Data()))
      delete can;
  }

  delete fDisplayList;
  fDisplayList = new TList();
  fDisplayList->SetOwner();

  TH2F* displayHisto = 0x0;

  AliMUONTriggerDisplay triggerDisplay;

  Int_t deType = ( perSlat ) ? kHslatEff : kHboardEff;
  AliMUONTriggerDisplay::EDisplayType displayType = ( perSlat ) ? 
    AliMUONTriggerDisplay::kDisplaySlats : AliMUONTriggerDisplay::kDisplayBoards;
  Int_t index = -1;

  TGraph* graph = 0x0;

  // Book histos
  for(Int_t ich=0; ich<AliMUONConstants::NTriggerCh(); ich++){
    Int_t currCh = 11 + ich;
    for(Int_t hType=0; hType<AliMUONTriggerEfficiencyCells::kNcounts - 1; hType++){
      index = GetIndex(deType, hType, ich);
      graph = (TGraph*)fEfficiencyObjects->At(index);
      if ( ! graph ) continue;
      histoName = graph->GetName();
      histoName += baseCanName;
      Int_t shift = 10*(index%((AliMUONTriggerEfficiencyCells::kNcounts - 1)*
			       AliMUONConstants::NTriggerCh()));
      can = new TCanvas(histoName.Data(), histoName.Data(), 100+shift, shift, 700, 700);
      can->SetRightMargin(0.14);
      can->SetLeftMargin(0.12);
      histoName.ReplaceAll(baseCanName.Data(), "Display");
      if ( show2Dhisto ) {
	displayHisto = 
	  (TH2F*)triggerDisplay.GetDisplayHistogram(graph, histoName,
						    displayType,
						    hType,currCh,histoName,
						    AliMUONTriggerDisplay::kShowZeroes);
	displayHisto->SetDirectory(0);
      }

      if ( show2Dhisto ){
	displayHisto->GetZaxis()->SetRangeUser(0.,1.);
	displayHisto->GetYaxis()->SetTitleOffset(1.4);
	displayHisto->SetStats(kFALSE);
	displayHisto->DrawCopy("COLZ");
	delete displayHisto;

	if ( deType == kHboardEff ){
	  histoName = Form("labels%iChamber%i", hType, currCh);
	  displayHisto = 
	    (TH2F*)triggerDisplay.GetBoardNumberHisto(histoName,currCh);
	  displayHisto->SetDirectory(0);
	  displayHisto->DrawCopy("textsame");
	  delete displayHisto;
	}
      }
      else {
	TGraphAsymmErrors* drawGraph = (TGraphAsymmErrors*)graph->Clone(histoName.Data());
	drawGraph->SetMarkerStyle(20);
	drawGraph->SetMarkerSize(0.7);
	drawGraph->SetMarkerColor(kRed);
	fDisplayList->Add(drawGraph);
	drawGraph->Draw("ap");
      } // loop on chamber
    } // loop on count type
  } // loop on chamber
}


//__________________________________________________________________________
Bool_t AliMUONTriggerChamberEfficiency::LowStatisticsSettings(Bool_t useMeanValues)
{
  //
  /// In case of low statistics, fill the local board efficiency with
  /// the average value of the RPC
  //

  if ( useMeanValues )
    AliInfo("Boards filled with the average efficiency of the RPC");

  FillFromList(useMeanValues);

  return kTRUE;
}


//__________________________________________________________________________
Int_t
AliMUONTriggerChamberEfficiency::GetIndex(Int_t histoType, Int_t countType, 
					  Int_t chamber) const
{
  //
  /// Return the index of the object in the array
  //

  const Int_t kNtypes = AliMUONTriggerEfficiencyCells::kNcounts - 1;
  const Int_t kNchambers = AliMUONConstants::NTriggerCh();
  return 
    histoType * kNtypes * kNchambers + 
    chamber * kNtypes +
    countType;
  //countType * kNchambers +
  //chamber;
}

//_____________________________________________________________________________
TGraphAsymmErrors* AliMUONTriggerChamberEfficiency::GetEfficiencyGraph(TH1* histoNum, TH1* histoDen)
{
  //
  /// Create the graph of efficiency from the numerator and denominator
  /// histogram in such a way to have a point set also for
  /// detection elements with efficiency = 0 or non calculated
  //

  TGraphAsymmErrors* auxGraph = 0x0;
  if ( histoDen ) auxGraph = new TGraphAsymmErrors(histoNum,histoDen);
  else auxGraph = new TGraphAsymmErrors(histoNum);

  Int_t npoints = histoNum->GetNbinsX();
  TGraphAsymmErrors* effGraph = new TGraphAsymmErrors(npoints);
  Double_t oldX, oldY;
  for ( Int_t ibin=0; ibin<npoints; ibin++ ) {
    Int_t foundPoint = -1;
    for (Int_t ipt=0; ipt<auxGraph->GetN(); ipt++) {
      auxGraph->GetPoint(ipt, oldX, oldY);
      if ( oldX > histoNum->GetBinLowEdge(ibin+1) &&
	   oldX < histoNum->GetBinLowEdge(ibin+2) ) {
	foundPoint = ipt;
	break;
      }
    }
    Double_t currX = ( foundPoint < 0 ) ? histoNum->GetBinCenter(ibin+1) : oldX; 
    Double_t currY = ( foundPoint < 0 ) ? 0. : oldY;
    Double_t eyl   = ( foundPoint < 0 ) ? 0. : auxGraph->GetErrorYlow(foundPoint);
    Double_t eyh   = ( foundPoint < 0 ) ? 0. : auxGraph->GetErrorYhigh(foundPoint);
    effGraph->SetPoint(ibin, currX, currY);
    effGraph->SetPointError(ibin, 0., 0., eyl, eyh);
  }

  delete auxGraph;

  return effGraph;
}
