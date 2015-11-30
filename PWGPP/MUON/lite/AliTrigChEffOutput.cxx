/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliTrigChEffOutput.h"

// ROOT includes
#include "TObjString.h"
#include "TObjArray.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"

/// \cond CLASSIMP
ClassImp(AliTrigChEffOutput) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliTrigChEffOutput::AliTrigChEffOutput ( TObjArray* outputList, const char* name ) :
AliMuonAnalysisOutput(outputList,name),
  fTrackSelKeys(0x0),
  fCountTypeKeys(0x0),
  fHistoTypeKeys(0x0),
  fEffMethodKeys(0x0),
  fMatchTrigKeys(0x0)
{
  /// Ctor.
  InitKeys();
}

//________________________________________________________________________
AliTrigChEffOutput::AliTrigChEffOutput ( const char *filename, const char* outputName ) :
  AliMuonAnalysisOutput(filename,outputName),
  fTrackSelKeys(0x0),
  fCountTypeKeys(0x0),
  fHistoTypeKeys(0x0),
  fEffMethodKeys(0x0),
  fMatchTrigKeys(0x0)
{
  /// Ctor.
  InitKeys();
}

//________________________________________________________________________
AliTrigChEffOutput::~AliTrigChEffOutput()
{
  //
  /// Destructor
  //
  delete fTrackSelKeys;
  delete fCountTypeKeys;
  delete fHistoTypeKeys;
  delete fEffMethodKeys;
  delete fMatchTrigKeys;
}

//________________________________________________________________________
void AliTrigChEffOutput::InitKeys()
{
  /// Initialize keys
  TString matchTrigNames = "Nopt Apt Lpt Hpt";
  fMatchTrigKeys = matchTrigNames.Tokenize(" ");

  TString countTypeNames = "bendPlaneCount nonBendPlaneCount bothPlanesCount allTracksCount";
  fCountTypeKeys = countTypeNames.Tokenize(" ");

  TString histoTypeKeys = "Chamber Slat Board checkRejectedBoard";
  fHistoTypeKeys = histoTypeKeys.Tokenize(" ");

  TString effMethodKeys = "FromTrk FromTrg";
  fEffMethodKeys = effMethodKeys.Tokenize(" ");

  TString trackSelNames = "Match NoSelMatch";
  fTrackSelKeys = trackSelNames.Tokenize(" ");
}

//___________________________________________________________________________
TList* AliTrigChEffOutput::GetEffHistoList ( TString physSel, TString trigClassNames, TString centrality, Int_t itrackSel, Int_t imatch, Int_t imethod )
{
  /// Get the list of objects for the efficiency calculation
  /// merging the splitted output of the fMergeableCollection
  /// The obtained list can be converted in the efficiency map used in simulations
  /// in a backward compatible way
  
  if ( ! GetMergeableCollection() ) return 0x0;
  TList* outList = new TList();
  outList->SetOwner();
  TString histoName = "";
  TH1* histo = 0x0;
  Bool_t isOk = kTRUE;
  for ( Int_t icount=0; icount<kNcounts; ++icount ) {
    histoName = GetHistoName(kHchamberEff, icount, -1, itrackSel, imatch, imethod);
    histo = static_cast<TH1*>(GetSum(physSel, trigClassNames, centrality, histoName));
    if ( histo ) {
      histoName = GetHistoName(kHchamberEff, icount, -1, -1, -1, -1);
      histo->SetName(histoName.Data());
      histo->SetTitle(histoName.Data());
    }
    else {
      histo = GetCountHisto(kHchamberEff, icount, -1, -1, -1, -1);
      isOk = kFALSE;
    }
    histo->SetDirectory(0);
    outList->Add(histo);
  }
  for ( Int_t icount=0; icount<kNcounts; ++icount ) {
    for ( Int_t ich=0; ich<4; ++ich ) {
      histoName = GetHistoName(kHslatEff, icount, ich, itrackSel, imatch, imethod);
      histo = static_cast<TH1*>(GetSum(physSel, trigClassNames, centrality, histoName));
      if ( histo ) {
        histoName = GetHistoName(kHslatEff, icount, ich, -1, -1, -1);
        histo->SetName(histoName.Data());
        histo->SetTitle(histoName.Data());
      }
      else {
        histo = GetCountHisto(kHslatEff, icount, ich, -1, -1, -1);
        isOk = kFALSE;
      }
      histo->SetDirectory(0);
      outList->Add(histo);
    }
  }
  for ( Int_t icount=0; icount<kNcounts; ++icount ) {
    for ( Int_t ich=0; ich<4; ++ich ) {
      histoName = GetHistoName(kHboardEff, icount, ich, itrackSel, imatch, imethod);
      histo = static_cast<TH1*>(GetSum(physSel, trigClassNames, centrality, histoName));
      if ( histo ) {
        histoName = GetHistoName(kHboardEff, icount, ich, -1, -1, -1);
        histo->SetName(histoName.Data());
        histo->SetTitle(histoName.Data());
      }
      else {
        histo = GetCountHisto(kHboardEff, icount, ich, -1, -1, -1);
        isOk = kFALSE;
      }
      histo->SetDirectory(0);
      outList->Add(histo);
    }
  }
  
  histoName = GetHistoName(kHcheckBoard, -1, -1, itrackSel, imatch, imethod);
  histo = static_cast<TH1*>(GetSum(physSel, trigClassNames, centrality, histoName));
  if ( histo ) {
    histoName = GetHistoName(kHcheckBoard, -1, -1, -1, -1, -1);
    histo->SetName(histoName.Data());
    histo->SetTitle(histoName.Data());
  }
  else {
    histo = GetCountHisto(kHcheckBoard, -1, -1, -1, -1, -1);
  }
  histo->SetDirectory(0);
  outList->Add(histo);
  
  return outList;
}

//___________________________________________________________________________
TH1* AliTrigChEffOutput::GetCountHisto ( Int_t itype, Int_t icount, Int_t ichamber, Int_t itrackSel, Int_t imatch, Int_t imethod )
{
  //
  /// Get histogram with counts for efficiency calculation
  //

  Int_t nBoardBins = 234;
  Float_t boardLow = 1.-0.5, boardHigh = (Float_t)nBoardBins+1.-0.5;
  const Char_t* boardName = "board";

  TString histoName = "";
  TH1* histo = 0x0;
  switch ( itype ) {
    case kHchamberEff:
      histoName = GetHistoName(kHchamberEff, icount, -1, itrackSel, imatch, imethod);
      histo = new TH1F(histoName, histoName, 4, 11.-0.5, 4.+11.-0.5);
      histo->GetXaxis()->SetTitle("chamber");
      histo->GetYaxis()->SetTitle("counts");
      break;
    case kHslatEff:
      histoName = GetHistoName(kHslatEff, icount, ichamber, itrackSel, imatch, imethod);
      histo = new TH1F(histoName, histoName, 18, 0.-0.5, 18.-0.5);
      histo->GetXaxis()->SetTitle("slat");
      histo->GetYaxis()->SetTitle("counts");
      break;
    case kHboardEff:
      histoName = GetHistoName(kHboardEff, icount, ichamber, itrackSel, imatch, imethod);
      histo = new TH1F(histoName, histoName, nBoardBins, boardLow, boardHigh);
      histo->GetXaxis()->SetTitle(boardName);
      histo->GetYaxis()->SetTitle("counts");
      break;
    case kHcheckBoard:
      histoName = GetHistoName(kHcheckBoard, -1, -1, itrackSel, imatch, imethod);
      histo = new TH2F(histoName.Data(), "Rejected tracks motivation", 5, 20.5, 25.5, nBoardBins, boardLow, boardHigh);
      histo->GetXaxis()->SetBinLabel(1,"Many pads");
      histo->GetXaxis()->SetBinLabel(2,"Few pads");
      histo->GetXaxis()->SetBinLabel(3,"Outside geom");
      histo->GetXaxis()->SetBinLabel(4,"Tracker track");
      histo->GetXaxis()->SetBinLabel(5,"Masked board");
      histo->GetYaxis()->SetTitle(boardName);
      break;
    default:
      return 0x0;
  }

  return histo;
}

//___________________________________________________________________________
TString AliTrigChEffOutput::GetHistoName(Int_t itype, Int_t icount, Int_t ichamber, Int_t itrackSel, Int_t imatch, Int_t imethod)
{
  /// Get histogram index
  TString histoName = "";
  if ( itype < kHcheckBoard && icount >= 0 ) histoName += static_cast<TObjString*>(fCountTypeKeys->At(icount))->GetString();
  if ( itype >= 0 ) histoName += static_cast<TObjString*>(fHistoTypeKeys->At(itype))->String();
  if ( itype != kHchamberEff && ichamber >= 0 ) histoName += Form("Ch%i", 11+ichamber);
  if ( itrackSel >= 0 ) histoName += static_cast<TObjString*>(fTrackSelKeys->At(itrackSel))->String();
  if ( imatch >= 0 ) histoName += static_cast<TObjString*>(fMatchTrigKeys->At(imatch))->String();
  if ( imethod >= 0 ) histoName += static_cast<TObjString*>(fEffMethodKeys->At(imethod))->String();
  return histoName;
}
