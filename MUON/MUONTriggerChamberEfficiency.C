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

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "Riostream.h"

// ROOT includes
#include "TGrid.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerChamberEfficiency.h"
#include "AliCDBManager.h"
#include "AliCDBRunRange.h"

#endif

/// \ingroup macros
/// \file MUONTriggerChamberEfficiency.C
/// \brief Macro to view and save the trigger chamber efficiency map 
/// calculated during reconstruction.
///
/// Efficiency map can be made available for next simulation.
///
/// \author Diego Stocco, Subatech, Nantes

void MUONTriggerChamberEfficiency(TString inputFile = "./MUON.TriggerEfficiencyMap.root",
				  TString outputCDB = "",
				  Int_t firstRun=0, Int_t lastRun = AliCDBRunRange::Infinity()
)
{
/// \param inputFile (default "./MUON.TriggerEfficiencyMaps.root")
///     File with the numerator and denominator histos for efficiency calculation
///     (It is the output of the PWG3/muon/AliAnalysisTaskTrigChEff analysis
/// \param outputCDB (default "")
///     add the map on the specified CDB
/// \param firstRun (default 0)
///     first run of validity for CDB object
/// \param lastRun (default AliCDBRunRange::Infinity())
///     last run of validity for CDB Object

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  AliMUONTriggerEfficiencyCells* effMap = new AliMUONTriggerEfficiencyCells(inputFile.Data());

  if ( outputCDB.IsNull() ){
    // Draw the efficiency and exit
    AliCDBManager::Instance()->SetRun(firstRun);
    AliMUONTriggerChamberEfficiency* trigChEff = new AliMUONTriggerChamberEfficiency(effMap);
 
    trigChEff->DisplayEfficiency(kFALSE,kFALSE);
    return;
  }


  // Write efficiency on OCDB

  AliCDBManager::Instance()->SetSpecificStorage("MUON/Calib/TriggerEfficiency", outputCDB.Data());
  
  AliMUONCDB::WriteToCDB(effMap, "MUON/Calib/TriggerEfficiency", firstRun, lastRun, "Measured efficiencies");
}

//____________________________________________________________
void ShowOCDBmap(Int_t runNumber = 0, TString specificCDB="", TString ocdbPath = "local://$ALICE_ROOT/OCDB", TString runType="Full")
{
/// \param runNumber (default 0)
///     run number
/// \param specificCDB (default "")
///     specific CDB for trigger efficiency
/// \param ocdbPath(default "local://$ALICE_ROOT/OCDB")
///     path to OCDB
  if ( ocdbPath.BeginsWith("alien://") || ocdbPath.BeginsWith("raw://"))
    TGrid::Connect("alien://");

  if (!ocdbPath.CompareTo("MC"))
    AliCDBManager::Instance()->SetDefaultStorage(ocdbPath.Data(),runType.Data());
  else
    AliCDBManager::Instance()->SetDefaultStorage(ocdbPath.Data());
  
  if ( !specificCDB.IsNull() )
    AliCDBManager::Instance()->SetSpecificStorage("MUON/Calib/TriggerEfficiency", specificCDB.Data());
  AliCDBManager::Instance()->SetRun(runNumber);
  AliMUONCalibrationData calib(runNumber);

  AliMUONTriggerChamberEfficiency* trigChEff = new AliMUONTriggerChamberEfficiency(calib.TriggerEfficiency());
  trigChEff->DisplayEfficiency(kFALSE,kFALSE);
}


//____________________________________________________________
void FillHisto(TH1* histo, Int_t nevents)
{
  /// Fill histogram with global value
  for ( Int_t ibin=1; ibin<=histo->GetXaxis()->GetNbins(); ++ibin ) {
    Double_t binCenter = histo->GetXaxis()->GetBinCenter(ibin);
    for ( Int_t ievent=0; ievent<nevents; ++ievent ) {
      histo->Fill(binCenter);
    }
  }
}

//____________________________________________________________
void BuildDefaultMap(TString outFilename="/tmp/defTrigChEff.root", Double_t globalValue = 1., Int_t nevents = 100000)
{
  /// Build default map (all boards with the same chosen value)
  
  // Create histograms
  enum { kBendingEff, kNonBendingEff, kBothPlanesEff, kAllTracks, kNcounts};
  TString countTypeName[kNcounts] = {"bendPlane", "nonBendPlane","bothPlanes", "allTracks"};

  const Char_t* yAxisTitle = "counts";

  const Int_t kNboards = 234; //AliMpConstants::NofLocalBoards();
  const Int_t kFirstTrigCh = 11;//AliMpConstants::NofTrackingChambers()+1;
  const Int_t kNchambers = 4;
  const Int_t kNslats = 18;

  Int_t chamberBins = kNchambers;
  Float_t chamberLow = kFirstTrigCh-0.5, chamberHigh = kFirstTrigCh+kNchambers-0.5;
  const Char_t* chamberName = "chamber";

  Int_t slatBins = kNslats;
  Float_t slatLow = 0-0.5, slatHigh = kNslats-0.5;
  const Char_t* slatName = "slat";

  Int_t boardBins = kNboards;
  Float_t boardLow = 1-0.5, boardHigh = kNboards+1.-0.5;
  const Char_t* boardName = "board";

  TString baseName, histoName, histoTitle;
  TList* histoList = new TList();
  histoList->SetOwner();
  
  TH1F* histo;

  for(Int_t icount=0; icount<kNcounts; icount++){
    histoName = Form("%sCountChamber", countTypeName[icount].Data());
    histo = new TH1F(histoName, histoName,
                     chamberBins, chamberLow, chamberHigh);
    histo->GetXaxis()->SetTitle(chamberName);
    histo->GetYaxis()->SetTitle(yAxisTitle);
    Double_t nfills = ( icount == kAllTracks ) ? nevents : globalValue * (Double_t)nevents;
    FillHisto(histo, (Int_t)nfills);
    histoList->AddLast(histo);
  } // loop on counts

  for(Int_t icount=0; icount<kNcounts; icount++){
    for(Int_t ch=0; ch<kNchambers; ch++){
      histoName = Form("%sCountSlatCh%i", countTypeName[icount].Data(), kFirstTrigCh+ch);
      histo = new TH1F(histoName, histoName,
                       slatBins, slatLow, slatHigh);
      histo->GetXaxis()->SetTitle(slatName);
      histo->GetYaxis()->SetTitle(yAxisTitle);
      Double_t nfills = ( icount == kAllTracks ) ? nevents : globalValue * (Double_t)nevents;
      FillHisto(histo, (Int_t)nfills);
      histoList->AddLast(histo);
    } // loop on chamber
  } // loop on counts

  for(Int_t icount=0; icount<kNcounts; icount++){
    for(Int_t ch=0; ch<kNchambers; ch++){
      histoName = Form("%sCountBoardCh%i", countTypeName[icount].Data(), kFirstTrigCh+ch);
      histo = new TH1F(histoName, histoName,
                       boardBins, boardLow, boardHigh);
      histo->GetXaxis()->SetTitle(boardName);
      histo->GetYaxis()->SetTitle(yAxisTitle);
      Double_t nfills = ( icount == kAllTracks ) ? nevents : globalValue * (Double_t)nevents;
      FillHisto(histo, (Int_t)nfills);
      histoList->AddLast(histo);
    } // loop on chamber
  } // loop on counts

  TFile* outFile = TFile::Open(outFilename,"create");
  histoList->Write("triggerChamberEff",TObject::kSingleKey);
  outFile->Close();
}
