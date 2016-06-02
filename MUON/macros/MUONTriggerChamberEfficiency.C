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
#include "TMap.h"
#include "TSystem.h"
#include "TMath.h"
#include "TArrayI.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerChamberEfficiency.h"
#include "AliCDBManager.h"
#include "AliCDBRunRange.h"
#include "AliCDBId.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"

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


//____________________________________________________________
void BuildSystematicMap ( TString runList, Double_t globalSyst = -0.02, TString outputCDB = "local://CDB", TString inputCDB = "raw://", TString systematicPerBoard = "" )
{
  /// Build map with systematic uncertainties
  /// The file systematicPerBoard is a text file containing the relative
  /// deltas w.r.t. original efficiency in the format:
  /// chamber(11-14) board(1-234) plane(0,1) delta
  /// where 0 = bending plane and 1 = non-bending plane 

  // Read run list
  TArrayI runListArray(1000);
  Int_t nruns=0;
  if ( gSystem->AccessPathName(runList.Data()) ) {
    runList.ReplaceAll(","," ");
    TObjArray* arr = runList.Tokenize(" ");
    TObjString* objStr = 0x0;
    TIter nextRun(arr);
    while ( (objStr = static_cast<TObjString*>(nextRun())) ) {
      if ( objStr->String().IsDigit() ) runListArray[nruns++] = objStr->String().Atoi();
    }
    delete arr;
  }
  else {
    ifstream inRunList(runList.Data());
    TString currLine = "";
    while ( ! inRunList.eof() ) {
      currLine.ReadLine(inRunList);
      currLine.ReplaceAll(","," ");
      TObjArray* arr = currLine.Tokenize(" ");
      TObjString* objStr = 0x0;
      TIter nextRun(arr);
      while ( (objStr = static_cast<TObjString*>(nextRun())) ) {
        if ( objStr->String().IsDigit() ) runListArray[nruns++] = objStr->String().Atoi();
      }
      delete arr;
    }
    inRunList.close();
  }
  runListArray.Set(nruns);

  // Read maps from input CDB
  AliCDBManager* mgr = AliCDBManager::Instance();
  mgr->SetDefaultStorage(inputCDB.Data());
  TMap effMapsList;
  AliCDBId* prevId = 0x0;
  for ( Int_t irun=0; irun<runListArray.GetSize(); irun++ ) {
    mgr->SetRun(runListArray[irun]);
    AliCDBEntry* entry = mgr->Get(AliCDBPath("MUON/Calib/TriggerEfficiency"));
    const AliCDBId cdbId = entry->GetId();
    if ( prevId && cdbId.IsEqual(prevId) ) continue;
    delete prevId;
    prevId = static_cast<AliCDBId*>(cdbId.Clone());

    AliMUONTriggerEfficiencyCells* effMap = static_cast<AliMUONTriggerEfficiencyCells*>(entry->GetObject());
    effMapsList.Add(new TObjString(Form("%i_%i",cdbId.GetFirstRun(),cdbId.GetLastRun())),effMap->Clone());
  }
  delete prevId;


  // Set global systematic variation
  Int_t nPoints = 234*4*3;
  TArrayD systDiff(nPoints);
  systDiff.Reset(globalSyst);
  if ( ! systematicPerBoard.IsNull() && gSystem->AccessPathName(systematicPerBoard.Data()) == 0 ) {
    // If file with specific chamber variation exists, read values
    TString currLine = "";
    ifstream inFile(systematicPerBoard.Data());
    while ( ! inFile.eof() ) {
      currLine.ReadLine(inFile);
      TObjArray* arr = currLine.Tokenize(" ");
      if ( arr->GetEntries() >= 4 && static_cast<TObjString*>(arr->UncheckedAt(0))->String().IsDigit() ) {
        Int_t ich = static_cast<TObjString*>(arr->UncheckedAt(0))->String().Atoi()%11;
        Int_t iboard = static_cast<TObjString*>(arr->UncheckedAt(1))->String().Atoi()-1;
        Int_t icount = static_cast<TObjString*>(arr->UncheckedAt(2))->String().Atoi();
        Double_t systErr = static_cast<TObjString*>(arr->UncheckedAt(3))->String().Atof();
        Int_t idx = ich*234*3+iboard*3+icount;
        if ( idx < nPoints ) systDiff[idx] = systErr;
      }
      delete arr;
    }
    inFile.close();

    // Set both planes uncertainty
    // by default take the maximum uncertainty between bending and non-bending
    for ( Int_t ich=0; ich<4; ich++ ) {
      for ( Int_t iboard=0; iboard<234; iboard++ ) {
        Int_t idx = ich*234*3+iboard*3;
        Double_t diffBend = systDiff[idx];
        Double_t diffNonBend = systDiff[idx+1];
        systDiff[idx+2] = ( TMath::Abs(diffBend) > TMath::Abs(diffNonBend) ) ? diffBend : diffNonBend;
      }
    }
  }

  mgr->SetDefaultStorage(outputCDB.Data());

  enum { kBendingEff, kNonBendingEff, kBothPlanesEff};
  TString countTypeName[3] = {"bendPlane", "nonBendPlane","bothPlanes"};

  TIter next(&effMapsList);
  TObject* key = 0x0;
  while ( ( key = next() ) ) {
    TString runRange = key->GetName();
    AliMUONTriggerEfficiencyCells* effMap = static_cast<AliMUONTriggerEfficiencyCells*>(effMapsList.GetValue(key));
    for ( Int_t ich=0; ich<4; ich++ ) {
      TString histoName = Form("allTracksCountBoardCh%i", 11+ich);
      TH1* auxHisto = static_cast<TH1*>(effMap->GetHistoList()->FindObject(histoName));
      for ( Int_t icount=0; icount<3; icount++ ) {
        histoName = Form("%sCountBoardCh%i", countTypeName[icount].Data(), 11+ich);
        TH1* histo = static_cast<TH1*>(effMap->GetHistoList()->FindObject(histoName));
        for ( Int_t iboard=0; iboard<234; iboard++ ) {
          Int_t idx = ich*234*3+iboard*3+icount;
          Int_t ibin = iboard+1;
          Double_t countDiff = TMath::Nint(systDiff[idx] * auxHisto->GetBinContent(ibin));
          Double_t countOrig = histo->GetBinContent(ibin);
          Double_t newCount = countOrig+countDiff;
          histo->SetBinContent(ibin,newCount);
          if ( histo->GetSumw2N() > 0 ) histo->SetBinError(ibin,TMath::Sqrt(newCount));
        }
      }
    }
    TObjArray* arr = runRange.Tokenize("_");
    Int_t firstRun = static_cast<TObjString*>(arr->UncheckedAt(0))->String().Atoi();
    Int_t lastRun = static_cast<TObjString*>(arr->UncheckedAt(1))->String().Atoi();
    AliMUONCDB::WriteToCDB(effMap, "MUON/Calib/TriggerEfficiency", firstRun, lastRun, "Measured efficiencies");
  }
}

//____________________________________________________________
void CompleteEfficiency(TString effFileWithHoles, TString effFileCompatible, TString outFilename)
{
  //
  /// When a local board or RPC is missing, the efficiency of other boards cannot be calculated
  /// If an efficiency file of the same period is available, it could be used to fill the missing information
  //
  
  TList* histoList[2] = {0x0, 0x0};
  TString filenames[2] = {effFileWithHoles, effFileCompatible};
  for ( Int_t ifile=0; ifile<2; ifile++ ) {
    TFile* file = TFile::Open(filenames[ifile].Data());
    if ( ! file ) {
      printf("Fatal: cannot find %s\n", filenames[ifile].Data());
      return;
    }
    histoList[ifile] = static_cast<TList*> (file->FindObjectAny("triggerChamberEff"));
    if ( ! histoList[ifile] ) {
      printf("Cannot find histo list in file %s\n", filenames[ifile].Data());
      return;
    }
  }
  
  TString detElemName[2] = {"Slat", "Board"};
  enum { kBendingEff, kNonBendingEff, kBothPlanesEff, kAllTracks, kNcounts};
  TString countTypeName[kNcounts] = {"bendPlane", "nonBendPlane","bothPlanes", "allTracks"};
  
  Bool_t isChanged = kFALSE;
  TString histoName = "";
  for ( Int_t idet=0; idet<2; idet++ ) {
    for ( Int_t ich=11; ich<=14; ich++ ) {
      histoName = Form("%sCount%sCh%i", countTypeName[kAllTracks].Data(), detElemName[idet].Data(), ich);
      TH1* allTracksHisto = static_cast<TH1*> (histoList[0]->FindObject(histoName.Data()));
      for ( Int_t ibin=1; ibin<=allTracksHisto->GetXaxis()->GetNbins(); ibin++ ) {
        if ( allTracksHisto->GetBinContent(ibin) > 0. ) continue;
        isChanged = kTRUE;
        printf("Modifying info for Ch %i %s %3i\n", ich, detElemName[idet].Data(), (Int_t)allTracksHisto->GetXaxis()->GetBinCenter(ibin));
        // If allTracks has no entries, it means that efficiency could not be calculated for this bin:
        // fill information from the compatible histogram
        
        // Check the statistics collected by the switched off detection element
        Double_t nTracks = 0;
        for ( Int_t jch=11; jch<=14; jch++ ) {
          histoName = Form("%sCount%sCh%i", countTypeName[kAllTracks].Data(), detElemName[idet].Data(), jch);
          TH1* allTracksOtherCh = static_cast<TH1*> (histoList[0]->FindObject(histoName.Data()));
          nTracks = allTracksOtherCh->GetBinContent(ibin);
          if ( nTracks > 0. ) {
            //printf("Statistics for %s : %g\n", histoName.Data(), nTracks); // REMEMBER TO CUT
            break;
          }
        }
        
        histoName = Form("%sCount%sCh%i", countTypeName[kAllTracks].Data(), detElemName[idet].Data(), ich);
        TH1* allTracksHistoAux = static_cast<TH1*> (histoList[1]->FindObject(histoName.Data()));
        Double_t nTracksNew = allTracksHistoAux->GetBinContent(ibin);
        if ( nTracksNew == 0.) {
          printf("Warning: new histogram has no entries for Ch %i %s %3i\n", ich, detElemName[idet].Data(), (Int_t)allTracksHisto->GetXaxis()->GetBinCenter(ibin));
          continue;
        }
        Double_t scaleFactor = TMath::Min(nTracksNew, nTracks) / nTracksNew;
        //printf("Statistics ineff %g  new %g  scaleFactor %g\n", nTracks, nTracksNew, scaleFactor); // REMEMBER TO CUT
        
        for ( Int_t icount=0; icount<kNcounts; icount++ ) {
          histoName = Form("%sCount%sCh%i", countTypeName[icount].Data(), detElemName[idet].Data(), ich);
          TH1* auxHisto = static_cast<TH1*> (histoList[1]->FindObject(histoName.Data()));
          TH1* histo = static_cast<TH1*> (histoList[0]->FindObject(histoName.Data()));
          histo->SetBinContent(ibin, auxHisto->GetBinContent(ibin) * scaleFactor);
          if ( histo->GetSumw2N() > 0 ) histo->SetBinError(ibin, auxHisto->GetBinError(ibin)*scaleFactor);
        } // loop on cont types
      } // loop on histogram bins
    } // loop on chamber
  } // loop on detection element (slat or board)
  
  if ( ! isChanged ) {
    printf("Input histograms not modified\n");
    return;
  }
  TFile* outFile = TFile::Open(outFilename,"create");
  histoList[0]->Write("triggerChamberEff",TObject::kSingleKey);
  outFile->Close();  
}


