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

/// \ingroup macros
/// \file CheckTriggerDCS
/// \brief Macro to check the MTR DCS aliases and
/// possibly recover missing information


#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TGraph.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGrid.h"
#include "TMap.h"
#include "TObjArray.h"
#include "THashList.h"
#include "TAxis.h"
#include "TTimeStamp.h"
#include "TObjString.h"

// STEER includes
#include "AliCDBManager.h"
#include "AliDCSValue.h"

// MUON includes
#include "AliMUONCalibrationData.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDCSNamer.h"
#include "AliMpCDB.h"
#include "AliMUONCDB.h"

#endif

enum {kDcsHV, kDcsI};

//________________________________________________
TGraph* GetGraph ( Int_t imeas, Int_t ich, Int_t irpc, TList* graphList )
{
  TString graphName = Form("%s_RPC%i_ch%i",imeas == kDcsI ? "I" : "HV", irpc, ich+1);
  TGraph* graph = static_cast<TGraph*>(graphList->FindObject(graphName.Data()));
  if ( ! graph ) {
    graph = new TGraph();
    graph->SetName(graphName.Data());
    graph->GetXaxis()->SetTitle("Time");
    graph->GetXaxis()->SetTimeDisplay(1);
    graph->GetYaxis()->SetTitle(imeas == kDcsI ? "RPC Currents (#muA)" : "RPC HV (V)");
    graphList->Add(graph);
  }
  return graph;
}


//________________________________________________
void PrintQueryOptions()
{
  printf("Missing points can be obtained in the correct format by querying the amanda server asking for the advanced options:\n");
  printf("Identification: alias; Time: Epoch\n");
  printf("Identification: alias; Time: Epoch\n");
}


//________________________________________________
Bool_t UpdateMap ( TMap* triggerDcsMap, TString recoveredPointsFilename, UInt_t startTime, UInt_t stopTime, THashList& missingAliases )
{
  /// Update the map with aliases with missing information
  if ( gSystem->AccessPathName(recoveredPointsFilename.Data()) ) {
    printf("Error: cannot open file %s\n",recoveredPointsFilename.Data());
    return kFALSE;
  }

  TString currLine = "";
  ifstream inFile(recoveredPointsFilename.Data());
  while ( ! inFile.eof() ) {
    currLine.ReadLine(inFile);
    if ( ! currLine.Contains(";") ) continue;
    TObjArray* arr = currLine.Tokenize(";");
    if ( arr->GetEntries() == 3 ) {
      Double_t timeStampFloat = static_cast<TObjString*>(arr->At(0))->String().Atof();
      UInt_t timeStamp = (UInt_t)timeStampFloat;
      TString alias = static_cast<TObjString*>(arr->At(1))->GetString();
      Float_t value = static_cast<TObjString*>(arr->At(2))->String().Atof();
      delete arr;
      if ( ! missingAliases.FindObject(alias.Data()) ) {
        printf("Warning: alias %s is not in the list of missing points\n",alias.Data());
        continue;
      }
      if ( timeStamp < startTime || timeStamp > stopTime ) {
        printf("Warning: the point %s %u is not in the time range (%u, %u)\n",alias.Data(),timeStamp,startTime,stopTime);
        continue;
      }
      TObjArray* values = static_cast<TObjArray*>(triggerDcsMap->GetValue(alias.Data()));
      if ( ! values ) {
        values = new TObjArray();
        triggerDcsMap->Add(new TObjString(alias),values);
      }
      AliDCSValue* val = new AliDCSValue(value,timeStamp);
      values->Add(val);
    }
    else {
      printf("Warning: ill formed point: %s\n",currLine.Data());
      PrintQueryOptions();
      delete arr;
    }
  }
  inFile.close();

  return kTRUE;
}


//________________________________________________
void CompareDCSvalues ( Int_t runNumber, const char* cdb1, const char* cdb2 )
{
  /// Compare trigger DCS Points
  TString cdb[2] = {cdb1, cdb2};
  TMap* maps[2] = {0x0,0x0};
  AliCDBManager* mgr = AliCDBManager::Instance();
  for ( Int_t icdb=0; icdb<2; icdb++ ) {
    mgr->SetDefaultStorage(cdb[icdb]);
    mgr->SetRun(runNumber);
    AliMUONCalibrationData calibrationData(runNumber);
    TMap* triggerDcsMap = calibrationData.TriggerDCS();
    if ( ! triggerDcsMap ) {
      printf("Error: cannot find TriggerDCS in %s",cdb[icdb].Data());
      return;
    }
    maps[icdb] = static_cast<TMap*>(triggerDcsMap->Clone());
  }

  printf("\n\n Checking for differences between %s and %s :\n\n",cdb[0].Data(),cdb[1].Data());

  TIter next(maps[0]);
  TObjString* key = 0x0;
  while ( (key =  static_cast<TObjString*>(next())) ) {
    TObjArray* values1 = static_cast<TObjArray*>(maps[0]->GetValue(key->GetName()));
    TObjArray* values2 = static_cast<TObjArray*>(maps[1]->GetValue(key->GetName()));
    if ( ! values2 ) {
      printf("Alias %s only in %s\n",key->GetName(),cdb[0].Data());
      continue;
    }
    if ( values1->GetEntriesFast() == values2->GetEntriesFast() ) {
      for ( Int_t ival=0; ival<values1->GetEntries(); ival++ ) {
        AliDCSValue* val1 = static_cast<AliDCSValue*>(values1->UncheckedAt(ival));
        AliDCSValue* val2 = static_cast<AliDCSValue*>(values2->UncheckedAt(ival));
        if ( val1->GetTimeStamp() != val2->GetTimeStamp() ) {
          printf("Alias %s: value %i. Timestamps %u != %u\n",key->GetName(),ival,val1->GetTimeStamp(), val2->GetTimeStamp());
        }
        if ( val1->GetFloat() != val2->GetFloat() ) {
          printf("Alias %s: value %i. Float %g != %g\n",key->GetName(),ival,val1->GetFloat(), val2->GetFloat());
        }
      }
    }
    else {
      printf("Alias %s. Number of points: %i != %i\n",key->GetName(),values1->GetEntriesFast(),values2->GetEntriesFast());
    }
    maps[1]->DeleteEntry(static_cast<TPair*>(maps[1]->FindObject(key->GetName()))->Key());
  }

  TIter next1(maps[1]);
  while ( (key =  static_cast<TObjString*>(next1())) ) {
    printf("Alias %s only in %s\n",key->GetName(),cdb[1].Data());
    maps[1]->GetValue(key->GetName())->Print();
  }

  for ( Int_t icdb=0; icdb<2; icdb++ ) delete maps[icdb
                                                   ];
}



//________________________________________________
void CheckTriggerDCS ( Int_t runNumber, TString ocdbPath="raw://", TString recoveredPointsFilename = "", TString outDirRecovery = "CDB" )
{
  /// Get HV and currents values for one trigger chamber
  
  THashList graphList;
  
  if (ocdbPath.Contains("alien")) TGrid::Connect("alien://");

  AliCDBManager* mgr = AliCDBManager::Instance();
  mgr->SetDefaultStorage(ocdbPath.Data());
  mgr->SetRun(runNumber);
  mgr->SetSpecificStorage("MUON/Calib/MappingData","local://$ALICE_ROOT/OCDB"); // SAVE TIME
  
  AliMpCDB::LoadDDLStore();
  
  AliMUONCalibrationData calibrationData(AliCDBManager::Instance()->GetRun());
  
  TMap* triggerDcsMap = calibrationData.TriggerDCS();
  
  if ( ! triggerDcsMap ) {
    printf("Cannot fill DCS histos, as triggerDcsMap is NULL\n");
    return;
  }
  
  AliMpDCSNamer triggerDcsNamer("TRIGGER");

  UInt_t minTimeStamp = 0xFFFFFFFF;

  UInt_t maxTimeStamp = 0;

  THashList missingPoints;
  missingPoints.SetOwner();
  
  AliMpDEIterator deIt;
  deIt.First();
  while ( !deIt.IsDone() )
  {
    Int_t detElemId = deIt.CurrentDEId();
    
    if ( AliMpDEManager::GetStationType(detElemId) == AliMp::kStationTrigger) {
      
      Int_t ich = AliMpDEManager::GetChamberId(detElemId);
      Int_t irpc = detElemId%100;

      for ( Int_t imeas=0; imeas<AliMpDCSNamer::kNDCSMeas; imeas++ ){
        TString currAlias = triggerDcsNamer.DCSAliasName(detElemId, 0, imeas);

        TObjArray* values = static_cast<TObjArray*>(triggerDcsMap->GetValue(currAlias.Data()));
        if ( values ) {
          TIter next(values);
          AliDCSValue* val = 0x0;
          Int_t ipoint = 0;
          while ( ( val = static_cast<AliDCSValue*>(next()) ) ) {
            Float_t hvi = val->GetFloat();
            UInt_t timeStamp = val->GetTimeStamp();
            if ( timeStamp < minTimeStamp ) minTimeStamp = timeStamp;
            if ( timeStamp > maxTimeStamp ) maxTimeStamp = timeStamp;
            TGraph* graph = GetGraph(imeas,ich,irpc,&graphList);
            if ( ipoint == 0 ) {
              TTimeStamp ts(timeStamp);
              graph->SetTitle(Form("Date: %u (yyyymmdd)\n",ts.GetDate()));
            }
            graph->SetPoint(ipoint++,timeStamp,hvi);
          } // loop on values
        }
        else {
          printf("Did not find expected alias (%s) for DE %d\n",currAlias.Data(),detElemId);
          missingPoints.Add(new TObjString(currAlias));
        }
      } // loop on measured types (HV and currents)
    } // if (stationType == kStationTrigger)
    
    deIt.Next();
  }

  if ( missingPoints.GetEntries() > 0 ) {
    if ( recoveredPointsFilename.IsNull() ) {
      printf("\n\nNo values found for the following aliases:\n");
      TIter next(&missingPoints);
      TObject* obj = 0x0;
      while ( (obj = next() ) ) printf("%s\n",obj->GetName());
      printf("\n");
      PrintQueryOptions();
      printf("\nIn the time range\n");
      UInt_t year, month, day, hour, min, sec;
      for ( Int_t itime=0; itime<2; itime++) {
        UInt_t timeStamp = ( itime == 0 ) ? minTimeStamp : maxTimeStamp;
        TTimeStamp ts(timeStamp);
        ts.GetDate(kTRUE,0,&year,&month,&day);
        ts.GetTime(kTRUE,0,&hour,&min,&sec);
        TString usTime = "AM";
        if ( hour >= 12 ) {
          hour = hour - 12;
          usTime = "PM";
        }
        if ( hour == 0 ) hour = 12;
        ts.Print();
        printf("  Amanda format: %i/%i/%i %i:%i:%i %s\n",month,day,year,hour,min,sec,usTime.Data());
      }
    }
    else {
      printf("\nRecover missing points with data in %s\n",recoveredPointsFilename.Data());

      TMap* clonedTriggerDCSmap = static_cast<TMap*>(triggerDcsMap->Clone());

      if ( UpdateMap(clonedTriggerDCSmap, recoveredPointsFilename, minTimeStamp, maxTimeStamp, missingPoints) ) {
        TString specStorageDir = "MUON/Calib/TriggerDCS";
        TString fullDir = Form("%s/%s",outDirRecovery.Data(),specStorageDir.Data());
        if ( gSystem->AccessPathName(fullDir) ) gSystem->Exec(Form("mkdir -p %s",fullDir.Data()));
        mgr->UnloadFromCache(specStorageDir.Data());
        mgr->SetSpecificStorage(specStorageDir.Data(),Form("local://%s",outDirRecovery.Data()));
        AliMUONCDB::WriteToCDB(clonedTriggerDCSmap, specStorageDir.Data(), runNumber, runNumber, "Recovered MTR DCS values");
      }
    }
  }

  
  for ( Int_t imeas=0; imeas<2; imeas++ ){
    for ( Int_t ich=10; ich<14; ich++ ){
      TString canName = Form("%s_ch%i",imeas == kDcsI ? "I" : "HV", ich+1);
      TCanvas* can = new TCanvas(canName.Data(), canName.Data(), 200, 10, 600, 600);
      TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
      leg->SetBorderSize(1);
      Int_t icolor=0;
      for ( Int_t irpc = 0; irpc<18; irpc++ ) {
        TGraph* graph = GetGraph(imeas,ich,irpc,&graphList);
        if ( graph->GetN() == 0 ) {
//          printf("No points found for %s\n",graph->GetName());
          continue;
        }
        TString drawOpt = "lp";
        if ( can->GetListOfPrimitives()->GetEntries() == 0 ) drawOpt.Append("a");
        if ( irpc<9 ) icolor = irpc+1;
        else icolor = (irpc-8)*10+1;
//        icolor++;
//        if ( icolor == 10 ) icolor++;
        graph->SetLineColor(icolor);
        graph->SetMarkerColor(icolor);
        graph->SetMarkerStyle(20+irpc%10);
        if ( imeas == kDcsHV ) graph->GetYaxis()->SetRangeUser(9500.,11000.);
        else graph->GetYaxis()->SetRangeUser(0,10);
        graph->GetXaxis()->SetTimeDisplay(1);
        graph->Draw(drawOpt.Data());
        leg->AddEntry(graph,Form("RPC %i",irpc),"lp");
      }
      leg->Draw();
    }
  }
}
