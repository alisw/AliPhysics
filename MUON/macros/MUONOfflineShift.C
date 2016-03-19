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
/// \file MUONOfflineShift.C
/// \brief Macro to be used to check raw data during MUON offline shifts.
///  
/// You NEED an access to the Grid to run this macro.
///
/// Basic usage is : 
///
/// MUONOfflineShift("path_to_raw_file","basename of output file"); > log
///
/// (the redirection to an output log file is recommended as the output from
/// this macro might be quite long...)
///
/// This will read the raw data and process it several times, varying what's done :
/// only decoding, decoding + zero-suppression, etc... (TBE)
///
/// Two outputs files will be created : 
///
/// - basename.root, containing AliMUONVTrackerData objects that can then 
///   be displayed using the mchview program (using mchview --use basename.root)
///
/// - basename.log, containing (for the moment) the occupancy numbers of the various detection elements
///
/// Unless you really know what you're doing, please use this macro together with 
/// the grid OCDB (i.e. don't set it to a local OCDB). There are now a lot of things
/// that are grabbed from the OCDB that are run dependent...
///
/// \author Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliMUONCDB.h"
#include "AliMUONPainterDataRegistry.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTrackerDataMaker.h"
#include "AliMUONVTrackerData.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDEIterator.h"
#include "AliRawReader.h"
#include <Riostream.h>
#include <TFile.h>
#include <TGrid.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TSystem.h>

#endif

//______________________________________________________________________________
Int_t DataMakerReading(const char* input,
                       TStopwatch& timer,
                       const char* cdbPath="",
                       Bool_t histogram=kFALSE,
                       Double_t xmin = 0.0,
                       Double_t xmax = 4096.0)
{
  /// Run over the data and calibrate it if so required (if cdbPath != "")
  
  TString fileName(gSystem->ExpandPathName(input));
  
  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  
  if (!rawReader) return 0;
  
  AliMUONVTrackerDataMaker* dm(0x0);
  
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return 0;
  
  if ( strlen(cdbPath) > 0 ) 
  {
    dm = new AliMUONTrackerDataMaker(recoParam,rawReader,cdbPath,histogram,xmin,xmax);
  }
  else  
  {
    dm = new AliMUONTrackerDataMaker(rawReader,histogram);
  }
  
  AliMUONPainterDataRegistry::Instance()->Register(dm);

  timer.Start(kTRUE);
  Int_t n(0);
  
  dm->SetRunning(kTRUE);
  
  while (dm->NextEvent())
  {
    ++n;
  }
  
  timer.Stop();
  
  return n;
}

//______________________________________________________________________________
void Print(const char* method, TStopwatch& timer, Int_t n)
{
  /// Get the timing for a given method
  
  cout << Form("%20s %10d events. Total CPU time %7.2f seconds.",
               method,n,timer.CpuTime());
  
  Double_t cpu = timer.CpuTime()/n;
  Double_t real = timer.RealTime()/n;
  
  cout << Form(" ms real/event = %7.2f ms CPU/event = %7.2f",real*1E3,cpu*1E3)
    << endl;
}

//______________________________________________________________________________
void Occupancy(ostream& outfile)
{
  /// Write occupancy numbers to output text file
  
  outfile << "-----------------------------------------------------" << endl;
  outfile << "Occupancy numbers" << endl;
  outfile << "-----------------------------------------------------" << endl;
  
  const Int_t occIndex = 2;
  
  AliMUONPainterDataRegistry* reg = AliMUONPainterDataRegistry::Instance();

  Int_t nofDataSources = reg->NumberOfDataSources();

  outfile << Form("%11s|"," ");
  
  for ( Int_t ids = 0; ids < nofDataSources; ++ ids ) 
  {
    AliMUONVTrackerData* data = reg->DataSource(ids);
    outfile << Form(" %13s |",data->GetName());
  }

  outfile << endl;

  for ( Int_t chamberId = 0; chamberId < AliMpConstants::NofTrackingChambers(); ++chamberId ) 
  {
    Bool_t nonZero(kFALSE);
    for ( Int_t ids = 0; ids < nofDataSources && nonZero == kFALSE; ++ ids ) 
    {
      if ( reg->DataSource(ids)->Chamber(chamberId,occIndex) ) nonZero = kTRUE;
    }
    
    if ( !nonZero ) continue;
    
    outfile << Form("Chamber %2d |",chamberId);
    for ( Int_t ids = 0; ids < nofDataSources; ++ ids ) 
    {
      AliMUONVTrackerData* data = reg->DataSource(ids);
      outfile << Form("   %7.2f %%   |",100.0*data->Chamber(chamberId,occIndex));
    }
    outfile << endl;
    
    AliMpDEIterator it;
    it.First(chamberId);
    while (!it.IsDone())
    {
      Int_t detElemId = it.CurrentDEId();
      Bool_t nonZero(kFALSE);
      for ( Int_t ids = 0; ids < nofDataSources && nonZero == kFALSE; ++ ids ) 
      {
        AliMUONVTrackerData* data = reg->DataSource(ids);
        if ( data->DetectionElement(detElemId,occIndex) > 0 ) 
        {
          nonZero = kTRUE;
        }
      }
      
      if ( nonZero ) 
      {
        outfile << Form("   DE %04d |",detElemId);
        for ( Int_t ids = 0; ids < nofDataSources; ++ ids ) 
        {
          AliMUONVTrackerData* data = reg->DataSource(ids);
          outfile << Form("   %7.2f %%   |",100.0*data->DetectionElement(detElemId,occIndex));
        }
        outfile << endl;
      }
      it.Next();
    }
  }
}

//______________________________________________________________________________
void MUONOfflineShift(const char* input="alien:///alice/data/2009/LHC09a/000067495/raw/09000067495031.10.root", 
                      const char* outputBase="67495031.10",
                      const char* ocdbPath="alien://folder=/alice/data/2009/OCDB")
{
  /// Entry point of the macro. 

  AliCodeTimer::Instance()->Reset();
  
  TGrid::Connect("alien://");
  
  AliRawReader* rawReader = AliRawReader::Create(input);
  
  rawReader->NextEvent();
  
  Int_t runNumber = rawReader->GetRunNumber();
    
  delete rawReader;
  
  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
  AliCDBManager::Instance()->SetRun(runNumber);

  if (!AliMUONCDB::LoadMapping()) return;
  
  TStopwatch timer1;
  TStopwatch timer2;
  TStopwatch timer3;
  TStopwatch timer4;
  TStopwatch timer5;
  
  Int_t n1 = DataMakerReading(input,timer1,"","",kTRUE,0,0);

  Int_t n2 = DataMakerReading(input,timer2,ocdbPath,"NOGAIN");

  Int_t n3 = DataMakerReading(input,timer3,ocdbPath,"GAINCONSTANTCAPA");

  Int_t n4 = DataMakerReading(input,timer4,ocdbPath,"GAIN");

  Int_t n5 = DataMakerReading(input,timer5,ocdbPath,"INJECTIONGAIN");

  Print("DataMakerReading(HRAW)",timer1,n1);  
  Print("DataMakerReading(HCALZ)",timer2,n2);
  Print("DataMakerReading(HCALG)",timer3,n3);
  Print("DataMakerReading(HCALC)",timer4,n4);
  Print("DataMakerReading(HCALE)",timer5,n5);
  
  AliMUONPainterDataRegistry* reg = AliMUONPainterDataRegistry::Instance();
  
  TFile f(gSystem->ExpandPathName(Form("%s.root",outputBase)),"RECREATE");
  ofstream out(gSystem->ExpandPathName(Form("%s.log",outputBase)));

  Occupancy(out);

  for ( Int_t i = 0; i < reg->NumberOfDataSources(); ++i )
  {
    AliMUONVTrackerData* data = reg->DataSource(i);
    data->Write();
  }
  
  f.Close();
  
  AliCodeTimer::Instance()->Print();

}
