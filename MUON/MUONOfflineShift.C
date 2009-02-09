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
/// \author Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliMUONPainterRegistry.h"
#include "AliMUONTrackerCalibratedDataMaker.h"
#include "AliMUONTrackerRawDataMaker.h"
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
                       const char* calibMode="",
                       Bool_t histogram=kFALSE,
                       Double_t xmin = 0.0,
                       Double_t xmax = 4096.0,
                       Bool_t fastDecoder=kFALSE)
{
  /// Run over the data and calibrate it if so required (if cdbPath != "")
  /// calibMode can be :
  /// - NOGAIN           : only zero-suppression will be done
  /// - GAINCONSTANTCAPA : zero-suppression + gain, but with a single capa value for all channels
  /// - GAIN             : zero-suppression + gain w/ individual capacitance per channel.
  
  TString fileName(gSystem->ExpandPathName(input));
  
  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  
  if (!rawReader) return 0;
  
  AliMUONVTrackerDataMaker* dm(0x0);
  
  if ( strlen(cdbPath) > 0 ) 
  {
    dm = new AliMUONTrackerCalibratedDataMaker(rawReader,cdbPath,calibMode,histogram,xmin,xmax);
  }
  else  
  {
    dm = new AliMUONTrackerRawDataMaker(rawReader,histogram,fastDecoder);
  }
  
  AliMUONPainterRegistry::Instance()->Register(dm);

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
  
  AliMUONPainterRegistry* reg = AliMUONPainterRegistry::Instance();

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
void MUONOfflineShift(const char* input="alien:///alice/data/2008/LHC08a/000021931/raw/08000021931001.50.root", 
                      const char* outputBase="21931.001.50",
                      const char* ocdbPath="alien://folder=/alice/data/2008/LHC08a/OCDB")
{
  /// Entry point of the macro. 
  /// Example of syntax for an input file (from alien)
  ///
  /// alien::///alice/data/2007/LHC07w/000014493/raw/07000014493001.10.root
  /// 
  /// and for an OCDB path : 
  ///
  /// alien://folder=/alice/data/2007/LHC07w/OCDB
  ///

  TGrid::Connect("alien://");
  
  AliCodeTimer::Instance()->Reset();
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0);
  AliMpCDB::LoadDDLStore();
  
  TStopwatch timer0;
  TStopwatch timer1;
  TStopwatch timer2;
  TStopwatch timer3;
  TStopwatch timer4;
  
  Int_t n0 = DataMakerReading(input,timer0,"","",kTRUE,0,0,kFALSE); // using old decoder

  Int_t n1 = DataMakerReading(input,timer1,"","",kTRUE,0,0,kTRUE); // using new decoder

  Int_t n2 = DataMakerReading(input,timer2,ocdbPath,"NOGAIN",kTRUE);

  Int_t n3 = DataMakerReading(input,timer3,ocdbPath,"GAINCONSTANTCAPA",kTRUE);

  Int_t n4 = DataMakerReading(input,timer4,ocdbPath,"GAIN",kTRUE);

  Print("DataMakerReading(HRAW)",timer0,n0);  
  Print("DataMakerReading(HRAW)",timer1,n1);  
  Print("DataMakerReading(HCALZ)",timer2,n2);
  Print("DataMakerReading(HCALG)",timer3,n3);
  Print("DataMakerReading(HCALC)",timer4,n4);
  
  AliMUONPainterRegistry* reg = AliMUONPainterRegistry::Instance();
  
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
