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
/// MUONOfflineShift("path_to_raw_file","basename of output file")
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
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include <Riostream.h>
#include <TFile.h>
#include <TGrid.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TSystem.h>

#endif

//______________________________________________________________________________
Int_t DataMakerReading(AliRawReader* rawReader, 
                       TStopwatch& timer,
                       const char* cdbPath="",
                       const char* calibMode="",
                       Bool_t histogram=kFALSE,
                       Double_t xmin = 0.0,
                       Double_t xmax = 4096.0)
{
  /// Run over the data and calibrate it if so required (if cdbPath != "")
  /// calibMode can be :
  /// - NOGAIN           : only zero-suppression will be done
  /// - GAINCONSTANTCAPA : zero-suppression + gain, but with a single capa value for all channels
  /// - GAIN             : zero-suppression + gain w/ individual capacitance per channel.
  
  rawReader->RewindEvents();
  
  AliMUONVTrackerDataMaker* dm(0x0);
  
  if ( strlen(cdbPath) > 0 ) 
  {
    dm = new AliMUONTrackerCalibratedDataMaker(rawReader,cdbPath,calibMode,histogram,xmin,xmax);
  }
  else  
  {
    dm = new AliMUONTrackerRawDataMaker(rawReader,kTRUE);
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
  
  delete dm;
  
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
void MUONOfflineShift(const char* input, const char* outputBase,
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
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetRun(0);
  AliMpCDB::LoadDDLStore();
  
  TString fileName(gSystem->ExpandPathName(input));
  
  AliRawReader* rawReader(0x0);
  
  // check extention to choose the rawdata file format
  if (fileName.EndsWith(".root")) 
  {
    rawReader = new AliRawReaderRoot(fileName);
  }
  else if (!fileName.IsNull()) 
  {
    rawReader = new AliRawReaderDate(fileName); // DATE file
  }
  
  if (!rawReader) return;
  
  TStopwatch timer2;
  TStopwatch timer3;
  
  Int_t n2 = DataMakerReading(rawReader,timer2);
  
  Int_t n3 = DataMakerReading(rawReader,timer3,ocdbPath);
  
  Print("DataMakerReading(RAW)",timer2,n2);
  
  Print("DataMakerReading(CAL)",timer3,n3);
  
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
  
}
