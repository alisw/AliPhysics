/*
 MCH DA for online occupancy
 
 Contact: Laurent Aphecetche <laurent.aphecetche@subatech.in2p3.fr>, Jean-Luc Charvet <jean-luc.charvet@cea.fr>, Alberto Baldisseri <alberto.baldisseri@cea.fr>
 Link: 
 Run Type: PHYSICS STANDALONE
 DA Type: MON
 Number of events needed: all (or at least as much as possible...)
 Input Files: 09000094301009.10.raw
 Output Files: mch.occupancy, to be exported to the DAQ FXS
 Trigger types used: PHYSICS_EVENT
*/

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

///
/// MUON TRACKER DA to compute the hit count at manu level.
///
/// In the end, this DA produces an ASCII file containing
/// the hit count of all the manus (that were seen in the data flow)
/// of the MUON Tracker (and the number of seen events, so we can
/// later on compute the occupancy)
///
/// $Id$

#include "AliMUON2DMap.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONRawStreamTrackerHP.h"
#include "AliMpConstants.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawReaderDate.h"
#include "Riostream.h"
#include "TPluginManager.h"
#include "TROOT.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TStopwatch.h"
#include "daqDA.h"
#include "event.h"
#include "monitor.h"
#include "signal.h"

#ifdef ALI_AMORE
#include <AmoreDA.h>
#include "TObjString.h"
#include "TSystem.h"
#include <sstream>
#endif

const char* OUTPUT_FILE = "mch.occupancy";
const char* DAVERSION = "MUONTRKOCCda v1.4 ($Id$)";

//______________________________________________________________________________
void Add(AliMUONVStore& destStore, const AliMUONVStore& srcStore)
{
  /// Add all elements from srcStore to destStore
  /// Each element of srcStore is supposed to be an AliMUONCalibParamNI,
  /// with ID0=busPatchId and ID1=manuId
  
  TIter next(srcStore.CreateIterator());
  AliMUONVCalibParam* source;
  
  while ( ( source = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    AliMUONCalibParamNI* dest = static_cast<AliMUONCalibParamNI*>(destStore.FindObject(source->ID0(),source->ID1()));
    if (!dest)
    {
      dest = static_cast<AliMUONCalibParamNI*>(source->Clone());
      destStore.Add(dest);
    }
    else
    {
      for ( Int_t i = 0; i < source->Size(); ++i ) 
      {
        for ( Int_t j = 0; j  < source->Dimension(); ++j ) 
        {
          dest->SetValueAsIntFast(i,j,dest->ValueAsIntFast(i,j)+source->ValueAsIntFast(i,j));
        }
      }
    }
  }
}

//______________________________________________________________________________
void GenerateOutputFile(const AliMUONVStore& store, ostream& out, 
                        Int_t runNumber, Int_t nevents)
{
  /// Write the channel hit count (grouped by manu) in the output file.
  
  TIter next(store.CreateIterator());
  AliMUONVCalibParam* manu;
  
  out << "//===========================================================================" << endl;
  out << "//  Hit counter file calculated by " << __FILE__ << endl;
  out << "//===========================================================================" << endl;
  out << "//" << endl;
  out << "//       * Run Number          : " << runNumber << endl;
  out << "//       * File Creation Date  : " << TTimeStamp().AsString("l") << endl;
  out << "//---------------------------------------------------------------------------" << endl;
  out << "//  BP   MANU  SUM_N  NEVENTS" << endl;
  out << "//---------------------------------------------------------------------------" << endl;
  
  while ( ( manu = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t sum(0);
//    Int_t nevents(0);
    
    for ( Int_t i = 0; i < manu->Size(); ++i ) 
    {
      sum += manu->ValueAsInt(i);
      //      nevents = TMath::Max(nevents,manu->ValueAsInt(i,1));
      // nevents = TMath::Max(nevents,manu->ValueAsInt(i,1));
    }
    
    out << Form("%5d %5d %10d %10d",manu->ID0(),manu->ID1(),sum,nevents) << endl;
  }
}
  
//______________________________________________________________________________
int main(int argc, char **argv) 
{
  /// Main method.
  ///
  /// We loop over all physics events.
  /// For each event we store the channels that were hit for that event.
  /// If the event is good, we then increment the list of channels hit for
  /// the whole run.
  /// We delete the store for a single event and move to next event.
  ///
  /// In the end we output an ASCII file with the necessary information 
  /// to compute the occupancy later on, i.e. the number of times channels
  /// were seen per manu, and the number of events.
  ///
  
  signal(SIGSEGV,SIG_DFL); // to be able to get core dumps...
  
  TStopwatch timers;
  timers.Start(kTRUE); 
  
  ios::sync_with_stdio();

  cout << "Running " << DAVERSION << endl;
  
  if ( argc < 2 ) 
  {
    cout << "Wrong number of arguments" << endl;
    cout << "Usage : " << argv[0] << " datasource1 [datasource2] ..." << endl;
    return -1;
  }
  
  // needed for streamer application
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");
  
  Int_t numberOfEvents(0);
  Int_t numberOfPhysicsEvent(0);
  Int_t numberOfBadEvents(0);
  Int_t numberOfUsedEvents(0);
  
  AliMUON2DMap oneEventData(kTRUE);
  AliMUON2DMap accumulatedData(kTRUE);
  
  UInt_t runNumber(0);

  for ( Int_t i = 1; i < argc; ++i ) 
  {    
    int status;
    AliRawReaderDate* rawReader(0x0);
    
// define data source : 
   status=monitorSetDataSource(argv[i]);
    if (status!=0) 
    {
      printf("MCH Occupancy DA ERROR: monitorSetDataSource() failed: %s\n", monitorDecodeError(status));
      return -1;
    }
    
// Declare monitoring program 
    status=monitorDeclareMp("MUON_TRK_OCC");
    if (status!=0) 
    {
      printf("MCH Occupancy DA ERROR: monitorDeclareMp() failed: %s\n", monitorDecodeError(status));
      return -1;
    }
// Define wait event timeout - 1s max
    monitorSetNowait();
    monitorSetNoWaitNetworkTimeout(1000);
    
    for(;;)
    {
      struct eventHeaderStruct *event(0x0);
      eventTypeType eventT;
      
      status=monitorGetEventDynamic((void **)&event);
      if (status!=0)
      {
        printf("MCH Occupancy DA ERROR: %s\n", monitorDecodeError(status));
        delete event;
        break;
      }

      /* check shutdown condition */
      if (daqDA_checkShutdown())
      {
        delete event;
        break;
      }
      
      /* retry if got no event */
      if (event==NULL) continue;

      ++numberOfEvents;

      eventT=event->eventType;
      if ((eventT == END_OF_RUN)||(eventT == END_OF_RUN_FILES)) 
      {
        delete event;
        break;
      }
      if (eventT != PHYSICS_EVENT) 
      {
        delete event;
        continue;
      }
                 
      ++numberOfPhysicsEvent;
      
      rawReader = new AliRawReaderDate((void*)event);
      
      if ( rawReader->GetRunNumber() != runNumber )
      {
        if ( runNumber != 0 ) 
        {
          cout << "Uh oh. That's bad... Changing of run number ???" << endl;
          delete event;
          delete rawReader;
          return -9999;
        }
        runNumber = rawReader->GetRunNumber();
      }
      
      AliMUONRawStreamTrackerHP stream(rawReader);
      
      stream.DisableWarnings();
      
      oneEventData.Clear();
      
      Int_t buspatchId;
      UShort_t  manuId;
      UChar_t manuChannel;
      UShort_t adc;
      
      stream.First();
            
      while ( stream.Next(buspatchId,manuId,manuChannel,adc,kTRUE) )
      {    
        AliMUONVCalibParam* one = static_cast<AliMUONVCalibParam*>(oneEventData.FindObject(buspatchId,manuId));
        
        if (!one)
        {
          one = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),buspatchId,manuId);
          oneEventData.Add(one);
        }
        
        one->SetValueAsInt(manuChannel,0,one->ValueAsInt(manuChannel,0)+1);
      }
      
      Bool_t badEvent = stream.HasPaddingError() || stream.HasGlitchError();
      
      if ( !badEvent )
      {
        ++numberOfUsedEvents;
        Add(accumulatedData,oneEventData);
      }
      else
      {
        ++numberOfBadEvents;
      }

      delete event;
      
      delete rawReader;
    }    
  }
  
  
  cout << Form("%12d events processed : %12d physics %d used ones %d bad ones",
               numberOfEvents,numberOfPhysicsEvent,numberOfUsedEvents,numberOfBadEvents) << endl;
    
  ofstream fout(OUTPUT_FILE);
  
  GenerateOutputFile(accumulatedData,fout,runNumber,numberOfUsedEvents);
  
  fout.close();
  
#ifdef ALI_AMORE
  
  if ( numberOfUsedEvents ) // do not update the AMORE pool with an empty object...
  {
    // Send occupancy store (as a big string) to the AMORE DB
    amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
    
    ostringstream str;
    
    GenerateOutputFile(accumulatedData,str,runNumber,numberOfUsedEvents);
    
    TObjString occupancyAsString(str.str().c_str());
    
    Int_t status = amoreDA.Send("Occupancy",&occupancyAsString);
    if ( status )
    {
      cerr << "ERROR : Failed to write occupancies in the AMORE database : " << status << endl;
    } 
  }

#endif
  
  /* store the result file on FXS */  
  if (daqDA_FES_storeFile(OUTPUT_FILE,"OCCUPANCY")) return -9;

  timers.Stop();
  printf("\nExecution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

  return 0;
}
