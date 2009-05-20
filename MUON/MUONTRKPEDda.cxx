/*
Contact: Jean-Luc Charvet <jean-luc.charvet@cea.fr>
Link: http://aliceinfo.cern.ch/static/Offline/dimuon/muon_html/README_Mchda
Run Type: PEDESTAL
	DA Type: LDC
	Number of events needed: 400 events for pedestal run
	Input Files: Rawdata file (DATE format)
	Output Files: local dir (not persistent) -> MUONTRKPEDda_<run#>.ped 
	FXS -> run<#>_MCH_<ldc>_PEDESTALS
	Trigger types used:
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

/* $Id$ */

/*
	-------------------------------------------------------------------------
	2009-05-12 New version: MUONTRKPEDda.cxx,v 1.0
	-------------------------------------------------------------------------

	Version for MUONTRKPEDda MUON tracking
	(A. Baldisseri, J.-L. Charvet)


 Rem:  AliMUON2DMap stores all channels, even those which are not connected
 if pedMean == -1, channel not connected to a pad  


*/
extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <math.h> 

//AliRoot
#include "AliMUONRawStreamTrackerHP.h"
#include "AliRawReader.h"
#include "AliMUONVStore.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMpIntPair.h"
#include "AliMpConstants.h"
#include "AliRawDataErrorLog.h"
#include "AliMUONTrackerIO.h"

//ROOT
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TFitter.h"
#include "TObjString.h"
#include "THashTable.h"
#include <THashList.h>
//
//AMORE
//
#ifdef ALI_AMORE
#include <AmoreDA.h>
#endif

#include "AliMUONPedestal.h"
#include "AliMUONErrorCounter.h"


// main routine
int main(Int_t argc, Char_t **argv) 
{

  Int_t status=0;
  TStopwatch timers;
  timers.Start(kTRUE); 

  // needed for streamer application
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()"); 

  // needed for Minuit plugin
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer",
					"Minuit",
					"TMinuitMinimizer",
					"Minuit",
					"TMinuitMinimizer(const char*)");

  Char_t prefixDA[256]="MUONTRKPEDda"; // program prefix
//   cout << argv[0];
 
  Int_t skipEvents = 0;
  Int_t maxEvents  = 1000000;
  Int_t maxDateEvents  = 1000000;
  Char_t inputFile[256]="";

  Int_t  nDateEvents = 0;
  Int_t nGlitchErrors= 0;
  Int_t nParityErrors= 0;
  Int_t nPaddingErrors= 0;
  Int_t recoverParityErrors = 1;

  TString logOutputFile;

  Char_t flatFile[256]="";
  TString shuttleFile;

  Int_t nEventsRecovered = 0;
  Int_t nEvents = 0;
  UInt_t runNumber   = 0;
  Int_t  nChannel    = 0;
  ofstream filcout;
//   Int_t nIndex = -1; 

  // decode the input line
  for (Int_t i = 1; i < argc; i++) // argument 0 is the executable name
    {
      Char_t* arg;

      arg = argv[i];
      if (arg[0] != '-') 
	{
	  // If only one argument and no "-" => DA calling from ECS
	  if (argc == 2)
	    {
	      sprintf(inputFile,argv[i]);
	    }
	  continue;
	}
      switch (arg[1])
	{
	case 'f' : 
	  i++;
	  sprintf(inputFile,argv[i]);
	  break;
	case 'a' : 
	  i++;
	  shuttleFile = argv[i];
	  break;
	case 's' :
	  i++; 
	  skipEvents=atoi(argv[i]);
	  break;
	case 'm' :
	  i++; 
	  sscanf(argv[i],"%d",&maxDateEvents);
	  break;
	case 'n' :
	  i++; 
	  sscanf(argv[i],"%d",&maxEvents);
	  break;
	case 'p' : 
	  i++;
	  sscanf(argv[i],"%d",&recoverParityErrors);
	  break;
	case 'h' :
	  i++;
	  printf("\n******************* %s usage **********************",argv[0]);
	  printf("\nOnline (called from ECS) : %s <raw data file> (no inline options)\n",argv[0]);
	  printf("\n%s -options, the available options are :",argv[0]);
	  printf("\n-h help                    (this screen)");
	  printf("\n");
	  printf("\n Input");
	  printf("\n-f <raw data file>         (default = %s)",inputFile); 
	  printf("\n");
	  printf("\n Output");
	  printf("\n-a <Flat ASCII file>       (default = %s)",shuttleFile.Data()); 
	  printf("\n");
	  printf("\n Options");
	  printf("\n-m <max date events>       (default = %d)",maxDateEvents);
	  printf("\n-s <skip events>           (default = %d)",skipEvents);
	  printf("\n-n <max events>            (default = %d)",maxEvents);
	  printf("\n-p <Recover parity errors> (default = %d)",recoverParityErrors);

	  printf("\n\n");
	  exit(-1);
	default :
	  printf("%s : bad argument %s (please check %s -h)\n",argv[0],argv[i],argv[0]);
	  argc = 2; exit(-1); // exit if error
	} // end of switch  
    } // end of for i  

  // decoding the events

  UShort_t manuId;  
  UChar_t channelId;
  UShort_t charge;

  //Pedestal object
  AliMUONPedestal* muonPedestal = new AliMUONPedestal();
  muonPedestal->SetprefixDA(prefixDA);

  // Rawdeader, RawStreamHP
  AliRawReader* rawReader = AliRawReader::Create(inputFile);
  AliMUONRawStreamTrackerHP* rawStream  = new AliMUONRawStreamTrackerHP(rawReader);    
  rawStream->DisableWarnings();
  rawStream->EnabbleErrorLogger();

  cout << "\n" << prefixDA << " : Reading data from file " << inputFile  << endl;

  while (rawReader->NextEvent())
    {
      if (nDateEvents >= maxDateEvents) break;
      if (nEvents >= maxEvents) break;
      if (nDateEvents>0 &&  nDateEvents % 100 == 0) 	
	cout<<"Cumulated:  DATE events = " << nDateEvents << "   Used events = " << nEvents << endl;

      // check shutdown condition 
      if (daqDA_checkShutdown()) 
	break;

      //Skip events
      while (skipEvents)
	{
	  rawReader->NextEvent();
	  skipEvents--;
	}

      Int_t eventType = rawReader->GetType();
      runNumber = rawReader->GetRunNumber();

      // Output log file initialisations
      if(nDateEvents==0)
	{
	  sprintf(flatFile,"%s_%d.log",prefixDA,runNumber);
	  logOutputFile=flatFile;

	  filcout.open(logOutputFile.Data());
	  filcout<<"//=================================================" << endl;
	  filcout<<"//       " << prefixDA << " for run = " << runNumber << endl;
	  filcout<<"//=================================================" << endl;
	  filcout<<"//   * Date          : " << muonPedestal->GetDate()->AsString("l") << "\n" << endl;
	  cout<<"\n ********  " << prefixDA << " for run = " << runNumber << " ********\n" << endl;
	  cout<<" * Date          : " << muonPedestal->GetDate()->AsString("l") << "\n" << endl;

	}

      muonPedestal->SetAlifilcout(&filcout);

      nDateEvents++;
      if (eventType != PHYSICS_EVENT)
	continue; // for the moment

      // First lopp over DDL's to find good events
      // Error counters per event (counters in the decoding lib are for each DDL)
      Bool_t eventIsErrorMessage = kFALSE;
      int eventGlitchErrors = 0;
      int eventParityErrors = 0;
      int eventPaddingErrors = 0;
      rawStream->First();
      do
	{
	  if (rawStream->IsErrorMessage()) eventIsErrorMessage = kTRUE;
	  eventGlitchErrors += rawStream->GetGlitchErrors();
	  eventParityErrors += rawStream->GetParityErrors();
	  eventPaddingErrors += rawStream->GetPaddingErrors();
	} while(rawStream->NextDDL()); 

      AliMUONRawStreamTrackerHP::AliBusPatch* busPatch;
      if (!eventIsErrorMessage) 
	{
	  // Good events (no error) -> compute pedestal for all channels
	  rawStream->First(); 
	  while( (busPatch = (AliMUONRawStreamTrackerHP::AliBusPatch*) rawStream->Next())) 
	    {
	      for(int i = 0; i < busPatch->GetLength(); ++i)
		{
		  if (nEvents == 0) nChannel++;
		  busPatch->GetData(i, manuId, channelId, charge);
		  muonPedestal->MakePed(busPatch->GetBusPatchId(), (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
		}
	    }
	  nEvents++;
	}
      else
	{
	  // Events with errors
	  if (recoverParityErrors && eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
	    {
	      // Recover parity errors -> compute pedestal for all good buspatches
	      if ( TEST_SYSTEM_ATTRIBUTE( rawReader->GetAttributes(),
					  ATTR_ORBIT_BC )) 
		{
		  filcout <<"Event recovered -> Period:"<<EVENT_ID_GET_PERIOD( rawReader->GetEventId() )
			      <<" Orbit:"<<EVENT_ID_GET_ORBIT( rawReader->GetEventId() )
			      <<" BunchCrossing:"<<EVENT_ID_GET_BUNCH_CROSSING( rawReader->GetEventId() )<<endl;				
		} 
	      else 
		{
		  filcout <<"Event recovered -> nbInRun:"<<EVENT_ID_GET_NB_IN_RUN( rawReader->GetEventId() )
			      <<" burstNb:"<<EVENT_ID_GET_BURST_NB( rawReader->GetEventId() )
			      <<" nbInBurst:"<<EVENT_ID_GET_NB_IN_BURST( rawReader->GetEventId() )<<endl;
		}
	      rawStream->First();
	      while( (busPatch = (AliMUONRawStreamTrackerHP::AliBusPatch*) rawStream->Next())) 
		{
		  // Check the buspatch -> if error not use it in the pedestal calculation
		  int errorCount = 0;
		  for(int i = 0; i < busPatch->GetLength(); ++i)
		    {
		      if (!busPatch->IsParityOk(i)) errorCount++;
		    }
		  if (!errorCount) 
		    {
		      // Good buspatch
		      for(int i = 0; i < busPatch->GetLength(); ++i)
			{
			  if (nEvents == 0) nChannel++;
			  busPatch->GetData(i, manuId, channelId, charge);
			  // if (busPatch->GetBusPatchId()==1719 && manuId == 1 && channelId == 0) cout <<"Recovered charge "<<charge<<endl;
			  muonPedestal->MakePed(busPatch->GetBusPatchId(), (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
			}
		    }
		  else
		    {
		      char bpname[256];
		      AliMUONErrorCounter* errorCounter;
		      // Bad buspatch -> not used (just print)
		      filcout<<"bpId "<<busPatch->GetBusPatchId()<<" words "<<busPatch->GetLength()
				 <<" parity errors "<<errorCount<<endl;
		      // Number of events where this buspatch is missing
		      sprintf(bpname,"bp%d",busPatch->GetBusPatchId());						
		      if (!(errorCounter = (AliMUONErrorCounter*) (muonPedestal->GetErrorBuspatchTable()->FindObject(bpname))))
			{
			  // New buspatch
			  errorCounter = new AliMUONErrorCounter(busPatch->GetBusPatchId());
			  errorCounter->SetName(bpname);
			  muonPedestal->GetErrorBuspatchTable()->Add(errorCounter);
			}
		      else
			{
			  // Existing buspatch
			  errorCounter->Increment();
			}	
		      // errorCounter->Print();						
		    } // end of if (!errorCount)
		} // end of while( (busPatch = (AliMUONRawStreamTrackerHP ...
	      nEvents++;
	      nEventsRecovered++;
	    } //end of if (recoverParityErrors && eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
	  else
	    {
	      // Fatal errors reject the event
	      if ( TEST_SYSTEM_ATTRIBUTE( rawReader->GetAttributes(),
					  ATTR_ORBIT_BC )) 
		{
		  filcout <<"Event rejected -> Period:"<<EVENT_ID_GET_PERIOD( rawReader->GetEventId() )
			      <<" Orbit:"<<EVENT_ID_GET_ORBIT( rawReader->GetEventId() )
			      <<" BunchCrossing:"<<EVENT_ID_GET_BUNCH_CROSSING( rawReader->GetEventId() )<<endl;				
		} 
	      else 
		{
		  filcout <<"Event rejected -> nbInRun:"<<EVENT_ID_GET_NB_IN_RUN( rawReader->GetEventId() )
			      <<" burstNb:"<<EVENT_ID_GET_BURST_NB( rawReader->GetEventId() )
			      <<" nbInBurst:"<<EVENT_ID_GET_NB_IN_BURST( rawReader->GetEventId() )<<endl;

		}
	    } // end of if (!rawStream->GetGlitchErrors() && !rawStream->GetPaddingErrors() ...
	  filcout<<"Number of errors : Glitch "<<eventGlitchErrors
		     <<" Parity "<<eventParityErrors
		     <<" Padding "<<eventPaddingErrors<<endl;
	  filcout<<endl;			
	} // end of if (!rawStream->IsErrorMessage())

      if (eventGlitchErrors)  nGlitchErrors++;
      if (eventParityErrors)  nParityErrors++;
      if (eventPaddingErrors) nPaddingErrors++;

    } // while (rawReader->NextEvent())
  delete rawReader;
  delete rawStream;

  sprintf(flatFile,"%s_%d.ped",prefixDA,runNumber);
  if(shuttleFile.IsNull())shuttleFile=flatFile;
  muonPedestal->SetAliNEvents(nEvents);
  muonPedestal->SetAliRunNumber(runNumber);
  muonPedestal->SetAliNChannel(nChannel);
  muonPedestal->MakePedStore(shuttleFile);

  // writing some counters
  cout << endl;
  cout << prefixDA << " : Nb of DATE events           = " << nDateEvents    << endl;
  cout << prefixDA << " : Nb of Glitch errors         = "   << nGlitchErrors  << endl;
  cout << prefixDA << " : Nb of Parity errors         = "   << nParityErrors  << endl;
  cout << prefixDA << " : Nb of Padding errors        = "   << nPaddingErrors << endl;		
  cout << prefixDA << " : Nb of events recovered      = "   << nEventsRecovered<< endl;
  cout << prefixDA << " : Nb of events without errors = "   << nEvents-nEventsRecovered<< endl;
  cout << prefixDA << " : Nb of events used           = "   << nEvents        << endl;

  filcout << endl;
  filcout << prefixDA << " : Nb of DATE events           = " << nDateEvents    << endl;
  filcout << prefixDA << " : Nb of Glitch errors         = "   << nGlitchErrors << endl;
  filcout << prefixDA << " : Nb of Parity errors         = "   << nParityErrors << endl;
  filcout << prefixDA << " : Nb of Padding errors        = "   << nPaddingErrors << endl;
  filcout << prefixDA << " : Nb of events recovered      = "   << nEventsRecovered<< endl;	
  filcout << prefixDA << " : Nb of events without errors = "   << nEvents-nEventsRecovered<< endl;
  filcout << prefixDA << " : Nb of events used           = "   << nEvents        << endl;


  // Copying files to local DB folder defined by DAQ_DETDB_LOCAL
  Char_t *dir;
  dir= getenv("DAQ_DETDB_LOCAL");
  unsigned int nLastVersions = 2;
  cout << "\n *** Output files stored locally in " << dir << " (nb of previous versions = " << nLastVersions << ") ***" << endl;
  filcout << "\n *** Output files stored locally in " << dir << " (nb of previous versions = " << nLastVersions << ") ***" << endl;

  // ouput files
  cout << endl;
  cout << prefixDA << " : Output logfile         : " << logOutputFile  << endl;
  cout << prefixDA << " : Gain Histo file        : " << muonPedestal->GetHistoFileName() << endl;
  cout << prefixDA << " : Gain file (to SHUTTLE) : " << shuttleFile << endl;   

  filcout << endl;
  filcout << prefixDA << " : Output logfile         : " << logOutputFile  << endl;
  filcout << prefixDA << " : Gain Histo file        : " << muonPedestal->GetHistoFileName() << endl;
  filcout << prefixDA << " : Gain file (to SHUTTLE) : " << shuttleFile << endl;   

  status = daqDA_localDB_storeFile(muonPedestal->GetHistoFileName(),nLastVersions);
  status = daqDA_localDB_storeFile(shuttleFile.Data(),nLastVersions);
  status = daqDA_localDB_storeFile(logOutputFile.Data(),nLastVersions);

  // Transferring to OCDB via the SHUTTLE
  printf("\n *****  STORE FILE in FES ****** \n");

  // be sure that env variable DAQDALIB_PATH is set in script file
  //       gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/infoLogger");

  status = daqDA_FES_storeFile(shuttleFile.Data(),"PEDESTALS");
  if (status) 
    {
      printf(" Failed to export file : %d\n",status);
    }
  else printf(" %s successfully exported to FES  \n",shuttleFile.Data());

  filcout.close();
  timers.Stop();
  printf("\nExecution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());
  return status;
}