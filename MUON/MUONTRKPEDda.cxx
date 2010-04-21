/*
  Contact: Jean-Luc Charvet <jean-luc.charvet@cern.ch>  
  Link: http://aliceinfo.cern.ch/static/Offline/dimuon/muon_html/README_mchda.html
  Reference Run: 109302 (station 3 only)
  Run Type: PEDESTAL
  DA Type: LDC
  Number of events needed: 400 events for pedestal run
  Input Files: mutrkpedvalues and config_ldc-MTRK-S3-0 in path : /afs/cern.ch/user/j/jcharvet/public/DA_validation 
  Output Files: local dir (not persistent) -> MUONTRKPEDda.ped  FXS -> run<#>_MCH_<ldc>_PEDESTALS
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
        2010-04-18 New version: MUONTRKPEDda.cxx,v 1.6
	-------------------------------------------------------------------------

	Version for MUONTRKPEDda MUON tracking
	(A. Baldisseri, J.-L. Charvet)


 Rem:  AliMUON2DMap stores all channels, even those which are not connected
 if pedMean == -1, channel not connected to a pad  

&

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
#include "AliLog.h"
#include "AliMUONDspHeader.h"
#include "AliDAQ.h"

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
int main(Int_t argc, const char **argv) 
{
  Int_t status=0;
  TStopwatch timers;
  timers.Start(kTRUE); 

  const char* prefixDA = "MUONTRKPEDda"; // program prefix
  printf(" ######## Begin execution : %s ######## \n\n",prefixDA); 

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
//   cout << argv[0];
 
  Int_t skipEvents = 0;
  Int_t maxEvents  = 1000000;
  Int_t maxDateEvents  = 1000000;
  TString inputFile;

  Int_t  nDateEvents = 0;
  Int_t nGlitchErrors= 0;
  Int_t nParityErrors= 0;
  Int_t nPaddingErrors= 0;
  Int_t nTokenlostErrors= 0;
  Int_t recoverParityErrors = 1;

  TString logOutputFile;

  Char_t flatFile[256]="";
  TString shuttleFile;

  Int_t nEventsRecovered = 0;
  Int_t nEvents = 0;
  UInt_t runNumber   = 0;
  Int_t nConfig = 1;
  Int_t nEvthreshold = 10; //below this nb_evt pedestal are not calculated and forced to 4085 (sigma)
  ofstream filcout;

  // decode the input line
  for (Int_t i = 1; i < argc; i++) // argument 0 is the executable name
    {
      const char* arg = argv[i];

      arg = argv[i];
      if (arg[0] != '-') 
	{
	  // If only one argument and no "-" => DA calling from ECS
	  if (argc == 2)
	    {
        inputFile=argv[i];
	    }
	  continue;
	}
      switch (arg[1])
	{
	case 'f' : 
	  i++;
      inputFile=argv[i];
	  nConfig=0;
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
	case 'h' :
	  i++;
	  printf("\n******************* %s usage **********************",argv[0]);
	  printf("\nOnline (called from ECS) : %s <raw data file> (no inline options)\n",argv[0]);
	  printf("\n%s can be used locally only with options (without DiMuon configuration file)",argv[0]);
	  printf("\n%s -options, the available options are :",argv[0]);
	  printf("\n-h help                    (this screen)");
	  printf("\n");
	  printf("\n Input");
	  printf("\n-f <raw data file>         (default = %s)",inputFile.Data()); 
	  printf("\n");
	  printf("\n Output");
	  printf("\n-a <Flat ASCII file>       (default = %s)",shuttleFile.Data()); 
	  printf("\n");
	  printf("\n Options");
	  printf("\n-m <max date events>       (default = %d)",maxDateEvents);
	  printf("\n-s <skip events>           (default = %d)",skipEvents);
	  printf("\n-n <max events>            (default = %d)",maxEvents);

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

  Char_t dbfile[256]="";
  // nConfig=1 : Reading configuration (or not) status via "mutrkpedvalues" file located in DetDB
  if(nConfig)
    { 
      sprintf(dbfile,"mutrkpedvalues");
      status=daqDA_DB_getFile(dbfile,dbfile);
      if(status) {printf(" !!! Failed  : input file %s is missing, status = %d\n",dbfile,status); return -1; } 
      ifstream filein(dbfile,ios::in);
      filein >> nConfig;
      //      filein >> nEvthreshold;
    }
  else  printf(" ***  Config= %d: no configuration ascii file is used \n",nConfig); 
  muonPedestal->SetconfigDA(nConfig);
  muonPedestal->SetnEvthreshold(nEvthreshold);

  // nConfig=1: configuration ascii file config_$DATE_ROLE_NAME read from DetDB
  if(nConfig)
    {
      sprintf(dbfile,"config_%s",getenv("DATE_ROLE_NAME"));
      status=daqDA_DB_getFile(dbfile,dbfile);
      if(status) {printf(" !!! Failed  : Configuration file %s is missing, status = %d\n",dbfile,status); return -1; }
      //      else printf(" *** Copy ascii config file: %s from DetDB to working directory and reading ...*** \n",dbfile);
      muonPedestal->LoadConfig(dbfile);  
    } 

  // Rawdeader, RawStreamHP
  AliRawReader* rawReader = AliRawReader::Create(inputFile);
  AliMUONRawStreamTrackerHP* rawStream  = new AliMUONRawStreamTrackerHP(rawReader);    
//	rawStream->DisableWarnings();
  rawStream->EnabbleErrorLogger();
  //
  // kLowErrorDetail,     /// Logs minimal information in the error messages.
  // kMediumErrorDetail,  /// Logs a medium level of detail in the error messages.
  // kHighErrorDetail     /// Logs maximum information in the error messages.
  //  rawStream->SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kLowErrorDetail);
     rawStream->SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kMediumErrorDetail);
  //   rawStream->SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kHighErrorDetail);

  printf("\n %s : Reading data from file %s\n",prefixDA,inputFile.Data());

  Int_t tabTokenError[20][14];
  for ( Int_t i=0 ; i<20 ; i++) { for ( Int_t j=0 ; j<14 ; j++) { tabTokenError[i][j]=0;}	}

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
	  sprintf(flatFile,"%s.log",prefixDA);
	  logOutputFile=flatFile;
		AliLog::SetStreamOutput(&filcout); // Print details on logfile
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

      const char* detail = "";
      // First lopp over DDL's to find good events
      // Error counters per event (counters in the decoding lib are for each DDL)
      Bool_t eventIsErrorMessage = kFALSE;
      int eventGlitchErrors = 0;
      int eventParityErrors = 0;
      int eventPaddingErrors = 0;
      int eventTokenlostErrors = 0;
      rawStream->First();
      do
	{
	  if (rawStream->IsErrorMessage()) eventIsErrorMessage = kTRUE;
	  eventGlitchErrors += rawStream->GetGlitchErrors();
	  eventParityErrors += rawStream->GetParityErrors();
	  eventPaddingErrors += rawStream->GetPaddingErrors();
	  eventTokenlostErrors += rawStream->GetTokenLostErrors();
	  if (rawStream->GetTokenLostErrors())
	    {
	      nTokenlostErrors++;
	      const AliMUONRawStreamTrackerHP::AliBlockHeader*      blkHeader  = 0x0;
	      const AliMUONRawStreamTrackerHP::AliDspHeader*        dspHeader  = 0x0;
	      Int_t nBlock = rawStream->GetBlockCount();
	      for(Int_t iBlock = 0; iBlock < nBlock ;iBlock++)
		{
		  blkHeader = rawStream->GetBlockHeader(iBlock);
		  //		  printf("Block %d Total length %d\n",iBlock,blkHeader->GetTotalLength());
		  Int_t nDsp = rawStream->GetDspCount(iBlock);
		  //		  printf("Block %d DSP %d\n",iBlock,nDsp);		  
		  for(Int_t iDsp = 0; iDsp < nDsp ;iDsp++)
		    {
		      dspHeader =  blkHeader->GetDspHeader(iDsp);
		      //		      printf("Dsp %d Add %X\n",iDsp,dspHeader);
		      if (dspHeader->GetErrorWord())
			{
			  Int_t ddl = rawStream->GetDDL()  ; 
			  //	 Int_t ddl = AliDAQ::DdlID("MUONTRK", rawStream->GetDDL()) - 2560 ; // format 2560 + ddl
			  Int_t frt = (dspHeader->GetErrorWord() & 0xFFFF0000) >> 16 ; // 4*4bits right shift
			  tabTokenError[ddl][frt]++;
			  //	 printf(" DDL %d error word %X %d %d\n",ddl,dspHeader->GetErrorWord(),frt,tabTokenError[8][4]);
			}
		      
		    }
		}
	    }
	} while(rawStream->NextDDL()); 

      AliMUONRawStreamTrackerHP::AliBusPatch* busPatch;
      if (!eventIsErrorMessage) 
	{
	  // Good events (no error) -> compute pedestal for all channels
	  rawStream->First(); 
	  nEvents++;
	  muonPedestal->SetAliNCurrentEvents(nEvents);
	  while( (busPatch = (AliMUONRawStreamTrackerHP::AliBusPatch*) rawStream->Next())) 
	    {
	      for(int i = 0; i < busPatch->GetLength(); ++i)
		{
		  busPatch->GetData(i, manuId, channelId, charge);
		  muonPedestal->MakePed(busPatch->GetBusPatchId(), (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
		}
	    }
	}
      else
	{
	  // Events with errors
	  if (recoverParityErrors && eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
	    {
	      filcout << " ----------- Date Event recovered = " << nDateEvents <<  " ----------------" << endl;
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
	      nEvents++;
	      muonPedestal->SetAliNCurrentEvents(nEvents);
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
			  busPatch->GetData(i, manuId, channelId, charge);
			  muonPedestal->MakePed(busPatch->GetBusPatchId(), (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
			}
		    }
		  else
		    {
		      AliMUONErrorCounter* errorCounter;
		      // Bad buspatch -> not used (just print)
		      filcout<<"bpId "<<busPatch->GetBusPatchId()<<" words "<<busPatch->GetLength()
				 <<" parity errors "<<errorCount<<endl;
		      // Number of events where this buspatch is missing
		      if (!(errorCounter = (AliMUONErrorCounter*) (muonPedestal->GetErrorBuspatchTable()->FindObject(busPatch->GetBusPatchId()))))
			{
			  // New buspatch
			  errorCounter = new AliMUONErrorCounter(busPatch->GetBusPatchId());
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
	      //	      nEvents++;
	      nEventsRecovered++;
	    } //end of if (recoverParityErrors && eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
	  else
	    {
	      // Fatal errors reject the event
	      detail = Form(" ----------- Date Event rejected = %d  ----------------",nDateEvents);
	      filcout << detail << endl;
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
		     <<" Padding "<<eventPaddingErrors
		     <<" Token lost "<<eventTokenlostErrors<<endl;
	  filcout<<endl;			
	} // end of if (!rawStream->IsErrorMessage())

      if (eventGlitchErrors)  nGlitchErrors++;
      if (eventParityErrors)  nParityErrors++;
      if (eventPaddingErrors) nPaddingErrors++;
      //      if (eventTokenlostErrors) nTokenlostErrors++;
      //      muonPedestal->SetAliNCurrentEvents(nEvents);

    } // while (rawReader->NextEvent())
  delete rawReader;
  delete rawStream;

  sprintf(flatFile,"%s.ped",prefixDA);
  if(shuttleFile.IsNull())shuttleFile=flatFile;
  muonPedestal->SetAliNEvents(nEvents);
  muonPedestal->SetAliRunNumber(runNumber);
  
  muonPedestal->Finalize();  
  muonPedestal->MakeControlHistos();  
  if (!shuttleFile.IsNull())  
  {
    ofstream out(shuttleFile.Data());  
    muonPedestal->MakeASCIIoutput(out);
    out.close();
#ifdef ALI_AMORE
  //
  //Send objects to the AMORE DB
  //
    ostringstream stringout;
    muonPedestal->MakeASCIIoutput(stringout);
    
    amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
    TObjString peddata(stringout.str().c_str());
    Int_t amoreStatus = amoreDA.Send("Pedestals",&peddata);
    if ( amoreStatus )
      cout << "Warning: Failed to write Pedestals in the AMORE database : " << amoreStatus << endl;
    else 
      cout << "amoreDA.Send(Pedestals) ok" << endl;  
#else
    cout << "Warning: MCH DA not compiled with AMORE support" << endl;
#endif
    
  }

  // writing some counters
  cout << endl;
  cout << prefixDA << " : Nb of DATE events           = " << nDateEvents    << endl;
  cout << prefixDA << " : Nb of Glitch errors         = "   << nGlitchErrors  << endl;
  cout << prefixDA << " : Nb of Parity errors         = "   << nParityErrors  << endl;
  cout << prefixDA << " : Nb of Padding errors        = "   << nPaddingErrors << endl;	
  cout << prefixDA << " : Nb of Token lost errors     = "   << nTokenlostErrors << endl;
  cout << prefixDA << " : Nb of events recovered      = "   << nEventsRecovered<< endl;
  cout << prefixDA << " : Nb of events without errors = "   << nEvents-nEventsRecovered<< endl;
  cout << prefixDA << " : Nb of events used           = "   << nEvents        << endl;

  filcout << endl;
  filcout << prefixDA << " : Nb of DATE events           = " << nDateEvents    << endl;
  filcout << prefixDA << " : Nb of Glitch errors         = "   << nGlitchErrors << endl;
  filcout << prefixDA << " : Nb of Parity errors         = "   << nParityErrors << endl;
  filcout << prefixDA << " : Nb of Padding errors        = "   << nPaddingErrors << endl;
  filcout << prefixDA << " : Nb of Token lost errors     = "   << nTokenlostErrors << endl;
  filcout << prefixDA << " : Nb of events recovered      = "   << nEventsRecovered<< endl;	
  filcout << prefixDA << " : Nb of events without errors = "   << nEvents-nEventsRecovered<< endl;
  filcout << prefixDA << " : Nb of events used           = "   << nEvents        << endl;

  // ouput files
  cout << endl;
  cout << prefixDA << " : Output logfile         : " << logOutputFile  << endl;
  cout << prefixDA << " : Pedestal Histo file    : " << muonPedestal->GetHistoFileName() << endl;
  cout << prefixDA << " : Ped. file (to SHUTTLE) : " << shuttleFile << endl;   

  // Writing Token Error table
  if(nTokenlostErrors)
    {
      char* detail=Form("\nWarning: Token Lost occurence \n");
      printf("%s",detail);
      filcout <<  detail ;
      for ( Int_t i=0 ; i<20 ; i++) 
	{ 
	  for ( Int_t j=4 ; j<14 ; j++) 
	    { 
	      if(tabTokenError[i][j]>0)
		{
		  Int_t tab=tabTokenError[i][j];
		  Int_t frt=j/2-1;
		  Int_t station = i/4 +1;
		  if( j % 2 == 0)detail=Form(" in DDL= %d (station %d) and FRT%d ( Up ) => %d Token errors (address = 0x%X0000)",2560+i,station,frt,tab,j);
		  else detail=Form(" in DDL= %d (station %d) and FRT%d (Down) => %d Token errors (address = 0x%X0000)",2560+i,station,frt,tab,j);
		  printf("%s\n",detail);
		  filcout <<  detail << endl;
		}
	    }
	}
    }

  filcout << endl;
  filcout << prefixDA << " : Output logfile         : " << logOutputFile  << endl;
  filcout << prefixDA << " : Pedestal Histo file    : " << muonPedestal->GetHistoFileName() << endl;
  filcout << prefixDA << " : Ped. file (to SHUTTLE) : " << shuttleFile << endl;   

 // Copying files to local DB folder defined by DAQ_DETDB_LOCAL
  Char_t *dir;
  dir= getenv("DAQ_DETDB_LOCAL");
  unsigned int nLastVersions=50;
  cout << "\n ***  Local DataBase: " << dir << " (Max= " << nLastVersions << ") ***" << endl;
  status = daqDA_localDB_storeFile(muonPedestal->GetHistoFileName(),nLastVersions);
  if(status)printf(" Store file : %s   status = %d\n",muonPedestal->GetHistoFileName(),status);
  status = daqDA_localDB_storeFile(shuttleFile.Data(),nLastVersions);
  if(status)printf(" Store file : %s   status = %d\n",shuttleFile.Data(),status);
  status = daqDA_localDB_storeFile(logOutputFile.Data(),nLastVersions);
  if(status)printf(" Store file : %s   status = %d\n",logOutputFile.Data(),status);

  filcout.close();

  // Transferring to FES  (be sure that env variable DAQDALIB_PATH is set)
  printf("\n *****  STORE Pedestal FILE in FES ****** \n");
  status = daqDA_FES_storeFile(shuttleFile.Data(),"PEDESTALS");
  if (status) { printf(" !!! Failed to export file : %s , status = %d\n",shuttleFile.Data(),status); return -1; }
  //  else printf(" %s successfully exported to FES  \n",shuttleFile.Data());

  // Transferring to FES  (be sure that env variable DAQDALIB_PATH is set)
  printf("\n *****  STORE Configuration FILE in FES ****** \n");
  status = daqDA_FES_storeFile(dbfile,"CONFIG");
  if (status) { printf(" !!! Failed to export file : %s , status = %d\n",dbfile,status); return -1; }
  //  else printf(" %s successfully exported to FES  \n",dbfile);


  printf("\n ######## End execution : %s ######## \n",prefixDA); 
  timers.Stop();
  printf("\nExecution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());
  return status;
}
