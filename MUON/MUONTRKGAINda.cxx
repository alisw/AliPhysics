/*
 Contact: Jean-Luc Charvet <jean-luc.charvet@cea.fr>
 Link: http://aliceinfo.cern.ch/static/Offline/dimuon/muon_html/README_Mchda
 Run Type: CALIBRATION
 DA Type: LDC
 Number of events needed: 400 events for each calibration run
 Input Files: Rawdata file (DATE format)
 Output Files: local dir (not persistent) -> MUONTRKGAINda_<run#>.par
 FXS -> run<#>_MCH_<ldc>_GAINS
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
 2009-05-18 New version: MUONTRKGAINda.cxx,v 1.1
 -------------------------------------------------------------------------

 Version for MUONTRKGAINda MUON tracking
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
#include "AliMUONGain.h"
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

 Char_t prefixDA[256]="MUONTRKGAINda"; // program prefix 

 Int_t skipEvents = 0;
 Int_t maxEvents  = 1000000;
 Int_t maxDateEvents  = 1000000;
 Char_t inputFile[256]="";

 Int_t  nDateEvents = 0;
 Int_t nGlitchErrors= 0;
 Int_t nParityErrors= 0;
 Int_t nPaddingErrors= 0;

 TString logOutputFile;

 Char_t flatFile[256]="";
 TString shuttleFile;

 Int_t nEventsRecovered = 0;
 Int_t nEvents = 0;
 UInt_t runNumber   = 0;
 Int_t  nChannel    = 0;
 ofstream filcout;
 Int_t nIndex = -1; 

 // For DA Gain
 Int_t printLevel  = 0;  // printout (default=0, =1 =>.ped , => .peak & .param)
 Int_t plotLevel  = 1;  // plotout (default=1 => tree , =2 tree+Tgraph+fit)
 Int_t injCharge; 

 Int_t nEntries = daqDA_ECS_getTotalIteration(); // usually = 11 = Nb of calibration runs
 Int_t nInit=1;  // = 0 all DAC values ; = 1 DAC=0 excluded (default=1)
 Int_t nbpf1=6;  // nb of points for linear fit (default=6) 


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
// 	case 'd' :
// 	  i++; 
// 	  printLevel=atoi(argv[i]);
// 	  break;
// 	case 'g' :
// 	  i++; 
// 	  plotLevel=atoi(argv[i]);
// 	  break;
// 	case 'i' :
// 	  i++; 
// 	  nbpf1=atoi(argv[i]);
// 	  break;
// 	case 'j' :
// 	  i++; 
// 	  nInit=atoi(argv[i]);
// 	  break;
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
// 	  printf("\n-d <print level>           (default = %d)",printLevel);
// 	  printf("\n-g <plot level>            (default = %d)",plotLevel);
// 	  printf("\n-i <nb linear points>      (default = %d)",nbpf1);
// 	  printf("\n-j start point of fit      (default = %d)",nInit);
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


 //Gain object
 AliMUONGain* muonGain = new AliMUONGain();
 muonGain->SetprefixDA(prefixDA);
 muonGain->SetAliRootDataFileName(); // MUONTRKGAINda_data.root  

 // Reading current iteration
 nIndex = daqDA_ECS_getCurrentIteration();
 if(nIndex<0)nIndex=0; // compute gain directly from existing root data file




 if(nIndex==0) // compute gain from existing root data file: MUONTRKGAINda_data.root
   {
     sprintf(flatFile,"%s_data.log",prefixDA);
     logOutputFile=flatFile;
     filcout.open(logOutputFile.Data());
     muonGain->SetAlifilcout(&filcout);

     sprintf(flatFile,"%s_data.par",prefixDA);
     if(shuttleFile.IsNull())shuttleFile=flatFile;

     muonGain->SetAliPrintLevel(printLevel);
     muonGain->SetAliPlotLevel(plotLevel);
     muonGain->SetAliIndex(nIndex); // fIndex 
     muonGain->SetAliInit(nInit); // fnInit
     muonGain->SetAliNbpf1(nbpf1); // fnbpf1
     muonGain->MakeGainStore(shuttleFile);
#ifdef ALI_AMORE
     std::ifstream in(shuttleFile.Data());
     ostringstream stringout;
     char line[1024];
     while ( in.getline(line,1024) )
       stringout << line << "\n";  
     in.close();

     amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
     TObjString gaindata(stringout.str().c_str());
     status = amoreDA.Send("Gains",&gaindata);
     if ( status )
       cout << "Warning: Failed to write Pedestals in the AMORE database : " << status << endl;
     else 
       cout << "amoreDA.Send(Gains) ok" << endl;  
#else
     cout << "Warning: MCH DA not compiled with AMORE support" << endl;
#endif
   }



 if(nIndex>0) // normal case: reading calibration runs before computing gains
   {
     UShort_t manuId;  
     UChar_t channelId;
     UShort_t charge;
	  
     // DAC values stored in array vDAC (reading dbfile in DETDB)
     //   Int_t vDAC[11] = {0,200,400,800,1200,1600,2000,2500,3000,3500,4000}; // DAC values
     Int_t vDAC[11]; // DAC values
     Char_t *dbfile;
     dbfile="mutrkcalibvalues";
     cout << "\n *** Copy: " << dbfile << " from DetDB to working directory  *** \n" << endl;
     status=daqDA_DB_getFile(dbfile,dbfile);
     ifstream filein(dbfile,ios::in); Int_t k=0, kk;
     while (k<nEntries ) { filein >> kk >> vDAC[k] ; k++; }
     filein >> nInit; // = 0 all DAC values fitted ; = 1 DAC=0 excluded (default=1)
     filein >> nbpf1; // nb of points for linear fit (default=6) 
     filein >> printLevel;  // printout (default=0, =1 =>.ped /run, =2 => .peak & .param)
     filein >> plotLevel;  // plotout (default=1 => tree , =2 tree+Tgraph+fit)
     cout << "nInit=" << nInit << " Nb linear pts=" << nbpf1 << "    Print level=" << printLevel << "	  Plot Level=" << plotLevel << endl;
     muonGain->SetAliPrintLevel(printLevel);
     muonGain->SetAliPlotLevel(plotLevel);
		
     injCharge=vDAC[nIndex-1];

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
	      sprintf(flatFile,"%s.log",prefixDA);
	      logOutputFile=flatFile;

	      filcout.open(logOutputFile.Data());
	      filcout<<"//=================================================" << endl;
	      filcout<<"//       " << prefixDA << " for run = " << runNumber << "  (DAC=" << injCharge << ")" << endl;
	      filcout<<"//=================================================" << endl;
	      filcout<<"//   * Date          : " << muonGain->GetDate()->AsString("l") << "\n" << endl;

	      cout<<"\n ********  " << prefixDA << " for run = " << runNumber << "  (Index= " << nIndex << "/" << nEntries << "  DAC=" << injCharge << ") ********\n" << endl;
	      cout<<" * Date : " << muonGain->GetDate()->AsString("l") << "\n" << endl;
	    }

	  muonGain->SetAlifilcout(&filcout);

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
		      muonGain->MakePed(busPatch->GetBusPatchId(), (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
		    }
		}
	      nEvents++;
	    }
	  else
	    {
	      // Events with errors
	      if (eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
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
			      muonGain->MakePed(busPatch->GetBusPatchId(), (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
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
			  if (!(errorCounter = (AliMUONErrorCounter*) (muonGain->GetErrorBuspatchTable()->FindObject(bpname))))
			    {
			      // New buspatch
			      errorCounter = new AliMUONErrorCounter(busPatch->GetBusPatchId());
			      errorCounter->SetName(bpname);
			      muonGain->GetErrorBuspatchTable()->Add(errorCounter);
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
		} //end of if (eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
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

     // process and store mean peak values in .root file
    sprintf(flatFile,"%s.par",prefixDA);
     if(shuttleFile.IsNull())shuttleFile=flatFile;
     injCharge=vDAC[nIndex-1];
     muonGain->SetAliIndex(nIndex); // fIndex 
     muonGain->SetAliInjCharge(injCharge);
     muonGain->SetAliNEvents(nEvents);
     muonGain->SetAliRunNumber(runNumber);
     muonGain->SetAliNChannel(nChannel);
     muonGain->MakePedStoreForGain(shuttleFile);
     

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

     // Computing gain 
     if(nIndex==nEntries)
	{
	  muonGain->SetAliInit(nInit); // fnInit
	  muonGain->SetAliEntries(nEntries); // fnEntries
	  muonGain->SetAliNbpf1(nbpf1); // fnbpf1
     	  muonGain->MakeGainStore(shuttleFile);
#ifdef ALI_AMORE  
          std::ifstream in(shuttleFile.Data());
          ostringstream stringout;
          char line[1024];
          while ( in.getline(line,1024) )
  	    stringout << line << "\n";  
          in.close();
	  
          amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
          TObjString gaindata(stringout.str().c_str());
          status = amoreDA.Send("Gains",&gaindata);
          if ( status )
            cout << "Warning: Failed to write Pedestals in the AMORE database : " << status << endl;
          else 
            cout << "amoreDA.Send(Gains) ok" << endl;  
#else
          cout << "Warning: MCH DA not compiled with AMORE support" << endl;
#endif
	}

// ouput files
      filcout << endl;
      filcout << prefixDA << " : Root data file         : " << muonGain->GetRootDataFileName() << endl;
      filcout << prefixDA << " : Output logfile         : " << logOutputFile  << endl;
      filcout << prefixDA << " : Gain Histo file        : " << muonGain->GetHistoFileName() << endl;
      filcout << prefixDA << " : Gain file (to SHUTTLE) : " << shuttleFile << endl;

//	 Copying files to local DB folder defined by DAQ_DETDB_LOCAL
       Char_t *dir;
      dir= getenv("DAQ_DETDB_LOCAL");
      unsigned int nLastVersions = 99;
      cout << "\n ***  Local DataBase: " << dir << " (Max= " << nLastVersions << ") ***" << endl;
      status = daqDA_localDB_storeFile(logOutputFile.Data(),nLastVersions);
      printf(" Store file : %s   status = %d\n",logOutputFile.Data(),status);
      if(nIndex==nEntries)
       {
	 status = daqDA_localDB_storeFile(muonGain->GetRootDataFileName(),nLastVersions);
	 printf(" Store file : %s   status = %d\n",muonGain->GetRootDataFileName(),status);
	 status = daqDA_localDB_storeFile(muonGain->GetHistoFileName(),nLastVersions);
	 printf(" Store file : %s   status = %d\n",muonGain->GetHistoFileName(),status);
	 status = daqDA_localDB_storeFile(shuttleFile.Data(),nLastVersions);
	 printf(" Store file : %s   status = %d\n",shuttleFile.Data(),status);
       }      

    } // end (nIndex>0)

// ouput files
  cout << endl;
  cout << prefixDA << " : Root data file         : " << muonGain->GetRootDataFileName() << endl;
  cout << prefixDA << " : Output logfile         : " << logOutputFile  << endl;
  cout << prefixDA << " : Gain Histo file        : " << muonGain->GetHistoFileName() << endl;
  cout << prefixDA << " : Gain file (to SHUTTLE) : " << shuttleFile << endl;   
 

  // Transferring to OCDB via the SHUTTLE
  printf("\n *****  STORE FILE in FES ****** \n");

  // be sure that env variable DAQDALIB_PATH is set in script file
  //       gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/infoLogger");

  status = daqDA_FES_storeFile(shuttleFile.Data(),"GAINS");
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

