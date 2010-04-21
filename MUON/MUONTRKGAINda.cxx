/*
 Contact: Jean-Luc Charvet <jean-luc.charvet@cern.ch>
 Link: http://aliceinfo.cern.ch/static/Offline/dimuon/muon_html/README_mchda.html
 Reference Runs: (station 3)
 Index		Run
 1		109303
 2		109304
 3		109305
 4		109306
 5 		109307
 6 		109308
 7		109309
 8		109310
 9		109311
 10		109312
 11		109313
 Run Type: CALIBRATION
 DA Type: LDC
 Number of events needed: 400 events for each calibration run (11)
 Input Files:  mutrkcalibvalues and config_ldc-MTRK-S3-0 in path : /afs/cern.ch/user/j/jcharvet/public/DA_validation
 Output Files: local dir (not persistent) -> MUONTRKGAINda.par   FXS -> run<#>_MCH_<ldc>_GAINS
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
 2010-04-18 New version: MUONTRKGAINda.cxx,v 1.6
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
#include "AliMUONGain.h"
#include "AliMUONErrorCounter.h"

// main routine
int main(Int_t argc, const char** argv) 
{
  Int_t status=0;
  TStopwatch timers;
  timers.Start(kTRUE); 
  
  const char* prefixDA = "MUONTRKGAINda"; // program prefix
  printf(" ######## Begin execution : %s ######## \n\n",prefixDA); 
  
  TString inputFile;
  // decode the input line
  if (argc!=2) { printf("Wrong number of arguments\n");  return -1; }
  
  inputFile = argv[1];
  
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
  Int_t skipEvents = 0;
  Int_t maxEvents  = 1000000;
  Int_t maxDateEvents  = 1000000;
  
  Int_t nDateEvents = 0;
  Int_t nGlitchErrors= 0;
  Int_t nParityErrors= 0;
  Int_t nPaddingErrors= 0;
  Int_t nTokenlostErrors= 0;
  Int_t nEventsRecovered = 0;
  Int_t nEvents = 0;
  Int_t nEvthreshold = 10; //below this nb_evt the mean of the charge is not calculated and forced to 4085 (sigma)
  
  TString logOutputFile;
  Char_t flatFile[256]="";
  TString shuttleFile;
  
  UInt_t runNumber   = 0;
  ofstream filcout;
  Int_t nIndex = -1; 
  
  // For DA Gain
  Int_t nEntries = daqDA_ECS_getTotalIteration(); // usually = 11 = Nb of calibration runs
  Int_t nInit=1;  // = 0 all DAC values ; = 1 DAC=0 excluded (default=1)
  Int_t nbpf1=6;  // nb of points for linear fit (default=6) 
  Int_t printLevel  = 0;  // printout (default=0, =1 =>.ped , => .peak & .param)
  Int_t plotLevel  = 1;  // plotout (default=1 => tree , =2 tree+Tgraph+fit)
  Int_t nbev=0; 
  Int_t injCharge; 
  
  // Reading current iteration
  nIndex = daqDA_ECS_getCurrentIteration();
  if(nIndex<0 || nIndex>nEntries) {printf("\n Failed: nIndex = %d\n",nIndex); return -1 ;}
  
  //Gain object
  AliMUONGain* muonGain = new AliMUONGain();
  muonGain->SetprefixDA(prefixDA);
  muonGain->SetAliRootDataFileName(); // MUONTRKGAINda_data.root
  
  UShort_t manuId;  
  UChar_t channelId;
  UShort_t charge;
  
  // DAC values stored in array vDAC (reading dbfile in DETDB)
  //   Int_t vDAC[11] = {0,200,400,800,1200,1600,2000,2500,3000,3500,4000}
  Int_t nConfig = 1; // flag to read or not configuration ascii file in detDB
  Int_t vDAC[11]; // DAC values
  Char_t dbfile[256]="";
  sprintf(dbfile,"mutrkcalibvalues");
  status=daqDA_DB_getFile(dbfile,dbfile);
  if(status) {printf(" Failed  : input file %s is missing, status = %d\n",dbfile,status); return -1; } 
  
  ifstream filein(dbfile,ios::in); Int_t k=0, kk;
  while (k<nEntries ) { filein >> kk >> vDAC[k] ; k++; }
  injCharge=vDAC[nIndex-1];
  
  filein >> nInit; // = 0 all DAC values fitted ; = 1 DAC=0 excluded (default=1)
  filein >> nbpf1; // nb of points for linear fit (default=6) 
  filein >> printLevel;  // printout (default=0, =1 =>.ped /run, =2 => .peak & .param)
  filein >> plotLevel;  // plotout (default=1 => tree , =2 tree+Tgraph+fit)
  filein >> nConfig;  //nConfig (default=1 => read config in DetDB, otherwise =0)
  filein >> nbev;  // Nb of events to read  (default = 0 => reading all events)
  if(nbev>0)maxEvents=nbev;
  
  //  printf(" *** Copy: %s from DetDB to working directory  ***      Config= %d\n",dbfile,nConfig);
  printf(" Input parameters:  nInit= %d   Nb linear pts= %d   Print level= %d   Plot Level= %d    nConfig= %d",nInit,nbpf1,printLevel,plotLevel,nConfig);
  if(nbev==0)printf("\n");
  else printf("  Nb_max evt = %d\n",maxEvents);
  
  muonGain->SetAliPrintLevel(printLevel);
  muonGain->SetAliPlotLevel(plotLevel);
  muonGain->SetconfigDA(nConfig);
  muonGain->SetnEvthreshold(nEvthreshold);
  
  if(nConfig)
  {
    // Laod configuration ascii file from DetDB
    sprintf(dbfile,"config_%s",getenv("DATE_ROLE_NAME"));
    status=daqDA_DB_getFile(dbfile,dbfile);
    if(status) {printf(" Failed  : Configuration file %s is missing, status = %d\n",dbfile,status); return -1; }
    //    else printf(" *** Copy ascii config file: %s from DetDB to working directory and reading ...*** \n",dbfile);
    muonGain->LoadConfig(dbfile);  
  } 
  
  // Rawdeader, RawStreamHP
  AliRawReader* rawReader = AliRawReader::Create(inputFile.Data());
  AliMUONRawStreamTrackerHP* rawStream  = new AliMUONRawStreamTrackerHP(rawReader);    
  //  rawStream->DisableWarnings();
  rawStream->EnabbleErrorLogger();
  //
  // kLowErrorDetail,     /// Logs minimal information in the error messages.
  // kMediumErrorDetail,  /// Logs a medium level of detail in the error messages.
  // kHighErrorDetail     /// Logs maximum information in the error messages.
  //  rawStream->SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kLowErrorDetail);
     rawStream->SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kMediumErrorDetail);
  //   rawStream->SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kHighErrorDetail);
  
  cout << "\n" << prefixDA << " : Reading data from file " << inputFile.Data()  << endl;

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
      muonGain->SetAliNCurrentEvents(nEvents);
      while( (busPatch = (AliMUONRawStreamTrackerHP::AliBusPatch*) rawStream->Next())) 
	    {
	      for(int i = 0; i < busPatch->GetLength(); ++i)
        {
          busPatch->GetData(i, manuId, channelId, charge);
          muonGain->MakePed(busPatch->GetBusPatchId(), (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
        }
	    }
    }
    else
      {
	// Events with errors
	if (eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
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
	    muonGain->SetAliNCurrentEvents(nEvents);
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
	    //	      nEvents++;
	    nEventsRecovered++;
	  } //end of if (eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
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
    //if (eventTokenlostErrors) nTokenlostErrors++;

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
  muonGain->MakePedStoreForGain(shuttleFile);
  
  
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
		  Int_t frt = j/2-1;
		  Int_t station = i/4 +1;
		  if( j % 2 == 0)detail=Form(" in DDL= %d (station %d) and FRT%d ( Up ) => %d Token errors (address = 0x%X0000)",2560+i,station,frt,tab,j);
		  else detail=Form(" in DDL= %d (station %d) and FRT%d (Down) => %d Token errors (address = 0x%X0000)",2560+i,station,frt,tab,j);
		  printf("%s\n",detail);
		  filcout <<  detail << endl;
		}
	    }
	}
    }

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
  unsigned int nLastVersions = 50;
  printf("\n ***  Local DataBase: %s  (Max= %d) ***\n",dir,nLastVersions);
  status = daqDA_localDB_storeFile(logOutputFile.Data(),nLastVersions);
  if(status)printf(" Store file : %s   status = %d\n",logOutputFile.Data(),status);
  if(nIndex==nEntries)
  {
    status = daqDA_localDB_storeFile(muonGain->GetRootDataFileName(),nLastVersions);
    if(status)printf(" Store file : %s   status = %d\n",muonGain->GetRootDataFileName(),status);
    status = daqDA_localDB_storeFile(muonGain->GetHistoFileName(),nLastVersions);
    if(status)printf(" Store file : %s   status = %d\n",muonGain->GetHistoFileName(),status);
    status = daqDA_localDB_storeFile(shuttleFile.Data(),nLastVersions);
    if(status)printf(" Store file : %s   status = %d\n",shuttleFile.Data(),status);
  }      
  
  
  // ouput files
  cout << endl;
  cout << prefixDA << " : Root data file         : " << muonGain->GetRootDataFileName() << endl;
  cout << prefixDA << " : Output logfile         : " << logOutputFile  << endl;
  cout << prefixDA << " : Gain Histo file        : " << muonGain->GetHistoFileName() << endl;
  cout << prefixDA << " : Gain file (to SHUTTLE) : " << shuttleFile << endl;   
  
  filcout.close();
  
  // Transferring to calibration file to  FES
  // be sure that env variable DAQDALIB_PATH is set in script file
  //       gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/infoLogger");
  printf("\n *****  STORE calibration FILE to FES ****** \n");
  status = daqDA_FES_storeFile(shuttleFile.Data(),"GAINS");
  if (status) { printf(" Failed to export file : %s , status = %d\n",shuttleFile.Data(),status); return -1; }
  //  else printf(" %s successfully exported to FES  \n",shuttleFile.Data());
  
  printf("\n ######## End execution : %s ######## \n",prefixDA); 
  timers.Stop();
  printf("\nExecution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());
  return status;
}

