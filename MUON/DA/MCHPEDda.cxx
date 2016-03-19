/*
  Contact: Jean-Luc Charvet <jean-luc.charvet@cern.ch>  
  Link: http://aliceinfo.cern.ch/static/Offline/dimuon/muon_html/README_mchda.html
  Reference Run: 109302 (station 3 only)
  Run Type: PEDESTAL
  DA Type: LDC
  Number of events needed: 400 events for pedestal run
  Input Files: mutrkpedvalues and config_ldc-MTRK-S3-0 in path : /afs/cern.ch/user/j/jcharvet/public/DA_validation 
  Output Files: local dir (not persistent) -> MCHPEDda.ped  FXS -> run<#>_MCH_<ldc>_PEDESTALS
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
        2016-03-15 New version: MCHPEDda.cxx,v 1.10
	-------------------------------------------------------------------------

	Version for MCHPEDda MUON tracking
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
  Int_t status=0, status1=0;
  TStopwatch timers;
  timers.Start(kTRUE); 

  // Needed for streamer application
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
  Int_t errorDetail  = 1;
  Int_t checkTokenlost = 0;
  Int_t flag_histo=1;
  TString inputFile;

  Int_t  nDateEvents = 0;
  Int_t  nDateRejEvents = 0;
  Int_t nGlitchErrors= 0;
  Int_t nParityErrors= 0;
  Int_t nPaddingErrors= 0;
  Int_t nTokenlostErrors= 0;
  Int_t recoverParityErrors = 1;

  TString logOutputFile;

  Char_t flatFile[256]="";
  Char_t* detail;
  TString shuttleFile;

  Int_t nEventsRecovered = 0;
  Int_t nEvents = 0;
  UInt_t runNumber   = 0;
  Int_t nConfig = 1;
  Int_t nEvthreshold = 100; //below this nb_evt pedestal are not calculated and forced to 4085 (sigma)
  Int_t nSorting = 1 ; // pedestal sorting (ON by default)
  Int_t statusDA = 0 ; // DA return code 
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
      nConfig=0;
      //prefixLDC = "LDC"
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
	case 'v' :
	  i++; 
	  sscanf(argv[i],"%d",&errorDetail);
	  break;
	case 't' :
	  i++; 
	  sscanf(argv[i],"%d",&checkTokenlost);
	  break;
	case 'h' :
	  i++;
	  printf("\n******************* %s usage **********************",argv[0]);
	  printf("\nOnline (called from ECS) : %s <raw data file> (no inline options)\n",argv[0]);
	  printf("\n%s can be used locally only with options (without DiMuon configuration file: nConfig=%d)",argv[0],nConfig);
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
	  printf("\n-v <error detail>          (default = %d : 0=low 1=medium 2=high)",errorDetail);
	  printf("\n-t <token lost check>      (default = %d : 0=no , 1=yes)",checkTokenlost);

	  printf("\n\n");
	  exit(-1);
	default :
	  printf("%s : bad argument %s (please check %s -h)\n",argv[0],argv[i],argv[0]);
	  argc = 2; exit(-1); // exit if error
	} // end of switch  
    } // end of for i  

  UShort_t manuId;  
  UChar_t channelId;
  UShort_t charge;

  Char_t dum[256] ="";
  sprintf(dum,"%sPEDda",getenv("DATE_DETECTOR_CODE")); //  DATE_DETECTOR_CODE = MCH 
  const char* prefixDA = dum ; // program prefix
  const char* prefixLDC = getenv("DATE_ROLE_NAME"); // LDC name
  if(prefixLDC == NULL)  prefixLDC ="MCH" ;
  printf("%s : -------- Begin execution : %s --------  \n",prefixLDC,prefixDA); 

  //Pedestal object
  AliMUONPedestal* muonPedestal = new AliMUONPedestal();
  muonPedestal->SetprefixDA(prefixDA);
  muonPedestal->SetprefixLDC(prefixLDC);
  muonPedestal->SetStatusDA(statusDA);

  // Output log file initialisations
  sprintf(flatFile,"%s.log",prefixDA);
  logOutputFile=flatFile;
  AliLog::SetStreamOutput(&filcout); // Print details on logfile
  filcout.open(logOutputFile.Data());
  filcout<<"//=================================================" << endl;
  filcout<<"//" << prefixLDC << "       " << prefixDA  << endl;
  filcout<<"//=================================================" << endl;
  filcout<<"//  * Date  : " << muonPedestal->GetDate()->AsString("l") << "\n" << endl;

  muonPedestal->SetAlifilcout(&filcout);
  cout<<prefixLDC << " :  Date: " << muonPedestal->GetDate()->AsString("l") << "\n" << endl;


  Char_t dbfile[256]="";
  // nConfig=1 : Reading configuration (or not) status via "mutrkpedvalues" file located in DetDB
  if(nConfig)
    { 
      Int_t flag_hist, nEvthres,maxEvts;
      char line[80];
      sprintf(dbfile,"mutrkpedvalues");
      status=daqDA_DB_getFile(dbfile,dbfile);
      if(status) {fprintf(stderr," !!! ERROR  : input file %s is missing, status = %d\n",dbfile,status);
	fprintf(stderr,"\n%s : -------- %s ending in ERROR !!!! -------- (status= %d) \n",prefixLDC,prefixDA,-1);
	return -1; }
 
      ifstream filein(dbfile,ios::in);
      filein >> line >> nConfig ; cout << "mutrkpedvalues: " << line << nConfig << "  "  ;
      filein >> line >> nEvthres ; if(nEvthres !=0)nEvthreshold=nEvthres;  cout << line << nEvthreshold << "  " ; 
      filein >> line >> checkTokenlost ; cout << line << checkTokenlost << "  " ;
      filein >> line >> flag_histo ;  cout << line << flag_histo << "  " ; // =0 , 1-> standard , 2-> ntuple charge
      filein >> line >> maxEvts ;  if(maxEvts!=0){maxEvents=maxEvts;  cout << line << maxEvents;}
      //     filein >> line >> errorDetail ; cout << line << errorDetail << "  " ;
      //     filein >> line >> nSorting ; cout << line << nSorting << "  " ;
     cout << endl;
    }
  if(!nConfig)  printf("%s : ***  Config= %d: no configuration ascii file is used \n",prefixLDC,nConfig); 
  muonPedestal->SetconfigDA(nConfig);
  muonPedestal->SetnEvthreshold(nEvthreshold);
  muonPedestal->SetnSorting(nSorting);
  muonPedestal->SetHistos(flag_histo);

  // nConfig=1: configuration ascii file config_$DATE_ROLE_NAME read from DetDB
  if(nConfig)
    {
      sprintf(dbfile,"config_%s",getenv("DATE_ROLE_NAME"));
      status=daqDA_DB_getFile(dbfile,dbfile);
      if(status) {fprintf(stderr," !!! ERROR  : Configuration file %s is missing, status = %d\n",dbfile,status);
	fprintf(stderr,"\n%s : -------- %s ending in ERROR !!!! -------- (status= %d) \n",prefixLDC,prefixDA,-1);
	return -1; }
      muonPedestal->LoadConfig(dbfile);  
    } 

    if(flag_histo)muonPedestal->CreateControlHistos(); /// Generate pedestal histo rootfile

  // Rawdeader, RawStreamHP
  AliRawReader* rawReader = AliRawReader::Create(inputFile);
  AliMUONRawStreamTrackerHP* rawStream  = new AliMUONRawStreamTrackerHP(rawReader);    
//	rawStream->DisableWarnings();
  rawStream->EnabbleErrorLogger();
 
  switch (errorDetail)
    {
    case 0: rawStream->SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kLowErrorDetail); break;/// Logs minimal information in the error messages.
    case 1: rawStream->SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kMediumErrorDetail); break;/// Logs a medium level of detail in the error messages.
    case 2: rawStream->SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kHighErrorDetail); break;/// Logs a medium level of detail in the error messages.
    default: rawStream->SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kMediumErrorDetail); break;
    }

  printf("\n%s : Reading data from file %s\n",prefixLDC,inputFile.Data());

  Int_t tabTokenError[20][14];
  for ( Int_t i=0 ; i<20 ; i++) { for ( Int_t j=0 ; j<14 ; j++) { tabTokenError[i][j]=0;}	}

  while (rawReader->NextEvent())
    {
      Int_t eventType = rawReader->GetType();
      runNumber = rawReader->GetRunNumber();
      if(nDateEvents==0)  { filcout<<"//  ---->  RUN = " << runNumber << "\n" << endl;}
      if (nDateEvents >= maxDateEvents) break;
      if (nEvents >= maxEvents) break;
      if (nDateEvents>0 &&  nDateEvents % 100 == 0) 	
	cout<< prefixLDC << " :  DATE events = " << nDateEvents << "   Used events = " << nEvents << endl;

      // check shutdown condition 
      if (daqDA_checkShutdown())  break;
      //Skip events
      while (skipEvents)
	{
	  rawReader->NextEvent();
	  skipEvents--;
	}
      nDateEvents++;
      if (eventType != PHYSICS_EVENT)
	continue; // for the moment

      //      const char* detail = "";
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
//**	  if (rawStream->GetDDL()== 8 ) { 

	  if (rawStream->IsErrorMessage()) eventIsErrorMessage = kTRUE;
	  eventGlitchErrors += rawStream->GetGlitchErrors();
	  eventParityErrors += rawStream->GetParityErrors();
	  eventPaddingErrors += rawStream->GetPaddingErrors();
	  eventTokenlostErrors += rawStream->GetTokenLostErrors();
	  //	  cout << nDateEvents << "    "  << rawStream->GetDDL() << "   "  << eventTokenlostErrors << endl;
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
 //**	  }  //  if (rawStream->GetDDL()== 8 ) {
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
	  if (recoverParityErrors && eventParityErrors && !eventGlitchErrors && !eventPaddingErrors && !(eventTokenlostErrors && checkTokenlost))
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
	      nDateRejEvents++;
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

    } // while (rawReader->NextEvent())
  delete rawReader;
  delete rawStream;

  sprintf(flatFile,"%s.ped",prefixDA);
  if(shuttleFile.IsNull())shuttleFile=flatFile;
  muonPedestal->SetAliNEvents(nEvents);
  muonPedestal->SetAliRunNumber(runNumber);
  
  muonPedestal->Finalize();
  status = muonPedestal->GetStatusDA()  ; 
 
  // writing some counters
  detail=Form("\n%s : Nb of DATE events           = %d",prefixLDC,nDateEvents) ;                             cout << detail; filcout << detail ;
  detail=Form("\n%s : Nb of Glitch errors         = %d",prefixLDC,nGlitchErrors) ;                           cout << detail; filcout << detail ;
  detail=Form("\n%s : Nb of Parity errors         = %d",prefixLDC,nParityErrors) ;                           cout << detail; filcout << detail ;
  detail=Form("\n%s : Nb of Token lost errors     = %d",prefixLDC,nTokenlostErrors) ;                        cout << detail; filcout << detail ;
  detail=Form("\n%s : Nb of Rejected DATE events  = %d",prefixLDC,nDateRejEvents) ;                          cout << detail; filcout << detail ;
  detail=Form("\n%s : Nb of recovered events      = %d",prefixLDC,nEventsRecovered) ;                        cout << detail; filcout << detail ;
  detail=Form("\n%s : Nb of events without errors = %d",prefixLDC,nEvents-nEventsRecovered) ;                cout << detail; filcout << detail ;
  detail=Form("\n%s : Nb of used events           = %d (threshold= %d)\n\n",prefixLDC,nEvents,nEvthreshold); cout << detail; filcout << detail ;

  //2015-01-22  create filerror.txt file for checking (should be removed later) 
  //2015-01-22  FILE *fp;
  //2015-01-22  fp= fopen("filerror.txt","w+");
  // Writing Token Error table
  if(nTokenlostErrors)
    {
      detail=Form("%s : Error : Token Lost occurence \n",prefixLDC);
      filcout <<  detail ;
      //2015-01-22      printf("%s",detail);
      //2015-01-22      cerr << "ERROR :  " << detail;
      fprintf(stderr,"%s",detail);
      //2015-01-22      fprintf(fp,"%s",detail);
      for ( Int_t i=0 ; i<20 ; i++) 
	{ 
	  for ( Int_t j=4 ; j<14 ; j++) 
	    { 
	      if(tabTokenError[i][j]>0)
		{
		  Int_t tab=tabTokenError[i][j];
		  Int_t frt=j/2-1;
		  Int_t station = i/4 +1;
		  if( j % 2 == 0)detail=Form("%s : Error : in DDL= %d (station %d) and FRT%d ( Up ) => %d Token errors (address = 0x%X0000)",prefixLDC,2560+i,station,frt,tab,j);
		  else detail=Form("%s : Error : in DDL= %d (station %d) and FRT%d (Down) => %d Token errors (address = 0x%X0000)",prefixLDC,2560+i,station,frt,tab,j);
		  //2015-01-22		  printf("%s\n",detail);
		  fprintf(stderr,"%s\n",detail);
		  filcout <<  detail << endl;
		}
	    }
	}
    }
  //2015-01-22  fclose(fp);

  if (!shuttleFile.IsNull())  
    {
      ofstream out(shuttleFile.Data());  
      muonPedestal->MakeASCIIoutput(out); /// Generate pedestal output file
      out.close();
      detail=Form("%s : Pedestal file (to SHUTTLE) : %s\n",prefixLDC,shuttleFile.Data());
      cout << detail; filcout << detail  ;        
    }
  if(flag_histo) /// Generate pedestal histo rootfile
    {
      muonPedestal->MakeControlHistos(); 
      detail=Form("%s : Pedestal Histo file        : %s\n",prefixLDC,muonPedestal->GetHistoFileName());   
      cout << detail; filcout << detail ; 
     }       
  // .log files
  detail=Form("%s : Output logfile             : %s\n",prefixLDC,logOutputFile.Data());
  cout << detail; filcout << detail ;        


   // Transferring pedestal file to FES  (be sure that env variable DAQDALIB_PATH is set)
  cout << endl; 
  status1 = daqDA_FES_storeFile(shuttleFile.Data(),"PEDESTALS");
  if (status1) { detail=Form("%s: !!! ERROR: Failed to export pedestal file : %s to FES \n",prefixLDC,shuttleFile.Data()); 
    //2015-01-22    printf("%s",detail); filcout << detail ; status= -1; }
    fprintf(stderr,"%s",detail); filcout << detail ; status= -1; }

  // Transferring configuration file to FES  (be sure that env variable DAQDALIB_PATH is set)
  if(nConfig) 
    { cout << endl; 
      status1 = daqDA_FES_storeFile(dbfile,"CONFIG");
      if (status1) { detail=Form("%s: !!! ERROR: Failed to export configuration file : %s to FES \n",prefixLDC,dbfile); 
	//2015-01-22	printf("%s",detail); filcout << detail ; status=-1; }
	fprintf(stderr,"%s",detail); filcout << detail ; status=-1; }
    }

  filcout.close();

 // Copying files to local DB folder defined by DAQ_DETDB_LOCAL
  Char_t *dir;
  unsigned int nLastVersions=50;
  dir= getenv("DAQ_DETDB_LOCAL");
  if(dir != NULL)  {
    //    unsigned int nLastVersions=50;
    printf("\n%s : ---  Local DataBase: %s (Max= %d) ---\n",prefixLDC,dir,nLastVersions);
      if(!shuttleFile.IsNull())status1 = daqDA_localDB_storeFile(shuttleFile.Data(),nLastVersions);
      if(flag_histo)status1 = daqDA_localDB_storeFile(muonPedestal->GetHistoFileName(),nLastVersions);
      status1 = daqDA_localDB_storeFile(logOutputFile.Data(),nLastVersions);
  }

    cout << " " << endl; 
    
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
      {cout << prefixLDC << " :  !!! ERROR: Failed to write Pedestals in the AMORE database : " << amoreStatus << endl ; status=-1 ;}
    else 
      cout << prefixLDC << " : amoreDA.Send(Pedestals) ok" << endl;  
#else
    cout << prefixLDC << " : Warning: MCH DA not compiled with AMORE support" << endl;
#endif
    

  if(!status)printf("\n%s : -------- End execution : %s -------- (status= %d) \n",prefixLDC,prefixDA,status);
  else { fprintf(stderr,"\n%s : -------- %s ending in ERROR !!!! -------- (status= %d) \n",prefixLDC,prefixDA,status);}
  timers.Stop();
  printf("\n Execution time : R:%7.2fs C:%7.2fs\n",timers.RealTime(), timers.CpuTime());
  return status;
}
