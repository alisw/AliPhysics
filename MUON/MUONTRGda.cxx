/*

Version 1 for MUONTRGda MUON trigger
Working version for reading back raw data
(Ch. Finck)

*/
extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include "AliMUONRawStreamTrigger.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDDLTrigger.h"
// #include "AliMUONVStore.h"
// #include "AliMUON2DMap.h"
// #include "AliMUONCalibParamNF.h"
// #include "AliMpDDLStore.h"
// #include "AliMpIntPair.h"
#include "AliMpConstants.h"
#include "AliRawReaderDate.h"


//ROOT

#include "TString.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "TROOT.h"
#include "TPluginManager.h"


// global variables
TString command("pat");
UInt_t runNumber = 0;
Int_t nEvents = 0;


//*************************************************************//

// main routine
int main(Int_t argc, Char_t **argv) 
{
  
    // needed for streamer application
    gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					  "*",
					  "TStreamerInfo",
					  "RIO",
					  "TStreamerInfo()"); 


    Int_t printLevel = 0;  // Global variable defined as extern in the others .cxx files
    Int_t skipEvents = 0;
    Int_t maxEvents  = 1000000;
    Char_t inputFile[256];
    TString flatOutputFile;
    TString configFile;

// option handler

   // decode the input line
  for (Int_t i = 1; i < argc; i++) // argument 0 is the executable name
  {
      Char_t* arg;
      
      arg = argv[i];
      if (arg[0] != '-') continue;
      switch (arg[1])
	{
	case 'f' : 
	  i++;
	  sprintf(inputFile,argv[i]);
	  break;
	case 'a' : 
	  i++;
	  flatOutputFile = argv[i];
	  break;
	case 'c' : 
	  i++;
	  configFile = argv[i];
	  break;
	case 'e' : 
	  i++;
	  command = argv[i];
	  break;
	case 'd' :
	  i++; 
	  printLevel=atoi(argv[i]);
	  break;
	case 's' :
	  i++; 
	  skipEvents=atoi(argv[i]);
	  break;
	case 'n' :
	  i++; 
	  sscanf(argv[i],"%d",&maxEvents);
	  break;
	case 'h' :
	  i++;
	  printf("\n******************* %s usage **********************",argv[0]);
	  printf("\n%s -options, the available options are :",argv[0]);
	  printf("\n-h help                   (this screen)");
	  printf("\n");
	  printf("\n Input");
	  printf("\n-f <raw data file>        (default = %s)",inputFile); 
	  printf("\n");
	  printf("\n Output");
	  printf("\n-a <Flat ASCII file>      (default = %s)",flatOutputFile.Data()); 
	  printf("\n");
	  printf("\n Options");
	  printf("\n-d <print level>          (default = %d)",printLevel);
	  printf("\n-s <skip events>          (default = %d)",skipEvents);
	  printf("\n-n <max events>           (default = %d)",maxEvents);
	  printf("\n-e <execute pattern/scaler>     (default = %s)",command.Data());

	  printf("\n\n");
	  exit(-1);
	default :
	  printf("%s : bad argument %s (please check %s -h)\n",argv[0],argv[i],argv[0]);
	  argc = 2; exit(-1); // exit if error
	} // end of switch  
    } // end of for i  

  // set command to lower case
  command.ToLower();

  // decoding the events
  
  Int_t status;
  Int_t nDateEvents = 0;

  void* event;

  // containers
   AliMUONDDLTrigger*       ddlTrigger  = 0x0;
   AliMUONDarcHeader*       darcHeader  = 0x0;
   AliMUONRegHeader*        regHeader   = 0x0;
   AliMUONLocalStruct*      localStruct = 0x0;

  TStopwatch timers;

  timers.Start(kTRUE); 

  // once we have a configuration file in db
  // copy locally a file from daq detector config db 
  // The current detector is identified by detector code in variable
  // DATE_DETECTOR_CODE. It must be defined.
  // If environment variable DAQDA_TEST_DIR is defined, files are copied from DAQDA_TEST_DIR
  // instead of the database. The usual environment variables are not needed.
  if (!configFile.IsNull()) {
    status = daqDA_DB_getFile("myconfig", configFile.Data());
    if (status) {
      printf("Failed to get config file : %d\n",status);
      return -1;
    }
  }


  status = monitorSetDataSource(inputFile);
  if (status) {
    cerr << "ERROR : monitorSetDataSource status (hex) = " << hex << status
	      << " " << monitorDecodeError(status) << endl;
    return -1;
  }
  status = monitorDeclareMp("MUON Tracking monitoring");
  if (status) {
    cerr << "ERROR : monitorDeclareMp status (hex) = " << hex << status
	      << " " << monitorDecodeError(status) << endl;
    return -1;
  }

  cout << "MUONTRKda : Reading data from file " << inputFile <<endl;

  while(1) 
  {
    if (nEvents >= maxEvents) break;
    if (nEvents && nEvents % 100 == 0) 	
	cout<<"Cumulated events " << nEvents << endl;

    // check shutdown condition 
    if (daqDA_checkShutdown()) 
	break;

    // Skip Events if needed
    while (skipEvents) {
      status = monitorGetEventDynamic(&event);
      skipEvents--;
    }

    // starts reading
    status = monitorGetEventDynamic(&event);
    if (status < 0)  {
      cout<<"EOF found"<<endl;
      break;
    }

    nDateEvents++;

    // decoding rawdata headers
    AliRawReader *rawReader = new AliRawReaderDate(event);
 
    Int_t eventType = rawReader->GetType();
    runNumber = rawReader->GetRunNumber();
    

    if (eventType != PHYSICS_EVENT)
	continue; // for the moment

    nEvents++;
    if (printLevel) printf("\nEvent # %d\n",nEvents);

   // decoding MUON payload
    AliMUONRawStreamTrigger* rawStream  = new AliMUONRawStreamTrigger(rawReader);
    //rawStream->SetMaxReg(1);

   // loops over DDL 
    while((status = rawStream->NextDDL())) {

       if (printLevel) printf("iDDL %d\n", rawStream->GetDDL());

       ddlTrigger = rawStream->GetDDLTrigger();
       darcHeader = ddlTrigger->GetDarcHeader();

       if (printLevel) printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());

       // loop over regional structures
       Int_t nReg = darcHeader->GetRegHeaderEntries();
       for(Int_t iReg = 0; iReg < nReg; ++iReg){   //REG loop

	 if (printLevel) printf("RegionalId %d\n", iReg);

	 regHeader =  darcHeader->GetRegHeaderEntry(iReg);

	 // loop over local structures
	 Int_t nLocal = regHeader->GetLocalEntries();
	 for(Int_t iLocal = 0; iLocal < nLocal; ++iLocal) {  

	   localStruct = regHeader->GetLocalEntry(iLocal);

	   Int_t loStripX  = (Int_t)localStruct->GetXPos();
	   Int_t loStripY  = (Int_t)localStruct->GetYPos();
	   Int_t loDev     = (Int_t)localStruct->GetXDev();
	   
	   if (printLevel) printf("Index %d, XPos: %d, YPos: %d Dev: %d\n", 
				  localStruct->GetId(), loStripX, loStripY, loDev);

	   if (printLevel) printf("X pattern %x %x %x %x, Y pattern %x %x %x %x\n", 
				  localStruct->GetX1(), localStruct->GetX2(),localStruct->GetX3(),localStruct->GetX4(),
				  localStruct->GetY1(), localStruct->GetY2(),localStruct->GetY3(),localStruct->GetY4());
	   
	 } // iLocal
       } // iReg
     } // NextDDL

    delete rawReader;
    delete rawStream;

  } // while (1)


  timers.Stop();

  cout << "MUONTRKda : Run number                    : " << runNumber << endl;
  cout << "MUONTRKda : Flat ASCII file generated     : " << flatOutputFile << endl;
  // cout << "MUONTRKda : Histo file generated          : " << histoFileName  << endl;
  cout << "MUONTRKda : Nb of DATE events     = "         << nDateEvents    << endl;
  cout << "MUONTRKda : Nb of events used     = "         << nEvents        << endl;

  printf("Execution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

  return status;
}
