/*
  EMCAL DA for online calibration: for pedestal studies
  
  Contact: silvermy@ornl.gov
  Run Type: PEDESTAL 
  DA Type: LDC 
  Number of events needed: ~1000
  Frequency: expect to run this once per day initially 
             (maybe less frequent later)
  Input Files: argument list
  Output Files: RESULT_FILE=EMCALPED.root, to be exported to the DAQ FXS
  fileId:  FILE_ID=EMCALPED    

  Trigger types used: at least for now, all events in special PEDESTAL run
  (physics, calibration, systemSoftwareTrigger, detectorSoftwareTrigger)
  [When we have real data files later, we may restrict this further]
*/
/*
  This process reads RAW data from the files provided as command line arguments
  and save results (class itself) in a file (named from RESULT_FILE define 
  - see below).
*/

#define RESULT_FILE  "EMCALPED.root"
#define FILE_ID "EMCALPED"
#define AliDebugLevel() -1
#define FILE_PEDClassName "emcCalibPedestal"
const int kNRCU = 4;
/* LOCAL_DEBUG is used to bypass daq* calls, for local testing */
//#define LOCAL_DEBUG 1 // comment out to run normally

extern "C" {
#include <daqDA.h>
}
#include "event.h" /* in $DATE_COMMON_DEFS/; includes definition of event types */
#include "monitor.h" /* in $DATE_MONITOR_DIR/; monitor* interfaces */

#include "stdio.h"
#include "stdlib.h"

// ROOT includes
#include <TFile.h> 
#include <TROOT.h> 
#include <TPluginManager.h> 
#include <TSystem.h> 

//
//AliRoot includes
//
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawEventHeaderBase.h"
#include "AliCaloRawStreamV3.h"
#include "AliCaloAltroMapping.h"
#include "AliLog.h"

//
// EMC calibration-helper algorithm includes
//
#include "AliCaloCalibPedestal.h"

/*
  Main routine, EMC pedestal detector algorithm 
  Arguments: list of DATE raw data files
*/

int main(int argc, char **argv) {

  AliLog::SetClassDebugLevel("AliCaloRawStreamV3",-5);
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);

  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;  
  }

  /* magic line - for TStreamerInfo */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()"); 

  int i, status;

  /* log start of process */
  printf("EMCAL DA started - %s\n",__FILE__);

  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }
  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);

  /* Retrieve mapping files from DAQ DB */ 
  const char* mapFiles[kNRCU] = {"RCU0A.data","RCU1A.data","RCU0C.data","RCU1C.data"};

  for(Int_t iFile=0; iFile<kNRCU; iFile++) {
    int failed = daqDA_DB_getFile(mapFiles[iFile], mapFiles[iFile]);
    if(failed) { 
      printf("Cannot retrieve file %d : %s from DAQ DB. Exit now..\n",
	     iFile, mapFiles[iFile]);
#ifdef LOCAL_DEBUG
#else
            return -1;
#endif
    }
  }
  
  /* Open mapping files */
  AliCaloAltroMapping *mapping[kNRCU];
  TString path = "./";
  path += "RCU";
  TString path2;
 TString side[] = {"A","C"};//+ and - pseudarapidity supermodules
 for(Int_t j = 0; j < 2; j++){
   for(Int_t i = 0; i < 2; i++) {
     path2 = path;
     path2 += i;
     path2 +=side[j]; 
     path2 += ".data";
     mapping[j*2 + i] = new AliCaloAltroMapping(path2.Data());
   }
 }
  /* set up our analysis class */  
  AliCaloCalibPedestal * calibPedestal = new AliCaloCalibPedestal(AliCaloCalibPedestal::kEmCal); // pedestal and noise calibration
  calibPedestal->SetAltroMapping( mapping );

  AliRawReader *rawReader = NULL;
  int nevents=0;

  /* loop over RAW data files */
  for ( i=1; i<argc; i++ ) {

    /* define data source : this is argument i */
    printf("Processing file %s\n", argv[i]);
    status=monitorSetDataSource( argv[i] );
    if (status!=0) {
      printf("monitorSetDataSource() failed. Error=%s. Exiting ...\n", monitorDecodeError(status));
      return -1;
    }

    /* read until EOF */
    struct eventHeaderStruct *event;
    eventTypeType eventT;

    for ( ; ; ) { // infinite loop

      /* check shutdown condition */
      if (daqDA_checkShutdown()) {break;}

      /* get next event (blocking call until timeout) */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF) {
	printf ("End of File %d (%s) detected\n", i, argv[i]);
	break; /* end of monitoring file has been reached */
      }
      if (status!=0) {
	printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
	break;
      }

      /* retry if got no event */
      if (event==NULL) {
	continue;
      }
      eventT = event->eventType; /* just convenient shorthand */

      /* skip start/end of run events */
      if ( (eventT != physicsEvent) && (eventT != calibrationEvent) &&
	   (eventT != systemSoftwareTriggerEvent) && (eventT != detectorSoftwareTriggerEvent) ) {
	continue;
      }

      nevents++; // count how many acceptable events we have

      //  Pedestal calibration
      rawReader = new AliRawReaderDate((void*)event);
      calibPedestal->ProcessEvent(rawReader);
      delete rawReader;

      /* free resources */
      free(event);    

    } //until EOF
  } // loop over files

  //
  // write class to rootfile
  //

  printf ("%d physics/calibration events processed.\n",nevents);

  TFile f(RESULT_FILE, "recreate");
  if (!f.IsZombie()) { 
    f.cd();
    calibPedestal->Write(FILE_PEDClassName);
    f.Close();
    printf("Object saved to file \"%s\" as \"%s\".\n", 
	   RESULT_FILE, FILE_PEDClassName); 
  } 
  else {
    printf("Could not save the object to file \"%s\".\n", 
	   RESULT_FILE);
  }

  //
  // closing down; see if we can delete our analysis helper also
  //
  delete calibPedestal;
  for(Int_t iFile=0; iFile<kNRCU; iFile++) {
    if (mapping[iFile]) delete mapping[iFile];
  }

  /* store the result file on FES */
#ifdef LOCAL_DEBUG
#else
  status = daqDA_FES_storeFile(RESULT_FILE, FILE_ID);
#endif

  return status;
}
