/*
TPC DA for online calibration

Contact: Haavard.Helstrup@cern.ch
Link: 
Run Type: PHYSICS_RUN STANDALONE DAQ
DA Type: MON
Number of events needed: 500
Input Files: 
Output Files: tpcCE.root, to be exported to the DAQ FXS
Trigger types used: PHYSICS_EVENT

*/

/*

TPCCEda.cxx - calibration algorithm for TPC Central Electrode events

10/06/2007  sylvain.chapeland@cern.ch :  first version - clean skeleton based on DAQ DA case1
06/12/2007  haavard.helstrup@cern.ch  :  created CE DA based on pulser code

contact: marian.ivanov@cern.ch


This process reads RAW data from the files provided as command line arguments
and save results in a file (named from RESULT_FILE define - see below).

*/

#define RESULT_FILE "tpcCE.root"
#define MAPPING_FILE "tpcMapping.root"


#include <daqDA.h>
#include "event.h"
#include "monitor.h"
#include <stdio.h>
#include <stdlib.h>

//
//Root includes
//
#include <TFile.h>
#include "TROOT.h"
#include "TPluginManager.h"
//
//AliRoot includes
//
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliTPCmapper.h"
#include "AliTPCRawStream.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliMathBase.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include "TSystem.h"

//
// TPC calibration algorithm includes
//
#include "AliTPCCalibCE.h"




/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {

 gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                         "*",
                                         "TStreamerInfo",
                                         "RIO",
                                         "TStreamerInfo()");

  AliLog::SetClassDebugLevel("AliTPCRawStream",AliLog::kFatal);
  AliLog::SetClassDebugLevel("AliRawReaderDate",AliLog::kFatal);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",AliLog::kFatal);
  AliLog::SetModuleDebugLevel("TPC",AliLog::kFatal);
  AliLog::SetModuleDebugLevel("RAW",AliLog::kFatal);

  int i,status;
  AliTPCmapper *mapping = 0;   // The TPC mapping
  // if  test setup get parameters from $DAQDA_TEST_DIR 
  
  if (!mapping){
    /* copy locally the mapping file from daq detector config db */
    status = daqDA_DB_getFile(MAPPING_FILE,"./tpcMapping.root");
    if (status) {
      printf("Failed to get mapping file (%s) from DAQdetDB, status=%d\n", MAPPING_FILE, status);
      printf("Continue anyway ... maybe it works?\n");              // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return -1;   // temporarily uncommented for testing on pcald47 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    /* open the mapping file and retrieve mapping object */
    TFile *fileMapping = new TFile(MAPPING_FILE, "read");
    mapping = (AliTPCmapper*) fileMapping->Get("tpcMapping");
    delete fileMapping;
  }

  if (mapping == 0) {
    printf("Failed to get mapping object from %s.  ...\n", MAPPING_FILE);
    //return -1;
  } else {
    printf("Got mapping object from %s\n", MAPPING_FILE);
  }


  AliTPCCalibCE calibCE;   // pedestal and noise calibration

  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  /* log start of process */
  printf("TPC CE DA started - %s\n",__FILE__);


  /* set time bin range */
  calibCE.SetRangeTime(800,940);
  calibCE.SetAltroMapping(mapping->GetAltroMapping()); // Use altro mapping we got from daqDetDb

  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* loop over RAW data files */
  int nevents=0;
  for(i=1;i<argc;i++) {

    /* define data source : this is argument i */
    printf("Processing file %s\n", argv[i]);
    status=monitorSetDataSource( argv[i] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }

    /* read until EOF */
    while (true) {
      struct eventHeaderStruct *event;

      /* check shutdown condition */
      if (daqDA_checkShutdown()) {break;}

      /* get next event (blocking call until timeout) */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF) {
        printf ("End of File %d detected\n",i);
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
      nevents++;

      //  Pulser calibration

      AliRawReader *rawReader = new AliRawReaderDate((void*)event);
      calibCE.ProcessEvent(rawReader);
      delete rawReader;

      /* free resources */
      free(event);
    }
  }

  calibCE.Analyse(); 
  printf ("%d events processed\n",nevents);

  TFile * fileTPC = new TFile (RESULT_FILE,"recreate");
  calibCE.Write("calibCE");
  delete fileTPC;
  printf("Wrote %s\n",RESULT_FILE);

  /* store the result file on FES */

  status=daqDA_FES_storeFile(RESULT_FILE,RESULT_FILE);
  if (status) {
    status = -2;
  }

  return status;
}
