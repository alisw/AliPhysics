/*
TPC DA for online calibration

Contact: Haavard.Helstrup@cern.ch
Link:
Run Type: CALIBRATION_PULSER
DA Type: LDC
Number of events needed: 100
Input Files: 
Output Files: tpcPulser.root, to be exported to the DAQ FXS
fileId:   pulser
Trigger types used: CALIBRATION_EVENT

*/

/*

TPCda_pulser.cxx - calibration algorithm for TPC pulser events

10/06/2007  sylvain.chapeland@cern.ch :  first version - clean skeleton based on DAQ DA case1
30/09/2007  haavard.helstrup@cern.ch  :  created pulser DA based on pedestal code

contact: marian.ivanov@cern.ch


This process reads RAW data from the files provided as command line arguments
and save results in a file (named from RESULT_FILE define - see below).

*/

#define RESULT_FILE "tpcPulser.root"
#define FILE_ID "pulser"
#define MAPPING_FILE "tpcMapping.root"
#define AliDebugLevel() -1


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
#include "AliTPCCalibPulser.h"




/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {
  AliLog::SetClassDebugLevel("AliTPCRawStream",-5);
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);

 gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                         "*",
                                         "TStreamerInfo",
                                         "RIO",
                                         "TStreamerInfo()");


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


  AliTPCCalibPulser calibPulser;   // pedestal and noise calibration

  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  /* log start of process */
  printf("TPC Pulser DA started - %s\n",__FILE__);

  /* set time bin range */
  calibPulser.SetRangeTime(400,500);
  calibPulser.SetAltroMapping(mapping->GetAltroMapping()); // Use altro mapping we got from daqDetDb
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
//      calibPulser.ProcessEvent(rawReader);
      calibPulser.ProcessEventFast(rawReader);
      delete rawReader;

      /* free resources */
      free(event);
    }
  }

  calibPulser.Analyse(); 
  printf ("%d events processed\n",nevents);

  TFile * fileTPC = new TFile (RESULT_FILE,"recreate");
  calibPulser.Write("tpcCalibPulser");
  delete fileTPC;
  printf("Wrote %s\n",RESULT_FILE);

  /* store the result file on FES */

  status=daqDA_FES_storeFile(RESULT_FILE,FILE_ID);
  if (status) {
    status = -2;
  }

  return status;
}
