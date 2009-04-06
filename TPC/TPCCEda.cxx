/*
TPC DA for online calibration

Contact: Haavard.Helstrup@cern.ch
Link: 
Run Type: PHYSICS STANDALONE DAQ
DA Type: MON
Number of events needed: 500
Input Files: 
Output Files: tpcCE.root, to be exported to the DAQ FXS
fileId:   CE
Trigger types used: PHYSICS_EVENT

*/

/*

TPCCEda.cxx - calibration algorithm for TPC Central Electrode events

10/06/2007  sylvain.chapeland@cern.ch :  first version - clean skeleton based on DAQ DA case1
06/12/2007  haavard.helstrup@cern.ch  :  created CE DA based on pulser code
19/09/2008  J.Wiechula@gsi.de:           Added support for configuration files.

contact: marian.ivanov@cern.ch


This process reads RAW data from the files provided as command line arguments
and save results in a file (named from RESULT_FILE define - see below).

*/

#define RESULT_FILE "tpcCE.root"
#define FILE_ID "CE"
#define MAPPING_FILE "tpcMapping.root"
#define CONFIG_FILE "TPCCEda.conf"


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
#include "TSystem.h"
#include "TString.h"
#include "TObjString.h"
#include "TDatime.h"
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
#include "AliTPCConfigDA.h"
//
//AMORE
//
#include <AmoreDA.h>
//
// TPC calibration algorithm includes
//
#include "AliTPCCalibCE.h"




/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {
  /* log start of process */
  printf("TPC CE DA started - %s\n",__FILE__);

  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
  
  AliLog::SetClassDebugLevel("AliTPCRawStream",-5);
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);

 gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                         "*",
                                         "TStreamerInfo",
                                         "RIO",
                                         "TStreamerInfo()");

 /* declare monitoring program */
 int i,status;
 status=monitorDeclareMp( __FILE__ );
 if (status!=0) {
   printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
   return -1;
 }
 monitorSetNowait();
 monitorSetNoWaitNetworkTimeout(1000);

  //variables 
  AliTPCmapper *mapping = 0;   // The TPC mapping
  char localfile[255];
  unsigned long32 runNb=0;      //run number
  // configuration options 
  Bool_t fastDecoding = kFALSE;
 
  // if  test setup get parameters from $DAQDA_TEST_DIR 
  if (!mapping){
    /* copy locally the mapping file from daq detector config db */
    sprintf(localfile,"./%s",MAPPING_FILE);
    status = daqDA_DB_getFile(MAPPING_FILE,localfile);
    if (status) {
      printf("Failed to get mapping file (%s) from DAQdetDB, status=%d\n", MAPPING_FILE, status);
      return -1;
    }

    /* open the mapping file and retrieve mapping object */
    TFile *fileMapping = new TFile(MAPPING_FILE, "read");
    mapping = (AliTPCmapper*) fileMapping->Get("tpcMapping");
    delete fileMapping;
  }

  if (mapping == 0) {
    printf("Failed to get mapping object from %s.  ...\n", MAPPING_FILE);
    return -1;
  } else {
    printf("Got mapping object from %s\n", MAPPING_FILE);
  }
  
  //
  // DA configuration from configuration file
  //
 //retrieve configuration file
  sprintf(localfile,"./%s",CONFIG_FILE);
  status = daqDA_DB_getFile(CONFIG_FILE,localfile);
  if (status) {
    printf("Failed to get configuration file (%s) from DAQdetDB, status=%d\n", CONFIG_FILE, status);
    return -1;
  }
  AliTPCConfigDA config(CONFIG_FILE);
  // check configuration
  if ( (Int_t)config.GetValue("UseFastDecoder") == 1 ){
    printf("Info: The fast decoder will be used for the processing.\n");
    fastDecoding=kTRUE;
  }

  //create calibration object 
  AliTPCCalibCE calibCE(config.GetConfigurationMap());   // central electrode calibration
  calibCE.SetAltroMapping(mapping->GetAltroMapping()); // Use altro mapping we got from daqDetDb

  //===========================//
  // loop over RAW data files //
  //==========================//
  int nevents=0;
  for ( i=1; i<argc; i++) {

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
      
      /* skip start/end of run events */
      if ( (event->eventType != physicsEvent) && (event->eventType != calibrationEvent) )
        continue;

      /* retry if got no event */
      if (event==NULL)
        continue;
      
      nevents++;
      // get the run number
      runNb = event->eventRunNb;
      //  CE calibration
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);
      if ( fastDecoding ) calibCE.ProcessEventFast(rawReader);
      else calibCE.ProcessEvent(rawReader);
      delete rawReader;

      /* free resources */
      free(event);
    }
  }

  //
  // Analyse CE data and write them to rootfile
  //
  calibCE.Analyse();
  printf ("%d events processed\n",nevents);

  TFile * fileTPC = new TFile (RESULT_FILE,"recreate");
  calibCE.Write("tpcCalibCE");
  delete fileTPC;
  printf("Wrote %s\n",RESULT_FILE);

  /* store the result file on FES */

  status=daqDA_FES_storeFile(RESULT_FILE,FILE_ID);
  if (status) {
    status = -2;
  }

  //
  //Send objects to the AMORE DB
  //
  printf ("AMORE part\n");
  const char *amoreDANameorig=gSystem->Getenv("AMORE_DA_NAME");
  //cheet a little -- temporary solution (hopefully)
  //
  //currently amoreDA uses the environment variable AMORE_DA_NAME to create the mysql
  //table in which the calib objects are stored. This table is dropped each time AmoreDA
  //is initialised. This of course makes a problem if we would like to store different
  //calibration entries in the AMORE DB. Therefore in each DA which writes to the AMORE DB
  //the AMORE_DA_NAME env variable is overwritten.  
  gSystem->Setenv("AMORE_DA_NAME",Form("TPC-%s",FILE_ID));
  //
  // end cheet
  TDatime time;
  TObjString info(Form("Run: %u; Date: %s",runNb,time.AsString()));
  amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
  Int_t statusDA=0;
  statusDA+=amoreDA.Send("CET0",calibCE.GetCalPadT0());
  statusDA+=amoreDA.Send("CEQ",calibCE.GetCalPadQ());
  statusDA+=amoreDA.Send("CERMS",calibCE.GetCalPadRMS());
  statusDA+=amoreDA.Send("Info",&info);
  if ( statusDA!=0 )
    printf("Waring: Failed to write one of the calib objects to the AMORE database\n");
  if (amoreDANameorig) gSystem->Setenv("AMORE_DA_NAME",amoreDANameorig);

  return status;
}
