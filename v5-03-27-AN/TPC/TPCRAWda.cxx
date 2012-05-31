/*
TPC DA for online calibration

Contact: Jens.Wiechula@cern.ch
Link:
Run Type: PHYSICS STANDALONE
DA Type: MON
Number of events needed: 200
Input Files: /castor/cern.ch/alice/raw/global/2009/08/22/11/09000080958023.30.root
Output Files: tpcCalibRaw.root, to be exported to the DAQ FXS
fileId:   tpcCalibRaw
Trigger types used: PHYSICS_EVENT

*/

/*

TPCRAWda.cxx - calibration algorithm for L1 phase monitoring and drift velocity from last time bin determination

30/07/2009  Jens.Wiechula@cern.ch:     First implementation.
10/09/2009  Jens.Wiechula@cern.ch:     Add configuration file support. Export object to AMOREdb
                                       after a defined update interval for QA
26/01/2010  Jens.Wiechula@cern.ch:     Exclude laser triggers when running in a global partition

This process reads RAW data from the files provided as command line arguments
and save results in a file (named from RESULT_FILE define - see below).

*/

#define RESULT_FILE "tpcCalibRaw.root"
#define FILE_ID "tpcCalibRaw"
#define MAPPING_FILE "tpcMapping.root"
#define CONFIG_FILE "TPCRAWda.conf"
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
#include <TROOT.h>
#include <TPluginManager.h>
#include <TString.h>
#include <TObjString.h>
#include <TDatime.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TMap.h>
//
//AliRoot includes
//
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliTPCmapper.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliMathBase.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include "TSystem.h"
#include "AliTPCConfigDA.h"
//
//AMORE
//
#include <AmoreDA.h>
//
// TPC calibration algorithm includes
//
#include "AliTPCCalibRaw.h"


//functions, implementation below
void SendToAmoreDB(TObject *o, unsigned long32 runNb);

/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {
  /* log start of process */
  printf("TPCRAWda: DA started - %s\n",__FILE__);

  if (argc<2) {
    printf("TPCRAWda: Wrong number of arguments\n");
    return -1;
  }
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");


  /* declare monitoring program */
  int status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("TPCRAWda: monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }
  //Set network timeout
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);
  
  // variables 
  AliTPCmapper *mapping = 0;   // The TPC mapping
  char localfile[255];
  unsigned long32 runNb=0;      //run number
  // if  test setup get parameters from $DAQDA_TEST_DIR
   
  if (!mapping){
    /* copy locally the mapping file from daq detector config db */
    sprintf(localfile,"./%s",MAPPING_FILE);
    status = daqDA_DB_getFile(MAPPING_FILE,localfile);
    if (status) {
      printf("TPCRAWda: Failed to get mapping file (%s) from DAQdetDB, status=%d\n", MAPPING_FILE, status);
      return -1;
    }

    /* open the mapping file and retrieve mapping object */
    TFile *fileMapping = new TFile(MAPPING_FILE, "read");
    mapping = (AliTPCmapper*) fileMapping->Get("tpcMapping");
    delete fileMapping;
  }

  if (mapping == 0) {
    printf("TPCRAWda: Failed to get mapping object from %s.  ...\n", MAPPING_FILE);
    return -1;
  } else {
    printf("TPCRAWda: Got mapping object from %s\n", MAPPING_FILE);
  }

  //
  // DA configuration from configuration file
  //
  //retrieve configuration file
  sprintf(localfile,"./%s",CONFIG_FILE);
  status = daqDA_DB_getFile(CONFIG_FILE,localfile);
  if (status) {
    printf("TPCRAWda: Failed to get configuration file (%s) from DAQdetDB, status=%d\n", CONFIG_FILE, status);
    return -1;
  }
  AliTPCConfigDA config(CONFIG_FILE);
  
  //set default configuration options
  Double_t updateInterval=30; //seconds
  TString laserTriggerName("C0LSR-ABCE-NOPF-CENT");
  TString forceLaserTriggerId("-1");
  TString monitorAttributes=Form("%d",ATTR_ORIGINAL_EVENT);
  Bool_t  skipAmore=kFALSE;
  
  //amore update interval
  Double_t valConf=config.GetValue("AmoreUpdateInterval");
  if ( valConf>0 ) updateInterval=valConf;

  //laser trigger class name
  if ( config.GetConfigurationMap()->GetValue("LaserTriggerName") ) {
    laserTriggerName=config.GetConfigurationMap()->GetValue("LaserTriggerName")->GetName();
    printf("TPCRAWda: Laser trigger class name set to: %s.\n",laserTriggerName.Data());
  }
  //force laser trigger id
  if ( config.GetConfigurationMap()->GetValue("ForceLaserTriggerId") ) {
    forceLaserTriggerId=config.GetConfigurationMap()->GetValue("ForceLaserTriggerId")->GetName();
    printf("TPCRAWda: Force laser trigger Id: %s.\n",forceLaserTriggerId.Data());
  }
  //skip the amore part
  if ( config.GetConfigurationMap()->GetValue("SkipAmore") ) {
    skipAmore=((TObjString*)config.GetConfigurationMap()->GetValue("SkipAmore"))->GetString().Atoi();
    printf("TPCRAWda: Skip Amore set in config\n");
  }
  //monitoring Attributes
  if ( config.GetConfigurationMap()->GetValue("MonitorAttributes") ) {
    monitorAttributes=config.GetConfigurationMap()->GetValue("MonitorAttributes")->GetName();
    printf("TPCRAWda: Monitor attributes set in config: %s\n",monitorAttributes.Data());
  }
  
  
  //reject laser triggers in a global partition if we have interleaved laser events
  unsigned char classId=0;
  int retClassId=daqDA_getClassIdFromName(laserTriggerName.Data(),&classId);
  //chek if we shall force the laser trigger id. Mainly for test purposes
  if (forceLaserTriggerId!="-1"){
    retClassId=0;
    classId=static_cast<unsigned char>(forceLaserTriggerId.Atoi());
  }
  //create trigger mask
  if (retClassId==0){
    //interleaved laser in physics runs
    //reject laser triggered events
    TString triggerClasses;
    //TODO
    //TODO: in the next release of daq put 49 back to 50!!!
    //TODO
    for (unsigned char iclassId=0; iclassId<50; ++iclassId){
      if (iclassId==classId) continue; //exclude laser trigger
      triggerClasses+=Form("%u|",(unsigned int)iclassId);
    }
    triggerClasses.Chop();
    char *table[5] = {"PHY",
                      "Y",
                      const_cast<char*>(monitorAttributes.Data()),
                      const_cast<char*>(triggerClasses.Data()),NULL};    
    monitorDeclareTableExtended(table);
    printf("TPCRAWda: Using laser trigger class Id: %u\n",(unsigned int)classId);
    printf("TPCRAWda: Accepted trigger class Ids: %s\n",triggerClasses.Data());
    printf("TPCRAWda: Monitor attributes used: %s\n",monitorAttributes.Data());
  }

  //
  // create calibration object
  //
  AliTPCCalibRaw calibRaw(config.GetConfigurationMap());   // raw calibration algorithm
  calibRaw.SetAltroMapping(mapping->GetAltroMapping()); // Use altro mapping we got from daqDetDb
  
  //timer
  TStopwatch stopWatch;
  
  //===========================//
  // loop over RAW data files //
  //==========================//
  int nevents=0;
  for(int i=1;i<argc;i++) {

    /* define data source : this is argument i */
    printf("TPCRAWda: Processing file %s\n", argv[i]);
    status=monitorSetDataSource( argv[i] );
    if (status!=0) {
      printf("TPCRAWda: monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
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
        printf ("TPCRAWda: End of File %d detected\n",i);
        break; /* end of monitoring file has been reached */
      }

      if (status!=0) {
        printf("TPCRAWda: monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        break;
      }

      /* retry if got no event */
      if (event==NULL) {
        continue;
      }
      nevents++;
      // get the run number
      runNb = event->eventRunNb;
      //  Raw calibration
      calibRaw.ProcessEvent(event);
      // sending to AMOREdb
      if (stopWatch.RealTime()>updateInterval){
        if (!skipAmore) SendToAmoreDB(&calibRaw,runNb);
        stopWatch.Start();
      } else {
        stopWatch.Continue();
      }
        
      /* free resources */
      free(event);
    }
  }

  //
  // Analyse pulser data and write them to rootfile
  //
  calibRaw.Analyse();
  printf ("TPCRAWda: %d events processed\n",nevents);

  TFile * fileTPC = new TFile (RESULT_FILE,"recreate");
  calibRaw.Write("tpcCalibRaw");
  delete fileTPC;
  printf("TPCRAWda: Wrote %s\n",RESULT_FILE);

  /* store the result file on FES */

  status=daqDA_FES_storeFile(RESULT_FILE,FILE_ID);
  if (status) {
    status = -2;
  }
  //
  //Send objects to the AMORE DB
  //
  if (!skipAmore) {
    printf ("TPCRAWda: AMORE part\n");
    SendToAmoreDB(&calibRaw, runNb);
  }
  
  return status;
}

void SendToAmoreDB(TObject *o, unsigned long32 runNb)
{
  //AMORE
  const char *amoreDANameorig=gSystem->Getenv("AMORE_DA_NAME");
  //cheet a little -- temporary solution (hopefully)
  //
  //currently amoreDA uses the environment variable AMORE_DA_NAME to create the mysql
  //table in which the calib objects are stored. This table is dropped each time AmoreDA
  //is initialised. This of course makes a problem if we would like to store different
  //calibration entries in the AMORE DB. Therefore in each DA which writes to the AMORE DB
  //the AMORE_DA_NAME env variable is overwritten.
  
  gSystem->Setenv("AMORE_DA_NAME","TPC-RAW");
  //
  // end cheet
  TDatime time;
  TObjString info(Form("Run: %u; Date: %s",runNb,time.AsSQLString()));
  
  amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
  Int_t statusDA=0;
  statusDA+=amoreDA.Send("CalibRaw",o);
  statusDA+=amoreDA.Send("Info",&info);
  if ( statusDA!=0 )
    printf("TPCRAWda: Waring: Failed to write one of the calib objects to the AMORE database\n");
  // reset env var
  if (amoreDANameorig) gSystem->Setenv("AMORE_DA_NAME",amoreDANameorig);
}
