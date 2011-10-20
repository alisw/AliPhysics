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
#include <vector>
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
#include "TStopwatch.h"
#include "TMap.h"
#include "TGraph.h"
#include "TMath.h"
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


//functios, implementation below
void SendToAmoreDB(AliTPCCalibCE *calibCE, unsigned long32 runNb);
//for threaded processing


/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {
  /* log start of process */
  printf("TPCCEda: DA started - %s\n",__FILE__);
  
  if (argc<2) {
    printf("TPCCEda: Wrong number of arguments\n");
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
    printf("TPCCEda: monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);

  //variables
  AliTPCmapper *mapping = 0;   // The TPC mapping
  char localfile[255];
  unsigned long32 runNb=0;      //run number
  
 
  //
  // DA configuration from configuration file
  //
  //retrieve configuration file
  sprintf(localfile,"./%s",CONFIG_FILE);
  status = daqDA_DB_getFile(CONFIG_FILE,localfile);
  if (status) {
    printf("TPCCEda: Failed to get configuration file (%s) from DAQdetDB, status=%d\n", CONFIG_FILE, status);
    return -1;
  }
  AliTPCConfigDA config(CONFIG_FILE);
  // check configuration options
  TString laserTriggerName("C0LSR-ABCE-NOPF-CENT");
  TString monitoringType("YES");
  Int_t   forceTriggerId=-1;
  Int_t   saveOption=2; // how to store the object. See AliTPCCalibCE::DumpToFile
  Bool_t  skipAmore=kFALSE;
  
  if ( config.GetConfigurationMap()->GetValue("LaserTriggerName") ) {
    laserTriggerName=config.GetConfigurationMap()->GetValue("LaserTriggerName")->GetName();
    printf("TPCCEda: Laser trigger class name set to: %s.\n",laserTriggerName.Data());
  }

  if ( config.GetConfigurationMap()->GetValue("MonitoringType") ) {
    monitoringType=config.GetConfigurationMap()->GetValue("MonitoringType")->GetName();
    printf("TPCCEda: Monitoring type set to: %s.\n",monitoringType.Data());
  }

  if ( config.GetConfigurationMap()->GetValue("ForceLaserTriggerId") ) {
    forceTriggerId=TMath::Nint(config.GetValue("ForceLaserTriggerId"));
    printf("TPCCEda: Only processing triggers with Id: %d.\n",forceTriggerId);
  }
  
  if ( config.GetConfigurationMap()->GetValue("SaveOption") ) {
    saveOption=TMath::Nint(config.GetValue("SaveOption"));
    printf("TPCCEda: Saving option set to: %d.\n",saveOption);
  }

  if ( config.GetConfigurationMap()->GetValue("SkipAmore") ) {
    skipAmore=((TObjString*)config.GetConfigurationMap()->GetValue("SkipAmore"))->GetString().Atoi();
    printf("TPCCEda: Skip Amore set in config\n");
  }
  
  //subsribe to laser triggers only in physics partition
  //if the trigger class is not available the return value is -1
  //in this case we are most probably running as a standalone
  //  laser run and should request all events
  unsigned char classIdptr=0;
  int retClassId=daqDA_getClassIdFromName(laserTriggerName.Data(),&classIdptr);
  if (retClassId==0){
    //interleaved laser in physics runs
    //select proper trigger class id
    char c[5];
    snprintf(c,sizeof(c),"%u",(unsigned int)classIdptr);
    char *table[5] = {"PHY",const_cast<char*>(monitoringType.Data()),"*",c,NULL};
    monitorDeclareTableExtended(table);
    printf("TPCCEda: Using monitoring table: (PHY, %s, *, %s)\n",monitoringType.Data(),c);
  } else if (retClassId==-1){
    //global partition without laser triggered events
    //the DA should exit properly without processing
    printf("TPCCEda: Laser trigger class '%s' was not found among trigger class names. Will stop processing.\n",laserTriggerName.Data());
    return 0;
  } else if (retClassId==-2){
    //standalone case, accept all physics events
    char *table[5] = {"PHY","Y","*","*",NULL};
    monitorDeclareTableExtended(table);
    printf("TPCCEda: Using all trigger class Ids\n");
  } else {
    printf("TPCCEda: Unknown return value of 'daqDA_getClassIdFromName': %d\n",retClassId);
    return -2;
  }

  //see if we should force the trigger id
  if (forceTriggerId>-1){
    char c[5];
    sprintf(c,"%d",forceTriggerId);
    char *table[5] = {"PHY","Y","*",c,NULL};
    monitorDeclareTableExtended(table);
  }
  
  
  // if  test setup get parameters from $DAQDA_TEST_DIR
  if (!mapping){
    /* copy locally the mapping file from daq detector config db */
    sprintf(localfile,"./%s",MAPPING_FILE);
    status = daqDA_DB_getFile(MAPPING_FILE,localfile);
    if (status) {
      printf("TPCCEda: Failed to get mapping file (%s) from DAQdetDB, status=%d\n", MAPPING_FILE, status);
      return -1;
    }
    
    /* open the mapping file and retrieve mapping object */
    TFile *fileMapping = new TFile(MAPPING_FILE, "read");
    mapping = (AliTPCmapper*) fileMapping->Get("tpcMapping");
    delete fileMapping;
  }
  
  if (mapping == 0) {
    printf("TPCCEda: Failed to get mapping object from %s.  ...\n", MAPPING_FILE);
    return -1;
  } else {
    printf("TPCCEda: Got mapping object from %s\n", MAPPING_FILE);
  }
  
    
  //create calibration object
  AliTPCCalibCE *calibCE=new AliTPCCalibCE(config.GetConfigurationMap());   // central electrode calibration
  calibCE->SetAltroMapping(mapping->GetAltroMapping()); // Use altro mapping we got from daqDetDb

  //amore update interval
  Double_t updateInterval=300; //seconds
  Double_t valConf=config.GetValue("AmoreUpdateInterval");
  if ( valConf>0 ) updateInterval=valConf;
  //timer
  TStopwatch stopWatch;
  
  //===========================//
  // loop over RAW data files //
  //==========================//
  int nevents=0;
  int neventsOld=0;
  size_t counter=0;
  for ( i=1; i<argc; i++) {
    
    /* define data source : this is argument i */
    printf("TPCCEda: Processing file %s\n", argv[i]);
    status=monitorSetDataSource( argv[i] );
    if (status!=0) {
      printf("TPCCEda: monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
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
        printf ("TPCCEda: End of File %d detected\n",i);
        break; /* end of monitoring file has been reached */
      }
      
      if (status!=0) {
        printf("TPCCEda: monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        break;
      }
      
        /* retry if got no event */
      if (event==NULL){
        //use time in between bursts to
        // send the data to AMOREdb
        if (stopWatch.RealTime()>updateInterval){
          calibCE->Analyse();
          if (!skipAmore) SendToAmoreDB(calibCE,runNb);
          stopWatch.Start();
        } else {
          stopWatch.Continue();
        }
        //debug output
        if (nevents>neventsOld){
          printf ("TPCCEda: %d events processed, %d used\n",nevents,calibCE->GetNeventsProcessed());
          neventsOld=nevents;
        }
        
        continue;
      }
      
      /* skip start/end of run events */
      if ( (event->eventType != physicsEvent) && (event->eventType != calibrationEvent) ){
        free(event);
        continue;
      }
      
      
      // get the run number
      runNb = event->eventRunNb;
      
      // CE calibration
      calibCE->ProcessEvent(event);
      
      /* free resources */
      free(event);
      ++nevents;
    }
  }
  
  //
  // Analyse CE data and write them to rootfile
  //
  printf ("TPCCEda: %d events processed, %d used\n",nevents,calibCE->GetNeventsProcessed());

  //save data to file
  calibCE->DumpToFile(RESULT_FILE,Form("name=tpcCalibCE,type=%d",saveOption));
  printf("TPCCEda: Wrote %s\n",RESULT_FILE);
  
  /* store the result file on FES */
  status=daqDA_FES_storeFile(RESULT_FILE,FILE_ID);
  if (status) {
    status = -2;
  }

  if (!skipAmore){
    printf("TPCCEda: Amore part\n");
    calibCE->Analyse();
    SendToAmoreDB(calibCE,runNb);
  }
  
  delete calibCE;
  return status;
}


void SendToAmoreDB(AliTPCCalibCE *calibCE, unsigned long32 runNb)
{
  //AMORE
//   printf ("AMORE part\n");
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
  TGraph *grA=calibCE->MakeGraphTimeCE(-1,0,2);
  TGraph *grC=calibCE->MakeGraphTimeCE(-2,0,2);
  TDatime time;
  TObjString info(Form("Run: %u; Date: %s",runNb,time.AsSQLString()));
  amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
  Int_t statusDA=0;
  statusDA+=amoreDA.Send("CET0",calibCE->GetCalPadT0());
  statusDA+=amoreDA.Send("CEQ",calibCE->GetCalPadQ());
  statusDA+=amoreDA.Send("CERMS",calibCE->GetCalPadRMS());
  statusDA+=amoreDA.Send("DriftA",grA);
  statusDA+=amoreDA.Send("DriftC",grC);
  statusDA+=amoreDA.Send("Info",&info);
  if ( statusDA!=0 )
    printf("TPCCEda: Waring: Failed to write one of the calib objects to the AMORE database\n");
  // reset env var
  if (amoreDANameorig) gSystem->Setenv("AMORE_DA_NAME",amoreDANameorig);
  if (grA) delete grA;
  if (grC) delete grC;
}


