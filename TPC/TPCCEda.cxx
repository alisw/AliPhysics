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
#include <pthread.h>
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
void SendToAmoreDB(AliTPCCalibCE &calibCE, unsigned long32 runNb);
//for threaded processing
void *processEventBuffer(void *arg);

//common event processing variables for threaded processing
std::vector<eventHeaderStruct*> eventBuffer;
volatile int bStop = false;
struct timespec duree_nanosleep;
Int_t  forceNevents=-1;
Bool_t forceBufferEnds=kFALSE;


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
  
  //
  // thread wait time
  //
  duree_nanosleep.tv_sec=0;
  duree_nanosleep.tv_nsec=1000000; //1 ms
  
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
  // check configuration options
  TString laserTriggerName("C0LSR-ABCE-NOPF-CENT");
  size_t  bufferSize=30;
  Int_t   forceTriggerId=-1;
  Int_t   forceNeventsStandalone=-1;
  Bool_t  forceBufferEndsGlobal=kTRUE;
//   Bool_t  forceBufferEndsGlobalDummy=kFALSE;  
  if ( config.GetConfigurationMap()->GetValue("LaserTriggerName") ) {
    laserTriggerName=config.GetConfigurationMap()->GetValue("LaserTriggerName")->GetName();
    printf("Laser trigger class name set to: %s.\n",laserTriggerName.Data());
  }
  
  if ( config.GetConfigurationMap()->GetValue("BufferSize") ) {
    bufferSize=(size_t)config.GetValue("BufferSize");
    printf("Setting event buffer size to: %d.\n",bufferSize);
  }
  
  if ( config.GetConfigurationMap()->GetValue("ForceTriggerId") ) {
    forceTriggerId=TMath::Nint(config.GetValue("ForceTriggerId"));
    printf("Only processing triggers with Id: %d.\n",forceTriggerId);
  }
  
  if ( config.GetConfigurationMap()->GetValue("ForceBufferEndsGlobal") ) {
    forceBufferEndsGlobal=config.GetValue("ForceBufferEndsGlobal")!=0.;
    printf("Process all buffered events in global partition: %s.\n",forceBufferEndsGlobal?"yes":"no");
  }
  
  if ( config.GetConfigurationMap()->GetValue("ForceNMaxEvents") ) {
    forceNevents=TMath::Nint(config.GetValue("ForceNMaxEvents"));
    printf("Forcing maximum number of %d events.\n",forceNeventsStandalone);
  }
  
  //subsribe to laser triggers only in physics partition
  //if the trigger class is not available the return value is -1
  //in this case we are most probably running as a standalone
  //  laser run and should request all events
  unsigned char *classIdptr=0;
  int retClassId=daqDA_getClassIdFromName(laserTriggerName.Data(),classIdptr);
  if (retClassId==0){
    //interleaved laser in physics runs
    char c[5];
    sprintf(c,"%u",*classIdptr);
    char *table[5] = {"PHY","Y","*",c,NULL};
    monitorDeclareTableExtended(table);
    printf("Using trigger class Id: %s\n",c);
//     forceBufferEndsGlobal=forceBufferEndsGlobalDummy;
  } else if (retClassId==-1){
    //defaul case, accept all physics events
    char *table[3] = {"PHY","Y",NULL};
    monitorDeclareTableExtended(table);
    printf("Using all trigger class Ids\n");
//     forceNevents=forceNeventsStandalone;
  } else {
    printf("Unknown return value of 'daqDA_getClassIdFromName': %d\n",retClassId);
    return -2;
  }

  //see if we should force the trigger id
  if (forceTriggerId>-1){
    char c[5];
    sprintf(c,"%d",forceTriggerId);
    char *table[5] = {"PHY","Y","*",c,NULL};
    monitorDeclareTableExtended(table);
//     forceBufferEndsGlobal=forceBufferEndsGlobalDummy;
  }
  
  
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
  
    
  //create calibration object
  AliTPCCalibCE calibCE(config.GetConfigurationMap());   // central electrode calibration
  calibCE.SetAltroMapping(mapping->GetAltroMapping()); // Use altro mapping we got from daqDetDb

  //
  // start thread
  //
//   sleep(5);
  pthread_t threadId=0;
  int threadStatus=0;
  threadStatus = pthread_create( &threadId, NULL, processEventBuffer, (void*)(&calibCE));
  eventBuffer.resize(bufferSize); 
  struct timespec duree_out;
  //===========================//
  // loop over RAW data files //
  //==========================//
  int nevents=0;
  size_t counter=0;
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

      //check for predefined number of events
      if (forceNevents>0 && calibCE.GetNeventsProcessed()>=forceNevents) {
        printf("Requested number of events reached (%d).\n",forceNevents);
        break;
      }
      
      //buffer events, only read them if a buffer position is free
      if (eventBuffer[counter]==0) {
        
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
        if (event==NULL)
          continue;
      
      /* skip start/end of run events */
        if ( (event->eventType != physicsEvent) && (event->eventType != calibrationEvent) ){
          free(event);
          continue;
        }

      
        // get the run number
        runNb = event->eventRunNb;

//         printf(" trigger (%05d-%03d) = %8.8x %8.8x - %02u\n",nevents, calibCE.GetNeventsProcessed(),
//                event->eventTriggerPattern[1], event->eventTriggerPattern[0],event->eventType);

//         printf("filling buffer %d\n",counter);
        eventBuffer[counter]=event;
        
        ++nevents;
        ++counter;
        if (counter >= eventBuffer.size()) {counter=0;}
      }else{
//         printf("buffer already used: %d\n",counter);
        nanosleep(&duree_nanosleep,&duree_out);
      }
    }
  }

  //
  // wait for thread to end
  //
  if (!forceBufferEndsGlobal) bStop = true;
  else forceBufferEnds=forceBufferEndsGlobal;
  
  pthread_join( threadId, NULL);
//   printf("Event Processing Thread ended with: %d\n",threadStatus);
  
  //
  // free unprocessed events
  //
  for (size_t i=0;i<eventBuffer.size();++i){
    if (eventBuffer[i]) {
      free(eventBuffer[i]);
      eventBuffer[i]=0;
//       printf("freeing buffer %d\n",i);
    }
  }

  //
  // Analyse CE data and write them to rootfile
  //
  calibCE.Analyse();
  printf ("%d events processed, %d used\n",nevents,calibCE.GetNeventsProcessed());
  
  TFile * fileTPC = new TFile (RESULT_FILE,"recreate");
  calibCE.Write("tpcCalibCE");
  delete fileTPC;
  printf("Wrote %s\n",RESULT_FILE);
  
  /* store the result file on FES */
  
  status=daqDA_FES_storeFile(RESULT_FILE,FILE_ID);
  if (status) {
    status = -2;
  }
  
  SendToAmoreDB(calibCE,runNb);
  
  return status;
}

void *processEventBuffer(void *arg)
{
  //
  // event procssing thread functio
  //

  //cast argument
  AliTPCCalibCE *ce=(AliTPCCalibCE*)arg;
  AliTPCCalibCE &calibCE=*ce;

  size_t counter=0;
  unsigned long32 runNb=0;
  Bool_t published=kTRUE;
  struct timespec duree_out;
  
  struct eventHeaderStruct *event;

  //wait for the first buffer to be filled
  while (!eventBuffer[0]) nanosleep(&duree_nanosleep,&duree_out);
  //loop over buffer
  while (!bStop){
//     printf("testing buffer: %d\n",counter);
    if (eventBuffer[counter]) {
      event=eventBuffer[counter];
      runNb = event->eventRunNb;
//       printf("processing buffer: %d\n",counter);
      eventBuffer[counter]=0;
      calibCE.ProcessEvent(event);
      free(event);
      published=kFALSE;
    } else {
      //in case of empty buffer publish the results it this was not done
      if (!published) {
        SendToAmoreDB(calibCE,runNb);
        published=kTRUE;
      }
      nanosleep(&duree_nanosleep,&duree_out);
    }
    ++counter;
    if (counter >= eventBuffer.size()) {
      counter=0;
      if (forceBufferEnds) break;
    }
  }
}

void SendToAmoreDB(AliTPCCalibCE &calibCE, unsigned long32 runNb)
{
  //AMORE
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
  TGraph *grA=calibCE.MakeGraphTimeCE(-1,0,2);
  TGraph *grC=calibCE.MakeGraphTimeCE(-2,0,2);
  TDatime time;
  TObjString info(Form("Run: %u; Date: %s",runNb,time.AsSQLString()));
  amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
  Int_t statusDA=0;
  statusDA+=amoreDA.Send("CET0",calibCE.GetCalPadT0());
  statusDA+=amoreDA.Send("CEQ",calibCE.GetCalPadQ());
  statusDA+=amoreDA.Send("CERMS",calibCE.GetCalPadRMS());
  statusDA+=amoreDA.Send("DriftA",grA);
  statusDA+=amoreDA.Send("DriftC",grC);
  statusDA+=amoreDA.Send("Info",&info);
  if ( statusDA!=0 )
    printf("Waring: Failed to write one of the calib objects to the AMORE database\n");
  // reset env var
  if (amoreDANameorig) gSystem->Setenv("AMORE_DA_NAME",amoreDANameorig);
  if (grA) delete grA;
  if (grC) delete grC;
}


