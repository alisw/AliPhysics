/*


TRDda.cxx - calibration algorithm to be run on monitoring server
DAcase2.c

AliTRDCalibraFillHisto - average pulse height/ vdrift calibration
AliTRDCalibPadStatus - pad status calibration


This program connects to the DAQ data source passed as argument
and populates local "./result.txt" file with the ids of events received
during the run.

The program exits when being asked to shut down (daqDA_checkshutdown)
or End of Run event.

Messages on stdout are exported to DAQ log system.

contact: alice-datesupport@cern.ch

*/

#define RESULT_FILE "trdCalibrationv.root"


extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"
#include <stdio.h>
#include <stdlib.h>


//
//Root includes
//
#include <TObjArray.h>
#include <TString.h>
#include <TVectorF.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TFile.h>
//
//AliRoot includes
//
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliTRDRawStreamV2.h"
#include "AliCDBManager.h"

//
// TRD calibration algorithm includes
//
#include "AliTRDCalibraFillHisto.h"




/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {

  int status;

  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  /* define data source : this is argument 1 */  
  status=monitorSetDataSource( argv[1] );
  if (status!=0) {
    printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);
  

  /* log start of process */
  printf("TRD DA VDRIFT started\n");  


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  //Instance of AliCDBManager: needed by AliTRDRawStreamV2
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);
  //Instance of AliTRDCalibraFillHisto
  AliTRDCalibraFillHisto *calibra      = AliTRDCalibraFillHisto::Instance();
  // everythings are okey for vdrift
  Bool_t passvdrift  = kTRUE;    // if timebin okey
  Int_t  nbvdrift    = 0;     // number of events with entries for vdrift

  
  /* main loop (infinite) */
  for(;;) {
  //for(Int_t k = 0; k < 20; k++) {
    struct eventHeaderStruct *event;
    eventTypeType eventT;

    /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}
    
    /* get next event (blocking call until timeout) */
    status=monitorGetEventDynamic((void **)&event);
    if (status==MON_ERR_EOF) {
      printf ("End of File detected\n");
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

    if( ((Int_t)nevents_total)%100 == 0 ) printf(" event number %d (physic event number %d) will be processed\n",(Int_t) nevents_total,(Int_t) nevents_physics);  


    /* use event - here, just write event id to result file */
    eventT=event->eventType;
    
    //
    // vdrift calibration: run only for physics events
    //
    if ((eventT==PHYSICS_EVENT) && (passvdrift)) {
      //if (eventT==PHYSICS_EVENT) {
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);
      AliTRDRawStreamV2 *trdRawStream = new AliTRDRawStreamV2((AliRawReader *) rawReader);
      Int_t result = calibra->ProcessEventDAQ(trdRawStream,(Bool_t)nevents_physics);
      if(!result) passvdrift = kFALSE;
      else nbvdrift += (Int_t) result/2;
             
      
      delete trdRawStream;
      delete rawReader;
    } 
   
    if(eventT==PHYSICS_EVENT){
 
      nevents_physics++;
     
    }
    
    nevents_total++;


    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }
  }


  /* report progress */
  printf("%d events processed and %d used\n",nevents_total,nbvdrift);
  
  //
  // Write a local file and put it on the FX 
  //
  TFile  *fileTRD     = new TFile(RESULT_FILE,"recreate");
 
  if((nbvdrift > 0) && passvdrift){
    //Double_t *stat = calibra->StatH((TH2 *)(calibra->GetPH2d()),1);
    //if((stat[6] < 0.20) && (stat[5] > 3000.0)) {
    // write the histogram in any case to see if problems
    calibra->GetPH2d()->Write();
    //}
  }
  fileTRD->Close();
  status=0;
  // Export the file in any case to see if problems
  if(daqDA_FES_storeFile(RESULT_FILE,RESULT_FILE)) status = -2;
  
  delete fileTRD;  

  return status;
}
