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
#include "AliTRDRawStream.h"
#include "AliCDBManager.h"

//
// TRD calibration algorithm includes
//
#include "AliTRDCalibPadStatus.h"
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
  printf("DA example case2 monitoring program started\n");  


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  //Instance of AliCDBManager: needed by AliTRDRawStream
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);
  //Instance of AliTRDCalibraFillHisto
  AliTRDCalibraFillHisto *calibra      = AliTRDCalibraFillHisto::Instance();
  //AliTRDCalibPadStatus
  AliTRDCalibPadStatus      calibpad   = AliTRDCalibPadStatus();
  // pad status on: no zero suppression (special runs)
  Bool_t passpadstatus = kTRUE;
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

    printf(" event number %d (physic event number %d) will be processed\n",(Int_t) nevents_total,(Int_t) nevents_physics);  


    /* use event - here, just write event id to result file */
    eventT=event->eventType;
    //
    // pad status calibration: we try a first time to see if zero suppressed or not
    //
    if(passpadstatus){
      printf("pad status calibration\n");
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);
      AliTRDRawStream *trdRawStream = new AliTRDRawStream((AliRawReader *) rawReader);
      if(!calibpad.ProcessEvent(trdRawStream,(Bool_t)nevents_total)) passpadstatus = kFALSE;
      
      delete trdRawStream;
      delete rawReader;
    }
    //
    // vdrift calibration: run only for physics events
    //
    if ((eventT==PHYSICS_EVENT) && (!passpadstatus)  && (passvdrift)) {
      //if (eventT==PHYSICS_EVENT) {
      printf("vdrift calibration\n");
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);
      AliTRDRawStream *trdRawStream = new AliTRDRawStream((AliRawReader *) rawReader);
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

  //
  // Write a local file and put it on the FX 
  //
  TFile  *fileTRD     = new TFile("trdCalibration.root","recreate");
  Bool_t passwrite   = kFALSE;
  //
  // pad status
  //
  if(passpadstatus){
    // We do this at the shuttle maybe!
    //calibpad.AnalyseHisto();
    //AliTRDCalPadStatus *calPadStatus = calibpad.CreateCalPadStatus();
    calibpad.Write("calibpadstatus");
    passwrite  = kTRUE; 
  }
  //
  // vdrift
  //
  if((nbvdrift > 0) && passvdrift){
    Double_t *stat = calibra->StatH((TH2 *)(calibra->GetPH2d()),1);
    // write only of enough statistics
    if((stat[6] < 0.20) && (stat[5] > 3000.0)) {
      calibra->GetPH2d()->Write();
      passwrite = kTRUE;
    }
    
  }
  fileTRD->Close();
  status=0;
  if(passwrite) {
    if(daqDA_FES_storeFile("trdCalibration.root","trdCalibration.root")) status = -2;
  }
  delete fileTRD;  

  return status;
}
