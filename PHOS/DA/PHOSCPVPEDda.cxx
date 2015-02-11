/*
CPV PED DA for processing pedestal runs and producing pedestal tables for loading to CPV FEE.

Contact: Sergey Evdokimov <sevdokim@cern.ch>
Link: https://twiki.cern.ch/twiki/bin/view/ALICE/CPVda
Reference run: 211758 (/afs/cern.ch/user/s/sevdokim/public/CPV_run211758_pedestal.raw)
Run Type:  PHYSICS
DA Type: PED
Number of events needed: 2000 events
Input files: raw data file
Output files: thr?_??.dat CpvPeds.root
Trigger types used: PHYSICS_EVENT
*/

#include "event.h"
#include "monitor.h"
#include "daqDA.h"

#include <Riostream.h>
#include <stdlib.h>
#include <fstream>
#include <string>

//AliRoot
#include "AliPHOSCpvRawStream.h"
#include "AliPHOSCpvPedProducer.h"
#include "AliPHOSCpvParam.h"
#include "AliRawReaderDate.h"
#include "AliBitPacking.h"
#include "TMath.h"

//ROOT
#include "TROOT.h"
#include "TPluginManager.h"
#include "TSAXParser.h"
#include "TTree.h"

#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"
#include "TH2S.h"
#include "TH2F.h"
#include "TObject.h"
#include "TBenchmark.h"
#include "TMath.h"
#include "TRandom.h"


int main( int argc, char **argv )
{
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");
  Int_t status,print;
  Int_t sigcut=3;
  Bool_t turbo = kTRUE;

  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // log start of process
  printf("Cpv DA program started\n");

  /* report progress */
  daqDA_progressReport(0);

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

    /* report progress */
  daqDA_progressReport(5);


  // init event counter
  Int_t iPhysEvnt=0;
  Int_t iTotEvnt =0;

  AliPHOSCpvPedProducer * pedProducer = new AliPHOSCpvPedProducer(sigcut);

  // Reader
  AliRawReader * reader;

  /* report progress */
  daqDA_progressReport(10);

  /* main loop (infinite) */
  for(;;) { // infinite loop
    struct eventHeaderStruct *event;
    eventTypeType eventT;

    /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}

    // get next event
    status=monitorGetEventDynamic((void **)&event);
    if (status==MON_ERR_EOF) { // end of monitoring file has been reached
      printf("End of monitoring file has been reached! \n");
      break;
    }

    if (status!=0) {
      printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
      break;
    }

    // retry if got no event
    if (event==NULL) continue;

    // use event - here, just write event id to result file
    eventT=event->eventType;
    if (eventT==PHYSICS_EVENT) { //we use PHYSICS_EVENT for pedestal not CALIBRATION_EVENT???
      iTotEvnt++;
      reader = new AliRawReaderDate((void*)event);
      pedProducer->LoadNewEvent(reader);
      if(pedProducer->FillPedestal()) iPhysEvnt++;
      delete reader;
    } // if PHYSICS_EVENT

    free(event);

    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }
  }

  Printf(" Received %d events, %d good events",iTotEvnt,iPhysEvnt);
  /* report progress */
  daqDA_progressReport(80);

  for(int iDDL = 0; iDDL<AliPHOSCpvParam::kNDDL; iDDL++){
    if(pedProducer -> CalcPedestal(iDDL)){
      pedProducer -> WritePedFiles(iDDL);
      //for (int iCC = 0; iCC<AliPHOSCpvParam::kNRows){
      //	status=daqDA_DB_storeFile(Form("thr%d_%02d.dat", iDDL, iCC));
      //	if(status) printf("Failed to store thr%d_%02d.dat in DAQ DB!\n",iDDL, iCC);
      //	status=daqDA_FES_storeFile(Form("thr%d_%02d.dat", iDDL, iCC));
      //	if(status) printf("Failed to export thr%d_%02d.dat to DAQ FES!\n",iDDL, iCC);
      //}
    }
  }

  pedProducer->WriteAllHistsToFile("CpvPeds.root");
  status = daqDA_DB_storeFile("CpvPeds.root","CpvPeds.root");
  if(status) printf("Failed to store CpvPeds.root in DAQ DB!\n");

  /* report progress */
  daqDA_progressReport(95);


  return status;
}
