/*

HMPID PHYSICS DA for noise and DDL monitoring

Contact:                  Levente.Molnar@cern.ch, Giacomo.Volpe@ba.infn.it
Link:                     https://twiki.cern.ch/twiki/bin/view/ALICE/DAInstructions 
Run Type:                 PHSYICS and STANDALONE
DA Type:                  MON
Reference Run:            104157
Number of events needed:  ~ 2000 events
Input Files:              none
Output Files:             HmpPhysicsDaNoiseMap.root
Trigger types used:       PHSYICS and STANDALONE

*/

#define NUM_PHYSICS_EVENTS 2000
#define HMP_OUTPUT_FILE "HmpPhysicsDaNoiseMap.root"
#include "event.h"
#include "monitor.h"
//DAQ
extern "C" {
#include <daqDA.h>
}
//
#include <stdio.h>
#include <stdlib.h>
//ROOT
#include <TPluginManager.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH2.h>
//AliRoot
#include "AliRawReaderDate.h"
#include "AliHMPIDRawStream.h"
#include "AliHMPIDParam.h"


int main(int argc, char **argv) {

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()"); 
  int status;  
  if (argc<2) { printf("HMP PhysicsDa: Wrong number of arguments\n"); return -1; }

  /* define data source : this is argument 1 */  
  status=monitorSetDataSource( argv[1] );
  if (status!=0) {
    printf("HMP PhysicsDa: monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("HMP PhysicsDa: monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);
  

  /* log start of process */
  printf("HMP PhysicsDa: HMPID PHYSICS Monitoring DA program has started\n");  


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;
  
  struct eventHeaderStruct *event;
  eventTypeType eventT;
  
  Int_t runNumber=0;
  AliHMPIDDigit dig;
  TH2F *hHmpNoiseMaps=new TH2F("hHmpNoiseMaps","Noise Maps Ch: 0-7;Ch 0-7: pad X;Ch0, Ch1, Ch2, Ch3, Ch4, Ch5, Ch6 pad Y ;Noise level (\%)",160,0,160,144*7,0,144*7); //In y we have all 7 chambers
 
  /* main loop (infinite) */
  for(;;) {
    /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}
    
    /* get next event (blocking call until timeout) */
    status=monitorGetEventDynamic((void **)&event);
    if (status==MON_ERR_EOF) {
      printf ("HMP PhysicsDa: End of File detected\n");
      break; /* end of monitoring file has been reached */
    }
    
    if (status!=0) {
      printf("HMP PhysicsDa: monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
      break;
    }
    
    /* retry if got no event */
    if (event==NULL) {
      continue;
    }

    /* use event - here, just write event id to result file */
    eventT=event->eventType;
    
    if (eventT==PHYSICS_EVENT) {
      runNumber=(unsigned long)event->eventRunNb;   
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      AliHMPIDRawStream stream(reader);stream.SetTurbo(kTRUE);                  //raw data decoding without error checks SetTurbo(kTRUE)
      while(stream.Next())
         {
            for(Int_t iPad=0;iPad<stream.GetNPads();iPad++) {
              dig.SetPad(stream.GetPadArray()[iPad]);
              hHmpNoiseMaps->Fill( dig.PadChX(), dig.Ch()*144+dig.PadChY() );
            } //pads            
          }//while -- loop on det load in one event
          
      nevents_physics++;
    }//physics event
    nevents_total++;
  
    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("HMP PhysicsDa: EOR event detected\n");
      break;
    }
    /* Exit after the NUM_PHYSICS_EVENTS physics events */
    if( nevents_physics == NUM_PHYSICS_EVENTS ) 
    {
     printf("HMP PhysicsDa: The number of required physics events (%d) is reached!\n",NUM_PHYSICS_EVENTS);
     break;  
    }    
  }//main loop


  /* write report */
  printf("HMP PhysicsDa: Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* set histogram properties */

  hHmpNoiseMaps->SetTitle(Form("Run number: %d Tested events: %d",runNumber,nevents_physics)); 
  if(nevents_physics!=0) hHmpNoiseMaps->Scale(1.0/nevents_physics);


  /* write output file */
  TFile *fout = new TFile(HMP_OUTPUT_FILE,"recreate");
  hHmpNoiseMaps->Write();
  fout->Close();
  
  /* send file to DAQ FXS */
  status=0;
  status=daqDA_FES_storeFile(HMP_OUTPUT_FILE,HMP_OUTPUT_FILE);
  if(status) return -1;
    
  return status;
}
