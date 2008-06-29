/*********************************************************************************
- Contact:    Brigitte Cheynis     b.cheynis@ipnl.in2p3.fr
- Link:       http
- Raw data test file :          
- Reference run number : 	      
- Run Type:   PHYSICS
- DA Type:    MON
- Number of events needed: >=500
- Input Files:  argument list
- Output Files: local files  VZERO_Histos.root, V0_Ped_Width_Gain.dat
                FXS file     V0_Ped_Width_Gain.dat
- Trigger types used: PHYSICS_EVENT
**********************************************************************************/

/**********************************************************************************
*                                                                                 *
* VZERO Detector Algorithm used for extracting calibration parameters             *
*                                                                                 *
* This program connects to the DAQ data source passed as argument.                *
* It computes calibration parameters, populates local "./V0_Ped_Width_Gain.dat"   *            
* file and exports it to the FES.                                                 *
* The program exits when being asked to shut down (daqDA_checkshutdown)           *
* or End of Run event.                                                            *
* We have 128 channels instead of 64 as expected for V0 due to the two sets of    *
* charge integrators which are used by the FEE ...                                *
*                                                                                 *
***********************************************************************************/

// DATE
#include "event.h"
#include "monitor.h"
#include "daqDA.h"

//AliRoot
#include <AliVZERORawStream.h>
#include <AliRawReaderDate.h>
#include <AliRawReader.h>
#include <AliDAQ.h>

// standard
#include <stdio.h>
#include <stdlib.h>

//ROOT
#include "TROOT.h"
#include "TPluginManager.h"
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>

/* Main routine --- Arguments: monitoring data source */
      
int main(int argc, char **argv) {

/* magic line from Cvetan */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                    "*",
                    "TStreamerInfo",
                    "RIO",
                    "TStreamerInfo()");
  int status;
  if (argc!=2) {
     printf("Wrong number of arguments\n");
     return -1;
  }

//  printf(" argc = %d, argv = %s \n",argc, &(**argv));

  Int_t    kHighCut = 50; // high cut on pedestal distribution - to be tuned
  Int_t    kLowCut  = 30; // low cut on signal distribution - to be tuned
  Double_t ADCmean[128];
  Double_t ADCsigma[128];
  Double_t PEDmean[128];
  Double_t PEDsigma[128];
  
//___________________________________________________
// Book HISTOGRAMS - dynamics of p-p collisions -
      
  char     ADCname[6]; 
  char     PEDname[6]; 
  TH1F     *hADCname[128];
  TH1F     *hPEDname[128];  
  
  char  texte[12];
  for (Int_t i=0; i<128; i++) {
       sprintf(ADCname,"hADC%d",i);
       sprintf(texte,"ADC cell%d",i);
       hADCname[i]  = new TH1F(ADCname,texte,1024,0,1023);
       sprintf(PEDname,"hPED%d",i);
       sprintf(texte,"PED cell%d",i);
       hPEDname[i]  = new TH1F(PEDname,texte,1024,0,1023);
  }
//___________________________________________________ 
  
 
  /* open result file to be exported to FES */
  FILE *fp=NULL;
  fp=fopen("./V0_Ped_Width_Gain.dat","a");
  if (fp==NULL) {
      printf("Failed to open result file\n");
      return -1;}

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
  
  /* init counters on events */
  int nevents_physics=0;
  int nevents_total=0;

  /* loop on events (infinite) */
  for(;;) {
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
      if (event==NULL) continue;
               
      /* decode event */
      eventT=event->eventType;
	
      switch (event->eventType){
      
      case START_OF_RUN:
           break;
      
      case END_OF_RUN:
           printf("End Of Run detected\n");
           break;
      
      case PHYSICS_EVENT:
           nevents_physics++;
 	     
//            fprintf(flog,"Run #%lu, event size: %lu, BC:%u, Orbit:%u, Period:%u\n",
//                  (unsigned long)event->eventRunNb,
//                  (unsigned long)event->eventSize,
//                  EVENT_ID_GET_BUNCH_CROSSING(event->eventId),
//                  EVENT_ID_GET_ORBIT(event->eventId),
//                  EVENT_ID_GET_PERIOD(event->eventId) );
		 
	   AliRawReader *rawReader = new AliRawReaderDate((void*)event);
  
	   AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader); 
	   rawStream->Next();	
           for(Int_t i=0; i<64; i++) {
	       if(!rawStream->GetIntegratorFlag(i,10))
                   hADCname[i]->Fill(float(rawStream->GetADC(i)));       // even integrator - fills 0 to 63
	       else 
		   hADCname[i+64]->Fill(float(rawStream->GetADC(i)));    // odd integrator  - fills 64 to 123
	       for(Int_t j=0; j<21; j++) {
		   if(j==10) continue;
		   if(!rawStream->GetIntegratorFlag(i,j))
                       { hPEDname[i]->Fill(float(rawStream->GetPedestal(i,j))); }     // even integrator
		   else 
		       { hPEDname[i+64]->Fill(float(rawStream->GetPedestal(i,j))); }  // odd integrator 
	       }
           }    
           delete rawStream;
           rawStream = 0x0;      
           delete rawReader;
           rawReader = 0x0;	     	 						         
      } // end of switch on event type 
	
      nevents_total++;
      /* free resources */
      free(event);
	
      /* exit when last event received, no need to wait for TERM signal */
      if (eventT==END_OF_RUN) {
           printf("End Of Run event detected\n");
      break;
    }

  }  // loop over events
    
//________________________________________________________________________
//  Computes mean values, dumps them into the output text file
	
  for(Int_t i=0; i<128; i++) {
      hPEDname[i]->GetXaxis()->SetRange(0,kHighCut);
      PEDmean[i]  = hPEDname[i]->GetMean(); 
      PEDsigma[i] = hPEDname[i]->GetRMS(); 
      hADCname[i]->GetXaxis()->SetRange(kLowCut,1023);
      ADCmean[i]  = hADCname[i]->GetMean();
      ADCsigma[i] = hADCname[i]->GetRMS(); 
      fprintf(fp," %.3f %.3f %.3f %.3f\n",PEDmean[i],PEDsigma[i],ADCmean[i],ADCsigma[i]);
  } 

//________________________________________________________________________
// Write root file with histos for users further check - just in case - 

  TFile *histoFile = new TFile("VZERO_histos.root","RECREATE");
    
  for (Int_t i=0; i<128; i++) {
       hADCname[i]->GetXaxis()->SetRange(0,1023);
       hADCname[i]->Write(); 
       hPEDname[i]->Write(); }

  histoFile->Close(); 
  delete histoFile;
  
//________________________________________________________________________
   

  /* close result file */
  fclose(fp);
 
  /* export result file to FES */
  status=daqDA_FES_storeFile("./V0_Ped_Width_Gain.dat","V00da_results");
  if (status)    {
      printf("Failed to export file : %d\n",status);
      return -1; }

  return status;
}
