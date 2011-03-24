/*********************************************************************************
- Contact:    Brigitte Cheynis     b.cheynis@ipnl.in2p3.fr
- Link:       http
- Raw data test file :          
- Reference run number : 	      
- Run Type:   PHYSICS
- DA Type:    MON
- Number of events needed: >=500
- Input Files:  argument list
- Output Files: local files  VZERO_Histos.root, V0_Pedestals.dat (Online mapping)
                FXS file     V0_Ped_Width_Gain.dat (Offline mapping)
- Trigger types used: PHYSICS_EVENT
**********************************************************************************/

/**********************************************************************************
*                                                                                 *
* VZERO Detector Algorithm used for extracting calibration parameters             *
*                                                                                 *
* This program connects to the DAQ data source passed as argument.                *
* It computes calibration parameters, populates local "./V0_Ped_Width_Gain.dat"   *            
* file, exports it to the FES, and stores it into DAQ DB                          *
* The program exits when being asked to shut down (daqDA_checkshutdown)           *
* or on End of Run event.                                                         *
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

Int_t GetOfflineChannel(Int_t channel);

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
  
  // Online values (using FEE channel numbering), 
  // stored into local V0_Pedestals.dat:
  Double_t adcMean[128];
  Double_t adcSigma[128];
  Double_t pedMean[128];
  Double_t pedSigma[128];
  // Offline values(same but ordered as in aliroot for offliners)
  // stored into V0_Ped_Width_Gain.dat: 
  Double_t adcMeanOff[128];
  Double_t adcSigmaOff[128];
  Double_t pedMeanOff[128];
  Double_t pedSigmaOff[128];
       
//___________________________________________________
// Get cuts from V00DA.config file

  Int_t    kClockMin;   // = 16;   LHC Clock Min for pedestal calculation
  Int_t    kClockMax;   // = 19;   LHC Clock Max for pedestal calculation
  Int_t    kLowCut;     // = 60;   low cut on signal distribution - to be tuned
  Int_t    kHighCut;    // = 50;   high cut on pedestal distribution - to be tuned
  
  status = daqDA_DB_getFile("V00DA.config","./V00DA.config");
  if (status) {
      printf("Failed to get Config file (V00DA.config) from DAQ DB, status=%d\n", status);
      printf("Take default values of parameters for pedestal calculation \n");
      kClockMin  =  16; 
      kClockMax  =  19; 
      kLowCut    =  60;   
      kHighCut   =  50;  
  } else {
      /* open the config file and retrieve cuts */
      FILE *fpConfig = fopen("V00DA.config","r");
      int res = fscanf(fpConfig,"%d %d %d %d ",&kClockMin,&kClockMax,&kLowCut,&kHighCut);
      if(res!=4) {
	    printf("Failed to get values from Config file (V00DA.config): wrong file format - 4 integers are expected - \n");
	    kClockMin  =  16; 
            kClockMax  =  19; 
            kLowCut    =  60;   
            kHighCut   =  50; 
      }
      fclose(fpConfig);
  }
  
  printf("LHC Clock Min for pedestal calculation = %d; LHC Clock Max for pedestal calculation = %d; LowCut on signal = %d ; HighCut on pedestal = %d\n",
          kClockMin, kClockMax, kLowCut, kHighCut);

//___________________________________________________
// Book HISTOGRAMS - dynamics of p-p collisions -
      
  char     adcName[6]; 
  char     pedName[6]; 
  TH1F     *hADCname[128];
  TH1F     *hPEDname[128];  
  
  char  texte[12];
  for (Int_t i=0; i<128; i++) {
       sprintf(adcName,"hADC%d",i);
       sprintf(texte,"ADC cell%d",i);
       hADCname[i]  = new TH1F(adcName,texte,1024,-0.5, 1023.5);
       sprintf(pedName,"hPED%d",i);
       sprintf(texte,"PED cell%d",i);
       hPEDname[i]  = new TH1F(pedName,texte,1024,-0.5, 1023.5);
  }
//___________________________________________________ 

  /* open result file to be exported to FES */
  FILE *fpLocal=NULL;
  fpLocal=fopen("./V0_Pedestals.dat","w");
  if (fpLocal==NULL) {
      printf("Failed to open local result file\n");
      return -1;}
   
  /* open result file to be exported to FES */
  FILE *fp=NULL;
  fp=fopen("./V0_Ped_Width_Gain.dat","w");
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
  int neventsPhysics=0;
  int neventsTotal=0;

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
           neventsPhysics++;
 	     		 
	   AliRawReader *rawReader = new AliRawReaderDate((void*)event);
  
	   AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader); 
	   if (rawStream->Next()) {	
           for(Int_t i=0; i<64; i++) {
	   	Int_t nFlag = 0;
		for(Int_t j=kClockMin; j <= kClockMax; j++) {  // Check flags on clock range used for pedestal calculation
		   if((rawStream->GetBBFlag(i,j)) || (rawStream->GetBGFlag(i,j))) nFlag++; 
		}
		if(nFlag == 0){       // Fill 64*2 pedestal histograms  - 2 integrators -
		   for(Int_t j=kClockMin;j <= kClockMax;j++){
		       Int_t integrator = rawStream->GetIntegratorFlag(i,j);
		       Float_t pedestal = (float)(rawStream->GetPedestal(i,j));
		       hPEDname[i + 64 * integrator]->Fill(pedestal);
		   }	
		} 
		if((rawStream->GetBBFlag(i,10)) || (rawStream->GetBGFlag(i,10))){ // Charge
		    Int_t integrator = rawStream->GetIntegratorFlag(i,10);
		    Float_t charge = (float)(rawStream->GetADC(i));   // Fill 64*2 ADCmax histograms 
		    hADCname[i + 64 * integrator]->Fill(charge);
		}   			   
             } 
	   }   
           delete rawStream;
           rawStream = 0x0;      
           delete rawReader;
           rawReader = 0x0;	     	 						         
      } // end of switch on event type 
	
      neventsTotal++;
      /* free resources */
      free(event);
	
      /* exit when last event received, no need to wait for TERM signal */
      if (eventT==END_OF_RUN) {
           printf("End Of Run event detected\n");
      break;
    }

  }  // loop over events
  
  printf("%d physics events processed\n",neventsPhysics);
    
//___________________________________________________________________________
//  Computes mean values, converts FEE channels into Offline AliRoot channels
//  and dumps the ordered values into the output text file for SHUTTLE
	
  for(Int_t i=0; i<128; i++) {
      hPEDname[i]->GetXaxis()->SetRange(0,kHighCut);
      pedMean[i]  = hPEDname[i]->GetMean(); 
      pedSigma[i] = hPEDname[i]->GetRMS(); 
      hADCname[i]->GetXaxis()->SetRange(kLowCut,1024);
      adcMean[i]  = hADCname[i]->GetMean();
      adcSigma[i] = hADCname[i]->GetRMS(); 
//      printf(" i = %d, %.3f %.3f %.3f %.3f\n",i,pedMean[i],pedSigma[i],adcMean[i],adcSigma[i]);
      fprintf(fpLocal," %.3f %.3f %.3f %.3f\n",pedMean[i],pedSigma[i],
                                               adcMean[i],adcSigma[i]);
      if (i < 64) {
          Int_t j = GetOfflineChannel(i);     
          pedMeanOff[j]  = pedMean[i];
          pedSigmaOff[j] = pedSigma[i];
          adcMeanOff[j]  = adcMean[i];
          adcSigmaOff[j] = adcSigma[i]; }
      else{
          Int_t j = GetOfflineChannel(i-64);     
          pedMeanOff[j+64]  = pedMean[i];
          pedSigmaOff[j+64] = pedSigma[i];
          adcMeanOff[j+64]  = adcMean[i];
          adcSigmaOff[j+64] = adcSigma[i];
      }
  }
  
  for(Int_t j=0; j<128; j++) {
//      printf(" j = %d, %.3f %.3f %.3f %.3f\n",j,pedMeanOff[j],pedSigmaOff[j],
//                                                adcMeanOff[j],adcSigmaOff[j]);
      fprintf(fp," %.3f %.3f %.3f %.3f\n",pedMeanOff[j],pedSigmaOff[j],
                                          adcMeanOff[j],adcSigmaOff[j]);				       
  }
   
//________________________________________________________________________
// Write root file with histos for users further check - just in case - 

  TFile *histoFile = new TFile("VZERO_histos.root","RECREATE");
    
  for (Int_t i=0; i<128; i++) {
       hADCname[i]->GetXaxis()->SetRange(0,1024);
       hADCname[i]->Write(); 
       hPEDname[i]->Write(); }

  histoFile->Close(); 
  delete histoFile;
  
//________________________________________________________________________
   
  /* close local result file and FXS result file*/
  fclose(fpLocal);
  fclose(fp);
 
  /* export result file to FES */
  status=daqDA_FES_storeFile("./V0_Ped_Width_Gain.dat","V00da_results");
  if (status)    {
      printf("Failed to export file : %d\n",status);
      return -1; }

  /* store result file into Online DB */
  status=daqDA_DB_storeFile("./V0_Pedestals.dat","V00da_results");
  if (status)    {
      printf("Failed to store file into Online DB: %d\n",status);
      return -1; }
      
  return status;
}

 Int_t GetOfflineChannel(Int_t channel) {

// Channel mapping Online - Offline:
 
 Int_t fOfflineChannel[64] = {39, 38, 37, 36, 35, 34, 33, 32, 
                              47, 46, 45, 44, 43, 42, 41, 40, 
			      55, 54, 53, 52, 51, 50, 49, 48, 
			      63, 62, 61, 60, 59, 58, 57, 56,
			       7,  6,  5,  4,  3,  2,  1,  0, 
			      15, 14, 13, 12, 11, 10,  9,  8,
			      23, 22, 21, 20, 19, 18, 17, 16, 
			      31, 30, 29, 28, 27, 26, 25, 24};
 return	fOfflineChannel[channel];		      
}			      
