/*********************************************************************************
- Contact:    Michal Broz     Michal.Broz@cern.ch
- Link:       http
- Raw data test file :          
- Reference run number : 	      
- Run Type:   PEDESTAL/PHYSICS
- DA Type:    LDC/MON
- Number of events needed: >=1000
- Input Files:  argument list
- Output Files: local files  AD0_Histos.root, AD0_Pedestals_On.dat (Online mapping)
                FXS file     AD0_Pedestals_Off.dat (Offline mapping)
- Trigger types used: PHYSICS_EVENT
**********************************************************************************/

/**********************************************************************************
*                                                                                 *
* AD Detector Algorithm used for extracting calibration parameters             *
*                                                                                 *
* This program connects to the DAQ data source passed as argument.                *
* It computes calibration parameters, populates local "./AD0_Pedestals_On.dat"   *            
* file, exports it to the FES, and stores it into DAQ DB                          *
* The program exits when being asked to shut down (daqDA_checkshutdown)           *
* or on End of Run event.                                                         *
* We have 32 channels instead of 16 as expected for AD0 due to the two sets of    *
* charge integrators which are used by the FEE ...                                *
*                                                                                 *
***********************************************************************************/

// DATE
#include "event.h"
#include "monitor.h"
#include "daqDA.h"

//AliRoot
#include <AliADConst.h>
#include <AliADRawStream.h>
#include <AliRawReaderDate.h>
#include <AliRawReader.h>
#include <AliDAQ.h>

// standard
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>

//ROOT
#include "TROOT.h"
#include "TPluginManager.h"
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

/* Main routine --- Arguments: monitoring data source */
      
int main(int argc, char **argv) {

/* magic line from Cvetan */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo","*","TStreamerInfo","RIO","TStreamerInfo()");
  
  int status;
  if (argc!=2) {
     printf("Wrong number of arguments\n");
     return -1;
  }
  
  // Online values (using FEE channel numbering), 
  // stored into local AD0_Pedestals_On.dat:
  Double_t adcMean[32];
  Double_t adcSigma[32];
  Double_t pedMean[32];
  Double_t pedSigma[32];
  // Offline values(same but ordered as in aliroot for offliners)
  // stored into AD0_Pedestals_Off.dat: 
  Double_t adcMeanOff[32];
  Double_t adcSigmaOff[32];
  Double_t pedMeanOff[32];
  Double_t pedSigmaOff[32];
       
//___________________________________________________
// Get cuts from AD0_Pedestal_DA.config file

  Int_t    kClockMin;   // = 16;   LHC Clock Min for pedestal calculation
  Int_t    kClockMax;   // = 19;   LHC Clock Max for pedestal calculation
  Int_t    kLowCut;     // = 60;   low cut on signal distribution - to be tuned
  Int_t    kHighCut;    // = 50;   high cut on pedestal distribution - to be tuned
  Int_t    kClockMinRef;   // = 16;   LHC Clock Min for Flag checking
  Int_t    kClockMaxRef;   // = 20;   LHC Clock Max for Flag checking
  Float_t  kChi2Max; 		// = 1.  Maximum chi2 
  std::string runType;
  std::string configFile;
  Bool_t runFound = kFALSE;


  if(getenv("DATE_RUN_TYPE") == NULL){
  	printf("DATE_RUN_TYPE enviroment variable undefined");
	return -1;
	}
  

  status = daqDA_DB_getFile("AD0_Runtypes.config","AD0_Runtypes.config");
  
  if (status) printf("Failed to get Config file (AD0_Runtypes.config) from DAQ DB, status=%d\n", status);
  else{
  	ifstream infile("AD0_Runtypes.config");
  	while(infile>>runType>>configFile){ 
  		if(runType == getenv("DATE_RUN_TYPE")) {
			runFound = kTRUE; 
			break;
			}
		}
	}
  if(!runFound){
  	printf("Run type %s not found in config file",getenv("DATE_RUN_TYPE"));
	return -1;
	}
  
  status = daqDA_DB_getFile(configFile.c_str(),configFile.c_str());
  if (status) {
      printf("Failed to get Config file (%s) from DAQ DB, status=%d\n", configFile.c_str() , status);
      printf("Take default values of parameters for pedestal calculation \n");
      kClockMin  =  0; 
      kClockMax  =  6; 
      kLowCut    =  60;   
      kHighCut   =  50;  
      kClockMinRef  =  0; 
      kClockMaxRef  =  6; 
      kChi2Max		=  100.;
  } else {
      /* open the config file and retrieve cuts */
      FILE *fpConfig = fopen(configFile.c_str(),"r");
      int res = fscanf(fpConfig,"%d %d %d %d %d %d %f",&kClockMin,&kClockMax,&kLowCut,&kHighCut,&kClockMinRef,&kClockMaxRef,&kChi2Max);
      if(res!=7) {
	    printf("Failed to get values from Config file (%s): wrong file format - 7 integers are expected - \n",configFile.c_str());
	    kClockMin  =  0; 
        kClockMax  =  6; 
    	kLowCut    =  60;   
        kHighCut   =  50; 
      	kClockMinRef  =  0; 
      	kClockMaxRef  =  6; 
      	kChi2Max	  =  100.;
      }
      fclose(fpConfig);
  }
  
  printf("LHC Clock Min for pedestal calculation = %d; LHC Clock Max for pedestal calculation = %d; LowCut on signal = %d ; HighCut on pedestal = %d\n",
          kClockMin, kClockMax, kLowCut, kHighCut);

//___________________________________________________
// Book HISTOGRAMS - dynamics of p-p collisions -
      
  char     adcName[10]; 
  char     pedName[10]; 
  TH1F     *hADCname[32];
  TH1F     *hPEDname[32]; 
  
  char  texte[20];
  for (Int_t i=0; i<16; i++) {
  	for(Int_t j=0; j<2; j++){
       		sprintf(adcName,"hADC%d%d",i,j);
       		sprintf(texte,"ADC cell%d int%d",i,j);
       		hADCname[i + 16*j]  = new TH1F(adcName,texte,1024,-0.5, 1023.5);
       		sprintf(pedName,"hPED%d%d",i,j);
       		sprintf(texte,"PED cell%d int%d",i,j);
       		hPEDname[i + 16*j]  = new TH1F(pedName,texte,1024,-0.5, 1023.5);
       		}
       }
       
  char Chi2Name[20];
  TH2F *hChi2PerEvent[2];
  for (Int_t i=0; i<2; i++) {
  	sprintf(Chi2Name,"hChi2PerEventInt%d",i);
   	hChi2PerEvent[i] = new TH2F(Chi2Name,Chi2Name,16,0,16,100,0,10);
	}
   
//___________________________________________________ 

  /* open result file to be exported to FES */
  FILE *fpLocal=NULL;
  fpLocal=fopen("./AD0_Pedestals_On.dat","w");
  if (fpLocal==NULL) {
      printf("Failed to open local result file\n");
      return -1;}
   
  /* open result file to be exported to FES */
  FILE *fp=NULL;
  fp=fopen("./AD0_Pedestals_Off.dat","w");
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
  
	   AliADRawStream* rawStream  = new AliADRawStream(rawReader); 
	   if (rawStream->Next()) {	
           for(Int_t i=0; i<16; i++) {
	   		Int_t nFlag = 0;
			for(Int_t j=kClockMinRef; j <= kClockMaxRef; j++) {  // Check flags on clock range used for pedestal calculation
		   		if((rawStream->GetBBFlag(i,j)) || (rawStream->GetBGFlag(i,j))) nFlag++; 
			}
			if(nFlag == 0){       // Fill 16*2 pedestal histograms  - 2 integrators -
				Float_t sum[2] = {0.,0.};
				Float_t sumwi[2] = {0.,0.};
		   		for(Int_t j=kClockMin;j <= kClockMax;j++){
		       		Int_t integrator = rawStream->GetIntegratorFlag(i,j);
		       		Float_t pedestal = (float)(rawStream->GetPedestal(i,j));
		       		sum[integrator] += pedestal;
		   			sumwi[integrator] += 1.;
		   		}	
		   		Float_t mean[2] =  {0.,0.};
		   		Float_t chi2[2] =  {0.,0.};
	   			
	   			for(int ii=0;ii<2;ii++) if(sumwi[ii]>1.e-6) mean[ii] = sum[ii] / sumwi[ii];

				for(Int_t j=kClockMin;j <= kClockMax;j++){
		       		Int_t integrator = rawStream->GetIntegratorFlag(i,j);
		       		Float_t pedestal = (float)(rawStream->GetPedestal(i,j));
		   			chi2[integrator] += (mean[integrator] - pedestal) * (mean[integrator] - pedestal);
	   			}
				chi2[0] = chi2[0]/(kClockMax-kClockMin+1);
				chi2[1] = chi2[1]/(kClockMax-kClockMin+1);
				
				for(int ii=0;ii<2;ii++)hChi2PerEvent[ii]->Fill(i,chi2[ii]);
					
	   			if(chi2[0]<kChi2Max && chi2[1]<kChi2Max) {
					for(Int_t j=kClockMin;j <= kClockMax;j++){
		       				Int_t integrator = rawStream->GetIntegratorFlag(i,j);
		       				Float_t pedestal = (float)(rawStream->GetPedestal(i,j));
						hPEDname[i + 16*integrator]->Fill(pedestal);
	   				}
				}

			} 
			if((rawStream->GetBBFlag(i,10)) || (rawStream->GetBGFlag(i,10))){ // Charge
		    	Int_t integrator = rawStream->GetIntegratorFlag(i,10);
		    	Float_t charge = (float)(rawStream->GetADC(i));   // Fill 16*2 ADCmax histograms 
		    	hADCname[i + 16 * integrator]->Fill(charge);
			}   			   
    	} // End loop over channels
	   }   // End : if rawstream
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
	
  for(Int_t i=0; i<32; i++) {
      hPEDname[i]->GetXaxis()->SetRangeUser(0,kHighCut);
      pedMean[i]  = hPEDname[i]->GetMean(); 
      pedSigma[i] = hPEDname[i]->GetRMS(); 
      hADCname[i]->GetXaxis()->SetRangeUser(kLowCut,1024);
      adcMean[i]  = hADCname[i]->GetMean();
      adcSigma[i] = hADCname[i]->GetRMS(); 
//      printf(" i = %d, %.3f %.3f %.3f %.3f\n",i,pedMean[i],pedSigma[i],adcMean[i],adcSigma[i]);
      fprintf(fpLocal," %.3f %.3f %.3f %.3f\n",pedMean[i],pedSigma[i],adcMean[i],adcSigma[i]);
      
      if (i < 16) {
          Int_t j = kOfflineChannel[i];     
          pedMeanOff[j]  = pedMean[i];
          pedSigmaOff[j] = pedSigma[i];
          adcMeanOff[j]  = adcMean[i];
          adcSigmaOff[j] = adcSigma[i]; 
	  }
      else{
          Int_t j = kOfflineChannel[i-16];     
          pedMeanOff[j+16]  = pedMean[i];
          pedSigmaOff[j+16] = pedSigma[i];
          adcMeanOff[j+16]  = adcMean[i];
          adcSigmaOff[j+16] = adcSigma[i];
      }
  }
  
  for(Int_t j=0; j<32; j++) {
//      printf(" j = %d, %.3f %.3f %.3f %.3f\n",j,pedMeanOff[j],pedSigmaOff[j],
//                                                adcMeanOff[j],adcSigmaOff[j]);
      fprintf(fp," %.3f %.3f %.3f %.3f\n",pedMeanOff[j],pedSigmaOff[j],
                                          adcMeanOff[j],adcSigmaOff[j]);				       
  }
   
//________________________________________________________________________
// Write root file with histos for users further check - just in case - 

  TFile *histoFile = new TFile("AD0_Pedestals.root","RECREATE");
    
  for (Int_t i=0; i<32; i++) {
       hADCname[i]->GetXaxis()->SetRange(0,1024);
       hADCname[i]->Write(); 
       hPEDname[i]->Write(); 
       }
  hChi2PerEvent[0]->Write();
  hChi2PerEvent[1]->Write();
  histoFile->Close(); 
  delete histoFile;
  
//________________________________________________________________________
   
  /* close local result file and FXS result file*/
  fclose(fpLocal);
  fclose(fp);
 
  /* export result file to FES */
  status=daqDA_FES_storeFile("AD0_Pedestals_Off.dat","AD0da_results");
  if (status)    {
      printf("Failed to export file : %d\n",status);
      return -1; }

  /* store result file into Online DB in case of LDC mode*/
  if(runType == "PEDESTAL"){
  	status=daqDA_DB_storeFile("AD0_Pedestals_On.dat","AD0_Pedestals_On.dat");
  	if (status)    {
      		printf("Failed to store file into Online DB: %d\n",status);
      		return -1; }
	}
      
  return status;
}
    
