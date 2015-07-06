/*********************************************************************************
- Contact:    Michal Broz     Michal.Broz@cern.ch
- Link:       http
- Raw data test file :          
- Reference run number : 	      
- Run Type:   PHYSICS
- DA Type:    MON
- Number of events needed: >=50000
- Input Files:  argument list
- Output Files: FXS file  AD0_SlewingSplines.root (Offline mapping)
- Trigger types used: PHYSICS_EVENT
**********************************************************************************/

/**********************************************************************************
*                                                                                 *
* AD Detector Algorithm used for extracting time slewing splines                  *
*                                                                                 *
* This program connects to the DAQ data source passed as argument.                *
* It reads calibration parameters, from DAQ_DB "AD0_Pedestals_On.dat"   	  *            
* file, Compute 16 TSplines, place it to a root file in TList. 			  *
* File is stored in FXS.                          				  *
* The program exits when being asked to shut down (daqDA_checkshutdown)           *
* or on End of Run event.                                                         *
* Runs in monitoring mode, the splines needs solid statistics, so it asks 	  *
* for rather high number of events to be processed.				  *
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
#include <TString.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TSpline.h>
#include <TF1.h>
#include <TList.h>
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
         
//___________________________________________________
// Get cuts from AD0_Pedestal_DA.config file
  Int_t    kMinEvents; //Minimal number of events requested to make splines
  Int_t    kNPreClocks; //Pre-clocks for charge integration  
  Int_t    kNPostClocks; //Post-clocks for charge integration
  Int_t    kNPedSigma; //Cut on pedestal number of sigmas
  Int_t	   kNBinsCharge; //Number of bins in charge distribution
  Int_t	   kMinTime;
  Int_t    kMaxTime;
  
  status = daqDA_DB_getFile("AD0_TimeSlewing_DA.config","AD0_TimeSlewing_DA.config");
  if (status) {
      printf("Failed to get Config file AD0_TimeSlewing_DA.config from DAQ DB, status=%d\n", status);
      printf("Take default values of parameters for pedestal calculation \n");
       kMinEvents = 50000;
       kNPreClocks = 1;  
       kNPostClocks = 10;
       kNPedSigma = 3;
       kNBinsCharge = 1000;
       kMinTime = 1740;
       kMaxTime = 2252;
  } else {
      /* open the config file and retrieve cuts */
      FILE *fpConfig = fopen("AD0_TimeSlewing_DA.config","r");
      int res = fscanf(fpConfig,"%d %d %d %d %d %d %d",&kMinEvents,&kNPreClocks,&kNPostClocks,&kNPedSigma,&kNBinsCharge,&kMinTime,&kMaxTime);
      if(res!=7) {
	    printf("Failed to get values from Config file (AD0_TimeSlewing_DA.config): wrong file format - 4 integers are expected - \n");
       	    kMinEvents = 50000;
            kNPreClocks = 1;  
            kNPostClocks = 10;
            kNPedSigma = 3;
	    kNBinsCharge = 1000;
	    kMinTime = 1740;
	    kMaxTime = 2252;
      }
      fclose(fpConfig);
  }
  printf("Minimal number of events requested = %d; Integrating charge in [-%d + %d] LHC clocks around maximum; Number of pedestal sigma cut = %d; Number of bins on charge axis = %d; Time range [%d , %d]\n",
          kMinEvents, kNPreClocks, kNPostClocks, kNPedSigma, kNBinsCharge,kMinTime,kMaxTime);
  

  // Online pedestals (using FEE channel numbering), 
  // from DAQ_DB AD0_Pedestals_On.dat:
  Float_t adcMean[32];
  Float_t adcSigma[32];
  Float_t pedMean[32];
  Float_t pedSigma[32];
  
  status = daqDA_DB_getFile("AD0_Pedestals_On.dat","AD0_Pedestals_On.dat");
  if (status) {
      printf("Failed to get AD0_Pedestals_On.dat from DAQ DB, status=%d\n", status);
      printf("Take default pedestals \n");
      for(Int_t j=0; j<32; j++) {pedMean[j] = 0; pedSigma[j] = 0.1;}
      
  } else {
      /* open the pedestal file and retrieve values */
      FILE *fpPedestals = fopen("AD0_Pedestals_On.dat","r");
      for(Int_t j=0; j<32; j++) {
	  Int_t res = fscanf(fpPedestals,"%f %f %f %f",&pedMean[j], &pedSigma[j], &adcMean[j], &adcSigma[j]);
	  if (res != 4) printf("Failed to get values from file: wrong file format - 4 floats are expected - \n");
		}
      fclose(fpPedestals);
  }
  printf("Used pedestals:\n");
  for(Int_t iChannel = 0; iChannel < 16; ++iChannel) {
    	for(Int_t integrator = 0; integrator < 2; ++integrator){
    		if(integrator == 0)printf("ChOn = %d, Int = %d, Pedestal = %.3f, Width = %3f,", iChannel,integrator, pedMean[iChannel+16*integrator],pedSigma[iChannel+16*integrator]);
		else printf(" Int = %d, Pedestal = %.3f, Width = %3f\n", integrator, pedMean[iChannel+16*integrator],pedSigma[iChannel+16*integrator]);	
			}
		}
  
//___________________________________________________
// Book HISTOGRAMS - dynamics of p-p collisions -
      
  char     TimeVsChargeName[20]; 
  TH2F     *hTimeVsCharge[16];
  
  char  texte[20];
  for (Int_t i=0; i<16; i++) {
       	sprintf(TimeVsChargeName,"hTimeVsCharge%d",i);
       	sprintf(texte,"hTimeVsCharge ch%d",i);
       	hTimeVsCharge[i]  = new TH2F(TimeVsChargeName,texte,kNBinsCharge,-4,0,kMaxTime-kMinTime,kMinTime,kMaxTime);
        }
       
   
//___________________________________________________ 

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
	   
	   Float_t adc[16], time[16];	
           for(Int_t iChannel=0; iChannel<16; iChannel++) {
	   	adc[iChannel] = 0.0;
	   	time[iChannel] = rawStream->GetTime(iChannel);
	   	Float_t maxadc = 0;
      		Int_t imax = -1;
     		Float_t adcPedSub[21];
      		for(Int_t iClock=0; iClock<21; iClock++){
			Bool_t iIntegrator = rawStream->GetIntegratorFlag(iChannel,iClock);
			Int_t k = iChannel+16*iIntegrator;

			adcPedSub[iClock] = rawStream->GetPedestal(iChannel,iClock) - pedMean[k];
			if(adcPedSub[iClock] <= kNPedSigma*pedSigma[k]) {
	  			adcPedSub[iClock] = 0;
	  			continue;
				}
			if(adcPedSub[iClock] > maxadc) {
	 			maxadc = adcPedSub[iClock];
	  			imax   = iClock;
				}
      			}
        	if (imax != -1) {
			Int_t start = imax - kNPreClocks;
			if (start < 0) start = 0;
			Int_t end = imax + kNPostClocks;
			if (end > 20) end = 20;
			for(Int_t iClock = start; iClock <= end; iClock++) adc[iChannel] += adcPedSub[iClock];
			}
		else adc[iChannel] = 0.1;
		hTimeVsCharge[iChannel]->Fill(TMath::Log10(1/adc[iChannel]),time[iChannel]);
		
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
  if(neventsPhysics < kMinEvents){
  	printf("Number of events too low to evaluate time slewing splines, exiting\n");
	return status;
  	}  
//___________________________________________________________________________
//  Computes mean values, converts FEE channels into Offline AliRoot channels
//  and dumps the ordered values into the output text file for SHUTTLE
	
  TH1D *hTimeSlice = 0x0;
  TH1F *hMeanTimeVsCharge[16];

  for(Int_t i=0; i<16; i++){
	TString MeanTimeVsChargeName = "hMeanTimeVsCharge";
	MeanTimeVsChargeName += kOfflineChannel[i];
	hMeanTimeVsCharge[i] = new TH1F(MeanTimeVsChargeName.Data(),MeanTimeVsChargeName.Data(),kNBinsCharge,-4,0);
	hMeanTimeVsCharge[i]->GetXaxis()->SetTitle("Log10(1/charge) [ADC counts]");
	hMeanTimeVsCharge[i]->GetYaxis()->SetTitle("Mean leading time [TDC counts]");
	
	for(Int_t j=0; j<kNBinsCharge;j++){
  		hTimeSlice = hTimeVsCharge[i]->ProjectionY("hTimeSlice",j+1,j+1);
		if(hTimeSlice->GetEntries()<100 && j>kNBinsCharge/2){
			hMeanTimeVsCharge[i]->SetBinContent(j+1,hMeanTimeVsCharge[i]->GetBinContent(j));
			hMeanTimeVsCharge[i]->SetBinError(j+1,hMeanTimeVsCharge[i]->GetBinError(j));
			}
		else{
			hMeanTimeVsCharge[i]->SetBinContent(j+1,hTimeSlice->GetMean());
        		if(hTimeSlice->GetEntries()>0)hMeanTimeVsCharge[i]->SetBinError(j+1,1/TMath::Power(hTimeSlice->GetEntries(),2));
			}
		
  		}
	}
	
  TSpline3 *fTimeSlewingSpline[16];
  for(Int_t i=0; i<16; i++){
  	Int_t offlineChannel = kOfflineChannel[i];
	TString TimeSlewingSplineName = "hTimeSlewingSpline";
	TimeSlewingSplineName += offlineChannel;
	fTimeSlewingSpline[offlineChannel] = new TSpline3(hMeanTimeVsCharge[i]);
	fTimeSlewingSpline[offlineChannel]->SetName(TimeSlewingSplineName.Data());	
	}
  /*/	
  TF1 *fTimeSlewingFit[16];
  for(Int_t i=0; i<16; i++){
	TString TimeSlewingFitName = "hTimeSlewingFit";
	TimeSlewingFitName += i;
	fTimeSlewingFit[i] = new TF1(TimeSlewingFitName.Data(),"[0]+[1]*x^[2]+((x<0.04)?0:[3]*(x-0.04)+[4]*(x-0.04)^2+[5]*(x-0.04)^3)",1,4000);
	fTimeSlewingFit[i]->SetParameter(0,650);
	fTimeSlewingFit[i]->SetParameter(1,40);
	fTimeSlewingFit[i]->SetParameter(2,-0.5);
	fTimeSlewingFit[i]->SetLineColor(kMagenta);
	hMeanTimeVsCharge[i]->Fit(TimeSlewingFitName.Data(),"R0");
	}
  /*/
   
//________________________________________________________________________
// Write root file with splines

  TList *fListSplines = new TList();

  TFile *splineFile = new TFile("AD0_SlewingSplines.root","RECREATE");
  for(Int_t i=0; i<16; i++) fListSplines->Add(fTimeSlewingSpline[i]); 
  //for(Int_t i=0; i<16; i++) fListSplines->Add(hTimeVsCharge[i]);
  //for(Int_t i=0; i<16; i++) fListSplines->Add(hMeanTimeVsCharge[i]);
  fListSplines->Write("fListSplines",1);
  splineFile->Close();
  
//________________________________________________________________________
   
  /* export result file to FES */
  status=daqDA_FES_storeFile("AD0_SlewingSplines.root","AD0da_slewing");
  if (status)    {
      printf("Failed to export file : %d\n",status);
      return -1; }

      
  return status;
}
    
