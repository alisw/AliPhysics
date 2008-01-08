/********************************************************************************
*                                                                               *
* VZERO Detector Algorithm used for extracting calibration parameters           *
*                                                                               *
* This program reads the DDL data file passed as argument using the monitoring  *
* library.                                                                      *
* It computes calibration parameters, populates local "./V0_Ped_Width_Gain.dat" *           
* file and exports it to the FES                                                *
*                                                                               *
* The program reports about its processing progress.                            *
*                                                                               *
*********************************************************************************/

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

#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>


/* Main routine --- Arguments: list of DATE raw data files */
      
int main(int argc, char **argv) {

  int status;

  printf(" argc = %d, argv = %s \n",argc, &(**argv));

  Int_t    kHighCut = 50; // high cut on pedestal distribution - to be tuned
  Int_t    kLowCut  = 30; // low cut on signal distribution - to be tuned
  Double_t ADCmean[64];
  Double_t PEDmean[64];
  Double_t PEDsigma[64];
  
//___________________________________________________
// Book HISTOGRAMS - dynamics of p-p collisions -
      
  char     ADCname[6]; 
  char     PEDname[6]; 
  TH1F     *hADCname[64];
  TH1F     *hPEDname[64];  
  
  char  texte[12];
  for (Int_t i=0; i<64; i++) {
       sprintf(ADCname,"hADC%d",i);
       sprintf(texte,"ADC cell%d",i);
       hADCname[i]  = new TH1F(ADCname,texte,1024,0,2047);
       sprintf(PEDname,"hPED%d",i);
       sprintf(texte,"PED cell%d",i);
       hPEDname[i]  = new TH1F(PEDname,texte,1024,0,2047);}
//___________________________________________________ 
  
  
  /* log start of process */
  printf("VZERO DA program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2)   {
      printf("Wrong number of arguments\n");
      return -1;}

  /* open result file to be exported to FES */
  FILE *fp=NULL;
  fp=fopen("./V0_Ped_Width_Gain.dat","a");
  if (fp==NULL) {
      printf("Failed to open result file\n");
      return -1;}

  /* open log file to inform user */
  FILE *flog=NULL;
  flog=fopen("./V00log.txt","a");
  if (flog==NULL) {
      printf("Failed to open log file\n");
      return -1;  }
    
  /* report progress */
  daqDA_progressReport(10);


  /* init counters on events */
  int nevents_physics=0;
  int nevents_total=0;

  /* read the data files, considering n files */
  
  int n;
  
  for (n=1;n<argc;n++) {
   
    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
        printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
        return -1; }

    /* report progress - here indexed on the number of files */
    daqDA_progressReport(10+80*n/argc);

    /* read the data file */
    for(;;) {
        struct eventHeaderStruct *event;
        eventTypeType eventT;

        /* get next event */
        status=monitorGetEventDynamic((void **)&event);
        if (status==MON_ERR_EOF) break; /* end of monitoring file has been reached */
        if (status!=0) {
            printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
            return -1; }

        /* retry if got no event */
        if (event==NULL) break;
        
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
        
             fprintf(flog,"Run #%lu, event size: %lu, BC:%u, Orbit:%u, Period:%u\n",
                 (unsigned long)event->eventRunNb,
                 (unsigned long)event->eventSize,
                 EVENT_ID_GET_BUNCH_CROSSING(event->eventId),
                 EVENT_ID_GET_ORBIT(event->eventId),
                 EVENT_ID_GET_PERIOD(event->eventId) );
		 
	     AliRawReader *rawReader = new AliRawReaderDate((void*)event);
             rawReader->RequireHeader(kFALSE);	 
	     AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader);    
	     rawStream->Next();	
	     for(Int_t i=0; i<64; i++) {
                hADCname[i]->Fill(rawStream->GetADC(i)); 
		for(Int_t j=0; j<21; j++) {
		    if(j==10) continue;
                    hPEDname[i]->Fill(rawStream->GetPedestal(i,j)); }
             }    
             delete rawStream;
             rawStream = 0x0;      
             delete rawReader;
             rawReader = 0x0;	     	 						         
        } // end of switch on event type 
	
        nevents_total++;
        /* free resources */
        free(event);

    }  // loop over events
    
   }  // loop over data files
//________________________________________________________________________
//  Computes mean values, dumps them into the output text file
	
  for(Int_t i=0; i<64; i++) {
      hPEDname[i]->GetXaxis()->SetRange(0,kHighCut);
      PEDmean[i]  = hPEDname[i]->GetMean(); 
      PEDsigma[i] = hPEDname[i]->GetRMS(); 
      hADCname[i]->GetXaxis()->SetRange(kLowCut,2047);
      ADCmean[i] = hADCname[i]->GetMean() ; 
      fprintf(fp," %.3f %.3f %.3f\n",PEDmean[i],PEDsigma[i],ADCmean[i]);
  } 

//________________________________________________________________________
// Write root file with histos for users further check - just in case - 

  TFile *histoFile = new TFile("VZERO_histos.root","RECREATE");
    
  for (Int_t i=0; i<64; i++) {
       hADCname[i]->Write(); 
       hPEDname[i]->Write(); }

  histoFile->Close(); 
  
//________________________________________________________________________
   
  /* write report */
  fprintf(flog,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);
  printf("Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);
  
  /* close result and log files */
  fclose(fp);
  fclose(flog); 
  
  /* report progress */
  daqDA_progressReport(90);


  /* export result file to FES */
  status=daqDA_FES_storeFile("./V0_Ped_Width_Gain.dat","V00da_results");
  if (status)    {
      printf("Failed to export file : %d\n",status);
      return -1; }

  /* report progress */
  daqDA_progressReport(100);
  
  return status;
}
