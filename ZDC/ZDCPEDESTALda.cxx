/*

DAcase2.c

This program connects to the DAQ data source passed as argument
and populates local "./result.txt" file with the ids of events received
during the run.

The program exits when being asked to shut down (daqDA_checkshutdown)
or End of Run event.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalon pedestal runs
contact: Chiara.Oppedisano@cern.ch

*/

#include <stdio.h>
#include <Riostream.h>

// DATE
#include <daqDA.h>
#include <event.h>
#include <monitor.h>

//ROOT
#include <TRandom.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>

//AliRoot
#include <AliRawReaderDate.h>
#include <AliZDCRawStream.h>


/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {

  TH1F::AddDirectory(0);
  // --- Histograms for ADC pedestals 
  //     [22 signal channels x 2 gain chains + 2 reference PTMs]
  //
  TH1F *hPed[44], *hPedOutOfTime[44];
  TH2F *hPedCorr[44];
  //
  char namhist1[50], namhist2[50], namhist3[50];
  for(Int_t j=0; j<44; j++){
     if(j<10){
       sprintf(namhist1,"PedZN1_%d",j);
       sprintf(namhist2,"PedZN1OutOfTime_%d",j);
       sprintf(namhist3,"PedCorrZN1_%d",j);
     }
     else if(j>=10 && j<20){
       sprintf(namhist1,"PedZP1_%d",j-10);
       sprintf(namhist2,"PedZP1OutOfTime_%d",j-10);
       sprintf(namhist3,"PedCorrZP1_%d",j-10);
     }
     else if(j>=20 && j<24){
       sprintf(namhist1,"PedZEM_%d",j-20);
       sprintf(namhist2,"PedZEMOutOfTime_%d",j-20);
       sprintf(namhist3,"PedCorrZEM_%d",j-20);
     }
     else if(j>=24 && j<33){
       sprintf(namhist1,"PedZN2_%d",j-24);
       sprintf(namhist2,"PedZN2OutOfTime_%d",j-24);
       sprintf(namhist3,"PedCorrZN2_%d",j-24);
     }
     else if(j>=33 && j<43){
       sprintf(namhist1,"PedZP2_%d",j-33);
       sprintf(namhist2,"PedZP2OutOfTime_%d",j-33);
       sprintf(namhist3,"PedCorrZP2_%d",j-33);
     }
     hPed[j] = new TH1F(namhist1, namhist1, 100,0., 200.);
     hPedOutOfTime[j] = new TH1F(namhist2, namhist2, 100,0., 200.);
     hPedCorr[j] = new TH2F(namhist3,namhist3,100,0.,200.,100,0.,200.);
  }

  int status;
  
  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
  }
  

  /* define data source : this is argument 1 */  
  status = monitorSetDataSource( argv[1] );
  if(status!=0) {
    printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* declare monitoring program */
  status = monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);
  

  /* log start of process */
  printf("ZDC PEDESTAL monitoring program started\n");  
  
  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  struct equipmentStruct *equipment;
  int *eventEnd;
  int *eventData;
  int *equipmentEnd;
  int *equipmentData;
  int *equipmentID;

  struct eventHeaderStruct *event;
  eventTypeType eventT;
  Int_t iev=0;
  
  /* main loop (infinite) */
  for(;;) {
  
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

    iev++; 

    /* use event - here, just write event id to result file */
    eventT=event->eventType;
    
    if(eventT==PHYSICS_EVENT){
      //fprintf(fp,"Run #%lu, event size: %lu, BC:%u, Orbit:%u, Period:%u\n",
      //
      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      const AliRawDataHeader* header = reader->GetDataHeader();
      if(header) {
         UChar_t message = header->GetL1TriggerMessage();
	 if(message & 0x40000){ // DEDICATED PEDESTAL RUN
	    printf("\t L1 message -> PEDESTAL raw data\n");
	    continue;
	 }
	 else{
	    printf("\t L1 message -> NO PEDESTAL raw data found\n");
	    return -1;
	 }
      }
      //Commented until we won't have Raw Data Header...
      /*else{
         //printf("\t ERROR! No Raw Data Header found!!!\n");
	 //return -1;
      }*/
      //
    
      AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
      //
      if (!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
      Int_t counter=0;
      Int_t RawADC[44], RawADCoot[44];
      for(Int_t j=0; j<44; j++){
         RawADC[j]=0;
         RawADCoot[j]=0;
      }
      while(rawStreamZDC->Next()){
        Int_t index=-1;
        if(rawStreamZDC->IsADCDataWord()){
    	  if(rawStreamZDC->GetSector(0)==1 || rawStreamZDC->GetSector(0)==2){ // *** ZN1, ZP1
	    index = 10*(rawStreamZDC->GetSector(0)-1)+rawStreamZDC->GetSector(1)+5*rawStreamZDC->GetADCGain();
	  }
	  else if(rawStreamZDC->GetSector(0)==3){ // *** ZEM 
	    index = 10*(rawStreamZDC->GetSector(0)-1)+(rawStreamZDC->GetSector(1)-1)+2*rawStreamZDC->GetADCGain();
	  }
	  else if(rawStreamZDC->GetSector(0)==4 || rawStreamZDC->GetSector(0)==5){ // *** ZN2, ZP2
	    index = 10*(rawStreamZDC->GetSector(0)-2)+rawStreamZDC->GetSector(1)+5*rawStreamZDC->GetADCGain()+4;
	  }
	  if(counter<44){
	    hPed[index]->Fill(rawStreamZDC->GetADCValue()); 
	    RawADC[counter] = rawStreamZDC->GetADCValue();
	  }
	  else{ 
	    hPedOutOfTime[index]->Fill(rawStreamZDC->GetADCValue());
	    RawADCoot[counter-44] = rawStreamZDC->GetADCValue();
	  }
	  counter++;
         }//IsADCDataWord()
	 //
	 if(counter == 88){ // Last ADC channel
           for(Int_t k=0; k<44; k++){
              hPedCorr[k]->Fill(RawADCoot[k], RawADC[k]);
           }
	 }
       }
       //
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
  
  /* Analysis of the histograms */
  //
  FILE *fileShuttle;
  fileShuttle = fopen("ZDCPedestal.dat","w");
  //
  Float_t MeanPed[44], MeanPedWidth[44], 
   	MeanPedOOT[44], MeanPedWidthOOT[44],
   	CorrCoeff0[44], CorrCoeff1[44];
  // --- Out-of-time pedestals
  TF1 *ADCfunc[44];
  for(Int_t i=0; i<44; i++){
     hPed[i]->Fit("gaus","Q");
     ADCfunc[i] = hPed[i]->GetFunction("gaus");
     MeanPed[i] = ADCfunc[i]->GetParameter(1);
     MeanPedWidth[i] = ADCfunc[i]->GetParameter(2);
     fprintf(fileShuttle,"\t%f\t%f\n",MeanPed[i],MeanPedWidth[i]);
     //printf("\t MeanPed[%d] = %f\n",i, MeanPed[i]);
  }
  // --- Out-of-time pedestals
  TF1 *ADCootfunc[44];
  for(Int_t i=0; i<44; i++){
     hPedOutOfTime[i]->Fit("gaus","Q");
     ADCootfunc[i] = hPedOutOfTime[i]->GetFunction("gaus");
     MeanPedOOT[i] = ADCootfunc[i]->GetParameter(1);
     MeanPedWidthOOT[i] = ADCootfunc[i]->GetParameter(2);
     fprintf(fileShuttle,"\t%f\t%f\n",MeanPedOOT[i],MeanPedWidthOOT[i]);
     //printf("\t MeanPedOOT[%d] = %f\n",i, MeanPedOOT[i]);
  }
  //
  // --- Fit of correlations
  TProfile* hPedCorrProf[44];
  TF1 *ffunc[44];
  char namhist4[50];
  for(int i=0;i<44;i++) {
     sprintf(namhist4,"ADCvsOOT%d_Prof",i);
     hPedCorrProf[i] = hPedCorr[i]->ProfileX(namhist4,-1,-1,"S");
     hPedCorrProf[i]->SetName(namhist4);
     hPedCorrProf[i]->Fit("pol1","Q");
     ffunc[i] = hPedCorrProf[i]->GetFunction("pol1");
     CorrCoeff0[i] = ffunc[i]->GetParameter(0);
     CorrCoeff1[i] = ffunc[i]->GetParameter(1);
     fprintf(fileShuttle,"\t%f\t%f\n",CorrCoeff0[i],CorrCoeff1[i]);
     //printf("\t CorrCoeff0[%d] = %f, CorrCoeff1[%d] = %f\n",i, CorrCoeff0[i], i, CorrCoeff1[i]);
  }    
  //						       
  fclose(fileShuttle);
  

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);


  return status;
}
