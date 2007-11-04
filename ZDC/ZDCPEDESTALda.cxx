/*

DAcase2.c

This program connects to the DAQ data source passed as argument
and populates local "./result.txt" file with the ids of events received
during the run.

The program exits when being asked to shut down (daqDA_checkshutdown)
or End of Run event.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalone pedestal runs

Contact: Chiara.Oppedisano@to.infn.it
Link: /afs/cern.ch/user/c/chiarao/public/RawPed.date
Run Type: STANDALONE_PEDESTAL_RUN
DA Type: MON
Number of events needed: no constraint (tipically ~10^3)
Input Files: 
Output Files: ZDCPedestal.dat
Trigger Types Used: Standalone Trigger

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
  //     [22 signal channels + 2 reference PTMs]  x 2 gain chains
  //
  int const kNChannels = 24;
  TH1F *hPedhg[kNChannels], *hPedOutOfTimehg[kNChannels];
  TH2F *hPedCorrhg[kNChannels];
  TH1F *hPedlg[kNChannels], *hPedOutOfTimelg[kNChannels];
  TH2F *hPedCorrlg[kNChannels];
  //
  char namhist1hg[50], namhist2hg[50], namhist3hg[50];
  char namhist1lg[50], namhist2lg[50], namhist3lg[50];
  for(Int_t j=0; j<kNChannels; j++){
     if(j<5){ // ZN1
       sprintf(namhist1hg,"PedZN1hg_%d",j);
       sprintf(namhist2hg,"PedZN1hgOutOfTime_%d",j);
       sprintf(namhist3hg,"PedCorrZN1hg_%d",j);
       //
       sprintf(namhist1lg,"PedZN1lg_%d",j);
       sprintf(namhist2lg,"PedZN1lgOutOfTime_%d",j);
       sprintf(namhist3lg,"PedCorrZN1lg_%d",j);
     }
     else if(j>=5 && j<10){ // ZP1
       sprintf(namhist1hg,"PedZP1hg_%d",j-5);
       sprintf(namhist2hg,"PedZP1hgOutOfTime_%d",j-5);
       sprintf(namhist3hg,"PedCorrZP1hg_%d",j-5);
       //
       sprintf(namhist1lg,"PedZP1lg_%d",j-5);
       sprintf(namhist2lg,"PedZP1lgOutOfTime_%d",j-5);
       sprintf(namhist3lg,"PedCorrZP1lg_%d",j-5);       
     }
     else if(j>=10 && j<12){ // ZEM
       sprintf(namhist1hg,"PedZEMhg_%d",j-10);
       sprintf(namhist2hg,"PedZEMhgOutOfTime_%d",j-10);
       sprintf(namhist3hg,"PedCorrZEMhg_%d",j-10);
       //
       sprintf(namhist1lg,"PedZEMlg_%d",j-10);
       sprintf(namhist2lg,"PedZEMlgOutOfTime_%d",j-10);
       sprintf(namhist3lg,"PedCorrZEMlg_%d",j-10);
     }
     else if(j>=12 && j<17){ // ZN2
       sprintf(namhist1hg,"PedZN2hg_%d",j-12);
       sprintf(namhist2hg,"PedZN2hgOutOfTime_%d",j-12);
       sprintf(namhist3hg,"PedCorrZN2hg_%d",j-12);
       //
       sprintf(namhist1lg,"PedZN2lg_%d",j-12);
       sprintf(namhist2lg,"PedZN2lgOutOfTime_%d",j-12);
       sprintf(namhist3lg,"PedCorrZN2lg_%d",j-12);
     }
     else if(j>=17 && j<22){ // ZP2
       sprintf(namhist1hg,"PedZP2hg_%d",j-17);
       sprintf(namhist2hg,"PedZP2hgOutOfTime_%d",j-17);
       sprintf(namhist3hg,"PedCorrZP2hg_%d",j-17);
       //
       sprintf(namhist1lg,"PedZP2lg_%d",j-17);
       sprintf(namhist2lg,"PedZP2lgOutOfTime_%d",j-17);
       sprintf(namhist3lg,"PedCorrZP2lg_%d",j-17);
     }
     else if(j>=22 && j<24){ //Reference PMs
       sprintf(namhist1hg,"PedRefhg_%d",j-22);
       sprintf(namhist2hg,"PedRefhgOutOfTime_%d",j-22);
       sprintf(namhist3hg,"PedCorrRefhg_%d",j-22);
       //
       sprintf(namhist1lg,"PedReflg_%d",j-22);
       sprintf(namhist2lg,"PedReflgOutOfTime_%d",j-22);
       sprintf(namhist3lg,"PedCorrReflg_%d",j-22);
     }
     // --- High gain chain histos
     hPedhg[j] = new TH1F(namhist1hg, namhist1hg, 100,0., 200.);
     hPedOutOfTimehg[j] = new TH1F(namhist2hg, namhist2hg, 100,0., 200.);
     hPedCorrhg[j] = new TH2F(namhist3hg,namhist3hg,100,0.,200.,100,0.,200.);
     // --- Low gain chain histos
     hPedlg[j] = new TH1F(namhist1lg, namhist1lg, 100,0., 600.);
     hPedOutOfTimelg[j] = new TH1F(namhist2lg, namhist2lg, 100,0., 600.);
     hPedCorrlg[j] = new TH2F(namhist3lg,namhist3lg,100,0.,600.,100,0.,600.);
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
         UChar_t message = header->GetAttributes();
	 if(message & 0x20){ // DEDICATED PEDESTAL RUN
	    printf("\t STANDALONE_PEDESTAL_RUN raw data found\n");
	    continue;
	 }
	 else{
	    printf("\t NO STANDALONE_PEDESTAL_RUN raw data found\n");
	    return -1;
	 }
      }
      //Commented until we won't have true Raw Data Header...
      //else{
      //   printf("\t ATTENTION! No Raw Data Header found!!!\n");
	// return -1;
      //}
      //
      AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
      //
      if (!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
      Int_t counter=0;
      Int_t RawADChg[kNChannels], RawADCoothg[kNChannels];
      Int_t RawADClg[kNChannels], RawADCootlg[kNChannels];
      for(Int_t j=0; j<kNChannels; j++){
         RawADChg[j]=0; RawADCoothg[j]=0;
         RawADClg[j]=0; RawADCootlg[j]=0;
      }
      while(rawStreamZDC->Next()){
        Int_t index=-1;
        if(rawStreamZDC->IsADCDataWord()){
	  if(rawStreamZDC->GetSector(1)!=5){ // Physics signals
    	    if(rawStreamZDC->GetSector(0)==1) index = rawStreamZDC->GetSector(1); // *** ZN1
	    else if(rawStreamZDC->GetSector(0)==2) index = rawStreamZDC->GetSector(1)+5; // *** ZP1
	    else if(rawStreamZDC->GetSector(0)==3) index = rawStreamZDC->GetSector(1)+9; // *** ZEM
	    else if(rawStreamZDC->GetSector(0)==4) index = rawStreamZDC->GetSector(1)+12; // *** ZN2
	    else if(rawStreamZDC->GetSector(0)==5) index = rawStreamZDC->GetSector(1)+17; // *** ZP2
	  }
	  else{ // Reference PMs
	    index = (rawStreamZDC->GetSector(0)-1)/3+22;
	  }
	  //
	  /*printf("\t counter %d index %d det %d quad %d res %d ADC %d\n", counter, index,
	  	rawStreamZDC->GetSector(0), rawStreamZDC->GetSector(1), 
		rawStreamZDC->GetADCGain(), rawStreamZDC->GetADCValue());
	  */
	  //
	  if(counter<2*kNChannels){ // --- In-time pedestals (1st 48 raw data)
	    if(rawStreamZDC->GetADCGain()==0){ 
	      hPedhg[index]->Fill(rawStreamZDC->GetADCValue()); 
	      RawADChg[index] = rawStreamZDC->GetADCValue();
	    }
	    else{
	      hPedlg[index]->Fill(rawStreamZDC->GetADCValue()); 
	      RawADClg[index] = rawStreamZDC->GetADCValue();
	    }
	  }
	  else{  // --- Out-of-time pedestals
	    if(rawStreamZDC->GetADCGain()==0){
	      hPedOutOfTimehg[index]->Fill(rawStreamZDC->GetADCValue());
	      RawADCoothg[index] = rawStreamZDC->GetADCValue();
	    }
	    else{
	      hPedOutOfTimelg[index]->Fill(rawStreamZDC->GetADCValue());
	      RawADCootlg[index] = rawStreamZDC->GetADCValue();
	    }
	  }
	  counter++;
         }//IsADCDataWord()
	 //
	 if(counter == 4*kNChannels){ // Last ADC channel -> Filling correlation histos
           for(Int_t k=0; k<kNChannels; k++){
	    hPedCorrhg[k]->Fill(RawADCoothg[k], RawADChg[k]);
	    hPedCorrlg[k]->Fill(RawADCootlg[k], RawADClg[k]);
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
  Float_t MeanPed[2*kNChannels], MeanPedWidth[2*kNChannels], 
   	MeanPedOOT[2*kNChannels], MeanPedWidthOOT[2*kNChannels],
   	CorrCoeff0[2*kNChannels], CorrCoeff1[2*kNChannels];
  // --- In-time pedestals
  TF1 *ADCfunchg[kNChannels];
  for(Int_t i=0; i<kNChannels; i++){
     hPedhg[i]->Fit("gaus","Q");
     ADCfunchg[i] = hPedhg[i]->GetFunction("gaus");
     MeanPed[i] = ADCfunchg[i]->GetParameter(1);
     MeanPedWidth[i] = ADCfunchg[i]->GetParameter(2);
     fprintf(fileShuttle,"\t%f\t%f\n",MeanPed[i],MeanPedWidth[i]);
     //printf("\t MeanPed[%d] = %f\n",i, MeanPed[i]);
  }
  TF1 *ADCfunclg[kNChannels];
  for(Int_t i=0; i<kNChannels; i++){
     hPedlg[i]->Fit("gaus","Q");
     ADCfunclg[i] = hPedlg[i]->GetFunction("gaus");
     MeanPed[i+kNChannels] = ADCfunclg[i]->GetParameter(1);
     MeanPedWidth[i+kNChannels] = ADCfunclg[i]->GetParameter(2);
     fprintf(fileShuttle,"\t%f\t%f\n",MeanPed[i+kNChannels],MeanPedWidth[i+kNChannels]);
     //printf("\t MeanPed[%d] = %f\n",i+kNChannels, MeanPed[i+kNChannels]);
  }
  // --- Out-of-time pedestals
  TF1 *ADCootfunchg[kNChannels];
  for(Int_t i=0; i<kNChannels; i++){
     hPedOutOfTimehg[i]->Fit("gaus","Q");
     ADCootfunchg[i] = hPedOutOfTimehg[i]->GetFunction("gaus");
     MeanPedOOT[i] = ADCootfunchg[i]->GetParameter(1);
     MeanPedWidthOOT[i] = ADCootfunchg[i]->GetParameter(2);
     fprintf(fileShuttle,"\t%f\t%f\n",MeanPedOOT[i],MeanPedWidthOOT[i]);
     //printf("\t MeanPedOOT[%d] = %f\n",i, MeanPedOOT[i]);
  }
  TF1 *ADCootfunclg[kNChannels];
  for(Int_t i=0; i<kNChannels; i++){
     hPedOutOfTimelg[i]->Fit("gaus","Q");
     ADCootfunclg[i] = hPedOutOfTimelg[i]->GetFunction("gaus");
     MeanPedOOT[i+kNChannels] = ADCootfunclg[i]->GetParameter(1);
     MeanPedWidthOOT[i+kNChannels] = ADCootfunclg[i]->GetParameter(2);
     fprintf(fileShuttle,"\t%f\t%f\n",MeanPedOOT[i+kNChannels],MeanPedWidthOOT[i+kNChannels]);
     //printf("\t MeanPedOOT[%d] = %f\n",i+kNChannels, MeanPedOOT[i+kNChannels]);
  }
  //
  // --- Correlations
  TProfile *hPedCorrProfhg[kNChannels], *hPedCorrProflg[kNChannels];
  TF1 *ffunchg[kNChannels], *ffunclg[kNChannels];
  char namhist4[50];
  for(int i=0;i<kNChannels;i++) {
     sprintf(namhist4,"ADCHRvsOOT%d_Prof",i);
     hPedCorrProfhg[i] = hPedCorrhg[i]->ProfileX(namhist4,-1,-1,"S");
     hPedCorrProfhg[i]->SetName(namhist4);
     hPedCorrProfhg[i]->Fit("pol1","Q");
     ffunchg[i] = hPedCorrProfhg[i]->GetFunction("pol1");
     CorrCoeff0[i] = ffunchg[i]->GetParameter(0);
     CorrCoeff1[i] = ffunchg[i]->GetParameter(1);
     fprintf(fileShuttle,"\t%f\t%f\n",CorrCoeff0[i],CorrCoeff1[i]);
     //printf("\t CorrCoeff0[%d] = %f, CorrCoeff1[%d] = %f\n",i, CorrCoeff0[i], i, CorrCoeff1[i]);
  }    
  for(int i=0;i<kNChannels;i++) {
     sprintf(namhist4,"ADCLRvsOOT%d_Prof",i);
     hPedCorrProflg[i] = hPedCorrlg[i]->ProfileX(namhist4,-1,-1,"S");
     hPedCorrProflg[i]->SetName(namhist4);
     hPedCorrProflg[i]->Fit("pol1","Q");
     ffunclg[i] = hPedCorrProflg[i]->GetFunction("pol1");
     CorrCoeff0[i+kNChannels] = ffunclg[i]->GetParameter(0);
     CorrCoeff1[i+kNChannels] = ffunclg[i]->GetParameter(1);
     fprintf(fileShuttle,"\t%f\t%f\n",CorrCoeff0[i+kNChannels],CorrCoeff1[i+kNChannels]);
     //printf("\t CorrCoeff0[%d] = %f, CorrCoeff1[%d] = %f\n",
     //		i+kNChannels, CorrCoeff0[i+kNChannels], i+kNChannels, CorrCoeff1[i+kNChannels]);
  }    
  //						       
  fclose(fileShuttle);
  //
  

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);


  return status;
}
