/*

This program reads the DAQ data files passed as argument using the monitoring library.


The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalone pedestal runs

Contact: Chiara.Oppedisano@to.infn.it
Link: 
Run Type: STANDALONE_PEDESTAL_RUN
DA Type: LDC
Number of events needed: no constraint (tipically ~10^3)
Input Files:  
Output Files: ZDCPedestal.dat, ZDCChMapping.dat
Trigger Types Used: Standalone Trigger

*/
#define PEDDATA_FILE  "ZDCPedestal.dat"
#define MAPDATA_FILE  "ZDCChMapping.dat"

#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>

// DATE
#include <event.h>
#include <monitor.h>
#include <daqDA.h>

//ROOT
#include <TRandom.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitter.h>

//AliRoot
#include <AliRawReaderDate.h>
#include <AliRawEventHeaderBase.h>
#include <AliZDCRawStream.h>


/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {
  
//  TVirtualFitter::SetDefaultFitter("Minuit2");

  int status = 0;
  int const kNChannels = 24;

  /* log start of process */
  printf("\n ZDC PEDESTAL program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // --- Histograms for ADC pedestals 
  //     [22 signal channels + 2 reference PTMs]  x 2 gain chains
  //
  TH1F::AddDirectory(0);
  //
  TH1F *hPedhg[kNChannels], *hPedOutOfTimehg[kNChannels];
  TH2F *hPedCorrhg[kNChannels];
  TH1F *hPedlg[kNChannels], *hPedOutOfTimelg[kNChannels];
  TH2F *hPedCorrlg[kNChannels];
  //
  char namhist1hg[50], namhist2hg[50], namhist3hg[50];
  char namhist1lg[50], namhist2lg[50], namhist3lg[50];
  for(Int_t j=0; j<kNChannels; j++){
     if(j<=4){ // ZNC
       sprintf(namhist1hg,"PedZNChg_%d",j);
       sprintf(namhist2hg,"PedZNChgOutOfTime_%d",j);
       sprintf(namhist3hg,"PedCorrZNChg_%d",j);
       //
       sprintf(namhist1lg,"PedZNClg_%d",j);
       sprintf(namhist2lg,"PedZNClgOutOfTime_%d",j);
       sprintf(namhist3lg,"PedCorrZNClg_%d",j);
     }
     else if(j>=5 && j<=9){ // ZPC
       sprintf(namhist1hg,"PedZPChg_%d",j-5);
       sprintf(namhist2hg,"PedZPChgOutOfTime_%d",j-5);
       sprintf(namhist3hg,"PedCorrZPChg_%d",j-5);
       //
       sprintf(namhist1lg,"PedZPClg_%d",j-5);
       sprintf(namhist2lg,"PedZPClgOutOfTime_%d",j-5);
       sprintf(namhist3lg,"PedCorrZPClg_%d",j-5);       
     }
     else if(j==10 || j==11){ // ZEM
       sprintf(namhist1hg,"PedZEMhg_%d",j-10);
       sprintf(namhist2hg,"PedZEMhgOutOfTime_%d",j-10);
       sprintf(namhist3hg,"PedCorrZEMhg_%d",j-10);
       //
       sprintf(namhist1lg,"PedZEMlg_%d",j-10);
       sprintf(namhist2lg,"PedZEMlgOutOfTime_%d",j-10);
       sprintf(namhist3lg,"PedCorrZEMlg_%d",j-10);
     }
     else if(j>=12 && j<=16){ // ZNA
       sprintf(namhist1hg,"PedZNAhg_%d",j-12);
       sprintf(namhist2hg,"PedZNAhgOutOfTime_%d",j-12);
       sprintf(namhist3hg,"PedCorrZNAhg_%d",j-12);
       //
       sprintf(namhist1lg,"PedZNAlg_%d",j-12);
       sprintf(namhist2lg,"PedZNAlgOutOfTime_%d",j-12);
       sprintf(namhist3lg,"PedCorrZNAlg_%d",j-12);
     }
     else if(j>=17 && j<=21){ // ZPA
       sprintf(namhist1hg,"PedZPAhg_%d",j-17);
       sprintf(namhist2hg,"PedZPAhgOutOfTime_%d",j-17);
       sprintf(namhist3hg,"PedCorrZPAhg_%d",j-17);
       //
       sprintf(namhist1lg,"PedZPAlg_%d",j-17);
       sprintf(namhist2lg,"PedZPAlgOutOfTime_%d",j-17);
       sprintf(namhist3lg,"PedCorrZPAlg_%d",j-17);
     }
     else if(j==22 || j==24){ //Reference PMs
       sprintf(namhist1hg,"PedRefhg_%d",j-22);
       sprintf(namhist2hg,"PedRefhgOutOfTime_%d",j-22);
       sprintf(namhist3hg,"PedCorrRefhg_%d",j-22);
       //
       sprintf(namhist1lg,"PedReflg_%d",j-22);
       sprintf(namhist2lg,"PedReflgOutOfTime_%d",j-22);
       sprintf(namhist3lg,"PedCorrReflg_%d",j-22);
     }
     // --- High gain chain histos
     hPedhg[j] = new TH1F(namhist1hg, namhist1hg, 200, 0., 200.);
     hPedOutOfTimehg[j] = new TH1F(namhist2hg, namhist2hg, 200, 0., 200.);
     hPedCorrhg[j] = new TH2F(namhist3hg,namhist3hg,100,0.,200.,100,0.,200.);
     // --- Low gain chain histos
     hPedlg[j] = new TH1F(namhist1lg, namhist1lg, 100, 0., 1000.);
     hPedOutOfTimelg[j] = new TH1F(namhist2lg, namhist2lg, 100, 0., 1000.);
     hPedCorrlg[j] = new TH2F(namhist3lg,namhist3lg,100,0.,1000.,100,0.,1000.);
  }


  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","w");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
  }
  
  FILE *mapFile4Shuttle;

  /* report progress */
  daqDA_progressReport(10);


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  /* read the data files */
  int n;
  for(n=1;n<argc;n++){
   
    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }

    /* report progress */
    /* in this example, indexed on the number of files */
    daqDA_progressReport(10+80*n/argc);

    /* read the file */
    for(;;) {
      struct eventHeaderStruct *event;
      eventTypeType eventT;

      /* get next event */
      status=monitorGetEventDynamic((void **)&event);
      if(status==MON_ERR_EOF) break; /* end of monitoring file has been reached */
      if(status!=0) {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        return -1;
      }

      /* retry if got no event */
      if(event==NULL) {
        break;
      }
      
      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      reader->Select("ZDC");
      // --- Reading event header
      //UInt_t evtype = reader->GetType();
      //printf("\n\t ZDCPEDESTALda -> ev. type %d\n",evtype);
      //printf("\t ZDCPEDESTALda -> run # %d\n",reader->GetRunNumber());
      //
      AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
        

      /* use event - here, just write event id to result file */
      eventT=event->eventType;
      
      Int_t ich=0;
      Int_t adcMod[2*kNChannels], adcCh[2*kNChannels], sigCode[2*kNChannels];
      Int_t det[2*kNChannels], sec[2*kNChannels];
      
      if(eventT==START_OF_DATA){
	  	
	rawStreamZDC->SetSODReading(kTRUE);
	
	// --------------------------------------------------------
	// --- Writing ascii data file for the Shuttle preprocessor
        mapFile4Shuttle = fopen(MAPDATA_FILE,"w");
	if(!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
        else{
	  while((rawStreamZDC->Next()) && (ich<2*kNChannels)){
            if(rawStreamZDC->IsChMapping()){
	      adcMod[ich] = rawStreamZDC->GetADCModFromMap(ich);
	      adcCh[ich] = rawStreamZDC->GetADCChFromMap(ich);
	      sigCode[ich] = rawStreamZDC->GetADCSignFromMap(ich);
	      det[ich] = rawStreamZDC->GetDetectorFromMap(ich);
	      sec[ich] = rawStreamZDC->GetTowerFromMap(ich);
	      //
	      fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",
	        ich,adcMod[ich],adcCh[ich],sigCode[ich],det[ich],sec[ich]);
	      //
	      //printf("ZDCPEDESTALda.cxx -> %d mod %d ch %d, code %d det %d, sec %d\n",
	      //   ich,adcMod[ich],adcCh[ich],sigCode[ich],det[ich],sec[ich]);
	      //
	      ich++;
	    }
	  }
	}
        fclose(mapFile4Shuttle);
      }// SOD event
      
      if(eventT==PHYSICS_EVENT){
	// --- Reading data header
        reader->ReadHeader();
        const AliRawDataHeader* header = reader->GetDataHeader();
        if(header){
         UChar_t message = header->GetAttributes();
         if(message & 0x20){ // PEDESTAL RUN
            //printf("\t STANDALONE_PEDESTAL RUN raw data found\n");
         }
         else{
            printf("ZDCPEDESTALda.cxx -> NO STANDALONE_PEDESTAL RUN raw data found\n");
            return -1;
         }
        }
        else{
           printf("\t ATTENTION! No Raw Data Header found!!!\n");
           return -1;
        }
	
	rawStreamZDC->SetSODReading(kTRUE);

  	if(!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");	
	//
	// ----- Setting ch. mapping -----
	for(Int_t jk=0; jk<2*kNChannels; jk++){
	  //printf("ZDCPEDESTALDA.cxx ->  ch.%d mod %d, ch %d, code %d det %d, sec %d\n",
	  //    jk,adcMod[jk],adcCh[jk],sigCode[jk],det[jk],sec[jk]);
	  rawStreamZDC->SetMapADCMod(jk, adcMod[jk]);
	  rawStreamZDC->SetMapADCCh(jk, adcCh[jk]);
	  rawStreamZDC->SetMapADCSig(jk, sigCode[jk]);
	  rawStreamZDC->SetMapDet(jk, det[jk]);
	  rawStreamZDC->SetMapTow(jk, sec[jk]);
	}
	//
  	Int_t iraw=0;
  	Int_t RawADChg[kNChannels], RawADCoothg[kNChannels];
  	Int_t RawADClg[kNChannels], RawADCootlg[kNChannels];
  	for(Int_t j=0; j<kNChannels; j++){
  	   RawADChg[j]=0; RawADCoothg[j]=0;
  	   RawADClg[j]=0; RawADCootlg[j]=0;
  	}
	//
  	while(rawStreamZDC->Next()){
  	 Int_t index=-1;
	 Int_t detector = rawStreamZDC->GetSector(0);
	 Int_t sector = rawStreamZDC->GetSector(1);
	 
  	 if(rawStreamZDC->IsADCDataWord() && (detector!=-1)){
	  if(sector!=5){ // Physics signals
    	    if(detector==1) index = sector; // *** ZNC
	    else if(detector==2) index = sector+5; // *** ZPC
	    else if(detector==3) index = sector+9; // *** ZEM
	    else if(detector==4) index = sector+12; // *** ZNA
	    else if(detector==5) index = sector+17; // *** ZPA
	  }
	  else{ // Reference PMs
	    index = (detector-1)/3+22;
	  }
	  //
	  if(index==-1) printf("ERROR in ZDCPEDESTALda.cxx -> det %d quad %d res %d index %d \n", 
	    detector,sector,rawStreamZDC->GetADCGain(),index);
	  
	   //
	   if(iraw<2*kNChannels){ // --- In-time pedestals (1st 48 raw data)
	    if(rawStreamZDC->GetADCGain()==0){ 
	      hPedhg[index]->Fill(rawStreamZDC->GetADCValue()); 
	      RawADChg[index] = rawStreamZDC->GetADCValue();
	      //
	      //printf("\t filling histo hPedhg[%d]\n",index);
	    }
	    else{
	      hPedlg[index]->Fill(rawStreamZDC->GetADCValue()); 
	      RawADClg[index] = rawStreamZDC->GetADCValue();
	      //
	      //printf("\t filling histo hPedlg[%d]\n",index);
	    }
  	   }
  	   else{  // --- Out-of-time pedestals
  	    if(rawStreamZDC->GetADCGain()==0){
  	      hPedOutOfTimehg[index]->Fill(rawStreamZDC->GetADCValue());
  	      RawADCoothg[index] = rawStreamZDC->GetADCValue();
	      //
	      //printf("\t filling histo hPedOutOfTimehg[%d]\n",index);
  	    }
  	    else{
  	      hPedOutOfTimelg[index]->Fill(rawStreamZDC->GetADCValue());
  	      RawADCootlg[index] = rawStreamZDC->GetADCValue();
	      //
	      //printf("\t filling histo hPedOutOfTimelg[%d]\n",index);
  	    }
  	   }
  	    iraw++;
  	  }//IsADCDataWord()
  	  //
  	  if(iraw == 4*kNChannels){ // Last ADC channel -> Filling correlation histos
  	    for(Int_t k=0; k<kNChannels; k++){
  	      hPedCorrhg[k]->Fill(RawADCoothg[k], RawADChg[k]);
  	      hPedCorrlg[k]->Fill(RawADCootlg[k], RawADClg[k]);
  	    }
  	  }
         }
         //
         nevents_physics++;
         //
	 delete reader;
         delete rawStreamZDC;

      }//(if PHYSICS_EVENT) 
      nevents_total++;

      /* free resources */
      free(event);
    
    }
  }  
  
  /* Analysis of the histograms */
  //
  FILE *fileShuttle;
  fileShuttle = fopen(PEDDATA_FILE,"w");
  //
  Float_t MeanPed[2*kNChannels], MeanPedWidth[2*kNChannels], 
   	MeanPedOOT[2*kNChannels], MeanPedWidthOOT[2*kNChannels];
  // --- In-time pedestals
  TF1 *ADCfunchg[kNChannels];
  for(Int_t i=0; i<kNChannels; i++){
     hPedhg[i]->Fit("gaus","Q");
     ADCfunchg[i] = hPedhg[i]->GetFunction("gaus");
     MeanPed[i] = (Double_t) ADCfunchg[i]->GetParameter(1);
     MeanPedWidth[i] = (Double_t)  ADCfunchg[i]->GetParameter(2);
     fprintf(fileShuttle,"\t%f\t%f\n",MeanPed[i],MeanPedWidth[i]);
     //printf("\t MeanPedhg[%d] = %f\n",i, MeanPed[i]);
  }
  TF1 *ADCfunclg[kNChannels];
  for(Int_t i=0; i<kNChannels; i++){
     hPedlg[i]->Fit("gaus","Q");
     ADCfunclg[i] = hPedlg[i]->GetFunction("gaus");
     MeanPed[i+kNChannels] = (Double_t)  ADCfunclg[i]->GetParameter(1);
     MeanPedWidth[i+kNChannels] = (Double_t)  ADCfunclg[i]->GetParameter(2);
     fprintf(fileShuttle,"\t%f\t%f\n",MeanPed[i+kNChannels],MeanPedWidth[i+kNChannels]);
     //printf("\t MeanPedlg[%d] = %f\n",i+kNChannels, MeanPed[i+kNChannels]);
  }
  // --- Out-of-time pedestals
  TF1 *ADCootfunchg[kNChannels];
  for(Int_t i=0; i<kNChannels; i++){
     hPedOutOfTimehg[i]->Fit("gaus","Q");
     ADCootfunchg[i] = hPedOutOfTimehg[i]->GetFunction("gaus");
     MeanPedOOT[i] = (Double_t)  ADCootfunchg[i]->GetParameter(1);
     MeanPedWidthOOT[i] = (Double_t)  ADCootfunchg[i]->GetParameter(2);
     fprintf(fileShuttle,"\t%f\t%f\n",MeanPedOOT[i],MeanPedWidthOOT[i]);
     //printf("\t MeanPedOOThg[%d] = %f\n",i, MeanPedOOT[i]);
  }
  TF1 *ADCootfunclg[kNChannels];
  for(Int_t i=0; i<kNChannels; i++){
     hPedOutOfTimelg[i]->Fit("gaus","Q");
     ADCootfunclg[i] = hPedOutOfTimelg[i]->GetFunction("gaus");
     MeanPedOOT[i+kNChannels] = (Double_t)  ADCootfunclg[i]->GetParameter(1);
     MeanPedWidthOOT[i+kNChannels] = (Double_t)  ADCootfunclg[i]->GetParameter(2);
     fprintf(fileShuttle,"\t%f\t%f\n",MeanPedOOT[i+kNChannels],MeanPedWidthOOT[i+kNChannels]);
     //printf("\t MeanPedOOTlg[%d] = %f\n",i+kNChannels, MeanPedOOT[i+kNChannels]);
  }
  
  // --- Correlations

  Float_t CorrCoeff0[2*kNChannels], CorrCoeff1[2*kNChannels];
  TProfile *hPedCorrProfhg[kNChannels], *hPedCorrProflg[kNChannels];
  TF1 *ffunchg[kNChannels], *ffunclg[kNChannels];
  char namhist4[50];
  for(int i=0;i<kNChannels;i++) {
     sprintf(namhist4,"ADCHRvsOOT%d_Prof",i);
     hPedCorrProfhg[i] = hPedCorrhg[i]->ProfileX(namhist4,-1,-1,"S");
     hPedCorrProfhg[i]->SetName(namhist4);
     hPedCorrProfhg[i]->Fit("pol1","Q");
     ffunchg[i] = hPedCorrProfhg[i]->GetFunction("pol1");
     CorrCoeff0[i] = (Double_t)  ffunchg[i]->GetParameter(0);
     CorrCoeff1[i] = (Double_t) ffunchg[i]->GetParameter(1);
     fprintf(fileShuttle,"\t%f\t%f\n",CorrCoeff0[i],CorrCoeff1[i]);
     //printf("\t CorrCoeff0[%d] = %f, CorrCoeff1[%d] = %f\n",i, CorrCoeff0[i], i, CorrCoeff1[i]);
  }    
  for(int i=0;i<kNChannels;i++) {
     sprintf(namhist4,"ADCLRvsOOT%d_Prof",i);
     hPedCorrProflg[i] = hPedCorrlg[i]->ProfileX(namhist4,-1,-1,"S");
     hPedCorrProflg[i]->SetName(namhist4);
     hPedCorrProflg[i]->Fit("pol1","Q");
     ffunclg[i] = hPedCorrProflg[i]->GetFunction("pol1");
     CorrCoeff0[i+kNChannels] =  (Double_t) ffunclg[i]->GetParameter(0);
     CorrCoeff1[i+kNChannels] =  (Double_t) ffunclg[i]->GetParameter(1);
     fprintf(fileShuttle,"\t%f\t%f\n",CorrCoeff0[i+kNChannels],CorrCoeff1[i+kNChannels]);
     //printf("\t CorrCoeff0[%d] = %f, CorrCoeff1[%d] = %f\n",
     //		i+kNChannels, CorrCoeff0[i+kNChannels], i+kNChannels, CorrCoeff1[i+kNChannels]);
  }    

  //						       
  fclose(fileShuttle);
  //
  for(Int_t j=0; j<kNChannels; j++){
     delete hPedhg[j];
     delete hPedOutOfTimehg[j];
     delete hPedCorrhg[j];
     delete hPedlg[j];
     delete hPedOutOfTimelg[j];
     delete hPedCorrlg[j];
  }

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);
  
  /* report progress */
  daqDA_progressReport(90);

  /* store the result files on FES */
  // [1] File with mapping
  status = daqDA_FES_storeFile(MAPDATA_FILE,MAPDATA_FILE);
  if(status){
    printf("Failed to export file %d to DAQ FES\n",status);
    return -1;
  }
  // [2] File with pedestal data
  status = daqDA_FES_storeFile(PEDDATA_FILE,PEDDATA_FILE);
  if(status){
    printf("Failed to export file %d to DAQ FES\n",status);
    return -1;
  }
  
  /* store the result files on DB */
  status = daqDA_DB_storeFile(PEDDATA_FILE,PEDDATA_FILE);  
  if(status){
    printf("Failed to export file %d to DAQ DB\n",status);
    return -1;
  }

  /* report progress */
  daqDA_progressReport(100);

  return status;
}
