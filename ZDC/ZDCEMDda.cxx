/*

DAcase2.c

This program connects to the DAQ data source passed as argument
and populates local "./result.txt" file with the ids of events received
during the run.

The program exits when being asked to shut down (daqDA_checkshutdown)
or End of Run event.

Messages on stdout are exported to DAQ log system.

Contact: Chiara.Oppedisano@to.infn.it
Link: /afs/cern.ch/user/c/chiarao/public/RawEMD.date
Run Type: STANDALONE_EMD_RUN
DA Type: MON
Number of events needed: at least ~10^3
Input Files: ZDCPedestal.dat
Output Files: ZDCEMDCalib.dat
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
  
  //
  // --- Preparing histos for EM dissociation spectra
  //
  TH1F* histoEMDRaw[4];
  TH1F* histoEMDCorr[4];

  char namhistr[50], namhistc[50];
  for(Int_t i=0; i<4; i++) {
     if(i==0){
       sprintf(namhistr,"ZN%d-EMDRaw",i+1);
       sprintf(namhistc,"ZN%d-EMDCorr",i+1);
     }
     else if(i==1){
       sprintf(namhistr,"ZP%d-EMDRaw",i);
       sprintf(namhistc,"ZP%d-EMDCorr",i);
     }
     else if(i==2){
       sprintf(namhistr,"ZN%d-EMDRaw",i);
       sprintf(namhistc,"ZN%d-EMDCorr",i);
     }
     else if(i==3){
       sprintf(namhistr,"ZP%d-EMDRaw",i-1);
       sprintf(namhistc,"ZP%d-EMDCorr",i-1);
     }
     histoEMDRaw[i] = new TH1F(namhistr,namhistr,100,0.,4000.);
     histoEMDCorr[i] = new TH1F(namhistc,namhistc,100,0.,4000.);
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
      //
      // *** To analyze EMD events you MUST have a pedestal data file!!!
      // *** -> check if a pedestal run has been analyzied
      FILE *filePed=NULL;
      filePed=fopen("./ZDCPedestal.dat","r");
      if (filePed==NULL) {
        printf("\t ERROR!!! You MUST have a ZDCPedestal.dat file!!!\n");
        return -1;
      }
      // 144 = 48 in-time + 48 out-of-time + 48 correlations
      Float_t readValues[2][144], MeanPed[44], MeanPedWidth[44], 
   	      MeanPedOOT[44], MeanPedWidthOOT[44],
   	      CorrCoeff0[44], CorrCoeff1[44];
      //
      for(int jj=0; jj<144; jj++){
        for(int ii=0; ii<2; ii++){
           fscanf(filePed,"%f",&readValues[ii][jj]);
        }
	if(jj<48){
	  MeanPed[jj] = readValues[0][jj];
	  MeanPedWidth[jj] = readValues[1][jj];
	}
	else if(jj>48 && jj<96){
	  MeanPedOOT[jj-48] = readValues[0][jj];
	  MeanPedWidthOOT[jj-48] = readValues[1][jj];
	}
	else if(jj>144){
	  CorrCoeff0[jj-96] = readValues[0][jj]; 
	  CorrCoeff1[jj-96] = readValues[1][jj];;
	}
      }
      //
      //
      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      const AliRawDataHeader* header = reader->GetDataHeader();
      if(header){
         UChar_t message = header->GetAttributes();
	 if(message & 0x70){ // DEDICATED EMD RUN
	    printf("\t STANDALONE_EMD_RUN raw data found\n");
	    continue;
	 }
	 else{
	    printf("\t NO STANDALONE_EMD_RUN raw data found\n");
	    return -1;
	 }
      }
      //Commented until we won't have true Raw Data Header...
      //else{
      //   printf("\t ATTENTION! No Raw Data Header found!!!\n");
      //   return -1;
      //}
      //
      AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
      //
      if (!rawStreamZDC->Next()) printf(" \t No raw data found!! ");
      //
      Float_t ZDCRawADC[4], ZDCCorrADC[4], ZDCCorrADCSum[4];
      for(Int_t g=0; g<4; g++){
           ZDCCorrADCSum[g] = 0.;
	   ZDCRawADC[g] = 0.;
      }
      //
      while(rawStreamZDC->Next()){
        if(rawStreamZDC->IsADCDataWord()){
	  Int_t DetIndex=999, PedIndex=999;
	  if(rawStreamZDC->GetSector(0) == 1 || rawStreamZDC->GetSector(0) == 2){
	    DetIndex = rawStreamZDC->GetSector(0)-1;
	    PedIndex = (rawStreamZDC->GetSector(0)+1)+4*rawStreamZDC->GetSector(1);
	  }
	  else if(rawStreamZDC->GetSector(0) == 4 || rawStreamZDC->GetSector(0) == 5){
	    DetIndex = rawStreamZDC->GetSector(0)-2;
	    PedIndex = (rawStreamZDC->GetSector(0)-2)+4*rawStreamZDC->GetSector(1)+24;
	  }
          //
	  if(rawStreamZDC->GetADCGain() == 1){ //EMD -> LR ADC
	    //
	    ZDCRawADC[DetIndex] += (Float_t) rawStreamZDC->GetADCValue();
	    // Mean pedestal subtraction 
	    Float_t Pedestal = MeanPed[PedIndex];
	    // Pedestal subtraction from correlation with out-of-time signals
	    //Float_t Pedestal = CorrCoeff0[PedIndex]+CorrCoeff1[PedIndex]*MeanPedOOT[PedIndex];
	    //
	    ZDCCorrADC[DetIndex] = (rawStreamZDC->GetADCValue()) - Pedestal;
	    ZDCCorrADCSum[DetIndex] += ZDCCorrADC[DetIndex];
	    //
	    /*printf("\t det %d quad %d res %d pedInd %d ADCCorr %d ZDCCorrADCSum[%d] = %d\n", 
	       rawStreamZDC->GetSector(0),rawStreamZDC->GetSector(1),
	       rawStreamZDC->GetADCGain(),PedIndex,  
	       (Int_t) (rawStreamZDC->GetADCValue() - Pedestal), DetIndex, 
	       (Int_t) ZDCCorrADCSum[DetIndex]);
	    */
	  }   
	}//IsADCDataWord()
	 //
       }
       //
       nevents_physics++;
       //
       for(Int_t j=0; j<4; j++){
          histoEMDRaw[j]->Fill(ZDCRawADC[j]);
          histoEMDCorr[j]->Fill(ZDCCorrADCSum[j]);
       }
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
  fileShuttle = fopen("ZDCEMDCalib.dat","w");
  //
  Int_t BinMax[4];
  Float_t YMax[4];
  Int_t NBinsx[4];
  Float_t MeanFitVal[4];
  TF1 *fitfun[4];
  for(Int_t k=0; k<4; k++){
     BinMax[k] = histoEMDCorr[k]->GetMaximumBin();
     YMax[k] = (histoEMDCorr[k]->GetXaxis())->GetXmax();
     NBinsx[k] = (histoEMDCorr[k]->GetXaxis())->GetNbins();
//     printf("\n\t Det%d -> BinMax = %d, ChXMax = %f\n", k+1, BinMax[k], BinMax[k]*YMax[k]/NBinsx[k]);
     histoEMDCorr[k]->Fit("gaus","Q","",BinMax[k]*YMax[k]/NBinsx[k]*0.7,BinMax[k]*YMax[k]/NBinsx[k]*1.25);
     fitfun[k] = histoEMDCorr[k]->GetFunction("gaus");
     MeanFitVal[k] = (Float_t) (fitfun[k]->GetParameter(1));
     printf("\n\t Mean Value from gaussian fit = %f\n", MeanFitVal[k]);
  }
  //
   Float_t CalibCoeff[6];
  /*for(Int_t j=0; j<6; j++){
     if(j<4) CalibCoeff[j] = 2.76/MeanFitVal[j];
     else  CalibCoeff[j] = 1.;
     fprintf(fileShuttle,"\t%f\n",CalibCoeff[j]);
  }
  */
  // --- For the moment we have sim data only for ZN1!!!
  for(Int_t j=0; j<6; j++){
     if(j==0) CalibCoeff[j] = 2.76/MeanFitVal[j];
     else if(j>0 && j<4) CalibCoeff[j] = CalibCoeff[0];
     else  CalibCoeff[j] = 1.;
     fprintf(fileShuttle,"\t%f\n",CalibCoeff[j]);
  }
  //						       
  fclose(fileShuttle);
  

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);


  return status;
}
