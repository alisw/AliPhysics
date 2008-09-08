/*

This program reads the DAQ data files passed as argument using the monitoring library.

It computes the average event size and populates local "./result.txt" file with the 
result.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalone pedestal runs

Contact: Chiara.Oppedisano@to.infn.it
Link: 
Run Type: STANDALONE_EMD_RUN
DA Type: LDC
Number of events needed: at least ~5*10^3
Input Files: ZDCPedestal.dat
Output Files: ZDCEMDCalib.dat, ZDCChMapping.dat
Trigger Types Used: Standalone Trigger

*/
#define PEDDATA_FILE  "ZDCPedestal.dat"
#define MAPDATA_FILE  "ZDCChMapping.dat"
#define EMDDATA_FILE  "ZDCEMDCalib.dat"

#include <stdio.h>
#include <Riostream.h>
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
#include <TFitter.h>

//AliRoot
#include <AliRawReaderDate.h>
#include <AliRawEventHeaderBase.h>
#include <AliZDCRawStream.h>


/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {
  
  TFitter *minuitFit = new TFitter(4);
  TVirtualFitter::SetFitter(minuitFit);

  int status = 0;

  /* log start of process */
  printf("ZDC EMD program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
  
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

  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
  }
  
  FILE *mapFile4Shuttle;

  // *** To analyze LASER events you MUST have a pedestal data file!!!
  // *** -> check if a pedestal run has been analyzied
  int read = 0;
  read = daqDA_DB_getFile(PEDDATA_FILE,PEDDATA_FILE);
  if(read){
    printf("\t ERROR!!! ZDCPedestal.dat file NOT FOUND in DAQ db!!!\n");
    return -1;
  }
  else printf("\t ZDCPedestal.dat file retrieved from DAQ db\n");
  
  FILE *filePed = fopen(PEDDATA_FILE,"r");
  if (filePed==NULL) {
    printf("\t ERROR!!! Can't open ZDCPedestal.dat file!!!\n");
    return -1;
  }

  // 144 = 48 in-time + 48 out-of-time + 48 correlations
  Float_t readValues[2][144], MeanPed[44], MeanPedWidth[44], 
  	MeanPedOOT[44], MeanPedWidthOOT[44];
  // ***************************************************
  //   Unless we have a narrow correlation to fit we
  //	don't fit and store in-time vs. out-of-time
  //	histograms -> mean pedstal subtracted!!!!!!
  // ***************************************************
  //Float_t CorrCoeff0[44], CorrCoeff1[44];
  //
  for(int jj=0; jj<144; jj++){
    for(int ii=0; ii<2; ii++){
       fscanf(filePed,"%f",&readValues[ii][jj]);
    }
    if(jj<48){
      MeanPed[jj] = readValues[0][jj];
      MeanPedWidth[jj] = readValues[1][jj];
      //printf("\t MeanPed[%d] = %1.1f\n",jj, MeanPed[jj]);
    }
    else if(jj>48 && jj<96){
      MeanPedOOT[jj-48] = readValues[0][jj];
      MeanPedWidthOOT[jj-48] = readValues[1][jj];
    }
    /*else if(jj>144){
      CorrCoeff0[jj-96] = readValues[0][jj]; 
      CorrCoeff1[jj-96] = readValues[1][jj];;
    }
    */
  }

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
    for(;;){
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
      //printf("\n\t ZDCEMDda -> ev. type %d\n",evtype);
      //printf("\t ZDCEMDda -> run # %d\n",reader->GetRunNumber());
      //
      AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
        

      /* use event - here, just write event id to result file */
      eventT=event->eventType;
      
      Int_t ich=0, adcMod[48], adcCh[48], sigCode[48], det[48], sec[48];
      if(eventT==START_OF_DATA){
	  	
	if(!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
        else{
	  while(rawStreamZDC->Next()){
            if(rawStreamZDC->IsChMapping()){
	      adcMod[ich] = rawStreamZDC->GetADCModFromMap(ich);
	      adcCh[ich] = rawStreamZDC->GetADCChFromMap(ich);
	      sigCode[ich] = rawStreamZDC->GetADCSignFromMap(ich);
	      det[ich] = rawStreamZDC->GetDetectorFromMap(ich);
	      sec[ich] = rawStreamZDC->GetTowerFromMap(ich);
	      ich++;
	    }
	  }
	}
	// --------------------------------------------------------
	// --- Writing ascii data file for the Shuttle preprocessor
        mapFile4Shuttle = fopen(MAPDATA_FILE,"w");
        for(Int_t i=0; i<ich; i++){
	   fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",i,
	     adcMod[i],adcCh[i],sigCode[i],det[i],sec[i]);
	   //
	   //printf("ZDCPEDESTALDA.cxx ->  ch.%d mod %d, ch %d, code %d det %d, sec %d\n",
	   //	   i,adcMod[i],adcCh[i],sigCode[i],det[i],sec[i]);
        }
        fclose(mapFile4Shuttle);
      }
    
    if(eventT==PHYSICS_EVENT){
      // --- Reading data header
      reader->ReadHeader();
      const AliRawDataHeader* header = reader->GetDataHeader();
      if(header){
         UChar_t message = header->GetAttributes();
	 if(message & 0x70){ // DEDICATED EMD RUN
	    //printf("\t STANDALONE_EMD_RUN raw data found\n");
	    continue;
	 }
	 else{
	    printf("\t NO STANDALONE_EMD_RUN raw data found\n");
	    return -1;
	 }
      }
      else{
         printf("\t ATTENTION! No Raw Data Header found!!!\n");
         return -1;
      }

      if (!rawStreamZDC->Next()) printf(" \t No raw data found!! ");
      //
      // ----- Setting ch. mapping -----
      for(Int_t jk=0; jk<48; jk++){
        rawStreamZDC->SetMapADCMod(jk, adcMod[jk]);
        rawStreamZDC->SetMapADCCh(jk, adcCh[jk]);
        rawStreamZDC->SetMapADCSig(jk, sigCode[jk]);
        rawStreamZDC->SetMapDet(jk, det[jk]);
        rawStreamZDC->SetMapTow(jk, sec[jk]);
      }
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
  }
    
  /* Analysis of the histograms */
  //
  FILE *fileShuttle;
  fileShuttle = fopen(EMDDATA_FILE,"w");
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
     //printf("\n\t Mean Value from gaussian fit = %f\n", MeanFitVal[k]);
  }
  //
  Float_t CalibCoeff[6];     
  Float_t icoeff[5];
  //
  for(Int_t j=0; j<10; j++){
     if(j<4){
       CalibCoeff[j] = MeanFitVal[j];
       fprintf(fileShuttle,"\t%f\n",CalibCoeff[j]);
     }
     // ZEM calib. coeff. = 1
     else if(j==4 || j==5){
       CalibCoeff[j] = 1.; 
       fprintf(fileShuttle,"\t%f\n",CalibCoeff[j]);
     }
     // Note -> For the moment the inter-calibration
     //	     coefficients are set to 1 
     else if(j>5){
       for(Int_t k=0; k<5; k++){  
         icoeff[k] = 1.;
         fprintf(fileShuttle,"\t%f",icoeff[k]);
         if(k==4) fprintf(fileShuttle,"\n");
       }
     }
  }
  //						       
  fclose(fileShuttle);
  
  for(Int_t ij=0; ij<4; ij++){
    delete histoEMDRaw[ij];
    delete histoEMDCorr[ij];
  }
  
  //delete minuitFit;
  TVirtualFitter::SetFitter(0);

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);
  
  /* report progress */
  daqDA_progressReport(90);

  /* store the result file on FES */
  status = daqDA_FES_storeFile(MAPDATA_FILE, MAPDATA_FILE);
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }
  //
  status = daqDA_FES_storeFile(EMDDATA_FILE, EMDDATA_FILE);
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }

  /* report progress */
  daqDA_progressReport(100);

  return status;
}
