/*

This program reads the DAQ data files passed as argument using the monitoring library.

It computes the average event size and populates local "./result.txt" file with the 
result.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalone CALIBRATION_EMD runs

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
#define ENCALIBDATA_FILE   "ZDCEnergyCalib.dat"
#define TOWCALIBDATA_FILE  "ZDCTowerCalib.dat"

#include <stdio.h>
#include <Riostream.h>
#include <Riostream.h>

// DATE
#include <daqDA.h>
#include <event.h>
#include <monitor.h>

//ROOT
#include <TROOT.h>
#include <TPluginManager.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitter.h>
#include "TMinuitMinimizer.h"

//AliRoot
#include <AliRawReaderDate.h>
#include <AliRawEventHeaderBase.h>
#include <AliZDCRawStream.h>


/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {
  

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()"); 

  TMinuitMinimizer m; 
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", "Minuit","TMinuitMinimizer",
      "Minuit", "TMinuitMinimizer(const char *)");
  TVirtualFitter::SetDefaultFitter("Minuit");

  int status = 0;
  // No. of ZDC cabled ch.
  int const kNChannels = 24;
  int const kNScChannels = 32;

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

  // *** To analyze EMD events you MUST have a pedestal data file!!!
  // *** -> check if a pedestal run has been analyzed
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
  Float_t readValues[2][3*2*kNChannels];
  Float_t MeanPed[2*kNChannels];
  Float_t CorrCoeff0[2*kNChannels], CorrCoeff1[2*kNChannels];
  // ***************************************************
  //   Unless we have a narrow correlation to fit we
  //	don't fit and store in-time vs. out-of-time
  //	histograms -> mean pedstal subtracted!!!!!!
  // ***************************************************
  //
  for(int jj=0; jj<6*kNChannels; jj++){
    for(int ii=0; ii<2; ii++){
       fscanf(filePed,"%f",&readValues[ii][jj]);
    }
    if(jj<kNChannels && jj<2*kNChannels){
      MeanPed[jj] = readValues[0][jj];
      //printf("\t MeanPedhg[%d] = %1.1f\n",jj, MeanPedhg[jj]);
    }
    else if(jj>2*kNChannels && jj>4*kNChannels){
      CorrCoeff0[jj-4*kNChannels] = readValues[0][jj]; 
      CorrCoeff1[jj-4*kNChannels] = readValues[1][jj];;
    }
  }

  /* report progress */
  daqDA_progressReport(10);


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  struct eventHeaderStruct *event;
  eventTypeType eventT;

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
      
      Int_t ich=0;
      Int_t adcMod[2*kNChannels], adcCh[2*kNChannels], sigCode[2*kNChannels];
      Int_t det[2*kNChannels], sec[2*kNChannels];
      for(Int_t y=0; y<2*kNChannels; y++){
        adcMod[y]=adcCh[y]=sigCode[y]=det[y]=sec[y]=0;
      }
      
      Int_t iScCh=0;
      Int_t scMod[kNScChannels], scCh[kNScChannels], scSigCode[kNScChannels];
      Int_t scDet[kNScChannels], scSec[kNScChannels];
      for(Int_t y=0; y<kNScChannels; y++){
        scMod[y]=scCh[y]=scSigCode[y]=scDet[y]=scSec[y]=0;
      }
      //
      Int_t modNum=-1, modType=-1;
      
      if(eventT==START_OF_DATA){
	  	
	rawStreamZDC->SetSODReading(kTRUE);
	
	// --------------------------------------------------------
	// --- Writing ascii data file for the Shuttle preprocessor
        mapFile4Shuttle = fopen(MAPDATA_FILE,"w");
	if(!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
        else{
	  while((rawStreamZDC->Next())){
            if(rawStreamZDC->IsHeaderMapping()){ // mapping header
	       modNum = rawStreamZDC->GetADCModule();
	       modType = rawStreamZDC->GetModType();
	    }
            if(rawStreamZDC->IsChMapping()){ 
	      if(modType==1){ // ADC mapping ----------------------
	        adcMod[ich]  = rawStreamZDC->GetADCModFromMap(ich);
	        adcCh[ich]   = rawStreamZDC->GetADCChFromMap(ich);
	        sigCode[ich] = rawStreamZDC->GetADCSignFromMap(ich);
	        det[ich]     = rawStreamZDC->GetDetectorFromMap(ich);
	        sec[ich]     = rawStreamZDC->GetTowerFromMap(ich);
	        //
	        fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",
	          ich,adcMod[ich],adcCh[ich],sigCode[ich],det[ich],sec[ich]);
	        //
	        //printf("  Mapping in DA -> %d ADC: mod %d ch %d, code %d det %d, sec %d\n",
	        //  ich,adcMod[ich],adcCh[ich],sigCode[ich],det[ich],sec[ich]);
	        //
	        ich++;
	      }
	      else if(modType==2){ //VME scaler mapping --------------------
	        scMod[iScCh]     = rawStreamZDC->GetScalerModFromMap(iScCh);
	        scCh[iScCh]      = rawStreamZDC->GetScalerChFromMap(iScCh);
	        scSigCode[iScCh] = rawStreamZDC->GetScalerSignFromMap(iScCh);
	        scDet[iScCh]     = rawStreamZDC->GetScDetectorFromMap(iScCh);
	        scSec[iScCh]    = rawStreamZDC->GetScTowerFromMap(iScCh);
	        //
	        fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",
	          iScCh,scMod[iScCh],scCh[iScCh],scSigCode[iScCh],scDet[iScCh],scSec[iScCh]);
	        //
	        //printf("  Mapping in DA -> %d Scaler: mod %d ch %d, code %d det %d, sec %d\n",
	        //  iScCh,scMod[iScCh],scCh[iScCh],scSigCode[iScCh],scDet[iScCh],scSec[iScCh]);
	        //
	        iScCh++;
	      }
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
      
      rawStreamZDC->SetSODReading(kTRUE);

      if (!rawStreamZDC->Next()) printf(" \t No raw data found!! ");
      //
      // ----- Setting ch. mapping -----
      for(Int_t jk=0; jk<2*kNChannels; jk++){
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
	Int_t det = rawStreamZDC->GetSector(0);
	Int_t quad = rawStreamZDC->GetSector(1);
        
	if(rawStreamZDC->IsADCDataWord() && !(rawStreamZDC->IsUnderflow())
	     && !(rawStreamZDC->IsOverflow()) && det!=-1
	     && (rawStreamZDC->GetADCGain() == 1)){ // Selecting LOW RES ch.s

	  //printf("  IsADCWord %d, IsUnderflow %d, IsOverflow %d\n",
	  //  rawStreamZDC->IsADCDataWord(),rawStreamZDC->IsUnderflow(),rawStreamZDC->IsOverflow());
	  
	  // Taking LOW RES channels -> channel+kNChannels !!!!
	  Int_t DetIndex=999, PedIndex=999;
	  if(det != 3 && quad != 5){ // Not ZEM nor PMRef
	    if(det == 1){
	      DetIndex = det-1;
	      PedIndex = quad+kNChannels;
	    }
	    else if(det==2){
	      DetIndex = det-1;
	      PedIndex = quad+5+kNChannels;
	    }
	    else if(det == 4){
	      DetIndex = det-2;
	      PedIndex = quad+12+kNChannels;
	    }
	    else if(det == 5){
	      DetIndex = det-2;
	      PedIndex = quad+17+kNChannels;
	    }
            //EMD -> LR ADCs
	    if(rawStreamZDC->GetADCGain() == 1 && (DetIndex!=999 || PedIndex!=999)){ 
	      //
	      ZDCRawADC[DetIndex] += (Float_t) rawStreamZDC->GetADCValue();
	      //
	      // Mean pedestal subtraction 
	      Float_t Pedestal = MeanPed[PedIndex];
	      // Pedestal subtraction from correlation with out-of-time signals
	      //Float_t Pedestal = CorrCoeff0[PedIndex]+CorrCoeff1[PedIndex]*MeanPedOOT[PedIndex];
	      //
	      ZDCCorrADC[DetIndex] = (rawStreamZDC->GetADCValue()) - Pedestal;
	      ZDCCorrADCSum[DetIndex] += ZDCCorrADC[DetIndex];
	      //
	      /*printf("\t det %d quad %d res %d pedInd %d detInd %d:"
	         "ADCCorr = %d, ZDCCorrADCSum = %d\n", 
	         det,quad,rawStreamZDC->GetADCGain(),PedIndex,DetIndex, 
	         (Int_t) ZDCCorrADC[DetIndex],(Int_t) ZDCCorrADCSum[DetIndex]);
	      */
	    }
	    if(DetIndex==999 || PedIndex==999) 
	    	printf(" WARNING! Detector a/o pedestal index are WRONG!!!\n");
 
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
    }//(if PHYSICS_EVENT)
    
    nevents_total++;

   }

   /* free resources */
   free(event);
  }
    
  /* Analysis of the histograms */
  //
  FILE *fileShuttle1 = fopen(ENCALIBDATA_FILE,"w");
  //
  Int_t BinMax[4];
  Float_t YMax[4];
  Int_t NBinsx[4];
  Float_t MeanFitVal[4];
  TF1 *fitfun[4];
  for(Int_t k=0; k<4; k++){
     if(histoEMDCorr[k]->GetEntries() == 0){
       printf("\n WARNING! Empty histos -> ending DA WITHOUT writing output\n\n");
       return -1;
     } 
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
  for(Int_t j=0; j<6; j++){
     if(j<4){
       CalibCoeff[j] = MeanFitVal[j];
       fprintf(fileShuttle1,"\t%f\n",CalibCoeff[j]);
     }
     // ZEM energy calib. coeff. = 1
     else if(j==4 || j==5){
       CalibCoeff[j] = 1.; 
       fprintf(fileShuttle1,"\t%f\n",CalibCoeff[j]);
     }
  }
  fclose(fileShuttle1);
  //
  FILE *fileShuttle2 = fopen(TOWCALIBDATA_FILE,"w");
  for(Int_t j=0; j<4; j++){
     // Note -> For the moment the inter-calibration coeff. are set to 1 
     for(Int_t k=0; k<5; k++){  
       icoeff[k] = 1.;
       fprintf(fileShuttle2,"\t%f",icoeff[k]);
       if(k==4) fprintf(fileShuttle2,"\n");
     }
  }
  fclose(fileShuttle2);
  
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
  status = daqDA_FES_storeFile(MAPDATA_FILE, "MAPPING");
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }
  //
  status = daqDA_FES_storeFile(ENCALIBDATA_FILE, "EMDENERGYCALIB");
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }
  //
  status = daqDA_FES_storeFile(TOWCALIBDATA_FILE, "EMDTOWERCALIB");
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }

  /* report progress */
  daqDA_progressReport(100);

  return status;
}
