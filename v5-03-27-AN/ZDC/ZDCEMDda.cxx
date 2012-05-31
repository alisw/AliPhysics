/*

This program reads the DAQ data files passed as argument using the monitoring library.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalone CALIBRATION_EMD runs

Contact: Chiara.Oppedisano@to.infn.it
Link: 
Run Type: CALIBRATION_EMD_RUN
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
  int const kNModules = 9;
  int const kNChannels = 24;
  int const kNScChannels = 32;
  Int_t kFirstADCGeo=0, kLastADCGeo=1; // NO out-of-time signals!!!
            
  Int_t iMod=-1;
  Int_t modGeo[kNModules], modType[kNModules],modNCh[kNModules];
  for(Int_t kl=0; kl<kNModules; kl++){
     modGeo[kl]=modType[kl]=modNCh[kl]=0;
  }
  
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
      
  Int_t itdcCh=0;
  Int_t tdcMod[kNScChannels], tdcCh[kNScChannels], tdcSigCode[kNScChannels];
  Int_t tdcDet[kNScChannels], tdcSec[kNScChannels];
  for(Int_t y=0; y<kNScChannels; y++){
    tdcMod[y]=tdcCh[y]=tdcSigCode[y]=tdcDet[y]=tdcSec[y]=-1;
  }

  /* log start of process */
  printf("ZDC EMD program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
  
  // --- Preparing histos for EM dissociation spectra
  //
  TH1F::AddDirectory(0);
  TH1F* histoEMDCorr[4];
  //
  //char namhistr[50];
  char namhistc[50];
  for(Int_t i=0; i<4; i++) {
     if(i==0) sprintf(namhistc,"ZN%d-EMDCorr",i+1);
     else if(i==1) sprintf(namhistc,"ZP%d-EMDCorr",i);
     else if(i==2) sprintf(namhistc,"ZN%d-EMDCorr",i);
     else if(i==3) sprintf(namhistc,"ZP%d-EMDCorr",i-1);
     histoEMDCorr[i] = new TH1F(namhistc,namhistc,200,0.,2000.);
  }

   // --- Preparing histos for tower inter-calibration
  //
/*  TH1F* histZNCtow[4]; TH1F* histZPCtow[4];
  TH1F* histZNAtow[4]; TH1F* histZPAtow[4];
  //
  char namhistznc[50], namhistzpc[50];
  char namhistzna[50], namhistzpa[50];
  for(Int_t i=0; i<4; i++) {
     sprintf(namhistznc,"ZNC-tow%d",i+1);
     sprintf(namhistzpc,"ZPC-tow%d",i+1);
     sprintf(namhistzna,"ZNA-tow%d",i+1);
     sprintf(namhistzpa,"ZPA-tow%d",i+1);
     //
     histZNCtow[i] = new TH1F(namhistznc,namhistznc,100,0.,4000.);
     histZPCtow[i] = new TH1F(namhistzpc,namhistzpc,100,0.,4000.);
     histZNAtow[i] = new TH1F(namhistzna,namhistzna,100,0.,4000.);
     histZPAtow[i] = new TH1F(namhistzpa,namhistzpa,100,0.,4000.);
  }
*/  
  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
  }
  /* report progress */
  daqDA_progressReport(10);

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
  Float_t readValues[2][6*kNChannels];
  Float_t MeanPedhg[kNChannels], MeanPedlg[kNChannels];
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
    if(jj<kNChannels){
      MeanPedhg[jj] = readValues[0][jj];
      //printf("\t MeanPedhg[%d] = %1.1f\n",jj, MeanPedhg[jj]);
    }
    else if(jj>=kNChannels && jj<2*kNChannels){
      MeanPedlg[jj-kNChannels] = readValues[0][jj];
      //printf("\t MeanPedlg[%d] = %1.1f\n",jj-kNChannels, MeanPedlg[jj-kNChannels]);
    }
    else if(jj>4*kNChannels){
      CorrCoeff0[jj-4*kNChannels] = readValues[0][jj]; 
      CorrCoeff1[jj-4*kNChannels] = readValues[1][jj];;
    }
  }
  
  FILE *mapFile4Shuttle;

  /* report progress */
  daqDA_progressReport(20);


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
    daqDA_progressReport(20+70*n/argc);

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
      
      if(eventT==START_OF_DATA){

	iMod=-1; ich=0; iScCh=0; itdcCh=0;
	  			  		
	rawStreamZDC->SetSODReading(kTRUE);
	
	// --------------------------------------------------------
	// --- Writing ascii data file for the Shuttle preprocessor
        mapFile4Shuttle = fopen(MAPDATA_FILE,"w");
	if(!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
        else{
	  while((rawStreamZDC->Next())){
            if(rawStreamZDC->IsHeaderMapping()){ // mapping header
	       iMod++;
	       modGeo[iMod]  = rawStreamZDC->GetADCModule();
	       modType[iMod] = rawStreamZDC->GetModType();
	       modNCh[iMod]  = rawStreamZDC->GetADCNChannels();
	    }
            if(rawStreamZDC->IsChMapping()){ 
	      if(modType[iMod]==1){ // ADC mapping ----------------------
	        adcMod[ich]  = rawStreamZDC->GetADCModFromMap(ich);
	        adcCh[ich]   = rawStreamZDC->GetADCChFromMap(ich);
	        sigCode[ich] = rawStreamZDC->GetADCSignFromMap(ich);
	        det[ich]     = rawStreamZDC->GetDetectorFromMap(ich);
	        sec[ich]     = rawStreamZDC->GetTowerFromMap(ich);
	        ich++;
	      }
	      else if(modType[iMod]==2){ //VME scaler mapping --------------------
	        scMod[iScCh]     = rawStreamZDC->GetScalerModFromMap(iScCh);
	        scCh[iScCh]      = rawStreamZDC->GetScalerChFromMap(iScCh);
	        scSigCode[iScCh] = rawStreamZDC->GetScalerSignFromMap(iScCh);
	        scDet[iScCh]     = rawStreamZDC->GetScDetectorFromMap(iScCh);
	        scSec[iScCh]     = rawStreamZDC->GetScTowerFromMap(iScCh);
	        iScCh++;
	      }
	      else if(modType[iMod]==6 && modGeo[iMod]==4){ // ZDC TDC mapping --------------------
	        tdcMod[itdcCh]     = rawStreamZDC->GetTDCModFromMap(itdcCh);
	        tdcCh[itdcCh]      = rawStreamZDC->GetTDCChFromMap(itdcCh);
	        tdcSigCode[itdcCh] = rawStreamZDC->GetTDCSignFromMap(itdcCh);
	        itdcCh++;
	      }
	    }
 	  }
	  // Writing data on output FXS file
	  for(Int_t is=0; is<2*kNChannels; is++){
	     fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",
	       is,adcMod[is],adcCh[is],sigCode[is],det[is],sec[is]);
    	     //printf("  EMD DA -> %d ADC: mod %d ch %d, code %d det %d, sec %d\n",
	     //  is,adcMod[is],adcCh[is],sigCode[is],det[is],sec[is]);
	  }
	  for(Int_t is=0; is<kNScChannels; is++){
	     fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",
	       is,scMod[is],scCh[is],scSigCode[is],scDet[is],scSec[is]);
 	     //printf("  EMD DA -> %d Scaler: mod %d ch %d, code %d det %d, sec %d\n",
	     //  is,scMod[is],scCh[is],scSigCode[is],scDet[is],scSec[is]);
	  }
	  for(Int_t is=0; is<kNScChannels; is++){
	     fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\n",
	       is,tdcMod[is],tdcCh[is],tdcSigCode[is]);
 	     //if(tdcMod[is]!=-1) printf("  Mapping DA -> %d TDC: mod %d ch %d, code %d\n",
	     //  is,tdcMod[is],tdcCh[is],tdcSigCode[is]);
	  }
	  for(Int_t is=0; is<kNModules; is++){
	     fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\n",
	     modGeo[is],modType[is],modNCh[is]);
	     //printf("  EMD DA -> Module mapping: geo %d type %d #ch %d\n",
	     //  modGeo[is],modType[is],modNCh[is]);
	  }
	  
	}
        fclose(mapFile4Shuttle);
      }// SOD event
    
     else if(eventT==PHYSICS_EVENT){
      // --- Reading data header
      reader->ReadHeader();
      const AliRawDataHeader* header = reader->GetDataHeader();
      if(header){
         UChar_t message = header->GetAttributes();
	 if((message & 0x70) == 0x70){ // DEDICATED EMD RUN
	    //printf("\t CALIBRATION_EMD_RUN raw data found\n");
	 }
	 else{
	    // Temporary commented to validate DA with simulated data!!!
	    //printf("\t NO CALIBRATION_EMD_RUN raw data found\n");
	    //return -1;
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
      Float_t ZDCCorrADC[4], ZDCCorrADCSum[4];
      for(Int_t g=0; g<4; g++){
           ZDCCorrADCSum[g] = 0.;
	   ZDCCorrADC[g] = 0.;
      }

      //
      while(rawStreamZDC->Next()){
	Int_t detID = rawStreamZDC->GetSector(0);
	Int_t quad = rawStreamZDC->GetSector(1);
	
	if((rawStreamZDC->IsADCDataWord()) && !(rawStreamZDC->IsUnderflow())
	   && !(rawStreamZDC->IsOverflow()) && (detID!=-1) && (detID!=3) 
	   && (rawStreamZDC->GetADCGain() == 1) && // Selecting LOW RES ch.s
           (rawStreamZDC->GetADCModule()>=kFirstADCGeo && rawStreamZDC->GetADCModule()<=kLastADCGeo)){
	  
	  // Taking LOW RES channels -> ch.+kNChannels !!!!
	  Int_t DetIndex=999, PedIndex=999;
	  // Not PMRef
	  if(quad!=5){ 
	    if(detID == 1){
	      DetIndex = detID-1;
	      PedIndex = quad; 
	    }
	    else if(detID == 2){
	      DetIndex = detID-1;
	      PedIndex = quad+5;
	    }
	    else if(detID == 4){
	      DetIndex = detID-2;
	      PedIndex = quad+12;
	    }
	    else if(detID == 5){
	      DetIndex = detID-2;
	      PedIndex = quad+17;
	    }
	    // Mean pedestal subtraction 
	    Float_t Pedestal = MeanPedlg[PedIndex];
	    // Pedestal subtraction from correlation with out-of-time signals
	    //Float_t Pedestal = CorrCoeff0[PedIndex]+CorrCoeff1[PedIndex]*MeanPedOOT[PedIndex];
	    
	    // Run 2010 -> we decide to fit only PMCs
/*	    if(DetIndex!=999 || PedIndex!=999){ 
	      if(quad==0){ 
	        Float_t corrADCval = (rawStreamZDC->GetADCValue()) - Pedestal;
		if(detID==1 || detID==2) histoEMDCorr[detID-1]->Fill(corrADCval);
		else if(detID==4 || detID==5) histoEMDCorr[detID-2]->Fill(corrADCval);
	      }
	    }
*/
	    if(DetIndex!=999 || PedIndex!=999){ 
	      //
	      ZDCCorrADC[DetIndex] = (rawStreamZDC->GetADCValue()) - Pedestal;
	      ZDCCorrADCSum[DetIndex] += ZDCCorrADC[DetIndex];
	      //
	      /*printf("\t det %d quad %d res %d pedInd %d "
	         "Pedestal %1.0f -> ADCCorr = %d ZDCCorrADCSum = %d\n", 
	         detID,quad,rawStreamZDC->GetADCGain(),PedIndex,Pedestal, 
	         (Int_t) ZDCCorrADC[DetIndex],(Int_t) ZDCCorrADCSum[DetIndex]);
	      */	 
	    }
	    // Not common PM 
	    /*if(quad!=0){
	      Float_t corrADCval = (rawStreamZDC->GetADCValue()) - Pedestal;
	      if(detID==1)      histZNCtow[quad-1]->Fill(corrADCval);
	      else if(detID==2) histZPCtow[quad-1]->Fill(corrADCval);
	      else if(detID==4) histZNAtow[quad-1]->Fill(corrADCval);
	      else if(detID==5) histZPAtow[quad-1]->Fill(corrADCval);
	      //
	      //printf("\t det %d tow %d fill histo w. value %1.0f\n", 
	      //  detID,quad,corrADCval);
	    }*/

	    if(DetIndex==999 || PedIndex==999) 
	       printf(" WARNING! Detector a/o pedestal index are WRONG!!!\n");
 
	  }//quad!=5
	}//IsADCDataWord()
	
       }
       //
       nevents_physics++;
       //
       delete reader;
       delete rawStreamZDC;
       //
       for(Int_t j=0; j<4; j++) histoEMDCorr[j]->Fill(ZDCCorrADCSum[j]);

    }//(if PHYSICS_EVENT)
      
    /* exit when last event received, no need to wait for TERM signal */
    else if(eventT==END_OF_RUN) {
      printf(" -> EOR event detected\n");
      break;
    }
    
    nevents_total++;

   }

   /* free resources */
   free(event);
  }
    
  /* Analysis of the histograms */
  //
  FILE *fileShuttle1 = fopen(ENCALIBDATA_FILE,"w");
  //
  Int_t BinMax[4]={0,0,0,0};
  Float_t YMax[4]={0.,0.,0.,0.};
  Int_t NBinsx[4]={0,0,0,0};
  Float_t MeanFitVal[4]={1.,1.,1.,1.};
  TF1 *fitfun[4];
  for(Int_t k=0; k<4; k++){
     if(histoEMDCorr[k]->GetEntries()!=0 && histoEMDCorr[k]->GetMean()>0){
       BinMax[k] = histoEMDCorr[k]->GetMaximumBin();
       //printf("\n\t histoEMDCorr[%d]: BinMax = %d", k, BinMax[k]);
       if(BinMax[k]<=1){
         printf("\n WARNING! Something wrong with det %d histo -> calibration coeff. set to 1\n\n", k);
         continue;
       }
       // 
       YMax[k] = (histoEMDCorr[k]->GetXaxis())->GetXmax();
       NBinsx[k] = (histoEMDCorr[k]->GetXaxis())->GetNbins();
       //printf(" ChXMax = %f\n", BinMax[k]*YMax[k]/NBinsx[k]);
       histoEMDCorr[k]->Fit("gaus","Q","",BinMax[k]*YMax[k]/NBinsx[k]*0.7,BinMax[k]*YMax[k]/NBinsx[k]*1.25);
       fitfun[k] = histoEMDCorr[k]->GetFunction("gaus");
       MeanFitVal[k] = (Float_t) (fitfun[k]->GetParameter(1));
       //printf("\t Mean Value from gaussian fit = %f\n", MeanFitVal[k]);
     }
  }
  //
  Float_t CalibCoeff[6];     
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
 
  FILE *fileShuttle2 = fopen(TOWCALIBDATA_FILE,"w");
  //Float_t meanvalznc[4], meanvalzpc[4], meanvalzna[4], meanvalzpa[4];
  for(Int_t j=0; j<4; j++){
     /*if(histZNCtow[j]->GetEntries() == 0){
       printf("\n WARNING! Empty histos -> ending DA WITHOUT writing output\n\n");
       return -1;
     } 
     meanvalznc[j] = histZNCtow[j]->GetMean();
     meanvalzpc[j] = histZPCtow[j]->GetMean();
     meanvalzna[j] = histZNAtow[j]->GetMean();
     meanvalzpa[j] = histZPAtow[j]->GetMean();*/
     
     // Note -> For the moment the inter-calibration coeff. are set to 1 
     for(Int_t k=0; k<5; k++){  
       Float_t icoeff = 1.;
       fprintf(fileShuttle2,"\t%f",icoeff);
       if(k==5) fprintf(fileShuttle2,"\n");
     }
  }
  //
  /*if(meanvalznc[1]!=0 && meanvalznc[2]!=0 && meanvalznc[3]!=0 && 
     meanvalzpc[1]!=0 && meanvalzpc[2]!=0 && meanvalzpc[3]!=0 &&
     meanvalzna[1]!=0 && meanvalzna[2]!=0 && meanvalzna[3]!=0 &&
     meanvalzpa[1]!=0 && meanvalzpa[2]!=0 && meanvalzpa[3]!=0){
    fprintf(fileShuttle2,"\t%f\t%f\t%f\t%f\n",
  	1.0,meanvalznc[0]/meanvalznc[1],meanvalznc[0]/meanvalznc[2],meanvalznc[0]/meanvalznc[3]);
    fprintf(fileShuttle2,"\t%f\t%f\t%f\t%f\n",
  	1.0,meanvalzpc[0]/meanvalzpc[1],meanvalzpc[0]/meanvalzpc[2],meanvalzpc[0]/meanvalzpc[3]);
    fprintf(fileShuttle2,"\t%f\t%f\t%f\t%f\n",
  	1.0,meanvalzna[0]/meanvalzna[1],meanvalzpc[0]/meanvalzna[2],meanvalzpc[0]/meanvalzna[3]);
    fprintf(fileShuttle2,"\t%f\t%f\t%f\t%f\n",
  	1.0,meanvalzpa[0]/meanvalzpa[1],meanvalzpc[0]/meanvalzpa[2],meanvalzpc[0]/meanvalzpa[3]);
  }
  else{
    printf("\n Tower intercalib. coeff. CAN'T be calculated (some mean values are ZERO)!!!\n\n");
    return -1;
  }*/
  fclose(fileShuttle2);
  
  for(Int_t ij=0; ij<4; ij++){
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
