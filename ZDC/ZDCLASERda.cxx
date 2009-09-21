/*

This program reads the DAQ data files passed as argument using the monitoring library.

It computes the average event size and populates local "./result.txt" file with the 
result.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalone pedestal runs

Contact: Chiara.Oppedisano@to.infn.it
Link: 
Run Type: STANDALONE_LASER_RUN
DA Type: LDC
Number of events needed: no constraint (tipically ~10^3)
Input Files: 
Output Files: ZDCLaser.dat
Trigger Types Used: Standalone Trigger

*/
#define PEDDATA_FILE  "ZDCPedestal.dat"
#define MAPDATA_FILE  "ZDCChMapping.dat"
#define LASHISTO_FILE "ZDCLaserHisto.root"
#define LASDATA_FILE  "ZDCLaserCalib.dat"

#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>

// DATE
#include <event.h>
#include <monitor.h>
#include <daqDA.h>

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
      Arguments: list of DATE raw data files
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
  int const kNChannels = 24;
  int const kNScChannels = 32;

  /* log start of process */
  printf("\n ZDC LASER program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // --- Histograms for LASER runs 
  //     20 signal channels + 2 reference PTMs
  //
  TH1F::AddDirectory(0);
  // --- Histos for reference PMTs (high gain chains)
  TH1F *hPMRefChg = new TH1F("hPMRefChg","hPMRefChg", 100,0.,1400.);
  TH1F *hPMRefAhg = new TH1F("hPMRefAhg","hPMRefAhg", 100,0.,1400.);
  TH1F *hPMRefClg = new TH1F("hPMRefClg","hPMRefClg", 100,0.,4000.);
  TH1F *hPMRefAlg = new TH1F("hPMRefAlg","hPMRefAlg", 100,0.,4000.);
  //
  // --- Histos for detector PMTs 
  TH1F *hZNChg[5], *hZPChg[5], *hZNAhg[5], *hZPAhg[5], *hZEMhg[2];
  TH1F *hZNClg[5], *hZPClg[5], *hZNAlg[5], *hZPAlg[5], *hZEMlg[2];
  char hnamZNChg[20], hnamZPChg[20], hnamZNAhg[20], hnamZPAhg[20];
  char hnamZNClg[20], hnamZPClg[20], hnamZNAlg[20], hnamZPAlg[20];
  char hnamZEMhg[20], hnamZEMlg[20];
  for(Int_t j=0; j<5; j++){
    sprintf(hnamZNChg,"ZNChg-tow%d",j);
    sprintf(hnamZPChg,"ZPChg-tow%d",j);
    sprintf(hnamZNAhg,"ZNAhg-tow%d",j);
    sprintf(hnamZPAhg,"ZPAhg-tow%d",j);
    //
    hZNChg[j] = new TH1F(hnamZNChg, hnamZNChg, 100, 0., 1400.);
    hZPChg[j] = new TH1F(hnamZPChg, hnamZPChg, 100, 0., 1400.);
    hZNAhg[j] = new TH1F(hnamZNAhg, hnamZNAhg, 100, 0., 1400.);
    hZPAhg[j] = new TH1F(hnamZPAhg, hnamZPAhg, 100, 0., 1400.);
    //
    sprintf(hnamZNClg,"ZNClg-tow%d",j);
    sprintf(hnamZPClg,"ZPClg-tow%d",j);
    sprintf(hnamZNAlg,"ZNAlg-tow%d",j);
    sprintf(hnamZPAlg,"ZPAlg-tow%d",j);
    //
    hZNClg[j] = new TH1F(hnamZNClg, hnamZNClg, 100, 0., 4000.);
    hZPClg[j] = new TH1F(hnamZPClg, hnamZPClg, 100, 0., 4000.);
    hZNAlg[j] = new TH1F(hnamZNAlg, hnamZNAlg, 100, 0., 4000.);
    hZPAlg[j] = new TH1F(hnamZPAlg, hnamZPAlg, 100, 0., 4000.);
    //
    if(j<2){
      sprintf(hnamZEMhg,"ZEM%dhg",j);
      sprintf(hnamZEMlg,"ZEM%dlg",j);
      //
      hZEMhg[j] = new TH1F(hnamZEMhg, hnamZEMhg, 100, 0., 1400.);      
      hZEMlg[j] = new TH1F(hnamZEMlg, hnamZEMlg, 100, 0., 4000.);      
    }
  }

  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
  }
  /* report progress */
  daqDA_progressReport(10);
        
  // *** To analyze LASER events you MUST have a pedestal data file!!!
  // *** -> check if a pedestal run has been analyzed
  int read = 0;
  read = daqDA_DB_getFile(PEDDATA_FILE, PEDDATA_FILE);
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
  for(n=1;n<argc;n++) {
   
    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }

    /* report progress */
    /* in this example, indexed on the number of files */
    daqDA_progressReport(20+70*n/argc);

    /* read the file */
    for(;;) {

      /* get next event */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF) break; /* end of monitoring file has been reached */
      if (status!=0) {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        return -1;
      }

      /* retry if got no event */
      if (event==NULL) {
        break;
      }

      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      reader->Select("ZDC");
      // --- Reading event header
      //UInt_t evtype = reader->GetType();
      //printf("\n\t ZDCLASERda -> ev. type %d\n",evtype);
      //printf("\t ZDCLASERda -> run # %d\n",reader->GetRunNumber());
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

      else if(eventT==PHYSICS_EVENT){
 	// --- Reading data header
        reader->ReadHeader();
        const AliRawDataHeader* header = reader->GetDataHeader();
        if(header) {
         UChar_t message = header->GetAttributes();
	 if(message & 0x30){ // DEDICATED LASER RUN
	    //printf("\t STANDALONE_LASER_RUN raw data found\n");
	    continue;
	 }
	 else{
	    printf("ZDCLASERda.cxx -> NO STANDALONE_LASER_RUN raw data found\n");
	    return -1;
	 }
  	}
  	else{
  	   printf("\t ATTENTION! No Raw Data Header found!!!\n");
  	   return -1;
  	}

	rawStreamZDC->SetSODReading(kTRUE);

  	if (!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
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
	while(rawStreamZDC->Next()){
  	  Int_t index=-1;
	  Int_t detector = rawStreamZDC->GetSector(0);
	  Int_t sector = rawStreamZDC->GetSector(1);
	  
  	  if(rawStreamZDC->IsADCDataWord() && !(rawStreamZDC->IsUnderflow())
	     && !(rawStreamZDC->IsOverflow()) && detector!=-1){
	    
	    //printf("  IsADCWord %d, IsUnderflow %d, IsOverflow %d\n",
	    //  rawStreamZDC->IsADCDataWord(),rawStreamZDC->IsUnderflow(),rawStreamZDC->IsOverflow());
 
 	    if(sector!=5){ // Physics signals
    	      if(detector==1) index = sector;        // *** ZNC
	      else if(detector==2) index = sector+5; // *** ZPC
	      else if(detector==3) index = sector+9; // *** ZEM
	      else if(detector==4) index = sector+12;// *** ZNA
	      else if(detector==5) index = sector+17;// *** ZPA
	    }
	    else{ // Reference PMs
	      index = (detector-1)/3+22;
	    }
	    //
	    if(index==-1) printf("ERROR in ZDCLASERda.cxx -> det %d quad %d res %d index %d ADC %d\n", 
	      detector, sector, rawStreamZDC->GetADCGain(), index, rawStreamZDC->GetADCValue());
	    
	    Float_t Pedestal=0.;
	    if(rawStreamZDC->GetADCGain()==0) Pedestal = MeanPedhg[index];
	    else if(rawStreamZDC->GetADCGain()==1) Pedestal = MeanPedlg[index];
	    //
	    Float_t CorrADC = rawStreamZDC->GetADCValue() - Pedestal;
	    //
	    //printf("\tdet %d sec %d res %d index %d ped %1.0f ADCcorr %1.0f\n", 
	    //  detector, sector, rawStreamZDC->GetADCGain(), index, Pedestal,CorrADC);
	    
	    // **** Detector PMs
	    if(sector!=5){
	      if(rawStreamZDC->GetADCGain()==0){ // --- High gain chain ---
	        // ---- side C
	        if(detector==1) hZNChg[sector]->Fill(CorrADC);
	        else if(detector==2) hZPChg[sector]->Fill(CorrADC);
	        // ---- side A
	        else if(detector==4) hZNAhg[sector]->Fill(CorrADC);
	        else if(detector==5) hZPAhg[sector]->Fill(CorrADC);
	        // ---- ZEM
		else if(detector==3) hZEMhg[sector-1]->Fill(CorrADC);
	      }
	      else if(rawStreamZDC->GetADCGain()==1){ // --- Low gain chain ---
	        // ---- side C
	        if(detector==1) hZNClg[sector]->Fill(CorrADC);
	        else if(detector==2) hZPClg[sector]->Fill(CorrADC);
	        // ---- side A
	        else if(detector==4) hZNAlg[sector]->Fill(CorrADC);
	        else if(detector==5) hZPAlg[sector]->Fill(CorrADC);
	        // ---- ZEM
		else if(detector==3) hZEMlg[sector-1]->Fill(CorrADC);
	      }
	    }
	    // **** Reference PMs
	    else if(sector==5){
	      if(rawStreamZDC->GetADCGain()==0){ // --- High gain chain ---
	        // ---- PMRef chain side C
	        if(detector==1) hPMRefChg->Fill(CorrADC);
	        // ---- PMRef side A
	        else if(detector==4) hPMRefAhg->Fill(CorrADC);
	      }
	      else if(rawStreamZDC->GetADCGain()==1){ // --- Low gain chain ---
	        // ---- PMRef chain side C
	        if(detector==1) hPMRefClg->Fill(CorrADC);
	        // ---- PMRef side A
	        else if(detector==4) hPMRefAlg->Fill(CorrADC);
	      }
	    }
  	  }//IsADCDataWord()+NOunderflow+NOoverflow
  	  //
         }
         //
         nevents_physics++;
         //
	 delete reader;
         delete rawStreamZDC;

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
  Int_t det[2*kNChannels], quad[2*kNChannels];
  Int_t maxBin[2*kNChannels], nBin[2*kNChannels];
  Float_t xMax[2*kNChannels], maxXval[2*kNChannels], xlow[2*kNChannels]; 
  Float_t mean[2*kNChannels], sigma[2*kNChannels];
  TF1 *fun[2*kNChannels];
  
  // ******** High gain chain ********
  for(Int_t k=0; k<5; k++){
    // --- ZNC
    det[k] = 1;
    quad[k] = k;
    maxBin[k] = hZNChg[k]->GetMaximumBin();
    nBin[k] = (hZNChg[k]->GetXaxis())->GetNbins();
    xMax[k] = (hZNChg[k]->GetXaxis())->GetXmax();
    if(nBin[k]!=0) maxXval[k] = maxBin[k]*xMax[k]/nBin[k];
    if(maxXval[k]-150.<0.) xlow[k]=0.;
    else xlow[k] = maxXval[k]-150.;
    // checking if histos are empty
    if(hZNChg[k]->GetEntries() == 0){
      printf("\n WARNING! Empty LASER histos -> ending DA WITHOUT writing output\n\n");
      return -1;
    } 
    //
    hZNChg[k]->Fit("gaus","Q","",xlow[k],maxXval[k]+150.);
    fun[k] = hZNChg[k]->GetFunction("gaus");
    mean[k]  = (Float_t) (fun[k]->GetParameter(1));
    sigma[k] = (Float_t) (fun[k]->GetParameter(2));
    // --- ZPC
    det[k+5] = 2;
    quad[k+5] = k;
    maxBin[k+5] = hZPChg[k]->GetMaximumBin();
    nBin[k+5] = (hZPChg[k]->GetXaxis())->GetNbins();
    xMax[k+5] = (hZPChg[k]->GetXaxis())->GetXmax();
    if(nBin[k+5]!=0) maxXval[k+5] = maxBin[k+5]*xMax[k+5]/nBin[k+5];
    if(maxXval[k+5]-150.<0.) xlow[k+5]=0.;
    else xlow[k+5] = maxXval[k+5]-150.;
    hZPChg[k]->Fit("gaus","Q","",xlow[k+5],maxXval[k+5]+150.);
    fun[k+5] = hZPChg[k]->GetFunction("gaus");
    mean[k+5]  = (Float_t) (fun[k+5]->GetParameter(1));
    sigma[k+5] = (Float_t) (fun[k+5]->GetParameter(2));
    // --- ZEM1
    if(k<2){
      det[k+10] = 3;
      quad[k+10] = k+1;
      maxBin[k+10] = hZEMhg[k]->GetMaximumBin();
      nBin[k+10] = (hZEMhg[k]->GetXaxis())->GetNbins();
      xMax[k+10] = (hZEMhg[k]->GetXaxis())->GetXmax();
      if(nBin[k+10]!=0) maxXval[k+10] = maxBin[k+10]*xMax[k+10]/nBin[k+10];
      if(maxXval[k+10]-150.<0.) xlow[k+10]=0.;
      else xlow[k+10] = maxXval[k+10]-150.;
      hZEMhg[k]->Fit("gaus","Q","",xlow[k+10],maxXval[k+10]+150.);
      fun[k+10] = hZEMhg[k]->GetFunction("gaus");
      mean[k+10]  = (Float_t) (fun[k+10]->GetParameter(1));
      sigma[k+10] = (Float_t) (fun[k+10]->GetParameter(2));
    }
    // --- ZNA
    det[k+12] = 4;
    quad[k+12] = k;
    maxBin[k+12] = hZNAhg[k]->GetMaximumBin();
    nBin[k+12] = (hZNAhg[k]->GetXaxis())->GetNbins();
    xMax[k+12] = (hZNAhg[k]->GetXaxis())->GetXmax();
    if(nBin[k+12]!=0) maxXval[k+12] = maxBin[k+12]*xMax[k+12]/nBin[k+12];
    if(maxXval[k+12]-150.<0.) xlow[k+12]=0.;
    else xlow[k+12] = maxXval[k+12]-150.;
    hZNAhg[k]->Fit("gaus","Q","",xlow[k+12],maxXval[k+12]+150.);
    fun[k+12] = hZNAhg[k]->GetFunction("gaus");
    mean[k+12]  = (Float_t) (fun[k+12]->GetParameter(1));
    sigma[k+12] = (Float_t) (fun[k+12]->GetParameter(2));
    // --- ZPA
    det[k+17] = 4;
    quad[k+17] = 5;
    maxBin[k+17] = hZPAhg[k]->GetMaximumBin();
    nBin[k+17] = (hZPAhg[k]->GetXaxis())->GetNbins();
    xMax[k+17] = (hZPAhg[k]->GetXaxis())->GetXmax();
    if(nBin[k+17]!=0) maxXval[k+17] = maxBin[k+17]*xMax[k+17]/nBin[k+17];
    if(maxXval[k+17]-150.<0.) xlow[k+17]=0.;
    else xlow[k+17] = maxXval[k+17]-150.;
    hZPAhg[k]->Fit("gaus","Q","",xlow[k+17],maxXval[k+17]+150.);
    fun[k+17] = hZPAhg[k]->GetFunction("gaus");
    mean[k+17]  = (Float_t) (fun[k+17]->GetParameter(1));
    sigma[k+17] = (Float_t) (fun[k+17]->GetParameter(2));    
  }
  // ~~~~~~~~ PM Ref side C ~~~~~~~~
  det[22] = 1;
  quad[22] = 5;
  maxBin[22] = hPMRefChg->GetMaximumBin();
  nBin[22] = (hPMRefChg->GetXaxis())->GetNbins();
  xMax[22] = (hPMRefChg->GetXaxis())->GetXmax();
  if(nBin[22]!=0) maxXval[22] = maxBin[22]*xMax[22]/nBin[22];
  if(maxXval[22]-150.<0.) xlow[22]=0.;
  else xlow[22] = maxXval[22];
  hPMRefChg->Fit("gaus","Q","",xlow[22],maxXval[22]+150.);
  fun[22] = hPMRefChg->GetFunction("gaus");
  mean[22]  = (Float_t) (fun[22]->GetParameter(1));
  sigma[22] = (Float_t) (fun[22]->GetParameter(2));
  // ~~~~~~~~ PM Ref side A ~~~~~~~~
  det[23] = 4;
  quad[23] = 5;
  maxBin[23] = hPMRefAhg->GetMaximumBin();
  nBin[23] = (hPMRefAhg->GetXaxis())->GetNbins();
  xMax[23] = (hPMRefAhg->GetXaxis())->GetXmax();
  if(nBin[23]!=0) maxXval[23] = maxBin[23]*xMax[23]/nBin[23];
  if(maxXval[23]-100.<0.) xlow[23]=0.;
  else xlow[23] = maxXval[23];
  hPMRefAhg->Fit("gaus","Q","",xlow[23],maxXval[23]+100.);
  fun[23] = hPMRefAhg->GetFunction("gaus");
  mean[23]  = (Float_t) (fun[23]->GetParameter(1));
  sigma[23] = (Float_t) (fun[23]->GetParameter(2));
  
  // ******** Low gain chain ********
  Int_t kOffset = 24;
  for(Int_t k=0; k<5; k++){
    // --- ZNC
    det[k+kOffset] = 1;
    quad[k+kOffset] = k;
    maxBin[k+kOffset] = hZNClg[k]->GetMaximumBin();
    nBin[k+kOffset] = (hZNClg[k]->GetXaxis())->GetNbins();
    xMax[k+kOffset] = (hZNClg[k]->GetXaxis())->GetXmax();
    if(nBin[k+kOffset]!=0) maxXval[k+kOffset] = maxBin[k+kOffset]*xMax[k+kOffset]/nBin[k+kOffset];
    if(maxXval[k+kOffset]-150.<0.) xlow[k+kOffset]=0.;
    else xlow[k+kOffset] = maxXval[k+kOffset]-150.;
    hZNClg[k]->Fit("gaus","Q","",xlow[k+kOffset],maxXval[k+kOffset]+150.);
    fun[k+kOffset] = hZNClg[k]->GetFunction("gaus");
    mean[k+kOffset]  = (Float_t) (fun[k+kOffset]->GetParameter(1));
    sigma[k+kOffset] = (Float_t) (fun[k+kOffset]->GetParameter(2));
    // --- ZPC
    det[k+kOffset+5] = 2;
    quad[k+kOffset+5] = k;
    maxBin[k+kOffset+5] = hZPClg[k]->GetMaximumBin();
    nBin[k+kOffset+5] = (hZPClg[k]->GetXaxis())->GetNbins();
    xMax[k+kOffset+5] = (hZPClg[k]->GetXaxis())->GetXmax();
    if(nBin[k+kOffset+5]!=0) maxXval[k+kOffset+5] = maxBin[k+kOffset+5]*xMax[k+kOffset+5]/nBin[k+kOffset+5];
    if(maxXval[k+kOffset+5]-150.<0.) xlow[k+kOffset+5]=0.;
    else xlow[k+kOffset+5] = maxXval[k+kOffset+5]-150.;
    hZPClg[k]->Fit("gaus","Q","",xlow[k+kOffset+5],maxXval[k+kOffset+5]+150.);
    fun[k+kOffset+5] = hZPClg[k]->GetFunction("gaus");
    mean[k+kOffset+5]  = (Float_t) (fun[k+kOffset+5]->GetParameter(1));
    sigma[k+kOffset+5] = (Float_t) (fun[k+kOffset+5]->GetParameter(2));
    // --- ZEM1
    if(k+kOffset<2){
      det[k+kOffset+10] = 3;
      quad[k+kOffset+10] = k+1;
      maxBin[k+kOffset+10] = hZEMlg[k]->GetMaximumBin();
      nBin[k+kOffset+10] = (hZEMlg[k]->GetXaxis())->GetNbins();
      xMax[k+kOffset+10] = (hZEMlg[k]->GetXaxis())->GetXmax();
      if(nBin[k+kOffset+10]!=0) maxXval[k+kOffset+10] = maxBin[k+kOffset+10]*xMax[k+kOffset+10]/nBin[k+kOffset+10];
      if(maxXval[k+kOffset+10]-150.<0.) xlow[k+kOffset+10]=0.;
      else xlow[k+kOffset+10] = maxXval[k+kOffset+10]-150.;
      hZEMlg[k]->Fit("gaus","Q","",xlow[k+kOffset+10],maxXval[k+kOffset+10]+150.);
      fun[k+kOffset+10] = hZEMlg[k]->GetFunction("gaus");
      mean[k+kOffset+10]  = (Float_t) (fun[k+kOffset+10]->GetParameter(1));
      sigma[k+kOffset+10] = (Float_t) (fun[k+kOffset+10]->GetParameter(2));
    }
    // --- ZNA
    det[k+kOffset+12] = 4;
    quad[k+kOffset+12] = k;
    maxBin[k+kOffset+12] = hZNAlg[k]->GetMaximumBin();
    nBin[k+kOffset+12] = (hZNAlg[k]->GetXaxis())->GetNbins();
    xMax[k+kOffset+12] = (hZNAlg[k]->GetXaxis())->GetXmax();
    if(nBin[k+kOffset+12]!=0) maxXval[k+kOffset+12] = maxBin[k+kOffset+12]*xMax[k+kOffset+12]/nBin[k+kOffset+12];
    if(maxXval[k+kOffset+12]-150.<0.) xlow[k+kOffset+12]=0.;
    else xlow[k+kOffset+12] = maxXval[k+kOffset+12]-150.;
    hZNAlg[k]->Fit("gaus","Q","",xlow[k+kOffset+12],maxXval[k+kOffset+12]+150.);
    fun[k+kOffset+12] = hZNAlg[k]->GetFunction("gaus");
    mean[k+kOffset+12]  = (Float_t) (fun[k+kOffset+12]->GetParameter(1));
    sigma[k+kOffset+12] = (Float_t) (fun[k+kOffset+12]->GetParameter(2));
    // --- ZPA
    det[k+kOffset+17] = 5;
    quad[k+kOffset+17] = k;
    maxBin[k+kOffset+17] = hZPAlg[k]->GetMaximumBin();
    nBin[k+kOffset+17] = (hZPAlg[k]->GetXaxis())->GetNbins();
    xMax[k+kOffset+17] = (hZPAlg[k]->GetXaxis())->GetXmax();
    if(nBin[k+kOffset+17]!=0) maxXval[k+kOffset+17] = maxBin[k+kOffset+17]*xMax[k+kOffset+17]/nBin[k+kOffset+17];
    if(maxXval[k+kOffset+17]-150.<0.) xlow[k+kOffset+17]=0.;
    else xlow[k+kOffset+17] = maxXval[k+kOffset+17]-150.;
    hZPAlg[k]->Fit("gaus","Q","",xlow[k+kOffset+17],maxXval[k+kOffset+17]+150.);
    fun[k+kOffset+17] = hZPAlg[k]->GetFunction("gaus");
    mean[k+kOffset+17]  = (Float_t) (fun[k+kOffset+17]->GetParameter(1));
    sigma[k+kOffset+17] = (Float_t) (fun[k+kOffset+17]->GetParameter(2));    
  }
  // ~~~~~~~~ PM Ref side C ~~~~~~~~
  det[46] = 1;
  quad[46] = 5;
  maxBin[46] = hPMRefClg->GetMaximumBin();
  nBin[46] = (hPMRefClg->GetXaxis())->GetNbins();
  xMax[46] = (hPMRefClg->GetXaxis())->GetXmax();
  if(nBin[46]!=0) maxXval[46] = maxBin[46]*xMax[46]/nBin[46];
  if(maxXval[46]-150.<0.) xlow[46]=0.;
  else xlow[46] = maxXval[46];
  hPMRefClg->Fit("gaus","Q","",xlow[46],maxXval[46]+150.);
  fun[46] = hPMRefClg->GetFunction("gaus");
  mean[46]  = (Float_t) (fun[46]->GetParameter(1));
  sigma[46] = (Float_t) (fun[46]->GetParameter(2));
  // ~~~~~~~~ PM Ref side A ~~~~~~~~
  det[47] = 4;
  quad[47] = 5;
  maxBin[47] = hPMRefAlg->GetMaximumBin();
  nBin[47] = (hPMRefAlg->GetXaxis())->GetNbins();
  xMax[47] = (hPMRefAlg->GetXaxis())->GetXmax();
  if(nBin[47]!=0) maxXval[47] = maxBin[47]*xMax[47]/nBin[47];
  if(maxXval[47]-100.<0.) xlow[47]=0.;
  else xlow[47] = maxXval[47];
  hPMRefAlg->Fit("gaus","Q","",xlow[47],maxXval[47]+100.);
  fun[47] = hPMRefAlg->GetFunction("gaus");
  mean[47]  = (Float_t) (fun[47]->GetParameter(1));
  sigma[47] = (Float_t) (fun[47]->GetParameter(2));
    
  FILE *fileShuttle;
  fileShuttle = fopen(LASDATA_FILE,"w");
  for(Int_t i=0; i<2*kNChannels; i++){
    fprintf(fileShuttle,"\t%d\t%d\t%f\t%f\n",det[i],quad[i],mean[i], sigma[i]); 
  }
  //						       
  fclose(fileShuttle);
  /* report progress */
  daqDA_progressReport(80);
  //
  TFile *histofile = new TFile(LASHISTO_FILE,"RECREATE");
  histofile->cd();
  for(int j=0; j<5; j++){
     hZNChg[j]->Write();
     hZPChg[j]->Write();
     hZNAhg[j]->Write();
     hZPAhg[j]->Write();
     hZNClg[j]->Write();
     hZPClg[j]->Write();
     hZNAlg[j]->Write();
     hZPAlg[j]->Write();  
     if(j<2){
       hZEMhg[j]->Write();
       hZEMlg[j]->Write();
    }
  }
  hPMRefChg->Write();
  hPMRefAhg->Write();
  hPMRefClg->Write();
  hPMRefAlg->Write();  
  //
  histofile->Close();
  //
  for(Int_t j=0; j<5; j++){
    delete hZNChg[j];
    delete hZPChg[j];
    delete hZNAhg[j];
    delete hZPAhg[j];
    delete hZNClg[j];
    delete hZPClg[j];
    delete hZNAlg[j];
    delete hZPAlg[j];
    if(j<2){
      delete hZEMhg[j];
      delete hZEMlg[j];
    }
  }
  delete hPMRefChg;
  delete hPMRefAhg;
  delete hPMRefClg;
  delete hPMRefAlg;

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);
  
  /* report progress */
  daqDA_progressReport(90);
  
  /* store the result file on FES */
  // [1] File with mapping
  status = daqDA_FES_storeFile(MAPDATA_FILE, "MAPPING");
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }
  //
  // [2] File with laser data
  status = daqDA_FES_storeFile(LASDATA_FILE, "LASERDATA");
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }
  // [3] File with laser histos
  status = daqDA_FES_storeFile(LASHISTO_FILE, "LASERHISTOS");
  if(status){
    printf("Failed to export pedestal histos file to DAQ FES\n");
    return -1;
  }

  /* report progress */
  daqDA_progressReport(100);

  return status;
}
