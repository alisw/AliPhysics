/*

This program reads the DAQ data files passed as argument using the monitoring library.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalone pedestal runs

Contact: Chiara.Oppedisano@to.infn.it
Link: 
Run Type: STANDALONE_LASER_RUN
DA Type: LDC
Number of events needed: no constraint (tipically ~10^3)
Input Files: ZDCPedestal.dat
Output Files: ZDCLaserCalib.dat, ZDCLaserHisto.root, ZDCChMapping.dat
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
#include <TObjArray.h>

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
  
  /* log start of process */
  printf("\n ZDC LASER program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  //TMinuitMinimizer m; 
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
  					"Minuit", 
					"TMinuitMinimizer",
      					"Minuit", 
					"TMinuitMinimizer(const char *)");
  TVirtualFitter::SetDefaultFitter("Minuit");

  
  int status = 0;
  int const kNModules = 9;
  int const kNChannels = 24;
  int const kNScChannels = 32;
  Int_t kFirstADCGeo=0, kLastADCGeo=1; // NO out-of-time signals!!!
      
 // *** initializations ***************************
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

  // *** Creating a container for the histos ***************************
  TObjArray *hList = new TObjArray(0);
  hList->SetOwner(kTRUE);

  // --- Histograms for LASER runs 
  //     20 signal channels + 2 reference PTMs
  //
  TH1F::AddDirectory(0);
  // --- Histos for reference PMTs (high gain chains)
  TH1F *hPMRefChg = new TH1F("hPMRefChg","hPMRefChg", 100,-100.5,1100.5);
  TH1F *hPMRefAhg = new TH1F("hPMRefAhg","hPMRefAhg", 100,-100.5,1100.5);
  TH1F *hPMRefClg = new TH1F("hPMRefClg","hPMRefClg", 100,-100.5,4900.5);
  TH1F *hPMRefAlg = new TH1F("hPMRefAlg","hPMRefAlg", 100,-100.5,4900.5);
  hList->Add(hPMRefChg);
  hList->Add(hPMRefAhg);
  hList->Add(hPMRefClg);
  hList->Add(hPMRefAlg);
  //
  // --- Histos for detector PMTs 
  TH1F *hZNChg[5], *hZPChg[5], *hZNAhg[5], *hZPAhg[5], *hZEMhg[2];
  TH1F *hZNClg[5], *hZPClg[5], *hZNAlg[5], *hZPAlg[5], *hZEMlg[2];
  for(Int_t j=0; j<5; j++){
    TString hname1 = Form("ZNChg-tow%d",j);
    hZNChg[j] = new TH1F(hname1, hname1, 100,-100.5,1100.5);
    TString hname2 = Form("ZPChg-tow%d",j);
    hZPChg[j] = new TH1F(hname2, hname2, 100,-100.5,1100.5);
    TString hname3 = Form("ZNAhg-tow%d",j);
    hZNAhg[j] = new TH1F(hname3, hname3, 100,-100.5,1100.5);
    TString hname4 = Form("ZPAhg-tow%d",j);
    hZPAhg[j] = new TH1F(hname4, hname4, 100,-100.5,1100.5);
    //
    hList->Add(hZNChg[j]);
    hList->Add(hZPChg[j]);
    hList->Add(hZNAhg[j]);
    hList->Add(hZPAhg[j]);
    //
    TString hname5 = Form("ZNClg-tow%d",j);
    hZNClg[j] = new TH1F(hname5, hname5, 200,-100.5,4900.5);
    TString hname6 = Form("ZPClg-tow%d",j);
    hZPClg[j] = new TH1F(hname6, hname6, 200,-100.5,4900.5);
    TString hname7 = Form("ZNAlg-tow%d",j);
    hZNAlg[j] = new TH1F(hname7, hname7, 200,-100.5,4900.5);
    TString hname8 = Form("ZPAlg-tow%d",j);
    hZPAlg[j] = new TH1F(hname8, hname8, 200,-100.5,4900.5);
    //
    hList->Add(hZNClg[j]);
    hList->Add(hZPClg[j]);
    hList->Add(hZNAlg[j]);
    hList->Add(hZPAlg[j]);
    //
    if(j<2){
      TString hname9 = Form("ZEM%dhg",j);
      TString hname10 = Form("ZEM%dlg",j);
      //
      hZEMhg[j] = new TH1F(hname9, hname9, 100,-100.5,1100.5);      
      hZEMlg[j] = new TH1F(hname10, hname10, 100,-100.5,4900.5);      
      //
      hList->Add(hZEMhg[j]);
      hList->Add(hZEMlg[j]);
    }
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
    /*else if(jj>4*kNChannels){
      CorrCoeff0[jj-4*kNChannels] = readValues[0][jj]; 
      CorrCoeff1[jj-4*kNChannels] = readValues[1][jj];;
    }*/
  }
  
  /* report progress */
  daqDA_progressReport(20);

  /* open mapping file for Shuttle */
  FILE *mapFile4Shuttle=NULL;
  mapFile4Shuttle = fopen(MAPDATA_FILE,"w");
  if(mapFile4Shuttle==NULL) {
    printf("Failed to open mapFile4Shuttle file\n");
    return -1;
  }

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
    daqDA_progressReport(10+80*n/argc);

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
      if (event==NULL) break;

      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      reader->Reset();
      reader->Select("ZDC");
      //
      AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
        

      /* use event - here, just write event id to result file */
      eventT=event->eventType;
      
      if(eventT==START_OF_DATA){
	
	iMod=-1; 
	ich=0; 
	iScCh=0;
	
	rawStreamZDC->SetSODReading(kTRUE);
	
	// --------------------------------------------------------
	// --- Writing ascii data file for the Shuttle preprocessor
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
    	     //printf("  Laser DA -> %d ADC: mod %d ch %d, code %d det %d, sec %d\n",
	     //  is,adcMod[is],adcCh[is],sigCode[is],det[is],sec[is]);
	  }
	  for(Int_t is=0; is<kNScChannels; is++){
	     fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",
	       is,scMod[is],scCh[is],scSigCode[is],scDet[is],scSec[is]);
 	     //printf("  Laser DA -> %d Scaler: mod %d ch %d, code %d det %d, sec %d\n",
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
	     //printf("  Laser DA -> Module mapping: geo %d type %d #ch %d\n",
	     //  modGeo[is],modType[is],modNCh[is]);
	  }
	  
	}
        fclose(mapFile4Shuttle);
      }// SOD event

      else if(eventT==PHYSICS_EVENT){
 	// --- Reading data header
/*        reader->ReadHeader();
        const AliRawDataHeader* header = reader->GetDataHeader();
        if(header) {
         UChar_t message = header->GetAttributes();
	 if((message & 0x30) == 0x30){ // DEDICATED LASER RUN
	    //printf("\t STANDALONE_LASER_RUN raw data found\n");
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
*/
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
	  
  	  if(rawStreamZDC->IsADCDataWord() && !(rawStreamZDC->IsUnderflow()) && 
	     !(rawStreamZDC->IsOverflow()) && detector!=-1 &&
             rawStreamZDC->GetADCModule()>=kFirstADCGeo && rawStreamZDC->GetADCModule()<=kLastADCGeo){
	    
 	    if(sector!=5){ // Physics signals
    	      if(detector==1)      index = sector;   // *** ZNC
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
	    if(rawStreamZDC->GetADCGain()==0)      Pedestal = MeanPedhg[index];
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
	        if(detector==1)      hZNChg[sector]->Fill(CorrADC);
	        else if(detector==2) hZPChg[sector]->Fill(CorrADC);
	        // ---- side A
	        else if(detector==4) hZNAhg[sector]->Fill(CorrADC);
	        else if(detector==5) hZPAhg[sector]->Fill(CorrADC);
	        // ---- ZEM
		else if(detector==3){
		  hZEMhg[sector-1]->Fill(CorrADC);
		}
	      }
	      else if(rawStreamZDC->GetADCGain()==1){ // --- Low gain chain ---
	        // ---- side C
	        if(detector==1)      hZNClg[sector]->Fill(CorrADC);
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
         delete rawStreamZDC;
	 delete reader;

      }//(if PHYSICS_EVENT) 

      /* exit when last event received, no need to wait for TERM signal */
      else if(eventT==END_OF_RUN) {
        printf(" -> EOR event detected\n");
        break;
      }
      
      
      nevents_total++;

      /* free resources */
      free(event);
    
    }
  }  
  
  /* Analysis of the histograms */
  //
  Int_t detector[2*kNChannels], quad[2*kNChannels];
  Int_t maxBin[2*kNChannels], nBin[2*kNChannels];
  Float_t xMax[2*kNChannels], maxXval[2*kNChannels], xlow[2*kNChannels]; 
  Float_t mean[2*kNChannels], sigma[2*kNChannels];
  for(Int_t t=0; t<2*kNChannels; t++){
    detector[t] = quad[t] = 0;
    maxBin[t] = nBin[t] = 0;
    xMax[t] = maxXval[t] = xlow[t] = 0.;
    mean[t] = sigma[t] = 0.;
  }
  Int_t atLeastOneHisto=0;
  
  // ******** High gain chain ********
  for(Int_t k=0; k<5; k++){
   // --- ZNC
   if(hZNChg[k]->GetMean()>0.){
    detector[k] = 1;
    quad[k] = k;
    maxBin[k] = hZNChg[k]->GetMaximumBin();
    nBin[k] = (hZNChg[k]->GetXaxis())->GetNbins();
    xMax[k] = (hZNChg[k]->GetXaxis())->GetXmax();
    if(nBin[k]!=0) maxXval[k] = maxBin[k]*xMax[k]/nBin[k];
    if(maxXval[k]-150.<0.) xlow[k]=0.;
    else xlow[k] = maxXval[k]-150.;
    //
    hZNChg[k]->Fit("gaus","Q","",xlow[k],maxXval[k]+150.);
    TF1 *fun = hZNChg[k]->GetFunction("gaus");
    mean[k]  = (Float_t) (fun->GetParameter(1));
    sigma[k] = (Float_t) (fun->GetParameter(2));
    atLeastOneHisto=1;
   }
   // --- ZPC
   if(hZPChg[k]->GetMean()>0.){
    detector[k+5] = 2;
    quad[k+5] = k;
    maxBin[k+5] = hZPChg[k]->GetMaximumBin();
    nBin[k+5] = (hZPChg[k]->GetXaxis())->GetNbins();
    xMax[k+5] = (hZPChg[k]->GetXaxis())->GetXmax();
    if(nBin[k+5]!=0) maxXval[k+5] = maxBin[k+5]*xMax[k+5]/nBin[k+5];
    if(maxXval[k+5]-150.<0.) xlow[k+5]=0.;
    else xlow[k+5] = maxXval[k+5]-150.;
    //
    hZPChg[k]->Fit("gaus","Q","",xlow[k+5],maxXval[k+5]+150.);
    TF1 *fun = hZPChg[k]->GetFunction("gaus");
    mean[k+5]  = (Float_t) (fun->GetParameter(1));
    sigma[k+5] = (Float_t) (fun->GetParameter(2));
    atLeastOneHisto=1; 
   }
   // --- ZNA
   if(hZNAhg[k]->GetMean()>0.){
    detector[k+12] = 4;
    quad[k+12] = k;
    maxBin[k+12] = hZNAhg[k]->GetMaximumBin();
    nBin[k+12] = (hZNAhg[k]->GetXaxis())->GetNbins();
    xMax[k+12] = (hZNAhg[k]->GetXaxis())->GetXmax();
    if(nBin[k+12]!=0) maxXval[k+12] = maxBin[k+12]*xMax[k+12]/nBin[k+12];
    if(maxXval[k+12]-150.<0.) xlow[k+12]=0.;
    else xlow[k+12] = maxXval[k+12]-150.;
      //
      hZNAhg[k]->Fit("gaus","Q","",xlow[k+12],maxXval[k+12]+150.);
      TF1 *fun = hZNAhg[k]->GetFunction("gaus");
      mean[k+12]  = (Float_t) (fun->GetParameter(1));
      sigma[k+12] = (Float_t) (fun->GetParameter(2));
      atLeastOneHisto=1; 
   }
   // --- ZPA
   if(hZPAhg[k]->GetMean()>0.){
    detector[k+17] = 4;
    quad[k+17] = 5;
    maxBin[k+17] = hZPAhg[k]->GetMaximumBin();
    nBin[k+17] = (hZPAhg[k]->GetXaxis())->GetNbins();
    xMax[k+17] = (hZPAhg[k]->GetXaxis())->GetXmax();
    if(nBin[k+17]!=0) maxXval[k+17] = maxBin[k+17]*xMax[k+17]/nBin[k+17];
    if(maxXval[k+17]-150.<0.) xlow[k+17]=0.;
    else xlow[k+17] = maxXval[k+17]-150.;
    //
    hZPAhg[k]->Fit("gaus","Q","",xlow[k+17],maxXval[k+17]+150.);
    TF1 *fun = hZPAhg[k]->GetFunction("gaus");
    mean[k+17]  = (Float_t) (fun->GetParameter(1));
    sigma[k+17] = (Float_t) (fun->GetParameter(2));    
    atLeastOneHisto=1; 
   }
  } // loop over 5 PMTs
  // ~~~~~~~~ PM Ref side C ~~~~~~~~
  if(hPMRefChg->GetMean()>0.){
   detector[22] = 1;
   quad[22] = 5;
   maxBin[22] = hPMRefChg->GetMaximumBin();
   nBin[22] = (hPMRefChg->GetXaxis())->GetNbins();
   xMax[22] = (hPMRefChg->GetXaxis())->GetXmax();
   if(nBin[22]!=0) maxXval[22] = maxBin[22]*xMax[22]/nBin[22];
   if(maxXval[22]-150.<0.) xlow[22]=0.;
   else xlow[22] = maxXval[22]-150.;
   //
   hPMRefChg->Fit("gaus","Q","",xlow[22],maxXval[22]+150.);
   TF1 *fun = hPMRefChg->GetFunction("gaus");
   mean[22]  = (Float_t) (fun->GetParameter(1));
   sigma[22] = (Float_t) (fun->GetParameter(2));
   atLeastOneHisto=1; 
  }
  // ~~~~~~~~ PM Ref side A ~~~~~~~~
  if(hPMRefAhg->GetMean()>0.){
   detector[23] = 4;
   quad[23] = 5;
   maxBin[23] = hPMRefAhg->GetMaximumBin();
   nBin[23] = (hPMRefAhg->GetXaxis())->GetNbins();
   xMax[23] = (hPMRefAhg->GetXaxis())->GetXmax();
   if(nBin[23]!=0) maxXval[23] = maxBin[23]*xMax[23]/nBin[23];
   if(maxXval[23]-100.<0.) xlow[23]=0.;
   else xlow[23] = maxXval[23]-150.;
   //
   hPMRefAhg->Fit("gaus","Q","",xlow[23],maxXval[23]+100.);
   TF1 *fun = hPMRefAhg->GetFunction("gaus");
   mean[23]  = (Float_t) (fun->GetParameter(1));
   sigma[23] = (Float_t) (fun->GetParameter(2));
   atLeastOneHisto=1; 
  }
  // ******** Low gain chain ********

  if(atLeastOneHisto==0){
    printf("\n WARNING! Empty LASER histos -> ending DA WITHOUT writing output\n\n");
    return -1;
  }
  
  FILE *fileShuttle=NULL;
  fileShuttle = fopen(LASDATA_FILE,"w");
  if(fileShuttle==NULL) {
    printf("Failed to open LASDATA file\n");
    return -1;
  }
  for(Int_t i=0; i<2*kNChannels; i++){
    fprintf(fileShuttle,"\t%d\t%d\t%f\t%f\n",detector[i],quad[i],mean[i], sigma[i]); 
  }
  //						       
  fclose(fileShuttle);
    
  /* report progress */
  daqDA_progressReport(90);
  //
  TFile *histofile = new TFile(LASHISTO_FILE,"RECREATE");
  histofile->cd();
  hList->Write();
  histofile->Close();
  
  /* store the result file on FES */
  // [1] File with mapping
  status = daqDA_FES_storeFile(MAPDATA_FILE, "MAPPING");
  if(status){
    printf("Failed to export file : %d\n",status);
    return -2;
  }
  //
  // [2] File with laser data
  status = daqDA_FES_storeFile(LASDATA_FILE, "LASERDATA");
  if(status){
    printf("Failed to export file : %d\n",status);
    return -3;
  }
  // [3] File with laser histos
  status = daqDA_FES_storeFile(LASHISTO_FILE, "LASERHISTOS");
  if(status){
    printf("Failed to export pedestal histos file to DAQ FES\n");
    return -4;
  }

  /* report progress */
  daqDA_progressReport(100);

  return status;
}
