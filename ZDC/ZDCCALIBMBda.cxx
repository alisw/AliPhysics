/*

This program reads the DAQ data files passed as argument using the monitoring library.

It computes the average event size and populates local "./result.txt" file with the 
result.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalone CALIBRATION_MB runs

Contact: Chiara.Oppedisano@to.infn.it
Link: 
Run Type: CALIBRATION_MB_RUN
DA Type: LDC
Number of events needed: at least 10^4
Input Files: ZDCPedestal.dat
Output Files: ZDCMBCalib.root, ZDCChMapping.dat
Trigger Types Used: Standalone Trigger

*/
#define PEDDATA_FILE  "ZDCPedestal.dat"
#define MAPDATA_FILE  "ZDCChMapping.dat"
#define ENCALIBDATA_FILE   "ZDCEnergyCalib.dat"
#define MBCALIBDATA_FILE   "ZDCMBCalib.root"

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
  int const kNModules = 10;
  int const kNChannels = 24;
  int const kNScChannels = 32;
  Int_t kFirstADCGeo=0, kLastADCGeo=1;
      
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
  TH2F * hZDCvsZEM  = new TH2F("hZDCvsZEM","EZDCTot vs. EZEM",100,0.,8.,100,0.,800.);
  TH2F * hZDCCvsZEM = new TH2F("hZDCCvsZEM","EZDCC vs. EZEM",100,0.,8.,100,0.,400.);
  TH2F * hZDCAvsZEM = new TH2F("hZDCAvsZEM","EZDCA vs. EZEM",100,0.,8.,100,0.,400.);
  
  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
  }
  
  FILE *mapFile4Shuttle;

  // *** To analyze CALIBRATION_MB events a pedestal data file is needed!!!
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

  /* report progress */
  daqDA_progressReport(10);
  
  FILE *fileEnCal = fopen(ENCALIBDATA_FILE,"r");
  if(fileEnCal==NULL) {
    printf("\t ERROR!!! Can't open ZDCEnergyCalib.dat file!!!\n");
    return -1;
  }
  
  Float_t calibCoeff[6];
  for(Int_t irow=0; irow<6; irow++){
     fscanf(fileEnCal,"%f",&calibCoeff[irow]);
  }

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
      
      if(eventT==START_OF_DATA){
	  		
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
    	     //printf("  CalibMB DA -> %d ADC: mod %d ch %d, code %d det %d, sec %d\n",
	     //  is,adcMod[is],adcCh[is],sigCode[is],det[is],sec[is]);
	  }
	  for(Int_t is=0; is<kNScChannels; is++){
	     fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",
	       is,scMod[is],scCh[is],scSigCode[is],scDet[is],scSec[is]);
 	     //printf("  CalibMB DA -> %d Scaler: mod %d ch %d, code %d det %d, sec %d\n",
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
	     //printf("  CalibMB DA -> Module mapping: geo %d type %d #ch %d\n",
	     //  modGeo[is],modType[is],modNCh[is]);
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
	 if((message & 0xf0) == 0x60){ // DEDICATED CALIBRATION MB RUN
	    //printf("\t CALIBRATION_MB_RUN raw data found\n");
	    continue;
	 }
	 else{
	    //printf("\t NO CALIBRAION_MB_RUN raw data found\n");
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
      Float_t rawZNC=0., rawZPC=0., rawZNA=0., rawZPA=0., rawZEM1=0., rawZEM2=0.;
      Float_t corrZNC=0., corrZPC=0., corrZNA=0., corrZPA=0., corrZEM1=0.,corrZEM2=0.;
      Float_t calZNC=0., calZPC=0., calZNA=0., calZPA=0., calZEM1=0., calZEM2=0.;
      Float_t calZDCTot=0., calZDCC=0., calZDCA=0., calZEM=0.;
      //
      while(rawStreamZDC->Next()){
	Int_t detector = rawStreamZDC->GetSector(0);
	Int_t quad = rawStreamZDC->GetSector(1);
  	Int_t index=-1;
        
	if( (rawStreamZDC->IsADCDataWord()) && !(rawStreamZDC->IsUnderflow())
	     && !(rawStreamZDC->IsOverflow()) && (detector!=-1)
	     && ((rawStreamZDC->GetADCGain()) == 0) && // Selecting HIGH RES ch.s
             (rawStreamZDC->GetADCModule()>=kFirstADCGeo) && (rawStreamZDC->GetADCModule()<=kLastADCGeo)){

	    //printf("  IsADCWord %d, IsUnderflow %d, IsOverflow %d\n",
	    //  rawStreamZDC->IsADCDataWord(),rawStreamZDC->IsUnderflow(),rawStreamZDC->IsOverflow());
	  
 	    if(quad!=5){ // Physics signals
    	      if(detector==1)      index = quad;   // *** ZNC
	      else if(detector==2) index = quad+5; // *** ZPC
	      else if(detector==3) index = quad+9; // *** ZEM
	      else if(detector==4) index = quad+12;// *** ZNA
	      else if(detector==5) index = quad+17;// *** ZPA
	    }
	    else continue;
	    //
	    if(index==-1) printf("ERROR in ZDCCALIBMBda.cxx -> det %d quad %d index %d ADC %d\n", 
	      detector, quad, index, rawStreamZDC->GetADCValue());
	    
	    // Mean pedestal subtraction 
	    Float_t Pedestal = MeanPedhg[index];
	    // Pedestal subtraction from correlation with out-of-time signals
	    //Float_t Pedestal = CorrCoeff0[index]+CorrCoeff1[index]*MeanPedOOT[index];
            //
	    if(index!=-1 && quad!=5){ 
	      //
	      if(detector==1){
	       rawZNC  = (Float_t) rawStreamZDC->GetADCValue();
	       corrZNC = rawZNC - Pedestal;
	      }
	      else if(detector==2){
	       rawZPC  = (Float_t) rawStreamZDC->GetADCValue();
	       corrZPC = rawZPC - Pedestal;
	      }
	      else if(detector==3){
	       if(quad==1){
	         rawZEM1  = (Float_t) rawStreamZDC->GetADCValue();
	         corrZEM1 = (rawStreamZDC->GetADCValue()) - Pedestal;
	       }
	       else{
	         rawZEM2  = (Float_t) rawStreamZDC->GetADCValue();
	         corrZEM2 = (rawStreamZDC->GetADCValue()) - Pedestal;
	       }
	      }
	      else if(detector==4){
	       rawZNA  = (Float_t) rawStreamZDC->GetADCValue();
	       corrZNA = rawZNA - Pedestal;
	      }
	      else if(detector==5){
	       rawZPA  = (Float_t) rawStreamZDC->GetADCValue();
	       corrZPA = rawZPA - Pedestal;
	      }
	    }
	}//IsADCDataWord()
	
       } // Next()
       
       calZNC = calibCoeff[0]*corrZNC;
       calZPC = calibCoeff[1]*corrZPC;
       calZNA = calibCoeff[2]*corrZNA;
       calZPA = calibCoeff[3]*corrZPA;
       calZEM1 = calibCoeff[4]*corrZEM1;
       calZEM2 = calibCoeff[4]*corrZEM2;
       calZDCTot = calZNC+calZPC+calZPA+calZNA;
       calZDCC = calZNC+calZPC;
       calZDCA = calZNA+calZPA;
       calZEM = calZEM1+calZEM2;
       //
       hZDCvsZEM ->Fill(calZDCTot/1000, calZEM/1000);
       hZDCCvsZEM->Fill(calZDCC/1000, calZEM/1000);
       hZDCAvsZEM->Fill(calZDCA/1000, calZEM/1000);

       nevents_physics++;
    }//(if PHYSICS_EVENT)
    
    nevents_total++;

   }

   /* free resources */
   free(event);
  }
    
  //
  TFile *histofile = new TFile(MBCALIBDATA_FILE,"RECREATE");
  histofile->cd();
  hZDCvsZEM ->Write();
  hZDCCvsZEM->Write();
  hZDCAvsZEM->Write();
  histofile->Close();
  
  if(hZDCvsZEM)  delete hZDCvsZEM ;
  if(hZDCCvsZEM) delete hZDCCvsZEM;
  if(hZDCAvsZEM) delete hZDCAvsZEM;
  
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
  status = daqDA_FES_storeFile(MBCALIBDATA_FILE, "MBCALIB");
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }
 
  /* report progress */
  daqDA_progressReport(100);

  return status;
}
