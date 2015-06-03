/*

This program reads the DAQ data files passed as argument using the monitoring library.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA to write mapping for ADC modules and VME scaler

Contact: Chiara.Oppedisano@to.infn.it
Link: 
Run Type: PHYSICS
DA Type: MON
Number of events needed:  
Input Files:  ZDCTDCHistoLimits.dat
Output Files: ZDCChMapping.dat ZDCTDCCalib.dat ZDCTDCHisto.root
Trigger Types Used: different trigger types are used

*/

#define MAPDATA_FILE  "ZDCChMapping.dat"
#define TDCDATA_FILE  "ZDCTDCCalib.dat"
#define TDCHISTO_FILE "ZDCTDCHisto.root"
#define HISTLIM_FILE  "ZDCTDCHistoLimits.dat"

#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>
#include <signal.h>

// DATE
#include <event.h>
#include <monitor.h>
#include <daqDA.h>

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

int main(int argc, char **argv) {

  // needed for streamer application
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()"); 
  // needed for Minuit plugin
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer",
                                        "Minuit",
                                        "TMinuitMinimizer",
                                        "Minuit",
                                        "TMinuitMinimizer(const char*)");

  TVirtualFitter::SetDefaultFitter("Minuit");
  
  /* log start of process */
  printf("\n ZDC MAPPING program started\n");  
  signal(SIGSEGV, SIG_DFL);

  /* check that we got some arguments = list of files */
  if(argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  char *monitor_table[] = { "ALL", "no", "PHY", "yes", "SOD", "all", NULL };
  int err = monitorDeclareTable(monitor_table);
  if(err){
    printf("monitorDeclareTable() failed: %s\n", monitorDecodeError(err));
    return -1;
  } 
  
  // *** initializations ***************************
  int status=0, nphys=0;
  int const kNModules = 9;
  int const kNChannels = 24;
  int const kNScChannels = 32;  
//  int const kZDCTDCGeo = 4;
  
  int itdc=0, iprevtdc=-1, ihittdc=0;
  float tdcData[6], tdcL0=-999.;	
  for(int ij=0; ij<6; ij++) tdcData[ij]=-999.;
  bool hasL0arrived=kFALSE;
  
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

  // Module type codes
  enum ZDCModules{kV965=1, kV830=2, kTRG=3, kTRGI=4, kPU=5, KV1290=6, kV775N=7}; 
  
  // Module type codes
  enum ZDCGeoAddr{kFirstADCGeo=0, kLastADCGeo=3, kADDADCGeo=5,
       kTDCFakeGeo=8, kZDCTDCGeo=4, kADDTDCGeo=6,
       kScalerGeo=16, kPUGeo=29, kTrigScales=30, kTrigHistory=31};
  
  // Signal codes for ZDC 
  // Same codes used in DAQ configuration file
  // To be changed ONLY IF this file is changed!!! 
  // **** DO NOT CHANGE THE FOLLOWING LINES!!! ****
  enum ZDCSignal{
       kNotConnected=0, kVoid=1,
       kZNAC=2, kZNA1=3, kZNA2=4, kZNA3=5, kZNA4=6,
       kZPAC=7, kZPA1=8, kZPA2=9, kZPA3=10, kZPA4=11,
       kZNCC=12, kZNC1=13, kZNC2=14, kZNC3=15, kZNC4=16,
       kZPCC=17, kZPC1=18, kZPC2=19, kZPC3=20, kZPC4=21,
       kZEM1=22, kZEM2=23,
       kZDCAMon=24, kZDCCMon=25,
       kZNACoot=26, kZNA1oot=27, kZNA2oot=28, kZNA3oot=29, kZNA4oot=30,
       kZPACoot=31, kZPA1oot=32, kZPA2oot=33, kZPA3oot=34, kZPA4oot=35,
       kZNCCoot=36, kZNC1oot=37, kZNC2oot=38, kZNC3oot=39, kZNC4oot=40,
       kZPCCoot=41, kZPC1oot=42, kZPC2oot=43, kZPC3oot=44, kZPC4oot=45,
       kZEM1oot=46, kZEM2oot=47,
       kZDCAMonoot=48, kZDCCMonoot=49,
       kL1MBI=50, kL1CNI=51, kL1SCI=52, kL1EMDI=53, kL0I=54, 
       kL1MBO=55, kL1CNO=56, kL1SCO=57, kL1EMDO=58, 
       kHMBCN=59, kHSCEMD=60,
       kZNACD=61, kZNA1D=62, kZNA2D=63, kZNA3D=64, kZNA4D=65,
       kZPACD=66, kZPA1D=67, kZPA2D=68, kZPA3D=69, kZPA4D=70,
       kZNCCD=71, kZNC1D=72, kZNC2D=73, kZNC3D=74, kZNC4D=75,
       kZPCCD=76, kZPC1D=77, kZPC2D=78, kZPC3D=79, kZPC4D=80,
       kZEM1D=81, kZEM2D=82,
       kZDCAMonD=83, kZDCCMonD=84,
       kZNAD=85, kZPAD=86, kZNCD=87, kZPCD=88, kZEMD=89,
       kZNA0D=90, kZPA0D=91, kZNC0D=92, kZPC0D=93, k1kHzD=94, kGate=95, kAD=96, kCD=97, 
       kAorCD=98, kAandCD=99, kZEMORD=100, kAorCorZEMORD=101, kAorCorZEMD=102, kAD0=103, kAD1=104, kAD2=105, 
       kAD3=106, kAD4=107, kAD5=108, kAD6=109, kAD7=110, kAD8=111, kAD9=112, kAD10=113, 
       kAD11=114, kAD12=115, kAD13=116, kAD14=117, kAD15=118, kAD0D=119, kAD1D=120, kAD2D=121,
       kAD3D=122, kAD4D=123, kAD5D=124, kAD6D=125, kAD7D=126, kAD8D=127, kAD9D=128, kAD10D=129,
       kAD11D=130, kAD12D=131, kAD13D=132, kAD14D=133, kAD15D=134, kL0=135,
       k1ZAC=136, k1ZED=137, k1ZMD=138, k1ZMB=139
       };
  
  // *** read histo limits from data file ***************************
  float hlimit[2][6];
  for(int k=0; k<6; k++){
     for(int j=0; j<2; j++) hlimit[j][k]=0;
  }
  int read = 0;
  read = daqDA_DB_getFile(HISTLIM_FILE, HISTLIM_FILE);
  if(read){
    printf("\t ERROR!!! ZDCTDCHistoLimits.dat file NOT FOUND in DAQ db!  0> Exiting....\n");
    return -11;
  }
  
  FILE *fileHistLim = NULL;
  fileHistLim = fopen(HISTLIM_FILE,"r");
  if(fileHistLim==NULL) {
    printf("\t ERROR!!! Can't open ZDCTDCHistoLimits.dat file!  0> Exiting....\n");
    return -12;
  }
  float readValues[2][6];
  for(int jj=0; jj<6; jj++){
     for(int ii=0; ii<2; ii++){
  	fscanf(fileHistLim,"%f",&readValues[ii][jj]);
        hlimit[ii][jj] = readValues[ii][jj];
     }
     //printf(" ch. %d limits x histo: %f %f\n",jj,hlimit[0][jj],hlimit[1][jj]);
  }
  fclose(fileHistLim);
  // *****************************************************************

  // *** Creating a container for the histos ***************************
  TObjArray *hList = new TObjArray(0);
  hList->SetOwner(kTRUE);
 
  TH1F * hTDC[6]={0x0,0x0,0x0,0x0,0x0,0x0};
  TString prefix = "TDC";
  for(int it=0; it<6; it++){
     TString det;
     if(it==0) det ="ZEM1";
     else if(it==1) det ="ZEM2";
     else if(it==2) det ="ZNC";
     else if(it==3) det ="ZPC";
     else if(it==4) det ="ZNA";
     else if(it==5) det ="ZPA";
     TString hname = prefix+det;
     float limit[2] = {0., 0.};
     limit[0] = hlimit[0][it];
     limit[1] = hlimit[1][it];
     //
     hTDC[it] = new TH1F(hname.Data(), hname.Data(), int(limit[1]-limit[0]), limit[0], limit[1]);
     hList->Add(hTDC[it]);
  }
  
  FILE *mapFile4Shuttle=NULL;

  /* read the data files */
  for(int n=1;n<argc;n++){
   
    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }
    
    /* declare monitoring program */
    status=monitorDeclareMp( __FILE__ );
    if (status!=0) {
      printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
      return -1;
    }
    monitorSetNowait();
    monitorSetNoWaitNetworkTimeout(1000);


    int iev = 0;
    //Bool_t sodRead = kFALSE;
    //while(!sodRead){
    
    /* loop on events (infinite) */
    for(;;) {
      
      if(nphys > 50000) break;
      
      struct eventHeaderStruct *event;
      eventTypeType eventT;
 
      /* check shutdown condition */
      if (daqDA_checkShutdown()) {break;}

      /* get next event */
      status=monitorGetEventDynamic((void **)&event);
      if(status==MON_ERR_EOF){
        printf ("End of File detected\n");
        break; /* end of monitoring file has been reached */
      }
      if(status!=0) {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        return -1;
      }

      /* retry if got no event */
      if(event==NULL) continue;
      
      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      reader->Reset();
      reader->Select("ZDC");
      //
      AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
        
	
      /* use event - here, just write event id to result file */
      eventT=event->eventType;
     
      /*  READ MAPPING FROM SOD EVENT*/ 
      if(eventT==START_OF_DATA){
	  	
	rawStreamZDC->SetSODReading(kTRUE);
	
	// --------------------------------------------------------
	// --- Writing ascii data file for the Shuttle preprocessor
        mapFile4Shuttle = fopen(MAPDATA_FILE,"w");
        if(mapFile4Shuttle==NULL) {
          printf("Failed to open mapFile4Shuttle file\n");
          return -1;
        }
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
	      else if(modType[iMod]==2){ // VME scaler mapping -------------
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
	     fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",is,adcMod[is],adcCh[is],sigCode[is],det[is],sec[is]);
    	     //printf("  Mapping DA -> %d ADC: mod %d ch %d, code %d det %d, sec %d\n",is,adcMod[is],adcCh[is],sigCode[is],det[is],sec[is]);
	  }
	  for(Int_t is=0; is<kNScChannels; is++){
	     fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",is,scMod[is],scCh[is],scSigCode[is],scDet[is],scSec[is]);
 	     //if(scMod[is]!=-1) printf("  Mapping DA -> %d Scaler: mod %d ch %d, code %d det %d, sec %d\n",is,scMod[is],scCh[is],scSigCode[is],scDet[is],scSec[is]);
	  }
	  for(Int_t is=0; is<kNScChannels; is++){
	     fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\n",is,tdcMod[is],tdcCh[is],tdcSigCode[is]);
 	     //if(tdcMod[is]!=-1) printf("  Mapping DA -> %d TDC: mod %d ch %d, code %d\n", is,tdcMod[is],tdcCh[is],tdcSigCode[is]);
	  }
	  for(Int_t is=0; is<kNModules; is++){
	     fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\n",
	     modGeo[is],modType[is],modNCh[is]);
	     //printf("  Mapping DA -> Module mapping: geo %d type %d #ch %d\n",modGeo[is],modType[is],modNCh[is]);
	  }
	  
        fclose(mapFile4Shuttle);
      }// SOD event
      /*  READING PHYSICS EVENTS*/ 
      else if(eventT==PHYSICS_EVENT){ 

	rawStreamZDC->SetSODReading(kTRUE);

  	// ----- Setting ch. mapping read from SOD -----
	// ADC
	for(int jk=0; jk<2*kNChannels; jk++){
	  rawStreamZDC->SetMapADCMod(jk, adcMod[jk]);
	  rawStreamZDC->SetMapADCCh(jk, adcCh[jk]);
	  rawStreamZDC->SetMapADCSig(jk, sigCode[jk]);
	  rawStreamZDC->SetMapDet(jk, det[jk]);
	  rawStreamZDC->SetMapTow(jk, sec[jk]);
	}
	// TDC
	for(int l=0; l<32; l++){
	   rawStreamZDC->SetTDCModFromMap(l, tdcMod[l]);
	   rawStreamZDC->SetTDCChFromMap(l, tdcCh[l]);
	   rawStreamZDC->SetTDCSignFromMap(l, tdcSigCode[l]);
//printf(" %d setting TDC ch. %d to signal %d\n",l, tdcCh[l], tdcSigCode[l]);
	}
	
  	while(rawStreamZDC->Next()){
         if(rawStreamZDC->GetADCModule()!=kZDCTDCGeo) continue; //skipping ADCs, scalers and trigger cards

          if(rawStreamZDC->GetADCModule()==kZDCTDCGeo && rawStreamZDC->IsZDCTDCDatum()==kTRUE){
             //
	     itdc = rawStreamZDC->GetChannel(); 
             if(itdc==iprevtdc) ihittdc++;
             else ihittdc=0;
             iprevtdc=itdc;
	     //
	     // NB -> THE TDC cabling has been varied in 12/2013 
	     // PLEASE CROSS-CHECK https://twiki.cern.ch/twiki/bin/view/ALICE/ZDCTDCCabl
	     int firstintch = 16;
	     int lastintch=23;
	     int tdcmapch=-1;
	     if((itdc>=firstintch && itdc<=lastintch) && ihittdc==0){
//printf(" **** TDC ch. %d cabled signal: %d\n",itdc, rawStreamZDC->GetTDCSignFromMap(itdc));
               int sigcode = rawStreamZDC->GetTDCSignFromMap(itdc);
	       //
	       if(sigcode==kZEM1D) tdcmapch=0;
	       else if(sigcode==kZEM2D) tdcmapch=1;
	       else if(sigcode==kZNCD) tdcmapch=2;
	       else if(sigcode==kZPCD) tdcmapch=3;
	       else if(sigcode==kZNAD) tdcmapch=4;
	       else if(sigcode==kZPAD) tdcmapch=5;
	       if(tdcmapch>=0) tdcData[tdcmapch] = 0.025*rawStreamZDC->GetZDCTDCDatum();
	       //
//if(tdcmapch>=0) printf("   *** TDC %d -> tdcData[%d] = %f  \n",itdc,  tdcmapch, tdcData[tdcmapch]);
	       //
	       else if(sigcode==kL0){
		  hasL0arrived=kTRUE;
	          tdcL0 = 0.025*rawStreamZDC->GetZDCTDCDatum();
	          //
//printf("    -> L0 = %f \n",tdcL0);
		  for(int ic=0; ic<6; ic++){
		    if(tdcData[ic]!=-999.){
//printf("   \t    reading tdcData[%d] = %f  \n",ic, tdcData[ic]);
		      hTDC[ic]->Fill(tdcData[ic]-tdcL0);
		      //printf(" ev.%d -> Filling histo%d:  %f ns\n",nphys,ic, tdcData[ic]-tdcL0);
		    }
	            // Resetting TDC values after L0 reading
  	            tdcData[ic]=-999.;
		  }
               }//L0
	     }//Loop on TDC
	  }
        }
	
        // Resetting TDC values after event reading
        if(hasL0arrived==kFALSE) for(int ic=0; ic<6; ic++) tdcData[ic]=-999.;
	
 	nphys++;
      
        delete rawStreamZDC;
        rawStreamZDC = 0x0;	
	delete reader;

      }//(if PHYSICS_EVENT) 
      else if(eventT==END_OF_RUN){
        printf("End Of Run detected\n");
        break;
      }
      
      iev++; 
      hasL0arrived==kFALSE;

      /* free resources */
      free(event);
 
    } // event loop    
      
  }
  printf(" \n # of physics events analyzed = %d\n\n", nphys);

  /* Analysis of the histograms */
  //
  FILE *fileShuttle;
  fileShuttle = fopen(TDCDATA_FILE,"w");
  //
  Float_t xUp=0., xLow=0., deltaX=0;
  Int_t binMax=0, nBinsx=0;
  Float_t mean[6]={0.,0.,0.,0.,0.,0.}, sigma[6]={0.,0.,0.,0.,0.,0.};
  TF1 *fitfun[6]={0x0,0x0,0x0,0x0,0x0,0x0};
  for(int k=0; k<6; k++){
    if(hTDC[k]->GetEntries()>0){
       binMax = hTDC[k]->GetMaximumBin();
       if(binMax<=1){
         printf("\n WARNING! Something wrong with det %d histo \n\n", k);
         continue;
       }
       // 
       xUp = (hTDC[k]->GetXaxis())->GetXmax();
       xLow = (hTDC[k]->GetXaxis())->GetXmin();
       deltaX = xUp-xLow;
       nBinsx = (hTDC[k]->GetXaxis())->GetNbins();
       //printf(" xMax = %f\n", xLow+binMax*deltaX/nBinsx);
       hTDC[k]->Fit("gaus","Q","",xLow+binMax*deltaX/nBinsx*0.75,xLow+binMax*deltaX/nBinsx*1.25);
       fitfun[k] = hTDC[k]->GetFunction("gaus");
       mean[k] = (Float_t) (fitfun[k]->GetParameter(1));
       sigma[k] = (Float_t) (fitfun[k]->GetParameter(2));
       //printf("\t Mean value from fit = %1.2f ns\n", mean[k]);
       //
       fprintf(fileShuttle,"\t%f\t%f\n",mean[k], sigma[k]);
     }
     else printf("  WARNING! TDC histo %d has no entries and can't be processed!\n",k);
  }
  //						       
  fclose(fileShuttle);
  
  TFile *histofile = new TFile(TDCHISTO_FILE,"RECREATE");
  histofile->cd();
  for(int k=0; k<6; k++) hTDC[k]->Write();						       
  histofile->Close();
  //
  for(Int_t j=0; j<6; j++) delete hTDC[j];
  
  /* store the result files on FES */
  // [1] File with mapping
  status = daqDA_FES_storeFile(MAPDATA_FILE, "MAPPING");
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }
  // [2] File with TDC data
  status = daqDA_FES_storeFile(TDCDATA_FILE, "TDCDATA");
  if(status){
    printf("Failed to export TDC data file to DAQ FES\n");
    return -1;
  }
  // [3] File with TDC histos
  status = daqDA_FES_storeFile(TDCHISTO_FILE, "TDCHISTOS");
  if(status){
    printf("Failed to export TDC histos file to DAQ FES\n");
    return -1;
  }

  return status;
}
