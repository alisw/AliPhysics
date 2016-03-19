/*
contact: Boris.Polishchuk@cern.ch
reference run: /alice/data/2010/LHC10a_PHOS/000112189/raw/10000112189003.10.root
run type: PEDESTAL
DA type: MON 
number of events needed: 200
input files: Mod0RCU0.data Mod0RCU1.data Mod0RCU2.data Mod0RCU3.data Mod1RCU0.data Mod1RCU1.data Mod1RCU2.data Mod1RCU3.data Mod2RCU0.data Mod2RCU1.data Mod2RCU2.data Mod2RCU3.data Mod3RCU0.data Mod3RCU1.data Mod3RCU2.data Mod3RCU3.data Mod4RCU0.data Mod4RCU1.data Mod4RCU2.data Mod4RCU3.data
Output files: PHOS_PED.root
Trigger types used: 
*/

#include "event.h"
#include "monitor.h"

extern "C" {
#include "daqDA.h"
}

#include <stdio.h>
#include <stdlib.h>

#include <TSystem.h>
#include <TROOT.h>
#include <TPluginManager.h>

#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliCaloAltroMapping.h"
#include "AliCaloRawStreamV3.h"
#include "AliLog.h"

int main(int argc, char **argv) 
{
  // Read raw data, decode it to samples,
  // calculate pedestals from presamples, 
  // evaluate the signal amplitude as a maximum sample, 
  // and fill histograms with pedestals and amplitudes

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");
  
  AliLog::SetGlobalDebugLevel(0) ;
  AliLog::SetGlobalLogLevel(AliLog::kFatal);

  int status;
  
  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  /* Retrieve mapping files from DAQ DB */ 
  const char* mapFiles[20] = {
    "Mod0RCU0.data",
    "Mod0RCU1.data",
    "Mod0RCU2.data",
    "Mod0RCU3.data",
    "Mod1RCU0.data",
    "Mod1RCU1.data",
    "Mod1RCU2.data",
    "Mod1RCU3.data",
    "Mod2RCU0.data",
    "Mod2RCU1.data",
    "Mod2RCU2.data",
    "Mod2RCU3.data",
    "Mod3RCU0.data",
    "Mod3RCU1.data",
    "Mod3RCU2.data",
    "Mod3RCU3.data",
    "Mod4RCU0.data",
    "Mod4RCU1.data",
    "Mod4RCU2.data",
    "Mod4RCU3.data"
  };

  for(Int_t iFile=0; iFile<20; iFile++) {
    int failed = daqDA_DB_getFile(mapFiles[iFile], mapFiles[iFile]);
    if(failed) { 
      printf("Cannot retrieve file %s from DAQ DB. Exit.\n",mapFiles[iFile]);
      return -1;
    }
  }
  
  /* Open mapping files */
  AliAltroMapping *mapping[20];
  TString path = "./";

  path += "Mod";
  TString path2;
  TString path3;
  Int_t iMap = 0;

  for(Int_t iMod = 0; iMod < 5; iMod++) {
    path2 = path;
    path2 += iMod;
    path2 += "RCU";

    for(Int_t iRCU=0; iRCU<4; iRCU++) {
      path3 = path2;
      path3 += iRCU;
      path3 += ".data";
      mapping[iMap] = new AliCaloAltroMapping(path3.Data());
      iMap++;
    }
  }
  
  /* define data source : this is argument 1 */  
  status=monitorSetDataSource( argv[1] );
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
  
  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);
  
  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  AliRawReader * reader = NULL;
  AliCaloRawStreamV3* stream = NULL;
  
  TString baseNamePed ="hPed";
  TString baseTitlePed="Ped in cell (";
  const char* sgain[3]={"LG","HG", "TRU"};

  const Int_t caloFlagMax=3,modMax=5,cellXMax=64,cellZMax=56;
  TH1F *hPed[5][3][64][56] = {};

  TH1F *hPedHiMean1m1 = new TH1F("hPedHiMean1m1","Mean pedestals in module 1, high gain" ,100,0.,100.);
  TH1F *hPedHiRMS1m1  = new TH1F("hPedHiRMS1m1" ,"RMS pedestals in module 1, high gain"  ,100,0.,50.);
  TH1F *hPedLoMean1m1 = new TH1F("hPedLoMean1m1","Mean pedestals in module 1, low gain"  ,100,0.,100.);
  TH1F *hPedLoRMS1m1  = new TH1F("hPedLoRMS1m1" ,"RMS pedestals in module 1, low gain"   ,100,0.,50.);
  TH1F *hPedTRUMean1m1 = new TH1F("hPedTRUMean1m1","Mean pedestals in module 1, TRU"     ,1000,0.,1000.);
  TH1F *hPedTRURMS1m1  = new TH1F("hPedTRURMS1m1" ,"RMS pedestals in module 1, TRU"      ,100,0.,50.);

  TH1F *hPedHiMean1m2 = new TH1F("hPedHiMean1m2","Mean pedestals in module 2, high gain" ,100,0.,100.);
  TH1F *hPedHiRMS1m2  = new TH1F("hPedHiRMS1m2" ,"RMS pedestals in module 2, high gain"  ,100,0.,50.);
  TH1F *hPedLoMean1m2 = new TH1F("hPedLoMean1m2","Mean pedestals in module 2, low gain"  ,100,0.,100.);
  TH1F *hPedLoRMS1m2  = new TH1F("hPedLoRMS1m2" ,"RMS pedestals in module 2, low gain"   ,100,0.,50.);
  TH1F *hPedTRUMean1m2 = new TH1F("hPedTRUMean1m2","Mean pedestals in module 2, TRU"     ,1000,0.,1000.);
  TH1F *hPedTRURMS1m2  = new TH1F("hPedTRURMS1m2" ,"RMS pedestals in module 2, TRU"      ,100,0.,50.);

  TH1F *hPedHiMean1m3 = new TH1F("hPedHiMean1m3","Mean pedestals in module 3, high gain" ,100,0.,100.);
  TH1F *hPedHiRMS1m3  = new TH1F("hPedHiRMS1m3" ,"RMS pedestals in module 3, high gain"  ,100,0.,50.);
  TH1F *hPedLoMean1m3 = new TH1F("hPedLoMean1m3","Mean pedestals in module 3, low gain"  ,100,0.,100.);
  TH1F *hPedLoRMS1m3  = new TH1F("hPedLoRMS1m3" ,"RMS pedestals in module 3, low gain"   ,100,0.,50.);
  TH1F *hPedTRUMean1m3 = new TH1F("hPedTRUMean1m3","Mean pedestals in module 3, TRU"     ,1000,0.,1000.);
  TH1F *hPedTRURMS1m3  = new TH1F("hPedTRURMS1m3" ,"RMS pedestals in module 3, TRU"      ,100,0.,50.);

  TH1F *hPedHiMean1m4 = new TH1F("hPedHiMean1m4","Mean pedestals in module 4, high gain" ,100,0.,100.);
  TH1F *hPedHiRMS1m4  = new TH1F("hPedHiRMS1m4" ,"RMS pedestals in module 4, high gain"  ,100,0.,50.);
  TH1F *hPedLoMean1m4 = new TH1F("hPedLoMean1m4","Mean pedestals in module 4, low gain"  ,100,0.,100.);
  TH1F *hPedLoRMS1m4  = new TH1F("hPedLoRMS1m4" ,"RMS pedestals in module 4, low gain"   ,100,0.,50.);
  TH1F *hPedTRUMean1m4 = new TH1F("hPedTRUMean1m4","Mean pedestals in module 4, TRU"     ,1000,0.,1000.);
  TH1F *hPedTRURMS1m4  = new TH1F("hPedTRURMS1m4" ,"RMS pedestals in module 4, TRU"      ,100,0.,50.);

  hPedHiMean1m1->Sumw2();
  hPedHiRMS1m1 ->Sumw2();
  hPedLoMean1m1->Sumw2();
  hPedLoRMS1m1 ->Sumw2();
  hPedTRUMean1m1->Sumw2();
  hPedTRURMS1m1 ->Sumw2();
  hPedHiMean1m2->Sumw2();
  hPedHiRMS1m2 ->Sumw2();
  hPedLoMean1m2->Sumw2();
  hPedLoRMS1m2 ->Sumw2();
  hPedTRUMean1m2->Sumw2();
  hPedTRURMS1m2 ->Sumw2();
  hPedHiMean1m3->Sumw2();
  hPedHiRMS1m3 ->Sumw2();
  hPedLoMean1m3->Sumw2();
  hPedLoRMS1m3 ->Sumw2();
  hPedTRUMean1m3->Sumw2();
  hPedTRURMS1m3 ->Sumw2();
  hPedHiMean1m4->Sumw2();
  hPedHiRMS1m4 ->Sumw2();
  hPedLoMean1m4->Sumw2();
  hPedLoRMS1m4 ->Sumw2();
  hPedTRUMean1m4->Sumw2();
  hPedTRURMS1m4 ->Sumw2();

  TH2F *hPedHiMeanm1  = new TH2F("hPedHiMeanm1","Mean pedestals in module 1, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiRMSm1   = new TH2F("hPedHiRMSm1" ,"R.M.S. of pedestals in module 1, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiNumm1   = new TH2F("hPedHiNumm1" ,"Number of pedestals in module 1, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoMeanm1  = new TH2F("hPedLoMeanm1","Mean pedestals in module 1, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoRMSm1   = new TH2F("hPedLoRMSm1" ,"R.M.S. of pedestals in module 1, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoNumm1   = new TH2F("hPedLoNumm1" ,"Number of pedestals in module 1, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);

  TH2F *hPedHiMeanm2  = new TH2F("hPedHiMeanm2","Mean pedestals in module 2, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiRMSm2   = new TH2F("hPedHiRMSm2" ,"R.M.S. of pedestals in module 2, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiNumm2   = new TH2F("hPedHiNumm2" ,"Number of pedestals in module 2, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoMeanm2  = new TH2F("hPedLoMeanm2","Mean pedestals in module 2, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoRMSm2   = new TH2F("hPedLoRMSm2" ,"R.M.S. of pedestals in module 2, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoNumm2   = new TH2F("hPedLoNumm2" ,"Number of pedestals in module 2, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);

  TH2F *hPedHiMeanm3  = new TH2F("hPedHiMeanm3","Mean pedestals in module 3, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiRMSm3   = new TH2F("hPedHiRMSm3" ,"R.M.S. of pedestals in module 3, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiNumm3   = new TH2F("hPedHiNumm3" ,"Number of pedestals in module 3, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoMeanm3  = new TH2F("hPedLoMeanm3","Mean pedestals in module 3, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoRMSm3   = new TH2F("hPedLoRMSm3" ,"R.M.S. of pedestals in module 3, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoNumm3   = new TH2F("hPedLoNumm3" ,"Number of pedestals in module 3, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);

  TH2F *hPedHiMeanm4  = new TH2F("hPedHiMeanm4","Mean pedestals in module 4, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiRMSm4   = new TH2F("hPedHiRMSm4" ,"R.M.S. of pedestals in module 4, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiNumm4   = new TH2F("hPedHiNumm4" ,"Number of pedestals in module 4, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoMeanm4  = new TH2F("hPedLoMeanm4","Mean pedestals in module 4, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoRMSm4   = new TH2F("hPedLoRMSm4" ,"R.M.S. of pedestals in module 4, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoNumm4   = new TH2F("hPedLoNumm4" ,"Number of pedestals in module 4, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);

  TH1I *hNBunches = new TH1I("hNBunches","Number of bunches",10,0,10);
  TH2I *hHWaddr   = new TH2I("hHWaddr","DDL is vs HW addr",216,0,216,4096,0,4096);
  TH1I *hModule   = new TH1I("hModule" ,"Module number", 5,0.,5);

  Int_t runNum=0;
  Int_t module,cellX,cellZ,caloFlag;

  /* main loop (infinite) */
  for(;;) {
    struct eventHeaderStruct *event;
    eventTypeType eventT;
    
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

    /* use event - here, just write event id to result file */
    eventT=event->eventType;
    
    if (eventT==PHYSICS_EVENT) {
      
      reader = new AliRawReaderDate((void*)event);
      stream =  new AliCaloRawStreamV3(reader,"PHOS",mapping);

      runNum = reader->GetRunNumber();
      while (stream->NextDDL()) {
	while (stream->NextChannel()) {
	  module   = stream->GetModule();
	  cellX    = stream->GetCellX();
	  cellZ    = stream->GetCellZ();
	  caloFlag = stream->GetCaloFlag();
	  if (caloFlag!=0 && caloFlag!=1) continue;

	  hHWaddr->Fill(stream->GetDDLNumber(),stream->GetHWAddress());
	  hModule->Fill(module);
	  if (!hPed[module][caloFlag][cellX][cellZ]) {
	    TString name  = baseNamePed;
	    TString title = baseTitlePed;
	    name +="_g"; name +=caloFlag;
	    name +="_m"; name +=module;
	    name +="_x"; name +=cellX;
	    name +="_z"; name +=cellZ;
	    
	    title +=module; title +=",";
	    title +=cellX; title +=",";
	    title +=cellZ; title +="), ";
	    title +=sgain[caloFlag];
	    
	    hPed[module][caloFlag][cellX][cellZ] = new TH1F(name,title,100,0.,100.);
	    hPed[module][caloFlag][cellX][cellZ]->Sumw2();
	    hPed[module][caloFlag][cellX][cellZ]->SetMarkerStyle(20);
	    hPed[module][caloFlag][cellX][cellZ]->SetOption("eph");
	  }
	    
	  Int_t nBunches = 0;
	  while (stream->NextBunch()) {
	    nBunches++;
	    const UShort_t *sig = stream->GetSignals();
	    Int_t sigLength = stream->GetBunchLength();
	    for (Int_t i = 0; i < sigLength; i++) {
	      hPed[module][caloFlag][cellX][cellZ]->Fill(sig[i]);
	    }
	  }
	  hNBunches->Fill(nBunches);
	} // end of NextChannel()
      } // end of NextDDL()
      
      delete reader;     
      delete stream;
      nevents_physics++;
    } // end of if(eventT==PHYSICS_EVENT)
    
    nevents_total++;
    
    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }

  } // end of inf. loop over events
  
  for(Int_t i = 0; i < 20; i++) delete mapping[i];

  // Fill 2-dim histograms for mean, rms and n pedestals
  
  for (Int_t mod=1; mod<=4; mod++) {
    for (Int_t caloFlag=0; caloFlag<2; caloFlag++) {
      for (Int_t cellX=0; cellX<cellXMax; cellX++) {
	for (Int_t cellZ=0; cellZ<cellZMax; cellZ++) {
	  if (hPed[mod][caloFlag][cellX][cellZ] != 0) {
	    if      (caloFlag == 0) {
	      if (mod==1) {
		hPedLoMean1m1->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMS1m1 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedLoMeanm1 ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMSm1  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedLoNumm1  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetEntries());
	      }
	      else if (mod==2) {
		hPedLoMean1m2->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMS1m2 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedLoMeanm2 ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMSm2  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedLoNumm2  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetEntries());
	      }
	      else if (mod==3) {
		hPedLoMean1m3->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMS1m3 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedLoMeanm3 ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMSm3  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedLoNumm3  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetEntries());
	      }
	      else if (mod==4) {
		hPedLoMean1m4->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMS1m4 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedLoMeanm4 ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMSm4  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedLoNumm4  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetEntries());
	      }
	    }
	    else if (caloFlag == 1) {
	      if (mod==1) {
		hPedHiMean1m1->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMS1m1 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS() );
		hPedHiMeanm1 ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMSm1  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedHiNumm1  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetEntries());
	      }
	      if (mod==2) {
		hPedHiMean1m2->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMS1m2 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS() );
		hPedHiMeanm2 ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMSm2  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedHiNumm2  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetEntries());
	      }
	      if (mod==3) {
		hPedHiMean1m3->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMS1m3 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS() );
		hPedHiMeanm3 ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMSm3  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedHiNumm3  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetEntries());
	      }
	      if (mod==4) {
		hPedHiMean1m4->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMS1m4 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedHiMeanm4 ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMSm4  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetRMS());
		hPedHiNumm4  ->Fill( cellX, cellZ, hPed[mod][caloFlag][cellX][cellZ]->GetEntries());
	      }
	    }
	    else if (caloFlag == 2) {
	      if (mod==1) {
		hPedTRUMean1m1->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedTRURMS1m1 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS() );
	      }
	      if (mod==2) {
		hPedTRUMean1m2->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedTRURMS1m2 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS() );
	      }
	      if (mod==3) {
		hPedTRUMean1m3->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedTRURMS1m3 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS() );
	      }
	      if (mod==4) {
		hPedTRUMean1m4->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetMean());
		hPedTRURMS1m4 ->Fill( hPed[mod][caloFlag][cellX][cellZ]->GetRMS() );
	      }
	    }
	  }
	}
      }
    }
  }
  
  // Write existing histograms to a root file
  
  TString fileName = "PHOS_PED.root";
  TFile *file = new TFile(fileName,"RECREATE");
  
  for (Int_t caloFlag=0; caloFlag<caloFlagMax; caloFlag++) {
    for (Int_t mod=0; mod<modMax; mod++) {
      for (Int_t cellX=0; cellX<cellXMax; cellX++) {
	for (Int_t cellZ=0; cellZ<cellZMax; cellZ++) {
	  if (hPed[mod][caloFlag][cellX][cellZ] != 0)
	    hPed[mod][caloFlag][cellX][cellZ]->Write();
	}
      }
    }
  }
  
  hPedHiMean1m1->Write();
  hPedHiRMS1m1 ->Write();
  hPedLoMean1m1->Write();
  hPedLoRMS1m1 ->Write();
  hPedHiMeanm1 ->Write();
  hPedHiRMSm1  ->Write();
  hPedHiNumm1  ->Write();
  hPedLoMeanm1 ->Write();
  hPedLoRMSm1  ->Write();
  hPedLoNumm1  ->Write();
  hPedTRUMean1m1->Write();
  hPedTRURMS1m1 ->Write();

  hPedHiMean1m2->Write();
  hPedHiRMS1m2 ->Write();
  hPedLoMean1m2->Write();
  hPedLoRMS1m2 ->Write();
  hPedHiMeanm2 ->Write();
  hPedHiRMSm2  ->Write();
  hPedHiNumm2  ->Write();
  hPedLoMeanm2 ->Write();
  hPedLoRMSm2  ->Write();
  hPedLoNumm2  ->Write();
  hPedTRUMean1m2->Write();
  hPedTRURMS1m2 ->Write();
  
  hPedHiMean1m3->Write();
  hPedHiRMS1m3 ->Write();
  hPedLoMean1m3->Write();
  hPedLoRMS1m3 ->Write();
  hPedHiMeanm3 ->Write();
  hPedHiRMSm3  ->Write();
  hPedHiNumm3  ->Write();
  hPedLoMeanm3 ->Write();
  hPedLoRMSm3  ->Write();
  hPedLoNumm3  ->Write();
  hPedTRUMean1m3->Write();
  hPedTRURMS1m3 ->Write();

  hPedHiMean1m4->Write();
  hPedHiRMS1m4 ->Write();
  hPedLoMean1m4->Write();
  hPedLoRMS1m4 ->Write();
  hPedHiMeanm4 ->Write();
  hPedHiRMSm4  ->Write();
  hPedHiNumm4  ->Write();
  hPedLoMeanm4 ->Write();
  hPedLoRMSm4  ->Write();
  hPedLoNumm4  ->Write();
  hPedTRUMean1m4->Write();
  hPedTRURMS1m4 ->Write();

  hNBunches  ->Write();
  hHWaddr    ->Write();
  hModule    ->Write();
  
  file->Close();
  
  /* Store output files to the File Exchange Server */
  daqDA_FES_storeFile(fileName.Data(),"PED");
  
  printf("%d physics events of %d total processed.\n",nevents_physics,nevents_total);
  return status;
}
