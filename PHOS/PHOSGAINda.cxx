/*
contact: Boris.Polishchuk@cern.ch
link: see comments in the $ALICE_ROOT/PHOS/AliPHOSRcuDA1.cxx
reference run: /castor/cern.ch/alice/phos/2007/10/04/18/07000008249001.1000.root
run type: PHYSICS
DA type: MON 
number of events needed: 1000
input files: RCU0.data  RCU1.data  RCU2.data  RCU3.data
Output files: PHOS_Module2_Calib.root
Trigger types used: PHYSICS
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

#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliPHOSRcuDA1.h"
#include "AliPHOSRawFitterv0.h"
#include "AliCaloAltroMapping.h"
#include "AliCaloRawStreamV3.h"


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

  int status;
  
  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
  
  /* Retrieve ZS parameters from DAQ DB */
  const char* zsfile = "zs.txt";
  int failZS = daqDA_DB_getFile(zsfile, zsfile);
  
  Int_t offset,threshold;
  
  if(!failZS) {
    FILE *f = fopen(zsfile,"r");
    int scan = fscanf(f,"%d %d",&offset,&threshold);
  }
  
  /* Retrieve mapping files from DAQ DB */ 
  const char* mapFiles[4] = {"RCU0.data","RCU1.data","RCU2.data","RCU3.data"};
  
  for(Int_t iFile=0; iFile<4; iFile++) {
    int failed = daqDA_DB_getFile(mapFiles[iFile], mapFiles[iFile]);
    if(failed) { 
      printf("Cannot retrieve file %s from DAQ DB. Exit.\n",mapFiles[iFile]);
      return -1;
    }
  }
  
  /* Open mapping files */
  AliAltroMapping *mapping[4];
  TString path = "./";
  path += "RCU";
  TString path2;
  for(Int_t i = 0; i < 4; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    mapping[i] = new AliCaloAltroMapping(path2.Data());
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

  AliRawReader *rawReader = NULL;
  AliPHOSRcuDA1* dAs[5];

  for(Int_t iMod=0; iMod<5; iMod++) {
    dAs[iMod] = 0;
  }
  
  Float_t e[64][56][2];
  Float_t t[64][56][2];

  for(Int_t iX=0; iX<64; iX++) {
    for(Int_t iZ=0; iZ<56; iZ++) {
      for(Int_t iGain=0; iGain<2; iGain++) {
	e[iX][iZ][iGain] = 0.;
	t[iX][iZ][iGain] = 0.;
      }
    }
  }
  
  Int_t cellX    = -1;
  Int_t cellZ    = -1;
  Int_t nBunches =  0;
  Int_t sigStart, sigLength;
  Int_t caloFlag;

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
      
      rawReader = new AliRawReaderDate((void*)event);
      AliCaloRawStreamV3 stream(rawReader,"PHOS",mapping);
      AliPHOSRawFitterv0 fitter;
      fitter.SubtractPedestals(kTRUE); // assume that data is non-ZS
      
      if(!failZS) {
	fitter.SubtractPedestals(kFALSE);
	fitter.SetAmpOffset(offset);
	fitter.SetAmpThreshold(threshold);
      }
      
      while (stream.NextDDL()) {
	while (stream.NextChannel()) {

	  cellX    = stream.GetCellX();
	  cellZ    = stream.GetCellZ();
	  caloFlag = stream.GetCaloFlag();  // 0=LG, 1=HG, 2=TRU
	  
	  if(caloFlag!=0 && caloFlag!=1) continue; //TRU data!
	  
	  // In case of oscillating signals with ZS, 
	  //a channel can have several bunches.
	  
	  nBunches = 0;
	  while (stream.NextBunch()) {
	    nBunches++;
	    if (nBunches > 1) continue;
	    sigStart  = stream.GetStartTimeBin();
	    sigLength = stream.GetBunchLength();
	    fitter.SetChannelGeo(stream.GetModule(),cellX,cellZ,caloFlag);
	    fitter.Eval(stream.GetSignals(),sigStart,sigLength);
	  } // End of NextBunch()
	  
	  if (nBunches>1) continue;
	  
	  e[cellX][cellZ][caloFlag] = fitter.GetEnergy();
	  t[cellX][cellZ][caloFlag] = fitter.GetTime();
	}
	
	if(stream.GetModule()<0 || stream.GetModule()>4) continue;
	
	if(dAs[stream.GetModule()])
	  dAs[stream.GetModule()]->FillHistograms(e,t);
	else {
	  dAs[stream.GetModule()] = new AliPHOSRcuDA1(stream.GetModule(),-1);
	  dAs[stream.GetModule()]->FillHistograms(e,t);
	}
	
	for(Int_t iX=0; iX<64; iX++) {
	  for(Int_t iZ=0; iZ<56; iZ++) {
	    for(Int_t iGain=0; iGain<2; iGain++) {
	      e[iX][iZ][iGain] = 0.;
	      t[iX][iZ][iGain] = 0.;
	    }
	  }
	}

      }

//       da1.FillHistograms(e,t);
//     //da1.UpdateHistoFile();
      
	delete rawReader;     
	nevents_physics++;
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
  
  for(Int_t i = 0; i < 4; i++) delete mapping[i];  
  
  /* Be sure that all histograms are saved */
  
  char localfile[128];
  
  for(Int_t iMod=0; iMod<5; iMod++) {
    if(!dAs[iMod]) continue;
    
    dAs[iMod]->UpdateHistoFile();
    dAs[iMod]->SetWriteToFile(kFALSE);    
    
    /* Store output files to the File Exchange Server */
    sprintf(localfile,"PHOS_Module%d_Calib.root",iMod);
    daqDA_FES_storeFile(localfile,"AMPLITUDES");
  }
  
  return status;
}
