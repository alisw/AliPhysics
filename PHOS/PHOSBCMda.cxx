/*
contact: Boris.Polishchuk@cern.ch
link: http://aliceinfo.cern.ch/static/phpBB3/viewtopic.php?f=4&t=17
reference run: /castor/cern.ch/alice/phos/2007/10/02/13/07000008232001.10.root
run type: LED
DA type: MON
number of events needed: 1000
number of events needed: 1000
input files: Mod0RCU0.data Mod0RCU1.data Mod0RCU2.data Mod0RCU3.data Mod1RCU0.data Mod1RCU1.data Mod1RCU2.data Mod1RCU3.data Mod2RCU0.data Mod2RCU1.data Mod2RCU2.data Mod2RCU3.data Mod3RCU0.data Mod3RCU1.data Mod3RCU2.data Mod3RCU3.data Mod4RCU0.data Mod4RCU1.data Mod4RCU2.data Mod4RCU3.data 
Output files: PHOS_Module2_BCM.root
Trigger types used: CALIBRATION_EVENT
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
#include "AliPHOSDA2.h"
#include "AliPHOSRawFitterv1.h"
#include "AliCaloAltroMapping.h"
#include "AliCaloRawStreamV3.h"
#include "AliLog.h"


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

  AliLog::SetGlobalDebugLevel(0) ;
  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  
  int status;
  
  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
  
  short offset, threshold;

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

  
  /* log start of process */
  printf("DA2 (bad channels search) started.\n");  


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  AliRawReader *rawReader = NULL;

  AliPHOSDA2* da2 = new AliPHOSDA2(2); // DA2 ("Checking for bad channels") for module2
  
  Float_t q[64][56][2];

  Int_t cellX    = -1;
  Int_t cellZ    = -1;
  Int_t nBunches =  0;
  Int_t nFired   = -1;
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
    
    if (eventT==PHYSICS_EVENT || eventT==CALIBRATION_EVENT) {
      
      for(Int_t iX=0; iX<64; iX++) {
	for(Int_t iZ=0; iZ<56; iZ++) {
	  for(Int_t iGain=0; iGain<2; iGain++) {
	    q[iX][iZ][iGain] = 0.;
	  }
	}
      }

      nFired = 0;

      rawReader = new AliRawReaderDate((void*)event);
      AliCaloRawStreamV3 stream(rawReader,"PHOS",mapping);
      AliPHOSRawFitterv1 fitter;
      fitter.SubtractPedestals(kTRUE); // assume that data is non-ZS
      
      while (stream.NextDDL()) {
	while (stream.NextChannel()) {

	  /* Retrieve ZS parameters from data*/
	  short value = stream.GetAltroCFG1();
	  bool ZeroSuppressionEnabled = (value >> 15) & 0x1;
	  bool AutomaticBaselineSubtraction = (value >> 14) & 0x1;
	  if(ZeroSuppressionEnabled) {
	    offset = (value >> 10) & 0xf;
	    threshold = value & 0x3ff;
	    fitter.SubtractPedestals(kFALSE);
	    fitter.SetAmpOffset(offset);
	    fitter.SetAmpThreshold(threshold);
	  }
	  
	  cellX    = stream.GetCellX();
	  cellZ    = stream.GetCellZ();
	  caloFlag = stream.GetCaloFlag();  // 0=LG, 1=HG, 2=TRU
	  
	  if(caloFlag!=0 && caloFlag!=1) continue; //TRU data!
	  
	  // In case of oscillating signals with ZS, a channel can have several bunches
	  nBunches = 0;
	  while (stream.NextBunch()) {
	    nBunches++;
	    sigStart  = stream.GetStartTimeBin();
	    sigLength = stream.GetBunchLength();
	    fitter.SetChannelGeo(stream.GetModule(),cellX,cellZ,caloFlag);
	    fitter.Eval(stream.GetSignals(),sigStart,sigLength);
	    q[cellX][cellZ][caloFlag] = fitter.GetSignalQuality();
	    printf("q[%d][%d][%d] = %.3f\n",cellX,cellZ,caloFlag,q[cellX][cellZ][caloFlag]);
	  } // End of NextBunch()
	  
	  if(caloFlag==1 && fitter.GetEnergy()>40)
	    nFired++;
	}
      }
      
      da2->FillQualityHistograms(q);
      da2->FillFiredCellsHistogram(nFired);
      //da1.UpdateHistoFile();
      
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
  
  for(Int_t i = 0; i < 20; i++) delete mapping[i];  

  /* Be sure that all histograms are saved */
  delete da2;
  
  /* Store output files to the File Exchange Server */
  daqDA_FES_storeFile("PHOS_Module2_BCM.root","BAD_CHANNELS");

  return status;
}
