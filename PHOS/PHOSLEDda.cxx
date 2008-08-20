/*
contact: Boris.Polishchuk@cern.ch
link: see comments in the $ALICE_ROOT/PHOS/AliPHOSRcuDA1.cxx
reference run: /castor/cern.ch/alice/phos/2007/10/04/18/07000008249001.1000.root
run type: STANDALONE
DA type: MON 
number of events needed: 1000
input files: RCU0.data  RCU1.data  RCU2.data  RCU3.data
Output files: PHOS_Module2_LED.root
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
#include "AliPHOSRcuDA1.h"
#include "AliPHOSRawDecoder.h"
#include "AliCaloAltroMapping.h"


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

  Int_t fMod = 2; // module 2!
  AliRawReader *rawReader = NULL;

  AliPHOSRcuDA1 da1(fMod,-1,0); // DA1 for module2, no input/output file
  
  Float_t e[64][56][2];
  Float_t t[64][56][2];

  Int_t gain = -1;
  Int_t X = -1;
  Int_t Z = -1;

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
      
      for(Int_t iX=0; iX<64; iX++) {
	for(Int_t iZ=0; iZ<56; iZ++) {
	  for(Int_t iGain=0; iGain<2; iGain++) {
	    e[iX][iZ][iGain] = 0.;
	    t[iX][iZ][iGain] = 0.;
	  }
	}
      }

      rawReader = new AliRawReaderDate((void*)event);
//       AliPHOSRawDecoderv1 dc(rawReader,mapping);
      AliPHOSRawDecoder dc(rawReader,mapping);
      dc.SubtractPedestals(kTRUE);
      
      while(dc.NextDigit()) {

	X = dc.GetRow() - 1;
	Z = dc.GetColumn() - 1;

	if(dc.IsLowGain()) gain = 0;
	else
	  gain = 1;
	
	e[X][Z][gain] = dc.GetEnergy();
	t[X][Z][gain] = dc.GetTime();
	
      }
      
      da1.FillHistograms(e,t);
      
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
  sprintf(localfile,"PHOS_Module%d_LED.root",fMod);
  TFile* f = new TFile(localfile,"recreate");
  
  const TH2F* h2=0;
  const TH1F* h1=0;

  Int_t nGood=0;    // >10 entries in peak
  Int_t nMax=-111;  // max. number of entries in peak
  Int_t iXmax=-1;
  Int_t iZmax=-1;

  for(Int_t iX=0; iX<64; iX++) {
    for(Int_t iZ=0; iZ<56; iZ++) {
      
      h1 = da1.GetHgLgRatioHistogram(iX,iZ); // High Gain/Low Gain ratio
      if(h1) {
	if(h1->GetMaximum()>10.) nGood++;
	if(h1->GetMaximum()>nMax) {nMax = (Int_t)h1->GetMaximum(); iXmax=iX; iZmax=iZ;}
	h1->Write(); 
      }
      
      for(Int_t iGain=0; iGain<2; iGain++) {
	h2 = da1.GetTimeEnergyHistogram(iX,iZ,iGain); // Time vs Energy
	if(h2) h2->Write();
      }
      
    }
  }
  
  f->Close();
  
  /* Store output files to the File Exchange Server */
  
  for(Int_t iMod=0; iMod<5; iMod++) {
    sprintf(localfile,"PHOS_Module%d_LED.root",iMod);
    daqDA_FES_storeFile(localfile,"LED");
  }
  
  printf("%d physics events of %d total processed.\n",nevents_physics,nevents_total);
  printf("%d histograms has >10 entries in maximum, max. is %d entries ",nGood,nMax);
  printf("(module,iX,iZ)=(%d,%d,%d)",fMod,iXmax,iZmax);
  printf("\n");

  return status;
}
