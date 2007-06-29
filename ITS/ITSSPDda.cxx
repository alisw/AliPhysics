////////////////////////////////////////////////////////////////////////////////
// This program can be run in two modes.                                      //
//                                                                            //
// 1. With the DAQ DA framework on. This is the default operating mode.       //
// Call this program with the name of the executable followed by the          //
// data files to process.                                                     //
//                                                                            //
// 2. Without the DAQ DA framework on. Define the SPD_DA_OFF environment var. //
// Call this program with the name of the executable followed by the          //
// runNr and the data files to process.                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

extern "C" {
#include "daqDA.h"
}
#include "event.h"
#include "monitor.h"
#include "AliRawReaderDate.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSOnlineSPDscan.h"
#include "AliITSOnlineSPDscanSingle.h"
#include "AliITSOnlineSPDscanMultiple.h"
#include "AliITSOnlineSPDscanMeanTh.h"
#include "AliITSOnlineSPDscanAnalyzer.h"
#include "AliLog.h"
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TPluginManager.h>


int main(int argc, char **argv) {
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // make directory structure (if not here already):
  system("mkdir ./calibResults >& /dev/null");
  system("mkdir ./calibResults/Noisy >& /dev/null");
  system("mkdir ./calibResults/NoisyToFXS >& /dev/null");
  system("mkdir ./calibResults/Parameters >& /dev/null");
  system("mkdir ./calibResults/Reference >& /dev/null");
  char *saveDirNoisy         = "./calibResults/Noisy";
#ifndef SPD_DA_OFF
  char *saveDirNoisyToFXS    = "./calibResults/NoisyToFXS";
#endif
  char *saveDirParameters    = "./calibResults/Parameters";
  char *saveDirRef           = "./calibResults/Reference";

  // This line is needed in case of a stand-alone application w/o
  // $ROOTSYS/etc/system.rootrc file
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");

  // turn off annoying warning messages
  new AliLog;
  AliLog::Instance()->SetGlobalDebugLevel(-20);

  // calib scan types
  enum calib_types{MINTH,MEANTH,DAC,UNIMA,NOISE,DELAY};


  // ********* STEP 1: Produce scan container files. ***********************************
  int startSeg = 1;

#ifndef SPD_DA_OFF
  int runNr = atoi(getenv("DATE_RUN_NUMBER"));
#else
  int runNr = atoi(argv[1]);
  startSeg = 2;
#endif

#ifndef SPD_DA_OFF
  UInt_t nrNoisyFilesProduced=0;
#endif

  // container objects
  AliITSOnlineSPDscan *scanObj[20];
  Bool_t bScanInit[20];
  for (UInt_t eqId=0; eqId<20; eqId++) {
    scanObj[eqId]=NULL;
    bScanInit[eqId]=kFALSE;
  }
  // header data variables
  UInt_t routerNr[20];
  Bool_t halfStaveScanned[20][6];
  UInt_t type[20];
  Bool_t dataFormat[20];
  UInt_t triggers[20];
  Bool_t chipPresent[20][6][10];
  UInt_t dacStart[20];  
  UInt_t dacEnd[20];
  UInt_t dacStep[20];
  UInt_t dacId[20];
  UInt_t rowStart[20];  
  UInt_t rowEnd[20];
  UInt_t rowValue[20];
  UInt_t dacValue[20];
  UInt_t dacHigh[20][6];
  UInt_t dacLow[20][6];
  UInt_t TPAmp[20][6];
  Bool_t minTHchipPresent[20][10];
  // current scan step flag
  Int_t currentStep[20];
  for (UInt_t eqId=0; eqId<20; eqId++) currentStep[eqId] = 9999;

  // loop over run segments
  for (int segNr=startSeg; segNr<argc; segNr++) {


    int status;

    /* define data source : this is argument 1 */  
    status=monitorSetDataSource( argv[segNr] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }
    /* declare monitoring program */
    status=monitorDeclareMp("ITS_SPD_CAL");
    if (status!=0) {
      printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
      return -1;
    }
    /* define wait event timeout - 1s max */
    monitorSetNowait();
    monitorSetNoWaitNetworkTimeout(1000);


    struct eventHeaderStruct *event;
    eventTypeType eventT;
    UInt_t eventNr=0;


    /* main loop (infinite) */
    for(;;) {

      /* check shutdown condition */
#ifndef SPD_DA_OFF
      if (daqDA_checkShutdown()) {break;}
#endif
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
	printf("event==NULL\n");
	continue;
      }

      eventT=event->eventType;
      if (eventT == PHYSICS_EVENT){


	eventNr++;
	//	printf("eventNr %d\n",eventNr);

	AliRawReader *reader = new AliRawReaderDate((void*)event);
	AliITSRawStreamSPD *str = new AliITSRawStreamSPD(reader);


	for (UInt_t eqId=0; eqId<20; eqId++) {

	  reader->Reset();
	  reader->Select("ITSSPD",eqId,eqId);

	  // Hit Event flags, specific for one event
	  Bool_t hitEventHSIncremented[6];
	  Bool_t hitEventChipIncremented[6][10];
	  for (UInt_t hs=0; hs<6; hs++) {
	    hitEventHSIncremented[hs] = kFALSE;
	    for (UInt_t chip=0; chip<10; chip++) {
	      hitEventChipIncremented[hs][chip] = kFALSE;
	    }
	  }

	  if (str->ReadCalibHeader()) {
	    // first check the type:
	    if (bScanInit[eqId] && type[eqId]!=str->GetHtype()) {
	      printf("Calib header problem. Type changed (%d -> %d)!\n",type[eqId],str->GetHtype());
	    }

	    // read calib values
	    routerNr[eqId]    = str->GetHrouterNr();
	    type[eqId]        = str->GetHtype();
	    dataFormat[eqId]  = str->GetHdataFormat();
	    triggers[eqId]    = str->GetHtriggers();
	    dacStart[eqId]    = str->GetHdacStart();
	    dacEnd[eqId]      = str->GetHdacEnd();
	    dacStep[eqId]     = str->GetHdacStep();
	    dacId[eqId]       = str->GetHdacId();
	    rowStart[eqId]    = str->GetHrowStart();
	    rowEnd[eqId]      = str->GetHrowEnd();
	    rowValue[eqId]    = str->GetHrowValue();
	    dacValue[eqId]    = str->GetHdacValue();
	    for (UInt_t hs=0; hs<6; hs++) {
	      halfStaveScanned[eqId][hs] = str->GetHhalfStaveScanned(hs);
	      dacHigh[eqId][hs]          = str->GetHdacHigh(hs);
	      dacLow[eqId][hs]           = str->GetHdacLow(hs);
	      TPAmp[eqId][hs]            = str->GetHTPAmp(hs);
	      for (UInt_t chip=0; chip<10; chip++) {
		chipPresent[eqId][hs][chip]      = str->GetHchipPresent(hs,chip);
	      }
	    }
	    for (UInt_t chip=0; chip<10; chip++) {
	      minTHchipPresent[eqId][chip] = str->GetHminTHchipPresent(chip);
	    }

	    currentStep[eqId] = (dacValue[eqId]-dacStart[eqId])/dacStep[eqId];
	    if (type[eqId]==DELAY) {
	      currentStep[eqId]=currentStep[eqId]*2;
	      dacValue[eqId]=dacValue[eqId]*2;
	      if (dacHigh[eqId][0]==128) { // misc_ctrl value
		currentStep[eqId]=currentStep[eqId]+1;
		dacValue[eqId]=dacValue[eqId]+1;
	      }
	    }

	    // router nr check:
	    if (routerNr[eqId]!=eqId) {
	      printf("Router nr problem? Router nr %d != EqID %d\n",routerNr[eqId],eqId);
	    }

	    if (!bScanInit[eqId]) {
	      // initialize container object
	      Char_t fileName[200];
	      sprintf(fileName,"%s/SPDcal_run_%d_eq_%d.root",saveDirRef,runNr,eqId);
	      switch (type[eqId]) {
	      case NOISE:
	      case UNIMA:
		scanObj[eqId] = new AliITSOnlineSPDscanSingle(fileName); 
		((AliITSOnlineSPDscanSingle*)scanObj[eqId])->ClearThis();
		bScanInit[eqId]=kTRUE;
		break;
	      case MINTH:
	      case DAC:
	      case DELAY:
		scanObj[eqId] = new AliITSOnlineSPDscanMultiple(fileName);
		scanObj[eqId]->ClearThis();
		bScanInit[eqId]=kTRUE;
		break;
	      case MEANTH: 
		scanObj[eqId] = new AliITSOnlineSPDscanMeanTh(fileName);
		scanObj[eqId]->ClearThis();
		bScanInit[eqId]=kTRUE;
		break;
	      default:
		printf("Unknown scan type: %d.\n",type[eqId]);
	      }
	      // some multiple scan data
	      if (type[eqId]==MINTH || type[eqId]==MEANTH || type[eqId]==DAC || type[eqId]==DELAY) {
		((AliITSOnlineSPDscanMultiple*)scanObj[eqId])->SetDacId(dacId[eqId]);
	      }
	      // some common data
	      scanObj[eqId]->SetRunNr((UInt_t)runNr);
	      scanObj[eqId]->SetRouterNr(routerNr[eqId]);
	      for (UInt_t hs=0; hs<6; hs++) {
		scanObj[eqId]->SetHalfStaveScanned(hs,halfStaveScanned[eqId][hs]);
	      }
	      scanObj[eqId]->SetType(type[eqId]);
	      scanObj[eqId]->SetDataFormat(dataFormat[eqId]);
	      for (Int_t hs=0; hs<6; hs++) {
	      
		// remove later when the chip present is set correctly !!!!!!!!!!!!!!!!!!!!!!!!!!!
		Bool_t halfStavePresent = str->GetHalfStavePresent(hs);
		// remove later when the chip present is set correctly !!!!!!!!!!!!!!!!!!!!!!!!!!!

		for (UInt_t chip=0; chip<10; chip++) {
		  scanObj[eqId]->SetChipPresent(hs,chip,chipPresent[eqId][hs][chip]);

		  // remove later when the chip present is set correctly !!!!!!!!!!!!!!!!!!!!!!!!!!!
		  if (halfStavePresent) scanObj[eqId]->SetChipPresent(hs,chip,kTRUE);
		  // remove later when the chip present is set correctly !!!!!!!!!!!!!!!!!!!!!!!!!!!

		}
	      }
	      scanObj[eqId]->SetRowStart(rowStart[eqId]);
	      scanObj[eqId]->SetRowEnd(rowEnd[eqId]);
	      scanObj[eqId]->SetDacStart(dacStart[eqId]);
	      scanObj[eqId]->SetDacEnd(dacEnd[eqId]);
	      scanObj[eqId]->SetDacStep(dacStep[eqId]);
	    }

	    if (type[eqId]==MINTH) {
	      scanObj[eqId]->SetTriggers(currentStep[eqId],triggers[eqId]);
	    }
	    if (type[eqId]==UNIMA || type[eqId]==NOISE) {
	      if (currentStep[eqId]==9999) printf("SPDcalibratorStep1 (eq %d): single step\n",eqId);
	      currentStep[eqId]=0;
	    }
	    if (type[eqId]==MINTH || type[eqId]==MEANTH || type[eqId]==DAC || type[eqId]==DELAY) {
	      ((AliITSOnlineSPDscanMultiple*)scanObj[eqId])->SetDacValue(currentStep[eqId],dacValue[eqId]);
	      if (type[eqId]==DELAY) {
		printf("SPDcalibratorStep1 (eq %d): DAC %d/%d , step %d\n",eqId,dacValue[eqId]/2,dacHigh[eqId][0],currentStep[eqId]);
	      }
	      else {
		printf("SPDcalibratorStep1 (eq %d): DAC %d , step %d\n",eqId,dacValue[eqId],currentStep[eqId]);
	      }
	    }
	    if (type[eqId]==MEANTH) {
	      for (Int_t hs=0; hs<6; hs++) {
		((AliITSOnlineSPDscanMeanTh*)scanObj[eqId])->SetDacLow(currentStep[eqId],hs,dacLow[eqId][hs]);
		((AliITSOnlineSPDscanMeanTh*)scanObj[eqId])->SetDacHigh(currentStep[eqId],hs,dacHigh[eqId][hs]);
		((AliITSOnlineSPDscanMeanTh*)scanObj[eqId])->SetTPAmp(currentStep[eqId],hs,TPAmp[eqId][hs]);
	      }
	    }


	  }
	
	  if (bScanInit[eqId]) {
	    while (str->Next()) {
	      UInt_t hs = str->GetHalfStaveNr();
	      UInt_t chip = str->GetChipAddr();
#ifndef SPD_DA_OFF
	      if (type[eqId]!=MINTH || minTHchipPresent[eqId][chip]) { 
#else
		if (type[eqId]!=MINTH || minTHchipPresent[eqId][chip] || runNr<=416900) { 
#endif
		scanObj[eqId]->IncrementHits(currentStep[eqId],hs,chip,str->GetChipCol(),str->GetChipRow());
	    
		if (!hitEventHSIncremented[hs]) {
		  scanObj[eqId]->IncrementHitEventsTot(currentStep[eqId],hs);
		  hitEventHSIncremented[hs]=kTRUE;
		}
	    
		if (!hitEventChipIncremented[hs][chip]) {
		  scanObj[eqId]->IncrementHitEvents(currentStep[eqId],hs,chip);
		  hitEventChipIncremented[hs][chip]=kTRUE;
		}
	      }

	    }

	    if (type[eqId]!=MINTH) { // for minth, triggers are set from header info
	      scanObj[eqId]->IncrementTriggers(currentStep[eqId]);
	    }

	  }

	}

	delete str;
	delete reader;

      }

      /* free resources */
      free(event);

    }

    status = monitorLogout();
    if (status != 0) {
      printf("monitorLogout() failed : %s\n",monitorDecodeError(status));
      return -1;
    }
    
#ifndef SPD_DA_OFF
    daqDA_progressReport((unsigned int)( ((Float_t)(segNr-startSeg+1))/(argc-startSeg)*50 ));
#else
    printf("progress: %d\n",(unsigned int)( ((Float_t)(segNr-startSeg+1))/(argc-startSeg)*50 ));
#endif

  }

  // clean up scan objects
  for (UInt_t eqId=0; eqId<20; eqId++) {
    if (scanObj[eqId]!=NULL) delete scanObj[eqId];
  }





  // ********* STEP 2: Analyze scan container files. ***********************************
#ifndef SPD_DA_OFF
  // clear noisyToFXS dir:
  Char_t command[200];
  sprintf(command,"cd %s; rm -f *",saveDirNoisyToFXS);
  system(command);
#endif

  AliITSOnlineSPDscanAnalyzer *analyzer;

  // *** *** *** start loop over equipments (eq_id)
  for (int eqId=0; eqId<20; eqId++) {

    Char_t fileName[200];
    sprintf(fileName,"%s/SPDcal_run_%d_eq_%d.root",saveDirRef,runNr,eqId);
    analyzer = new AliITSOnlineSPDscanAnalyzer(fileName);

    Int_t type  = analyzer->GetType();
    Int_t dacId = analyzer->GetDacId();
    if (type!=99) {
      if (type==DAC) {
	printf("SPD calibrator Step2: eqId %d, type %d, dacId %d\n",eqId,type,dacId);
      }
      else printf("SPD calibrator Step2: eqId %d type %d\n",eqId,type);  
    }




    if (type==UNIMA) {
    }
    else if (type==NOISE) {
      if (analyzer->ProcessNoisyPixels(saveDirNoisy)) {
	for (UInt_t module=0; module<240; module++) {
	  if (analyzer->SaveDeadNoisyPixels(module,saveDirNoisy)) {
#ifndef SPD_DA_OFF
	    nrNoisyFilesProduced++;
	    Char_t command[100];
	    sprintf(command,"cp %s/SPD_DeadNoisy_%d.root %s/.",saveDirNoisy,module,saveDirNoisyToFXS);
	    system(command);
#endif
	  }
	}
      }
    }
//    else if (type==DAC && dacId==42) {
//      Char_t ofileName[100];
//      sprintf(ofileName,"%s/delay_eq_%d.txt",saveDirParameters,eqId);
//      ofstream ofile;
//      ofile.open (ofileName);
//      for (UInt_t hs=0; hs<6; hs++) {
//	for (UInt_t chipNr=0; chipNr<10; chipNr++) {
//	  ofile << analyzer->GetDelay(hs,chipNr);
//	  ofile << "\t";
//	}
//	ofile << "\n";
//      }
//      ofile.close();
//#ifndef SPD_DA_OFF
//      Char_t id[20];
//      sprintf(id,"SPD_delay_%d",eqId);
//      Int_t status = daqDA_FES_storeFile(ofileName,id);
//      if (status) {
//	printf("Failed to export file %s , status %d\n",ofileName,status);
//      }
//#endif
//    }
    else if (type==MINTH || (type==DAC && dacId==39)) {
      Char_t ofileName[100];
      sprintf(ofileName,"%s/minth_eq_%d.txt",saveDirParameters,eqId);
      ofstream ofile;
      ofile.open (ofileName);
      for (UInt_t hs=0; hs<6; hs++) {
	for (UInt_t chipNr=0; chipNr<10; chipNr++) {
	  ofile << analyzer->GetMinTh(hs,chipNr);
	  ofile << "\t";
	}
	ofile << "\n";
      }
      ofile.close();
#ifndef SPD_DA_OFF
      Char_t id[20];
      sprintf(id,"SPD_minth_%d",eqId);
      Int_t status = daqDA_FES_storeFile(ofileName,id);
      if (status) {
	printf("Failed to export file %s , status %d\n",ofileName,status);
      }
#endif
    }
    else if (type==DELAY) {
      Char_t ofileName[100];
      sprintf(ofileName,"%s/delay_eq_%d.txt",saveDirParameters,eqId);
      ofstream ofile;
      ofile.open (ofileName);
      for (UInt_t hs=0; hs<6; hs++) {
	for (UInt_t chipNr=0; chipNr<10; chipNr++) {
	  UInt_t clockCycle = analyzer->GetDelay(hs,chipNr);
	  UInt_t delayCtrl = clockCycle/2;
	  UInt_t miscCtrl = 192;
	  if (clockCycle%2==1) miscCtrl = 128;
	  ofile << delayCtrl << "/" << miscCtrl;
	  ofile << "\t";
	}
	ofile << "\n";
      }
      ofile.close();
#ifndef SPD_DA_OFF
      Char_t id[20];
      sprintf(id,"SPD_delay_%d",eqId);
      Int_t status = daqDA_FES_storeFile(ofileName,id);
      if (status) {
	printf("Failed to export file %s , status %d\n",ofileName,status);
      }
#endif
    }


    delete analyzer;

#ifndef SPD_DA_OFF
    daqDA_progressReport((unsigned int)(50+(eqId+1)*2.5));
#else
    printf("progress: %d\n",(unsigned int)(50+(eqId+1)*2.5));
#endif

  }
  // *** *** *** end loop over equipments (eq_id)


#ifndef SPD_DA_OFF
  if (nrNoisyFilesProduced>0) {
    // send a tared file of all new noisy maps
    sprintf(command,"cd %s; tar -cf noisy.tar *",saveDirNoisyToFXS);
    printf("\n\n%s\n\n",command);
    system(command);
    Char_t fileName[200];
    sprintf(fileName,"%s/noisy.tar",saveDirNoisyToFXS);
    Char_t id[20];
    sprintf(id,"SPD_noisy");
    Int_t status = daqDA_FES_storeFile(fileName,id);
    if (status!=0) {
      printf("Failed to export file %s , status %d\n",fileName,status);
      return -1;
    }
  }
#endif


  return 0;
}
