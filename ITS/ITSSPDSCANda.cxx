/*
Contact: henrik.tydesjo@cern.ch
Link: tydes.home.cern.ch/tydes/doc/CalibrationOverview/CalibrationAlgorithms/
Run Type: DAQ_MIN_TH_SCAN,DAQ_MEAN_TH_SCAN,DAQ_UNIFORMITY_SCAN,DAQ_NOISY_PIX_SCAN,DAQ_PIX_DELAY_SCAN
DA Type: LDC
Number of events needed: Depending on scan type
Input Files: spd_standal_params,spd_perm_noisy ,  ./calibResults/ScanNoisy/* ,  raw data
Output Files: ./calibResults/ScanReference/* ,  ./calibResults/ScanDCSconfigToFXS/* ,  ./calibResults/ScanNoisyToFXS/* ,  ./calibResults/ScanNoisy/*
Trigger types used: PHYSICS
*/

////////////////////////////////////////////////////////////////////////////////
// This program can be compiled in two modes.                                 //
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

#ifndef SPD_DA_OFF
extern "C" {
#include "daqDA.h"
}
#endif
#include "event.h"
#include "monitor.h"
#include "AliRawReaderDate.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSOnlineSPDscan.h"
#include "AliITSOnlineSPDscanSingle.h"
#include "AliITSOnlineSPDscanMultiple.h"
#include "AliITSOnlineSPDscanMeanTh.h"
#include "AliITSOnlineSPDscanAnalyzer.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliLog.h"
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TPluginManager.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TFitter.h>

int main(int argc, char **argv) {
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // directory structure, hard coded
  char *saveDirNoisy         = "./calibResults/ScanNoisy";          // may NOT delete content
  char *saveDirNoisyToFXS    = "./calibResults/ScanNoisyToFXS";     //     may delete content
  char *saveDirDCSconfigToFXS= "./calibResults/ScanDCSconfigToFXS"; //     may delete content
  char *saveDirRef           = "./calibResults/ScanReference";      //     may delete content
  char *saveDirIdsToFXS      = "./calibResults/IdsToFXS";           //     may delete content
  char *configFilesDir       = "./configFiles";                     //     may delete content

  // make sure the directory structure is correct:
  system("mkdir ./calibResults >& /dev/null");
  system("mkdir ./calibResults/ScanNoisy >& /dev/null");
  system("mkdir ./calibResults/ScanNoisyToFXS >& /dev/null");
  system("mkdir ./calibResults/ScanDCSconfigToFXS >& /dev/null");
  system("mkdir ./calibResults/ScanReference >& /dev/null");
  system("mkdir ./calibResults/IdsToFXS >& /dev/null");
  system("mkdir ./configFiles >& /dev/null");
  // prameters config files
  TString paramsFileName = Form("%s/standal_params.txt",configFilesDir);
  TString permNoisyFileName = Form("%s/perm_noisy.txt",configFilesDir);

  TFitter *fitter = new TFitter(3);
  TVirtualFitter::SetFitter(fitter);

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


  // ********* STEP 0: Get configuration files from db (if there are any) , then read parameters*********
  UInt_t nrTuningParams = 0;
  TObjArray paramNames;  paramNames.SetOwner(kTRUE);
  TObjArray paramVals;  paramVals.SetOwner(kTRUE);
  
  // tuning parameters:
  Int_t status = 0;
#ifndef SPD_DA_OFF
  TString idp = "spd_standal_params";
  status=daqDA_DB_getFile(idp.Data(),paramsFileName.Data());
  if (status) {
    printf("Failed to get config file %s: status=%d. Using default tuning parameters.\n",idp.Data(),status);
    TString rmCmd = Form("rm -f %s",paramsFileName.Data());
    system(rmCmd.Data());
  }
#endif
  if (status==0) {
    ifstream paramsFile;
    paramsFile.open(paramsFileName.Data(), ifstream::in);
    if (paramsFile.fail()) {
      printf("No config file (%s) present. Using default tuning parameters.\n",paramsFileName.Data());
    }
    else {
      while(1) {
	Char_t paramN[50];
	Char_t paramV[50];
	paramsFile >> paramN;
	if (paramsFile.eof()) break;
	paramsFile >> paramV;
	paramNames.AddAtAndExpand(new TObjString(paramN),nrTuningParams);
	paramVals.AddAtAndExpand(new TObjString(paramV),nrTuningParams);
	nrTuningParams++;
	if (paramsFile.eof()) break;
      }
      paramsFile.close();
    }
  }
  //  for (UInt_t i=0; i<nrTuningParams; i++) {
  //  //    printf("Entry %d: N=%s , V=%s\n",i,((TString*)paramNames.At(i))->Data(),((TString*)paramVals.At(i))->Data());
  //    printf("Entry %d: N=%s , V=%s\n",i,((TObjString*)paramNames.At(i))->GetString().Data(),((TObjString*)paramVals.At(i))->GetString().Data());
  //  }

  // perm noisy list:
  Int_t permstatus = 0;
#ifndef SPD_DA_OFF
  TString idn = "spd_perm_noisy";
  permstatus=daqDA_DB_getFile(idn.Data(),permNoisyFileName.Data());
  if (permstatus) {
    printf("Failed to get config file %s: status=%d. No permanently noisy pixels will be added.\n",idn.Data(),permstatus);
    TString rmCmd = Form("rm -f %s",permNoisyFileName.Data());
    system(rmCmd.Data());
  }
#endif




  // ********* STEP 1: Produce scan container files (Reference Data). ***********************************
  int startSeg = 1;

#ifndef SPD_DA_OFF
  if (getenv("DATE_RUN_NUMBER")==0) {
    printf("DATE_RUN_NUMBER not properly set.\n");
    return -1;
  }
  int runNr = atoi(getenv("DATE_RUN_NUMBER"));
#else
  int runNr = atoi(argv[1]);
  startSeg = 2;
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
    
    
    UInt_t eventNr=0;

    /* main loop (infinite) */
    for(;;) {
      
      struct eventHeaderStruct *event;
      eventTypeType eventT;
      
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
	continue;
      }

      eventT=event->eventType;
      if (eventT == PHYSICS_EVENT) {
	
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
	  
	  if (str->ReadCalibHeader()>0) {
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
	    dacValue[eqId]    = str->GetHdacValue(); // this will change below for MEANTH scan

	    for (UInt_t hs=0; hs<6; hs++) {
	      halfStaveScanned[eqId][hs] = str->GetHhalfStaveScanned(hs);
	      dacHigh[eqId][hs]          = str->GetHdacHigh(hs);
	      dacLow[eqId][hs]           = str->GetHdacLow(hs);
	      TPAmp[eqId][hs]            = str->GetHTPAmp(hs);
	      for (UInt_t chip=0; chip<10; chip++) {
		chipPresent[eqId][hs][chip] = str->GetHchipPresent(hs,chip);
		if (type[eqId]==MEANTH && chipPresent[eqId][hs][chip]) dacValue[eqId] = str->GetHdacLow(hs);
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
	      TString fileName = Form("%s/SPDcal_run_%d_eq_%d.root",saveDirRef,runNr,eqId);
	      switch (type[eqId]) {
	      case NOISE:
	      case UNIMA:
		scanObj[eqId] = new AliITSOnlineSPDscanSingle(fileName.Data());
		((AliITSOnlineSPDscanSingle*)scanObj[eqId])->ClearThis();
		bScanInit[eqId]=kTRUE;
		break;
	      case MINTH:
	      case DAC:
	      case DELAY:
		scanObj[eqId] = new AliITSOnlineSPDscanMultiple(fileName.Data());
		scanObj[eqId]->ClearThis();
		bScanInit[eqId]=kTRUE;
		break;
	      case MEANTH: 
		scanObj[eqId] = new AliITSOnlineSPDscanMeanTh(fileName.Data());
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
		for (UInt_t chip=0; chip<10; chip++) {
		  scanObj[eqId]->SetChipPresent(hs,chip,chipPresent[eqId][hs][chip]);
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
	      //***remove last condition when minthpresent put correctly in calib header?
#ifndef SPD_DA_OFF
	      if (type[eqId]!=MINTH || minTHchipPresent[eqId][chip] || runNr<=416900) {
#else
	      if (type[eqId]!=MINTH || minTHchipPresent[eqId][chip] || runNr<=416900) {
#endif
		//*************************************************************************
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

#ifndef SPD_DA_OFF
    daqDA_progressReport((unsigned int)( ((Float_t)(segNr-startSeg+1))/(argc-startSeg)*50 ));
#else
    printf("progress: %d\n",(unsigned int)( ((Float_t)(segNr-startSeg+1))/(argc-startSeg)*50 ));
#endif

  }

  
  // clean up scan objects (which also saves them) , check if something happened...
    Bool_t somethingHappened = kFALSE;
  for (UInt_t eqId=0; eqId<20; eqId++) {
    if (scanObj[eqId]!=NULL) {
      delete scanObj[eqId];
      somethingHappened = kTRUE;
    }
  }
  if (!somethingHappened) {
    printf("WARNING: No data processed. Are the calibration headers missing?\n");
  }




  // ********* STEP 2: Analyze scan container files. ************************************************

  // clear noisyToFXS and DCSconfigToFXS dirs:
  TString command = Form("cd %s; rm -f *",saveDirNoisyToFXS);
  system(command.Data());
  TString command2 = Form("cd %s; rm -f *",saveDirDCSconfigToFXS);
  system(command2.Data());
  UInt_t nrNoisyFilesProduced=0;
  UInt_t nrDCSconfigFilesProduced=0;

  AliITSOnlineCalibrationSPDhandler* handler = new AliITSOnlineCalibrationSPDhandler();
  AliITSOnlineSPDscanAnalyzer *analyzer = NULL;
  AliITSOnlineCalibrationSPDhandler* handlerPermNoisy = NULL;
  // fill permanent noisy list to add later...
  if (permstatus==0) { 
    handlerPermNoisy = new AliITSOnlineCalibrationSPDhandler();
    UInt_t permNoisy = handlerPermNoisy->ReadNoisyFromText(permNoisyFileName.Data(),240); // 240 = read for all modules
    if (permNoisy>0) {
      printf("%d noisy pixels read from permanent list.\n",permNoisy);
    }
  }

  Bool_t reset_made = kFALSE;

  // *** *** *** start loop over equipments (eq_id)
  for (int eqId=0; eqId<20; eqId++) {

    // create analyzer for this eq
    TString fileName = Form("%s/SPDcal_run_%d_eq_%d.root",saveDirRef,runNr,eqId);
    analyzer = new AliITSOnlineSPDscanAnalyzer(fileName.Data(),handler);

    // configure analyzer with tuning parameters etc:
    for (UInt_t i=0; i<nrTuningParams; i++) {
      analyzer->SetParam(((TObjString*)paramNames.At(i))->GetString().Data(),((TObjString*)paramVals.At(i))->GetString().Data());
    }

    UInt_t type  = analyzer->GetType();
    UInt_t dacId = analyzer->GetDacId();
    UInt_t routerNr = analyzer->GetRouterNr();
    if (type!=99) {
      if (type==DAC) {
	printf("SPD scan calibrator Step2: eqId %d, type %d, dacId %d\n",eqId,type,dacId);
      }
      else printf("SPD scan calibrator Step2: eqId %d type %d\n",eqId,type);  
    }



    // algorithms for the different types of scans:

    if (type==UNIMA) {

    }

    else if (type==NOISE) {
      // read previous noisy list (clear if overwriting)
      handler->SetFileLocation(saveDirNoisy);
      if (analyzer->IsOverWriteSet() && !reset_made) {
	handler->ResetNoisy();
	handler->WriteToFilesAlways();
	reset_made=kTRUE;
      }
      else {
	handler->ReadFromFiles();
      }
      if (analyzer->ProcessNoisyPixels(/*saveDirNoisy*/)) {
	if (permstatus==0) {
	  handler->AddNoisyFrom(handlerPermNoisy);
	}
	// init dcs config text file
	TString dcsConfigFileName = Form("%s/dcsConfig_run_%d_eq_%d.txt",saveDirDCSconfigToFXS,runNr,eqId);
	ofstream dcsfile;
	dcsfile.open(dcsConfigFileName.Data());
	dcsfile << "[SPD SCAN]\n";
	dcsfile << "RunNumber=" << runNr << "\n";
	dcsfile << "Type=" << type << "\n";
	dcsfile << "Router=" << routerNr << "\n";
	dcsfile << "ActualDetConfiguration=" << "0,-1,-1\n"; // dummy values for now
	dcsfile << "[NOISY]\n";
	nrDCSconfigFilesProduced++;

	for (UInt_t hs=0; hs<6; hs++) {
	  for (UInt_t chip=0; chip<10; chip++) {
	    if (analyzer->IsChipPresent(hs,chip) || analyzer->IsOverWriteSet()) {
	      dcsfile << "-" << eqId << "," << hs << "," << chip << "\n";
	      UInt_t nrNoisy = handler->GetNrNoisyC(eqId,hs,chip);
	      for (UInt_t ind=0; ind<nrNoisy; ind++) {
		UInt_t col = handler->GetNoisyColAtC(eqId,hs,chip,ind);
		UInt_t row = handler->GetNoisyRowAtC(eqId,hs,chip,ind);
		dcsfile << col << "," << row << "\n";
	      }
	    }
	  }
	}
	handler->SetFileLocation(saveDirNoisy);
	handler->WriteNoisyToFile(eqId);
	handler->SetFileLocation(saveDirNoisyToFXS);
	handler->WriteNoisyToFile(eqId);
	nrNoisyFilesProduced++;
	
	dcsfile.close();
      }
    }

    else if (type==MINTH || (type==DAC && dacId==39)) {
      // init dcs config text file
      TString dcsConfigFileName = Form("%s/dcsConfig_run_%d_eq_%d.txt",saveDirDCSconfigToFXS,runNr,eqId);
      ofstream dcsfile;
      dcsfile.open(dcsConfigFileName.Data());
      dcsfile << "[SPD SCAN]\n";
      dcsfile << "RunNumber=" << runNr << "\n";
      dcsfile << "Type=" << type << "\n";
      dcsfile << "Router=" << routerNr << "\n";
      dcsfile << "ActualDetConfiguration=" << "0,-1,-1\n"; // dummy values for now
      dcsfile << "[DACvalues]\n";
      nrDCSconfigFilesProduced++;
      for (UInt_t hs=0; hs<6; hs++) {
	for (UInt_t chipNr=0; chipNr<10; chipNr++) {
	  Int_t minTh = -1;
	  if (analyzer->GetOnlineScan()->GetChipPresent(hs,chipNr)) {
	    minTh = analyzer->GetMinTh(hs,chipNr);
	    if (minTh!=-1) {
	      dcsfile << "39," << eqId << "," << hs << "," << chipNr << "=" << minTh << "\n";
	    }
	    else {
	      printf("MinTh failed for Eq %d , HS %d , Chip %d\n",eqId,hs,chipNr);
	    }
	  }
	}
      }
      dcsfile.close();
    }

    else if (type==DELAY) {
      // init dcs config text file
      TString dcsConfigFileName = Form("%s/dcsConfig_run_%d_eq_%d.txt",saveDirDCSconfigToFXS,runNr,eqId);
      ofstream dcsfile;
      dcsfile.open(dcsConfigFileName.Data());
      dcsfile << "[SPD SCAN]\n";
      dcsfile << "RunNumber=" << runNr << "\n";
      dcsfile << "Type=" << type << "\n";
      dcsfile << "Router=" << routerNr << "\n";
      dcsfile << "ActualDetCoonfiguration=" << "0,-1,-1\n"; // dummy values for now
      dcsfile << "[DACvalues]\n";
      nrDCSconfigFilesProduced++;
      for (UInt_t hs=0; hs<6; hs++) {
	for (UInt_t chipNr=0; chipNr<10; chipNr++) {
	  Int_t clockCycle = -1;
	  Int_t delayCtrl = -1;
	  Int_t miscCtrl = -1;
	  if (analyzer->GetOnlineScan()->GetChipPresent(hs,chipNr)) {
	    clockCycle = analyzer->GetDelay(hs,chipNr);
	    delayCtrl = clockCycle/2;
	    miscCtrl = 192;
	    if (clockCycle!=-1) {
	      if (clockCycle%2==1) miscCtrl = 128;
	      dcsfile << "42," << eqId << "," << hs << "," << chipNr << "=" << delayCtrl << "\n";
	      dcsfile << "43," << eqId << "," << hs << "," << chipNr << "=" << miscCtrl << "\n";
	    }
	    else {
	      printf("Delay failed for Eq %d , HS %d , Chip %d\n",eqId,hs,chipNr);
	    }
	  }
	}
      }
      dcsfile.close();
    }

    delete analyzer;

#ifndef SPD_DA_OFF
    daqDA_progressReport((unsigned int)(50+(eqId+1)*2.5));
#else
    printf("progress: %d\n",(unsigned int)(50+(eqId+1)*2.5));
#endif

  }
  // *** *** *** end loop over equipments (eq_id)


  delete handler;
  if (handlerPermNoisy!=NULL) {
    delete handlerPermNoisy;
  }

  printf("Opening id list file\n");
  TString idsFXSFileName = Form("%s/FXSids_run_%d.txt",saveDirIdsToFXS,runNr);
  ofstream idsFXSfile;
  idsFXSfile.open(idsFXSFileName.Data());



  // send noisy data to FXS
  if (nrNoisyFilesProduced>0) {
    printf("Preparing noisy files\n");
    // send a tared file of all new noisy maps
    TString command = Form("cd %s; tar -cf noisy_scan.tar *",saveDirNoisyToFXS);
    //    printf("\n\n%s\n\n",command.Data());
    system(command.Data());
    TString fileName = Form("%s/noisy_scan.tar",saveDirNoisyToFXS);
    TString id = "SPD_scan_noisy";
#ifndef SPD_DA_OFF
    Int_t status = daqDA_FES_storeFile(fileName.Data(),id.Data());
    if (status!=0) {
      printf("Failed to export file %s , status %d\n",fileName.Data(),status);
      return -1;
    }
#endif
    idsFXSfile << Form("%s\n",id.Data());
  }

  // send dcs config files to FXS
  if (nrDCSconfigFilesProduced>0) {
    printf("Preparing DCS config files\n");
    // send a tared file of all the dcsConfig text files
    TString command = Form("cd %s; tar -cf dcsConfig.tar *",saveDirDCSconfigToFXS);
    //    printf("\n\n%s\n\n",command.Data());
    system(command.Data());
    TString fileName = Form("%s/dcsConfig.tar",saveDirDCSconfigToFXS);
    TString id = "SPD_dcsConfig";
#ifndef SPD_DA_OFF
    Int_t status = daqDA_FES_storeFile(fileName.Data(),id.Data());
    if (status!=0) {
      printf("Failed to export file %s , status %d\n",fileName.Data(),status);
      return -1;
    }
#endif
    //idsFXSfile << Form("%s\n",id.Data()); // do NOT write this id (this is not for preprocessor)
  }

  // send reference data to FXS
  for (UInt_t eqId=0; eqId<20; eqId++) {
    if (bScanInit[eqId]) {
      printf("Preparing reference data for eq %d\n",eqId);
      TString fileName = Form("%s/SPDcal_run_%d_eq_%d.root",saveDirRef,runNr,eqId);
      TString id = Form("SPD_ref_scan_%d",eqId);
#ifndef SPD_DA_OFF
      Int_t status = daqDA_FES_storeFile(fileName.Data(),id.Data());
      if (status!=0) {
	printf("Failed to export file %s , status %d\n",fileName.Data(),status);
	return -1;
      }
#endif
      idsFXSfile << Form("%s\n",id.Data());
    }
  }


  printf("Preparing id list file\n");
  idsFXSfile.close();
  TString id = "SPD_id_list";
#ifndef SPD_DA_OFF
  status = daqDA_FES_storeFile(idsFXSFileName.Data(),id.Data());
  if (status!=0) {
    printf("Failed to export file %s , status %d\n",idsFXSFileName.Data(),status);
    return -1;
  }
#endif


  printf("DA finished.\n");

  return 0;
}
