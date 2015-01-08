/*
Contact: henrik.tydesjo@cern.ch
Link: tydes.home.cern.ch/tydes/doc/CalibrationOverview/CalibrationAlgorithms/
Run Type: PHYSICS
DA Type: MON
Number of events needed: Depending on muliplicity per event
Input Files: spd_physics_params ,  ./calibResults/DeadReferenceTmp/* ,  ./calibResults/DeadToFXS/*
Output Files: ./calibResults/NoisyReference/* ,  ./calibResults/DeadReference/* ,  ./calibResults/NoisyToFXS/* ,  ./calibResults/DeadReferenceTmp/*,./calibResults/DeadToFXS/*
Trigger types used: PHYSICS
*/

////////////////////////////////////////////////////////////////////////////////
// This program can be compiled in two modes.                                 //
//                                                                            //
// 1. Online. With the DAQ DA framework. This is the default operating mode.  //
//                                                                            //
// 2. Offline. Without the DAQ DA framework. Define the SPD_DA_OFF            //
//    environment var. Call this program with the name of the executable      //
//    followed by the runNr and the data files to process.                    //
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
#include "AliITSOnlineSPDphys.h"
#include "AliITSOnlineSPDphysAnalyzer.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliLog.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <TROOT.h>
#include <TPluginManager.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <cstdlib>


int main(int argc, char **argv) {
  if (argc<2) {
    printf("SPD DA ERROR: Wrong number of arguments\n");
    return -1;
  }

  UInt_t nrErrors = 0;

  // directory structure, hard coded                              //   FOR OFFLINE VERSION:
  char *saveDirDead          = "./calibResults/Dead";             // may NOT delete content
  char *saveDirDeadToFXS     = "./calibResults/DeadToFXS";        //     may delete content
  char *saveDirDeadRef       = "./calibResults/DeadReference";    //     may delete content
  char *saveDirDeadRefTmp    = "./calibResults/DeadReferenceTmp"; // may NOT delete content
  char *saveDirNoisyToFXS    = "./calibResults/NoisyToFXS";       //     may delete content
  char *saveDirNoisyRef      = "./calibResults/NoisyReference";   //     may delete content
  char *saveDirIdsToFXS      = "./calibResults/IdsToFXS";         //     may delete content
  char *configFilesDir       = "./configFiles";                   //     may delete content
  // clean up and make sure the directory structure is put up correctly:
#ifndef SPD_DA_OFF
  system("rm ./calibResults -rf >& /dev/null");
  system("rm ./ITS -rf >& /dev/null");
#else
  system("rm ./calibResults/DeadToFXS -rf >& /dev/null");
  system("rm ./calibResults/DeadReference -rf >& /dev/null");
  system("rm ./calibResults/NoisyToFXS -rf >& /dev/null");
  system("rm ./calibResults/NoisyReference -rf >& /dev/null");
  system("rm ./calibResults/IdsToFXS -rf >& /dev/null");
#endif
  system("mkdir ./calibResults >& /dev/null");
  system("mkdir ./calibResults/Dead >& /dev/null");
  system("mkdir ./calibResults/DeadToFXS >& /dev/null");
  system("mkdir ./calibResults/DeadReference >& /dev/null");
  system("mkdir ./calibResults/DeadReferenceTmp >& /dev/null");
  system("mkdir ./calibResults/NoisyToFXS >& /dev/null");
  system("mkdir ./calibResults/NoisyReference >& /dev/null");
  system("mkdir ./calibResults/IdsToFXS >& /dev/null");
  system("mkdir ./configFiles >& /dev/null");
#ifndef SPD_DA_OFF
  system("mkdir ./ITS >& /dev/null");
  system("mkdir ./ITS/Calib >& /dev/null");
  system("mkdir ./ITS/Calib/SPDNoisy >& /dev/null");
#endif
  // parameters config file
  TString paramsFileName = Form("%s/physics_params.txt",configFilesDir);

  // This line is needed in case of a stand-alone application w/o
  // $ROOTSYS/etc/system.rootrc file
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");

  // turn off annoying warning messages
  // NB: Should not be handled here
  AliLog* logger = AliLog::GetRootLogger();
  logger->SetGlobalDebugLevel(-20);

  // ********* STEP 0: Get configuration files from db (if there are any) , then read parameters*********
  UInt_t nrTuningParams = 0;
  TObjArray paramNames;  paramNames.SetOwner(kTRUE);
  TObjArray paramVals;  paramVals.SetOwner(kTRUE);
  
  // tuning parameters:
  Int_t par_status = 0;
#ifndef SPD_DA_OFF
  TString idp = "spd_physics_params";
  par_status=daqDA_DB_getFile(idp.Data(),paramsFileName.Data());
  if (par_status) {
    printf("SPD DA: Failed to get config file %s: status=%d. Using default tuning parameters.\n",idp.Data(),par_status);
    TString rmCmd = Form("rm -f %s",paramsFileName.Data());
    system(rmCmd.Data());
  }
#endif
  if (par_status==0) {
    ifstream paramsFile;
    paramsFile.open(paramsFileName.Data(), ifstream::in);
    if (paramsFile.fail()) {
      printf("SPD DA: No config file (%s) present. Using default tuning parameters.\n",paramsFileName.Data());
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





  // create calibration handler (needed already at step 1 in order to fill which eq,hs,chips are active
  AliITSOnlineCalibrationSPDhandler* handler = new AliITSOnlineCalibrationSPDhandler();

  // Read silent=dead+inactive info from previous calibrations
#ifndef SPD_DA_OFF
  for (UInt_t eq=0; eq<20; eq++) {
    Int_t getPreviousDead_status = 0;
    TString idpd = Form("spd_previous_dead_%d",eq);
    TString fileName = Form("%s/SPD_Dead_%d.root",saveDirDead,eq);
    getPreviousDead_status = daqDA_DB_getFile(idpd.Data(),fileName.Data());
    //        printf("daqDA_DB_getFile(%s,%s)\n",idpd.Data(),fileName.Data());
    if (getPreviousDead_status) {
      printf("SPD DA: Failed to get dead file %s: status=%d. Dead search will start from zero for this eq.\n",idpd.Data(),getPreviousDead_status);
      TString rmCmd = Form("rm -f %s",fileName.Data());
      system(rmCmd.Data());
    }
  }
#endif
  handler->SetFileLocation(saveDirDead);
  handler->ReadSilentFromFiles();
  printf("SPD DA: Number of single dead pixels from previous runs: %d\n",handler->GetNrDead());

  // Read hit-maps from previous runs
#ifndef SPD_DA_OFF
  for (UInt_t eq=0; eq<20; eq++) {
    Int_t getPreviousHitmap_status = 0;
    TString idph = Form("spd_previous_hitmap_%d",eq);
    TString fileName = Form("%s/SPDphys_dead_run_0_0_eq_%d.root",saveDirDeadRefTmp,eq);
    getPreviousHitmap_status = daqDA_DB_getFile(idph.Data(),fileName.Data());
    //    printf("daqDA_DB_getFile(%s,%s)\n",idph.Data(),fileName.Data());
    if (getPreviousHitmap_status) {
      printf("SPD DA: Failed to get previous hitmap file %s: status=%d. Dead search will start with empty hitmap for this eq.\n",idph.Data(),getPreviousHitmap_status);
      TString rmCmd = Form("rm -f %s",fileName.Data());
      system(rmCmd.Data());
    }
  }
#endif




  // ********* STEP 1: Produce phys container files (Reference Data). ***********************************

#ifndef SPD_DA_OFF
  if (getenv("DATE_RUN_NUMBER")==0) {
    printf("SPD DA ERROR: DATE_RUN_NUMBER not properly set.\n");
    return -1;
  }
  int runNr = atoi(getenv("DATE_RUN_NUMBER"));
#else
  int runNr = atoi(argv[1]);
  int startSeg = 2;
#endif


  // container objects
  AliITSOnlineSPDphys *physObj[20];
  Bool_t bPhysInit[20];
  for (UInt_t eq=0; eq<20; eq++) {
    physObj[eq]=NULL;
    bPhysInit[eq]=kFALSE;
  }


  UInt_t eventNr=0;

  // loop over run segments in case of offline mode
#ifdef SPD_DA_OFF
  for (int segNr=startSeg; segNr<argc; segNr++) {
#endif

    int status;

    /* define data source : */  
#ifndef SPD_DA_OFF
    status=monitorSetDataSource( argv[1] ); // should be "^SPD" in order to get full detector online
#else
    status=monitorSetDataSource( argv[segNr] );
#endif
    if (status!=0) {
      printf("SPD DA ERROR: monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }
    /* declare monitoring program */
    status=monitorDeclareMp("ITS_SPD_PHYS");
    if (status!=0) {
      printf("SPD DA ERROR: monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
      return -1;
    }
    /* define wait event timeout - 1s max */
    monitorSetNowait();
    monitorSetNoWaitNetworkTimeout(1000);


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
	printf ("SPD DA: End of File detected\n");
	break; /* end of monitoring file has been reached */
      }

      if (status!=0) {
	printf("SPD DA: monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
	break;
      }

      /* retry if got no event */
      if (event==NULL) {
	continue;
      }

      eventT=event->eventType;
      if (eventT == PHYSICS_EVENT){

	//	printf("eventNr %d\n",eventNr);

	AliRawReader *reader = new AliRawReaderDate((void*)event);
	AliITSRawStreamSPD *str = new AliITSRawStreamSPD(reader);

	while (str->Next()) {

	  Int_t eq = reader->GetDDLID();
	  // check that this hs is active in handler object
	  if (!(handler->IsActiveEq(eq))) {
	    printf("SPD DA: Found Eq (%d) , previously inactive\n",eq);
	    handler->ActivateEq(eq);
	  }
	  if (eq>=0 && eq<20) {
	    if (!bPhysInit[eq]) { // this code is duplicated for the moment... (see also below)
	      TString fileName = Form("%s/SPDphys_run_%d_eq_%d.root",saveDirNoisyRef,runNr,eq);
	      physObj[eq] = new AliITSOnlineSPDphys(fileName.Data());
	      physObj[eq]->AddRunNr(runNr);
	      physObj[eq]->SetEqNr(eq);
	      bPhysInit[eq]=kTRUE;
	    }

	    UInt_t hs = str->GetHalfStaveNr();
	    // check that this hs is active in handler object
	    if (!(handler->IsActiveHS(eq,hs))) {
	      printf("SPD DA: Found HS (%d,%d) , previously inactive\n",eq,hs);
	      handler->ActivateHS(eq,hs);
	    }
	    UInt_t chip = str->GetChipAddr();
	    // check that this chip is active in handler object
	    if (!(handler->IsActiveChip(eq,hs,chip))) {
	      printf("SPD DA: Found Chip (%d,%d,%d) , previously inactive\n",eq,hs,chip);
	      handler->ActivateChip(eq,hs,chip);
	    }
	    physObj[eq]->IncrementHits(hs,chip,str->GetChipCol(),str->GetChipRow());
	    
	  }
	}


	// check which eq and hs are active, for first event only
	if (eventNr==0) {
	  for (UInt_t eq=0; eq<20; eq++) {
	    // activate Eq and HSs in handler object
	    if (str->IsActiveEq(eq)) {
	      handler->ActivateEq(eq);
	      for (UInt_t hs=0; hs<6; hs++) {
		if (str->IsActiveHS(eq,hs)) {
		  handler->ActivateHS(eq,hs);
		  for (UInt_t chip=0; chip<10; chip++) {
		    if (str->IsActiveChip(eq,hs,chip)) {
		      handler->ActivateChip(eq,hs,chip);
		    }
		    else {
		      handler->ActivateChip(eq,hs,chip,kFALSE);
		    }
		  }
		}
		else {
		  handler->ActivateHS(eq,hs,kFALSE);
		}
	      }
	      if (!bPhysInit[eq]) { // this code is duplicated for the moment... (see also above)
		TString fileName = Form("%s/SPDphys_run_%d_eq_%d.root",saveDirNoisyRef,runNr,eq);
		physObj[eq] = new AliITSOnlineSPDphys(fileName.Data());
		physObj[eq]->AddRunNr(runNr);
		physObj[eq]->SetEqNr(eq);
		bPhysInit[eq]=kTRUE;
	      }
	    }
	    else {
	      handler->ActivateEq(eq,kFALSE);
	    }
	  }
	}


	for (UInt_t eq=0; eq<20; eq++) {
	  if (bPhysInit[eq]) {
	    physObj[eq]->IncrementNrEvents();
	  }
	}

	delete str;
	delete reader;

	eventNr++;

      }

      /* free resources */
      free(event);

    }
    

#ifdef SPD_DA_OFF
    printf("SPD DA: progress: %d\n",(unsigned int)( ((Float_t)(segNr-startSeg+1))/(argc-startSeg)*50 ));
  }
#endif

  // clean up phys objects (also saves them)
  for (UInt_t eq=0; eq<20; eq++) {
    if (physObj[eq]!=NULL) delete physObj[eq];
  }

  printf("SPD DA: %d events collected for this run.\n",eventNr);




  // ********* STEP 2: Analyze phys container files. ************************************************

  time_t timeStamp = time(NULL);
  //  printf("*** Start step2 , %d\n",time(NULL) - timeStamp);

  UInt_t firstRunNrDead = runNr;


  UInt_t nrEnoughStatNoisy = 0;
  UInt_t nrEqActiveNoisy = 0;
  Bool_t eqActiveNoisy[20];

  // *** *** *** start loop over equipments (eq_id)
  for (UInt_t eq=0; eq<20; eq++) {
    eqActiveNoisy[eq] = kFALSE;

    // create analyzer for this eq
    TString fileName = Form("%s/SPDphys_run_%d_eq_%d.root",saveDirNoisyRef,runNr,eq);
    AliITSOnlineSPDphysAnalyzer *noisyAnalyzer = new AliITSOnlineSPDphysAnalyzer(fileName.Data(),handler);

    // check data in container
    if (noisyAnalyzer->GetEqNr() != eq) {
      if (noisyAnalyzer->GetEqNr() != 999) {
	printf("SPD DA ERROR: Mismatching EqId in Container data and filename (%d!=%d). Skipping eq.\n",noisyAnalyzer->GetEqNr(),eq);
      }
      delete noisyAnalyzer;
      continue;
    }

    nrEqActiveNoisy++;
    eqActiveNoisy[eq] = kTRUE;

    // configure analyzer with tuning parameters etc:
    for (UInt_t i=0; i<nrTuningParams; i++) {
      noisyAnalyzer->SetParam(((TObjString*)paramNames.At(i))->GetString().Data(),((TObjString*)paramVals.At(i))->GetString().Data());
    }

    printf("SPD DA: SPD phys STEP 2: Noisy search for eq %d\n",eq);  

    // search for noisy pixels:
    nrEnoughStatNoisy += noisyAnalyzer->ProcessNoisyPixels();

    // copy this phys obj to temporary dead reference dir to process after noisy search
    TString fileNameDead = Form("%s/SPDphys_dead_run_0_0_eq_%d.root",saveDirDeadRefTmp,eq);
    AliITSOnlineSPDphys* physObj = new AliITSOnlineSPDphys(fileNameDead.Data());
    physObj->AddPhys(noisyAnalyzer->GetOnlinePhys());
    if (physObj->GetNrRuns()>0) {
      UInt_t firstRunNr = physObj->GetRunNr(0);
      if (firstRunNrDead>firstRunNr) {
	firstRunNrDead=firstRunNr;
      }
    }
    // remove noisy pixels from dead hitmap
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	for (UInt_t ind=0; ind<handler->GetNrNoisyC(eq,hs,chip); ind++) {
	  UInt_t col  = handler->GetNoisyColAtC(eq,hs,chip,ind);
	  UInt_t row  = handler->GetNoisyRowAtC(eq,hs,chip,ind);
	  physObj->AddHits(hs,chip,col,row,-noisyAnalyzer->GetOnlinePhys()->GetHits(hs,chip,col,row));
	}
      }
    }

    delete physObj;
    delete noisyAnalyzer;

#ifndef SPD_DA_OFF
    //    daqDA_progressReport((unsigned int)((eq+1)*2.5));
    printf("SPD DA: progress: %d , %ld\n",(unsigned int)(50+(eq+1)*2.5), time(NULL) - timeStamp);
#else
    //    printf("SPD DA: progress: %d\n",(unsigned int)(50+(eq+1)*1.25));
    printf("SPD DA: progress: %d , %ld\n",(unsigned int)(50+(eq+1)*1.25), time(NULL) - timeStamp);
#endif
  }
  // *** *** *** end loop over equipments (eq_id)

  printf("SPD DA: Noisy search finished. %d noisy pixels found. %d chips had enough statistics.\n",
	 handler->GetNrNoisy(),nrEnoughStatNoisy);
  handler->SetFileLocation(saveDirNoisyToFXS);
  UInt_t nrNoisyFilesProduced = handler->WriteNoisyToFiles();
#ifndef SPD_DA_OFF
  // save noisy to local OCDB
  handler->WriteNoisyToDB(runNr,9999999,".");
  // store local OCDB file in daq Det DB
  TString ocdb_fileName = Form("./ITS/Calib/SPDNoisy/Run%d_9999999_v0_s0.root",runNr);
  TString idod = "spd_noisy_ocdb";
  Int_t ocdb_store_status = daqDA_DB_storeFile(ocdb_fileName.Data(),idod.Data());
  //    printf("daqDA_DB_storeFile(%s,%s)\n",ocdb_fileName.Data(),idod.Data());
  if (ocdb_store_status) {
    printf("SPD DA ERROR: Failed to store noisy ocdb file %s in daqDetDB: status=%d.\n",idod.Data(),ocdb_store_status);
  }
#endif





  UInt_t nrEnoughStatChips = 0;
  UInt_t nrDeadChips = 0;
  UInt_t nrInefficientChips = 0;
  UInt_t nrEqActiveDead = 0;
  Bool_t eqActiveDead[20];

  // *** *** *** start loop over equipments (eq_id)
  for (UInt_t eq=0; eq<20; eq++) {
    eqActiveDead[eq] = kFALSE;

    // setup analyzer for dead search
    TString fileNameDead = Form("%s/SPDphys_dead_run_0_0_eq_%d.root",saveDirDeadRefTmp,eq);
    AliITSOnlineSPDphys* physObj = new AliITSOnlineSPDphys(fileNameDead.Data());
    AliITSOnlineSPDphysAnalyzer* deadAnalyzer = new AliITSOnlineSPDphysAnalyzer(physObj,handler);
    // check data in container
    if (deadAnalyzer->GetEqNr() != eq) {
      if (deadAnalyzer->GetEqNr() != 999) {
	printf("SPD DA ERROR: Mismatching EqId in Dead Container data and filename (%d!=%d). Skipping eq.\n",deadAnalyzer->GetEqNr(),eq);
      }
      delete deadAnalyzer;
      nrDeadChips+=60; // since this eq is inactive...
      continue;
    }

    nrEqActiveDead++;
    eqActiveDead[eq] = kTRUE;

    // configure analyzer with tuning parameters etc:
    for (UInt_t i=0; i<nrTuningParams; i++) {
      deadAnalyzer->SetParam(((TObjString*)paramNames.At(i))->GetString().Data(),((TObjString*)paramVals.At(i))->GetString().Data());
    }

    UInt_t nrEventsCollected = physObj->GetNrEvents();
    printf("SPD DA: SPD phys STEP 2: Dead search for eq %d (%d events)\n",eq,nrEventsCollected);  

    // search for dead pixels:
    nrEnoughStatChips += deadAnalyzer->ProcessDeadPixels();
    nrDeadChips += deadAnalyzer->GetNrDeadChips();
    nrInefficientChips += deadAnalyzer->GetNrInefficientChips();

    delete deadAnalyzer;


#ifndef SPD_DA_OFF
    //    daqDA_progressReport((unsigned int)(50+(eq+1)*2.5));
    printf("SPD DA: progress: %d , %ld\n",(unsigned int)(50+(eq+1)*2.5), time(NULL) - timeStamp);
#else
    //    printf("SPD DA: progress: %d\n",(unsigned int)(75+(eq+1)*1.25));
    printf("SPD DA: progress: %d , %ld\n",(unsigned int)(75+(eq+1)*1.25), time(NULL) - timeStamp);
#endif
  }
  // *** *** *** end loop over equipments (eq)



  
  printf("SPD DA: Dead search finished.\n");
  printf("SPD DA: %d single dead pixels , %d dead or inactive pixels in total.\n",handler->GetNrDead(),handler->GetNrSilent());
  printf("SPD DA: %d chips had enough statistics. %d chips are dead. %d chips are inefficient.\n",nrEnoughStatChips,nrDeadChips,nrInefficientChips);
  handler->SetFileLocation(saveDirDead);
  handler->WriteSilentToFilesAlways();
  handler->SetFileLocation(saveDirDeadToFXS);
  handler->WriteSilentToFilesAlways();


  //  printf("SPD DA: Opening id list file\n");
  printf("SPD DA: Opening id list file , %ld\n",time(NULL) - timeStamp);
  TString idsFXSFileName = Form("%s/FXSids_run_%d.txt",saveDirIdsToFXS,runNr);
  ofstream idsFXSfile;
  idsFXSfile.open(idsFXSFileName.Data());


  // store dead+inactive pixels in DB, used as starting point for later runs
#ifndef SPD_DA_OFF
  printf("store previous dead in DB , %ld\n",time(NULL) - timeStamp);
  for (UInt_t eq=0; eq<20; eq++) {
    Int_t storePreviousDead_status = 0;
    TString idpd = Form("spd_previous_dead_%d",eq);
    TString fileName = Form("%s/SPD_Dead_%d.root",saveDirDead,eq);
    storePreviousDead_status = daqDA_DB_storeFile(fileName.Data(),idpd.Data());
    //    printf("daqDA_DB_storeFile(%s,%s)\n",fileName.Data(),idpd.Data());
    if (storePreviousDead_status) {
      printf("SPD DA ERROR: Failed to store dead file %s in daqDetDB: status=%d.\n",idpd.Data(),storePreviousDead_status);
    }
  }
#endif


  // send (dead) reference data for this run to FXS - only if there is no chip in category "needsMoreStat"
  if (nrEnoughStatChips>0 && nrEnoughStatChips+nrDeadChips+nrInefficientChips == 1200) {
    printf("SPD DA: Dead calibration is complete.\n");
    printf("SPD DA: Preparing dead reference data\n");
    // send reference data for dead pixels to FXS
    TString tarFiles = "";
    for (UInt_t eq=0; eq<20; eq++) {
      if (eqActiveDead[eq]) {
	printf("SPD DA: Preparing dead pixels for eq %d\n",eq);
	// move file to ref dir
	TString fileName = Form("%s/SPDphys_dead_run_0_0_eq_%d.root",saveDirDeadRefTmp,eq);
	TString newFileName = Form("%s/SPDphys_dead_run_%d_%d_eq_%d.root",saveDirDeadRef,firstRunNrDead,runNr,eq);
	TString command = Form("mv -f %s %s",fileName.Data(),newFileName.Data());
	system(command.Data());
	tarFiles.Append(Form("SPDphys_dead_run_%d_%d_eq_%d.root ",firstRunNrDead,runNr,eq));
	// create empty hitmap file to send to DB later
	AliITSOnlineSPDphys* physObj = new AliITSOnlineSPDphys(fileName.Data());
	physObj->SetEqNr(eq);
	delete physObj;
      }
    }
    if  (tarFiles.Length() > 0) {   // make sure there are some files to send
      TString send_command = Form("cd %s; tar -cf ref_phys_dead.tar %s",saveDirDeadRef,tarFiles.Data());
      system(send_command.Data());
      TString fileName = Form("%s/ref_phys_dead.tar",saveDirDeadRef);
      TString id = "SPD_ref_phys_dead";
#ifndef SPD_DA_OFF
      status = daqDA_FES_storeFile(fileName.Data(),id.Data());
      if (status!=0) {
	printf("SPD DA ERROR: Failed to export file %s , status %d\n",fileName.Data(),status);
	nrErrors++;
      }
#endif
      idsFXSfile << Form("%s\n",id.Data());
    }
  }

  // store hitmaps (maybe empty) in DB, used as starting point for later runs:
#ifndef SPD_DA_OFF
  printf("store previous hitmaps in DB , %ld\n",time(NULL) - timeStamp);
  for (UInt_t eq=0; eq<20; eq++) {
    Int_t storePreviousHitmap_status = 0;
    TString idph = Form("spd_previous_hitmap_%d",eq);
    TString fileName = Form("%s/SPDphys_dead_run_0_0_eq_%d.root",saveDirDeadRefTmp,eq);
    storePreviousHitmap_status = daqDA_DB_storeFile(fileName.Data(),idph.Data());
    //    printf("daqDA_DB_storeFile(%s,%s)\n",fileName.Data(),idph.Data());
    if (storePreviousHitmap_status) {
      printf("SPD DA ERROR: Failed to store hitmap file %s in daqDetDB: status=%d.\n",idph.Data(),storePreviousHitmap_status);
    }
  }
#endif



  // send (noisy) reference data for this run to FXS
  //  printf("SPD DA: Preparing physics reference data for this run\n");
  printf("SPD DA: Preparing physics reference data for this run , %ld\n",time(NULL) - timeStamp);
  TString tarFiles = "";
  for (UInt_t eq=0; eq<20; eq++) {
    if (eqActiveNoisy[eq]) {
      tarFiles.Append(Form("SPDphys_run_%d_eq_%d.root ",runNr,eq));
    }
  }
  if  (tarFiles.Length() > 0) {   // make sure there are any files to send
    TString send_command = Form("cd %s; tar -cf ref_phys.tar %s",saveDirNoisyRef,tarFiles.Data());
    system(send_command.Data());
    TString fileName = Form("%s/ref_phys.tar",saveDirNoisyRef);
    TString id = "SPD_ref_phys";
#ifndef SPD_DA_OFF
    status = daqDA_FES_storeFile(fileName.Data(),id.Data());
    if (status!=0) {
      printf("SPD DA ERROR: Failed to export file %s , status %d\n",fileName.Data(),status);
      nrErrors++;
    }
#endif
    idsFXSfile << Form("%s\n",id.Data());
  }



  // send dead pixels to FXS
  //  printf("SPD DA: Preparing dead files for FXS\n");
  printf("SPD DA: Preparing dead files for FXS , %ld\n",time(NULL) - timeStamp);
  // send a tared file of all the dead files
  TString send_command = Form("cd %s; tar -cf dead_phys.tar *",saveDirDeadToFXS);
  //  printf("\n\n%s\n\n",command.Data());
  system(send_command.Data());
  TString fileName = Form("%s/dead_phys.tar",saveDirDeadToFXS);
  TString id = "SPD_phys_dead";
#ifndef SPD_DA_OFF
  Int_t send_status = daqDA_FES_storeFile(fileName.Data(),id.Data());
  if (send_status!=0) {
    printf("SPD DA ERROR: Failed to export file %s , status %d\n",fileName.Data(),send_status);
    nrErrors++;
  }
#endif
  idsFXSfile << Form("%s\n",id.Data());



  // send noisy pixels to FXS
  if (nrNoisyFilesProduced > 0) {   // make sure there is at least one file created
    //    printf("SPD DA: Preparing noisy files for FXS\n");
    printf("SPD DA: Preparing noisy files for FXS , %ld\n",time(NULL) - timeStamp);
    // send a tared file of all the noisy files
    TString command = Form("cd %s; tar -cf noisy_phys.tar *",saveDirNoisyToFXS);
    //    printf("\n\n%s\n\n",command.Data());
    system(command.Data());
    TString fileName = Form("%s/noisy_phys.tar",saveDirNoisyToFXS);
    TString id = "SPD_phys_noisy";
#ifndef SPD_DA_OFF
    status = daqDA_FES_storeFile(fileName.Data(),id.Data());
    if (status!=0) {
      printf("SPD DA ERROR: Failed to export file %s , status %d\n",fileName.Data(),status);
      nrErrors++;
    }
#endif
    idsFXSfile << Form("%s\n",id.Data());
  }



  // send ids file to FXS
  idsFXSfile.close();
  id = "SPD_id_list";
#ifndef SPD_DA_OFF
  status = daqDA_FES_storeFile(idsFXSFileName.Data(),id.Data());
  if (status!=0) {
    printf("SPD DA ERROR: Failed to export file %s , status %d\n",idsFXSFileName.Data(),status);
    nrErrors++;
  }
#endif




  delete handler;

  printf("SPD DA: End of all operations, %d errors, %ld\n", nrErrors, time(NULL) - timeStamp);

  if (nrErrors>0) return -1;

  return 0;

}
