/*
- "Contact:" - henrik.tydesjo@cern.ch
- "Link:" - 
- "Run Type:" - PHYSICS
- "DA Type:" - MON
- "Number of events needed:"
- "Input Files:" - daq db config files: spd_physics_params ,  previous dead ref: ./calibResults/DeadReferenceTmp/* , previous dead lists: ./calibResults/DeadToFXS/*
- "Output Files:" - Ref Data: ./calibResults/NoisyReference/* ,  Ref Data: ./calibResults/DeadReference/* ,  noisy lists: ./calibResults/NoisyToFXS/* ,  persistent files: ./calibResults/DeadReferenceTmp/*,./calibResults/DeadToFXS/*
- "Trigger types used:" PHYSICS
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

extern "C" {
#include "daqDA.h"
}
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
#include <TROOT.h>
#include <TPluginManager.h>
#include <TObjArray.h>
#include <TString.h>

int main(int argc, char **argv) {
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // directory structure, hard coded
  char *saveDirDead          = "./calibResults/Dead";             // may NOT delete content
  char *saveDirDeadToFXS     = "./calibResults/DeadToFXS";        //     may delete content
  char *saveDirDeadRef       = "./calibResults/DeadReference";    //     may delete content
  char *saveDirDeadRefTmp    = "./calibResults/DeadReferenceTmp"; // may NOT delete content
  char *saveDirNoisyToFXS    = "./calibResults/NoisyToFXS";       //     may delete content
  char *saveDirNoisyRef      = "./calibResults/NoisyReference";   //     may delete content
  char *saveDirIdsToFXS      = "./calibResults/IdsToFXS";         //     may delete content
  char *configFilesDir       = "./configFiles";                   //     may delete content
  // make sure the directory structure is put up correctly:
  system("mkdir ./calibResults >& /dev/null");
  system("mkdir ./calibResults/Dead >& /dev/null");
  system("mkdir ./calibResults/DeadToFXS >& /dev/null");
  system("mkdir ./calibResults/DeadReference >& /dev/null");
  system("mkdir ./calibResults/DeadReferenceTmp >& /dev/null");
  system("mkdir ./calibResults/NoisyToFXS >& /dev/null");
  system("mkdir ./calibResults/NoisyReference >& /dev/null");
  system("mkdir ./calibResults/IdsToFXS >& /dev/null");
  system("mkdir ./configFiles >& /dev/null");
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
  new AliLog;
  AliLog::Instance()->SetGlobalDebugLevel(-20);


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
    printf("Failed to get config file %s: status=%d. Using default tuning parameters.\n",idp.Data(),par_status);
  }
#endif
  if (par_status==0) {
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
	TString* paramNS = new TString(paramN);
	TString* paramVS = new TString(paramV);
	paramNames.AddAtAndExpand((TObject*)paramNS,nrTuningParams);
	paramVals.AddAtAndExpand((TObject*)paramVS,nrTuningParams);
	nrTuningParams++;
	if (paramsFile.eof()) break;
      }
      paramsFile.close();
    }
  }
  //  for (UInt_t i=0; i<nrTuningParams; i++) {
  //    printf("Entry %d: N=%s , V=%s\n",i,((TString*)paramNames.At(i))->Data(),((TString*)paramVals.At(i))->Data());
  //  }





  // ********* STEP 1: Produce phys container files (Reference Data). ***********************************

#ifndef SPD_DA_OFF
  int runNr = atoi(getenv("DATE_RUN_NUMBER"));
#else
  int runNr = atoi(argv[1]);
  int startSeg = 2;
#endif


  // container objects
  AliITSOnlineSPDphys *physObj[20];
  Bool_t bPhysInit[20];
  for (UInt_t eqId=0; eqId<20; eqId++) {
    physObj[eqId]=NULL;
    bPhysInit[eqId]=kFALSE;
  }


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
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }
    /* declare monitoring program */
    status=monitorDeclareMp("ITS_SPD_PHYS");
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
      if (eventT == PHYSICS_EVENT){

	eventNr++;
	//	printf("eventNr %d\n",eventNr);

	AliRawReader *reader = new AliRawReaderDate((void*)event);
	AliITSRawStreamSPD *str = new AliITSRawStreamSPD(reader);

	//	for (UInt_t eqId=0; eqId<20; eqId++) {
	//	  reader->Reset();
	//	  reader->Select("ITSSPD",eqId,eqId);

	while (str->Next()) {

	  Int_t eqId = reader->GetDDLID();
	  if (eqId>=0 && eqId<20) {
	    if (!bPhysInit[eqId]) {
	      TString fileName = Form("%s/SPDphys_run_%d_eq_%d.root",saveDirNoisyRef,runNr,eqId);
	      physObj[eqId] = new AliITSOnlineSPDphys(fileName.Data());
	      physObj[eqId]->AddRunNr(runNr);
	      physObj[eqId]->SetEqNr(eqId);
	      bPhysInit[eqId]=kTRUE;
	    }
	    
	    UInt_t hs = str->GetHalfStaveNr();
	    UInt_t chip = str->GetChipAddr();
	    physObj[eqId]->IncrementHits(hs,chip,str->GetChipCol(),str->GetChipRow());
	    
	  }
	}
	
	for (UInt_t eqId=0; eqId<20; eqId++) {
	  if (bPhysInit[eqId]) {
	    physObj[eqId]->IncrementNrEvents();
	  }
	}

	//	}

	delete str;
	delete reader;

      }

      /* free resources */
      free(event);

    }
    

#ifdef SPD_DA_OFF
    printf("progress: %d\n",(unsigned int)( ((Float_t)(segNr-startSeg+1))/(argc-startSeg)*50 ));
  }
#endif
  
  // clean up phys objects (also saves them)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    if (physObj[eqId]!=NULL) delete physObj[eqId];
  }





  // ********* STEP 2: Analyze phys container files. ************************************************

  // clear noisyToFXS dir:
  TString command;
  command = Form("cd %s; rm -f *",saveDirNoisyToFXS);
  system(command.Data());
  // clear deadToFXS dir:
  command = Form("cd %s; rm -f *",saveDirDeadToFXS);
  system(command.Data());


  // create calibration handler and read dead from previous calibrations
  AliITSOnlineCalibrationSPDhandler* handler = new AliITSOnlineCalibrationSPDhandler();
  handler->SetFileLocation(saveDirDead);
  handler->ReadDeadFromFiles();


  UInt_t firstRunNrDead = runNr;


  UInt_t nrEnoughStatNoisy = 0;
  UInt_t nrEqActiveNoisy = 0;
  Bool_t eqActiveNoisy[20];

  // *** *** *** start loop over equipments (eq_id)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    eqActiveNoisy[eqId] = kFALSE;

    // create analyzer for this eq
    TString fileName = Form("%s/SPDphys_run_%d_eq_%d.root",saveDirNoisyRef,runNr,eqId);
    AliITSOnlineSPDphysAnalyzer *noisyAnalyzer = new AliITSOnlineSPDphysAnalyzer(fileName.Data(),handler);

    // check data in container
    if (noisyAnalyzer->GetEqNr() != eqId) {
      if (noisyAnalyzer->GetEqNr() != 999) {
	printf("Error: Mismatching EqId in Container data and filename (%d!=%d). Skipping.\n",
	       noisyAnalyzer->GetEqNr(),eqId);
      }
      delete noisyAnalyzer;
      continue;
    }

    nrEqActiveNoisy++;
    eqActiveNoisy[eqId] = kTRUE;

    // configure analyzer with tuning parameters etc:
    for (UInt_t i=0; i<nrTuningParams; i++) {
      noisyAnalyzer->SetParam(((TString*)paramNames.At(i))->Data(),((TString*)paramVals.At(i))->Data());
    }

    printf("SPD phys STEP 2: Noisy search for eq %d\n",eqId);  

    // search for noisy pixels:
    nrEnoughStatNoisy += noisyAnalyzer->ProcessNoisyPixels();

    // copy this phys obj to temporary dead reference dir to process after noisy search
    TString fileNameDead = Form("%s/SPDphys_dead_run_0_0_eq_%d.root",saveDirDeadRefTmp,eqId);
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
	for (UInt_t ind=0; ind<handler->GetNrNoisyC(eqId,hs,chip); ind++) {
	  UInt_t col  = handler->GetNoisyColAtC(eqId,hs,chip,ind);
	  UInt_t row  = handler->GetNoisyRowAtC(eqId,hs,chip,ind);
	  physObj->AddHits(hs,chip,col,row,-noisyAnalyzer->GetOnlinePhys()->GetHits(hs,chip,col,row));
	}
      }
    }

    delete physObj;
    delete noisyAnalyzer;

#ifndef SPD_DA_OFF
    daqDA_progressReport((unsigned int)((eqId+1)*2.5));
#else
    printf("progress: %d\n",(unsigned int)(50+(eqId+1)*1.25));
#endif
  }
  // *** *** *** end loop over equipments (eq_id)

  printf("Noisy search finished. %d noisy pixels found. %d chips (%d) had enough statistics.\n",
	 handler->GetNrNoisy(),nrEnoughStatNoisy,nrEqActiveNoisy*60);
  handler->SetFileLocation(saveDirNoisyToFXS);
  handler->WriteNoisyToFiles();








  UInt_t nrEnoughStatChips = 0;
  UInt_t nrDeadChips = 0;
  UInt_t nrInefficientChips = 0;
  UInt_t nrEqActiveDead = 0;
  Bool_t eqActiveDead[20];

  // *** *** *** start loop over equipments (eq_id)
  for (UInt_t eqId=0; eqId<20; eqId++) {
    eqActiveDead[eqId] = kFALSE;

    // setup analyzer for dead search
    TString fileNameDead = Form("%s/SPDphys_dead_run_0_0_eq_%d.root",saveDirDeadRefTmp,eqId);
    AliITSOnlineSPDphys* physObj = new AliITSOnlineSPDphys(fileNameDead.Data());
    AliITSOnlineSPDphysAnalyzer* deadAnalyzer = new AliITSOnlineSPDphysAnalyzer(physObj,handler);
    // check data in container
    if (deadAnalyzer->GetEqNr() != eqId) {
      if (deadAnalyzer->GetEqNr() != 999) {
	printf("Error: Mismatching EqId in Dead Container data and filename (%d!=%d). Skipping.\n",
	       deadAnalyzer->GetEqNr(),eqId);
      }
      delete deadAnalyzer;
      continue;
    }

    nrEqActiveDead++;
    eqActiveDead[eqId] = kTRUE;

    // configure analyzer with tuning parameters etc:
    for (UInt_t i=0; i<nrTuningParams; i++) {
      deadAnalyzer->SetParam(((TString*)paramNames.At(i))->Data(),((TString*)paramVals.At(i))->Data());
    }

    printf("SPD phys STEP 2: Dead search for eq %d\n",eqId);  

    // search for dead pixels:
    nrEnoughStatChips += deadAnalyzer->ProcessDeadPixels();
    nrDeadChips += deadAnalyzer->GetNrDeadChips();
    nrInefficientChips += deadAnalyzer->GetNrInefficientChips();

    delete deadAnalyzer;

#ifndef SPD_DA_OFF
    daqDA_progressReport((unsigned int)(50+(eqId+1)*2.5));
#else
    printf("progress: %d\n",(unsigned int)(75+(eqId+1)*1.25));
#endif
  }
  // *** *** *** end loop over equipments (eq_id)

  
  printf("Dead search finished. %d dead pixels in total.\n%d chips (%d) had enough statistics. %d chips were dead. %d chips were inefficient.\n",handler->GetNrDead(),nrEnoughStatChips,nrEqActiveDead*60,nrDeadChips,nrInefficientChips);
  handler->SetFileLocation(saveDirDead);
  handler->WriteDeadToFilesAlways();

  handler->SetFileLocation(saveDirDeadToFXS);
  handler->WriteDeadToFilesAlways();
// *** old code (used if not all dead data should be uploaded)
//  UInt_t nrDeadFilesToTar = 0;
//  handler->SetFileLocation(saveDirDeadToFXS);
//  for (UInt_t eqId=0; eqId<20; eqId++) {
//    if (eqActiveDead[eqId]) {
//      handler->WriteDeadToFile(eqId);
//      nrDeadFilesToTar++;
//    }
//  }
//***  

  TString idsFXSFileName = Form("%s/FXSids_run_%d.txt",saveDirIdsToFXS,runNr);
  ofstream idsFXSfile;
  idsFXSfile.open(idsFXSFileName.Data());


  // if there is no chip in category "needsMoreStat"
  if (nrEnoughStatChips+nrDeadChips+nrInefficientChips == nrEqActiveDead*60) {
    // calibration is complete
    printf("Dead calibration is complete.\n");



    // send reference data for dead pixels to FXS
    TString tarFiles = "";
    for (UInt_t eqId=0; eqId<20; eqId++) {
      if (eqActiveDead[eqId]) {
	// move file to ref dir
	TString fileName = Form("%s/SPDphys_dead_run_0_0_eq_%d.root",saveDirDeadRefTmp,eqId);
	TString newFileName = Form("%s/SPDphys_dead_run_%d_%d_eq_%d.root",saveDirDeadRef,firstRunNrDead,runNr,eqId);
	TString command = Form("mv -f %s %s",fileName.Data(),newFileName.Data());
	system(command.Data());

	tarFiles.Append(Form("SPDphys_dead_run_%d_%d_eq_%d.root ",firstRunNrDead,runNr,eqId));
      }
    }
    TString send_command = Form("cd %s; tar -cf ref_phys_dead.tar %s",saveDirDeadRef,tarFiles.Data());
    system(send_command.Data());
    TString fileName = Form("%s/ref_phys_dead.tar",saveDirDeadRef);
    TString id = "SPD_ref_phys_dead";
#ifndef SPD_DA_OFF
    Int_t status = daqDA_FES_storeFile(fileName.Data(),id.Data());
    if (status!=0) {
      printf("Failed to export file %s , status %d\n",fileName.Data(),status);
      return -1;
    }
#endif
    idsFXSfile << Form("%s\n",id.Data());
//    for (UInt_t eqId=0; eqId<20; eqId++) { // OLD CODE NOT TARED
//      if (eqActiveDead[eqId]) {
//	TString fileName = Form("%s/SPDphys_dead_run_0_0_eq_%d.root",saveDirDeadRefTmp,eqId);
//	// move file to ref dir
//	TString newFileName = Form("%s/SPDphys_dead_run_%d_%d_eq_%d.root",saveDirDeadRef,firstRunNrDead,runNr,eqId);
//	TString command = Form("mv -f %s %s",fileName.Data(),newFileName.Data());
//	system(command.Data());
//	// send ref data to FXS
//	TString id = Form("SPD_ref_phys_dead_%d",eqId);
//#ifndef SPD_DA_OFF
//	Int_t status = daqDA_FES_storeFile(newFileName.Data(),id.Data());
//	if (status!=0) {
//	  printf("Failed to export file %s , status %d\n",newFileName.Data(),status);
//	  return -1;
//	}
//#endif
//	idsFXSfile << Form("%s\n",id.Data());
//      }
//      else {
//	TString command = Form("rm -f %s/SPDphys_dead_run_0_0_eq_%d.root",saveDirDeadRefTmp,eqId);
//	system(command.Data());
//      }
//    }
  }


  // send reference data for this run to FXS
  TString tarFiles = "";
  for (UInt_t eqId=0; eqId<20; eqId++) {
    if (eqActiveNoisy[eqId]) {
      tarFiles.Append(Form("SPDphys_run_%d_eq_%d.root ",runNr,eqId));
    }
  }
  TString send_command = Form("cd %s; tar -cf ref_phys.tar %s",saveDirNoisyRef,tarFiles.Data());
  system(send_command.Data());
  TString fileName = Form("%s/ref_phys.tar",saveDirNoisyRef);
  TString id = "SPD_ref_phys";
#ifndef SPD_DA_OFF
  Int_t status = daqDA_FES_storeFile(fileName.Data(),id.Data());
  if (status!=0) {
    printf("Failed to export file %s , status %d\n",fileName.Data(),status);
    return -1;
  }
#endif
  idsFXSfile << Form("%s\n",id.Data());


//  // send reference data for this run to FXS - OLD CODE NOT TARED
//  for (UInt_t eqId=0; eqId<20; eqId++) {
//    if (eqActiveNoisy[eqId]) {
//      TString fileName = Form("%s/SPDphys_run_%d_eq_%d.root",saveDirNoisyRef,runNr,eqId);
//      TString id = Form("SPD_ref_phys_%d",eqId);
//#ifndef SPD_DA_OFF
//      Int_t status = daqDA_FES_storeFile(fileName.Data(),id.Data());
//      if (status!=0) {
//	printf("Failed to export file %s , status %d\n",fileName.Data(),status);
//	return -1;
//      }
//      idsFXSfile << Form("%s\n",id.Data());
//    }
//  }
//#endif


  // send dead pixels to FXS
  //  if (nrDeadFilesToTar>0) { //*** old code (used if not all dead data should be uploaded)
  // send a tared file of all the dead files
  send_command = Form("cd %s; tar -cf dead_phys.tar *",saveDirDeadToFXS);
  //  printf("\n\n%s\n\n",command.Data());
  system(send_command.Data());
  fileName = Form("%s/dead_phys.tar",saveDirDeadToFXS);
  id = "SPD_phys_dead";
#ifndef SPD_DA_OFF
  Int_t send_status = daqDA_FES_storeFile(fileName.Data(),id.Data());
  if (send_status!=0) {
    printf("Failed to export file %s , status %d\n",fileName.Data(),send_status);
    return -1;
  }
#endif
  idsFXSfile << Form("%s\n",id.Data());
  //  } //*** old code (used if not all dead data should be uploaded)


  // send noisy pixels to FXS
  if (handler->GetNrNoisy()>0) { // there must be at least one file created
    // send a tared file of all the noisy files
    TString command = Form("cd %s; tar -cf noisy_phys.tar *",saveDirNoisyToFXS);
    //    printf("\n\n%s\n\n",command.Data());
    system(command.Data());
    TString fileName = Form("%s/noisy_phys.tar",saveDirNoisyToFXS);
    TString id = "SPD_phys_noisy";
#ifndef SPD_DA_OFF
    Int_t status = daqDA_FES_storeFile(fileName.Data(),id.Data());
    if (status!=0) {
      printf("Failed to export file %s , status %d\n",fileName.Data(),status);
      return -1;
    }
#endif
    idsFXSfile << Form("%s\n",id.Data());
  }


  // send ids file to FXS
  idsFXSfile.close();
  id = "SPD_id_list";
#ifndef SPD_DA_OFF
  Int_t status = daqDA_FES_storeFile(idsFXSFileName.Data(),id.Data());
  if (status!=0) {
    printf("Failed to export file %s , status %d\n",fileName.Data(),status);
    return -1;
  }
#endif




  delete handler;


  return 0;
}
