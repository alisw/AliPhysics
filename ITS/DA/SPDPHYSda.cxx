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
#include <time.h>
#include <TROOT.h>
#include <TPluginManager.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <cstdlib>


// directory structure, hard coded                              //   FOR OFFLINE VERSION:
char *saveDirDeadToFXS     = "./calibResults/DeadToFXS";        //     may delete content
char *saveDirNoisyToFXS    = "./calibResults/NoisyToFXS";       //     may delete content
char *saveDirNoisyRef      = "./calibResults/NoisyReference";   //     may delete content
char *configFilesDir       = "./configFiles";                   //     may delete content
char *saveDirIdsToFXS      = "./calibResults/IdsToFXS";         //     may delete content


// container objects
AliITSOnlineSPDphys *physObj[20];
AliITSOnlineSPDphys *physObjPartial[20];

TString idsFXSFileName;
ofstream idsFXSfile;
Int_t nrErrors=0;

// ********* STEP 0: Get configuration files from db (if there are any) , then read parameters*********
TObjArray paramNames; 
TObjArray paramVals; 
UInt_t nrTuningParams = 0;

AliITSOnlineCalibrationSPDhandler *handler = NULL;


//Int_t AnalyzeNoisy(AliITSOnlineSPDphys *physobj, UInt_t eq);
//Int_t AnalyzeDead(AliITSOnlineSPDphys *physobj,UInt_t eq);
Int_t AnalyzeNoisy(TString filename, UInt_t eq);
Int_t AnalyzeDead(TString filename,UInt_t eq);
void SendRefPhys(TString tarfiles);
void SendNoisy(TString tarfiles);
void SendDead(TString tarfiles);

int main(int argc, char **argv) {
	if (argc<2) {
		printf("SPD DA ERROR: Wrong number of arguments\n");
		return -1;
	}

	UInt_t nrErrors = 0;

	// clean up and make sure the directory structure is put up correctly:
	system("rm ./calibResults -rf >& /dev/null");
	system("mkdir ./calibResults >& /dev/null");
	system("mkdir ./calibResults/DeadToFXS >& /dev/null");
	system("mkdir ./calibResults/NoisyToFXS >& /dev/null");
	system("mkdir ./calibResults/NoisyReference >& /dev/null"); 
	system("mkdir ./configFiles >& /dev/null");
	system("mkdir ./calibResults/IdsToFXS >& /dev/null");

	// parameters config file
	TString paramsFileName = Form("%s/physics_params.txt",configFilesDir);
	paramNames.SetOwner();
	paramVals.SetOwner();

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


	// tuning parameters:
	Int_t par_status = 0;
	TString idp = "spd_physics_params";
	par_status=daqDA_DB_getFile(idp.Data(),paramsFileName.Data());
	if (par_status) {
		printf("SPD DA: Failed to get config file %s: status=%d. Using default tuning parameters.\n",idp.Data(),par_status);
		TString rmCmd = Form("rm -f %s",paramsFileName.Data());
		system(rmCmd.Data());
	}
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
				if (paramsFile.eof()) {
					break;
				}

				paramsFile >> paramV;
				paramNames.AddAtAndExpand(new TObjString(paramN),nrTuningParams);
				paramVals.AddAtAndExpand(new TObjString(paramV),nrTuningParams);
				nrTuningParams++;
				if (paramsFile.eof()) {
					break;
				}
			}
			paramsFile.close();
		}
	}


					// last param is the number of events for cyclic publishing 
					// nrTuningParams is the number for params for the analyzer classes
					// decrease of 1 unit is necessary
					nrTuningParams--;

	// create calibration handler (needed already at step 1 in order to fill which eq,hs,chips are active
	if(!handler) handler = new AliITSOnlineCalibrationSPDhandler();

	// ********* STEP 1: Produce phys container files (Reference Data). ***********************************

	if (getenv("DATE_RUN_NUMBER")==0) {
		printf("SPD DA ERROR: DATE_RUN_NUMBER not properly set.\n");
		return -1;
	}
	int runNr = atoi(getenv("DATE_RUN_NUMBER"));


	// container objects
	//AliITSOnlineSPDphys *physObjPartial[20];
	//Initialize the handler (deactivate everything first, then only active parts will be active)
	Bool_t backUpCalibInfo[20][6][10]; // needed in case few events are collected
	Bool_t bPhysInit[20];
	for (UInt_t eq=0; eq<20; eq++) {
		physObj[eq]=NULL;
		physObjPartial[eq]=NULL;
		bPhysInit[eq]=kFALSE;
		handler->ActivateEq(eq,kFALSE);
		for(Int_t hs=0; hs<6; hs++){
			handler->ActivateHS(eq,hs,kFALSE);
			for(Int_t chip=0; chip<10; chip++){
				handler->ActivateChip(eq,hs,chip,kFALSE);
				backUpCalibInfo[eq][hs][chip]=kFALSE;
			}
		}
	}


	UInt_t eventNr=0;
        TString lastPar = ((TObjString*)paramVals.At(nrTuningParams))->GetString();

	Int_t part = ((TObjString*)paramVals.At(nrTuningParams))->GetString().Atoi();
	printf("*** n Events for cyclic publishing in FXS %i ***\n",part);

	int status;

	/* define data source : */  
	status=monitorSetDataSource( argv[1] );
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

	idsFXSFileName = Form("%s/FXSids_run_%d.txt",saveDirIdsToFXS,runNr);			  

	/* main loop (infinite) */
	for(;;) {

		struct eventHeaderStruct *event; 
		eventTypeType eventT;

		/* check shutdown condition */
		if (daqDA_checkShutdown()) {break;} 

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
				if (eq>=0 && eq<20) {
					if (!bPhysInit[eq]) { // this code is duplicated for the moment... (see also below)
						TString fileName = Form("%s/SPDphys_run_%d_eq_%d.root",saveDirNoisyRef,runNr,eq);
						physObj[eq] = new AliITSOnlineSPDphys(fileName.Data());
						physObj[eq]->AddRunNr(runNr);
						physObj[eq]->SetEqNr(eq);
						bPhysInit[eq]=kTRUE;
					}
					if(str->IsActiveEq(eq)) handler->ActivateEq(eq,kTRUE);
					UInt_t hs = str->GetHalfStaveNr();
					handler->ActivateHS(eq,hs);
					UInt_t chip = str->GetChipAddr();
					handler->ActivateChip(eq,hs,chip);

					physObj[eq]->IncrementHits(hs,chip,str->GetChipCol(),str->GetChipRow());	    
				}
			}




			for (UInt_t eq=0; eq<20; eq++) {
				if(bPhysInit[eq]){
					physObj[eq]->IncrementNrEvents();
				}
			}


			if(eventNr==0){
				// in case few events are taken the handler 
				// is made aware of the hs included in daq
				for (UInt_t eq=0; eq<20; eq++) {
					for (UInt_t hs=0; hs<6; hs++) {
						for (UInt_t chip=0; chip<10; chip++) if(str->IsActiveChip(eq,hs,chip)) backUpCalibInfo[eq][hs][chip]=kTRUE;
					}

				}   
			}

			else if(eventNr%part==0){
				nrErrors=0;
				TString tarFiles="";
				TString tarFilesDead="";
				TString tarFilesNoisy="";
				idsFXSfile.open(idsFXSFileName.Data());
				idsFXSfile.flush();
   
	                       printf("publishing cycle %i : eventNr %i, part %i \n",eventNr/part,eventNr,part);

				for (UInt_t eq=0; eq<20; eq++) {
					if(bPhysInit[eq]){
						TString fileNamePart = Form("%s/SPDphys_run_%d_eq_%d_part_%d.root",saveDirNoisyRef,runNr,eq,eventNr/part);
						physObjPartial[eq] = new AliITSOnlineSPDphys(fileNamePart.Data());
						physObjPartial[eq]->AddPhys(physObj[eq]);
						delete physObjPartial[eq];
						tarFiles.Append(Form("SPDphys_run_%d_eq_%d_part_%d.root ",runNr,eq,eventNr/part));

						//noisy search
						AnalyzeNoisy(fileNamePart.Data(),eq);
						handler->SetFileLocation(saveDirNoisyToFXS);
						if(handler->GetNrNoisyEq(eq)>0){
							handler->WriteNoisyToFile(eq);
							tarFilesNoisy.Append(Form("SPD_Noisy_%d.root ",eq));
						}


						//dead search                                                                                 
						AnalyzeDead(fileNamePart.Data(),eq);
						/*handler->SetFileLocation(saveDirDeadToFXS);
						  handler->WriteSilentToFile(eq);
						  tarFilesDead.Append(Form("SPD_Dead_%d.root ",eq));
						  */
					}
					tarFilesDead.Append(Form("SPD_Dead_%d.root ",eq));
				}

				/*printf("Saving noisy\n");
				  handler->SetFileLocation(saveDirNoisyToFXS);
				  printf("File Location set\n");
				  handler->WriteNoisyToFiles();
				  printf("Saving dead\n");*/
				handler->SetFileLocation(saveDirDeadToFXS);
				//printf("File location set\n");
				handler->WriteSilentToFilesAlways();
				//handler->PrintEqSummary();
				//handler->PrintDead();
				//handler->PrintNoisy();

				//send partial outputs
				SendRefPhys(tarFiles.Data());
				SendDead(tarFilesDead.Data());
				SendNoisy(tarFilesNoisy.Data());
				// send ids file to FXS
				idsFXSfile.close();
				TString id = "SPD_id_list";
				Int_t send_status = daqDA_FES_storeFile(idsFXSFileName.Data(),id.Data());
				if (send_status!=0) {
					printf("SPD DA ERROR: Failed to export file %s , status %d\n",idsFXSFileName.Data(),send_status);
					nrErrors++;
				}
				Printf("Sending data for partition %d. N. of errors = %d\n",eventNr/part,nrErrors);
			}


			delete str;
			delete reader;


			eventNr++;
			//if(eventNr==5000) break;
		}

		/* free resources */
		free(event);

	}

	// In case of few number of collected events
	// the configuration is the one provided by the dcs
	if(eventNr<part){
		for(Int_t eq=0; eq<20; eq++) {
			Int_t countHsInEq=0;
			for(Int_t hs=0; hs<6; hs++) {
				Int_t countChipInHs=0; 
				for(Int_t chip=0; chip<10; chip++) {
					if(backUpCalibInfo[eq][hs][chip]) {
						handler->ActivateChip(eq,hs,chip);
						countChipInHs++;
					}
					if(countChipInHs==10) {
						handler->ActivateHS(eq,hs);
						countHsInEq++;
					}
				}
			}
			if(countHsInEq==6) handler->ActivateEq(eq);
		}
	}
	// clean up phys objects (also saves them)
	printf("\n");
	//Printf("SPD DA: global phys obj analysis\n");
	nrErrors=0;
	TString tarFiles="";
	TString tarFilesDead="";
	TString tarFilesNoisy="";
	idsFXSFileName = Form("%s/FXSids_run_%d.txt",saveDirIdsToFXS,runNr);
	idsFXSfile.open(idsFXSFileName.Data());
	idsFXSfile.flush();
	for (UInt_t eq=0; eq<20; eq++) {
		if (physObj[eq]!=NULL){
			TString fileName = Form("%s/SPDphys_run_%d_eq_%d.root",saveDirNoisyRef,runNr,eq);
			//Printf("NEVENTS = %d \t %d",physObj[eq]->GetNrEvents(), eventNr);
			delete physObj[eq];
			tarFiles.Append(Form("SPDphys_run_%d_eq_%d.root ",runNr,eq));
			//noisy search                                                                                                                   
			AnalyzeNoisy(fileName.Data(),eq);
			handler->SetFileLocation(saveDirNoisyToFXS);
			if(handler->GetNrNoisyEq(eq)>0){
				handler->WriteNoisyToFile(eq);
				tarFilesNoisy.Append(Form("SPD_Noisy_%d.root ",eq));
			}
			//dead search                                                                                                                    
			AnalyzeDead(fileName.Data(),eq);
			/*handler->SetFileLocation(saveDirDeadToFXS);
			  handler->WriteSilentToFile(eq);
			  tarFilesDead.Append(Form("SPD_Dead_%d.root ",eq));
			  */
		}
		tarFilesDead.Append(Form("SPD_Dead_%d.root ",eq));
	}
	/*handler->SetFileLocation(saveDirNoisyToFXS);
	  handler->WriteNoisyToFiles();
	  */
	handler->SetFileLocation(saveDirDeadToFXS);
	handler->WriteSilentToFilesAlways();
	//handler->PrintEqSummary();
	// handler->PrintDead();
	//handler->PrintNoisy();

	//send complete outputs
	SendRefPhys(tarFiles.Data());
	SendDead(tarFilesDead.Data());
	SendNoisy(tarFilesNoisy.Data());
	// send ids file to FXS
	idsFXSfile.close();
	TString id = "SPD_id_list";
	Int_t send_status = daqDA_FES_storeFile(idsFXSFileName.Data(),id.Data());
	if (send_status!=0) {
		printf("SPD DA ERROR: Failed to export file %s , status %d\n",idsFXSFileName.Data(),send_status);
		nrErrors++;
	}
	Printf("Sending complete data. N. errors =  %d\n",nrErrors);
	printf("SPD DA: %d events collected for this run. \n",eventNr);
	printf("SPD DA: %d  dead pixels \n",handler->GetNrBad()-handler->GetNrNoisy());
	printf("SPD DA: %d noisy pixels \n",handler->GetNrNoisy());


	return 0;
}
//_______________________________________________________________________________________________________________________________________
Int_t AnalyzeNoisy(TString filename, UInt_t eq){

//	printf("Analyze noisy\n");
//	Noisy pixels should be updated
         handler->ResetNoisyForEq(eq);
	AliITSOnlineSPDphysAnalyzer *noisyAnalyzer = new AliITSOnlineSPDphysAnalyzer(filename.Data(),handler);
	if(!noisyAnalyzer) {printf("SPD DA: noisy analizer not initialized equipment %i\n",eq); return 1;}


	// configure analyzer with tuning parameters etc:
	//printf("Configuring noisyAnalyzer\n");
	for (UInt_t i=0; i<nrTuningParams; i++) {
		noisyAnalyzer->SetParam(((TObjString*)paramNames.At(i))->GetString().Data(),((TObjString*)paramVals.At(i))->GetString().Data());
	}

	//printf("SPD DA: SPD phys STEP 2: Noisy search for eq %d\n",noisyAnalyzer->GetEqNr());

	// search for noisy pixels:
	noisyAnalyzer->ProcessNoisyPixelsFast();
	//printf("SPD DA: Noisy pixels processed\n");
	if(noisyAnalyzer!=NULL) delete noisyAnalyzer;
	//printf("SPD DA: noisyAnalyzer deleted\n");

	return 0;
}
//_______________________________________________________________________________________________________________________________________
Int_t AnalyzeDead(TString filename, UInt_t eq){

	//printf("Analyze dead\n");

	AliITSOnlineSPDphysAnalyzer *deadAnalyzer = new AliITSOnlineSPDphysAnalyzer(filename.Data(),handler);

	// configure analyzer with tuning parameters etc:
	for (UInt_t i=0; i<nrTuningParams; i++) {
		deadAnalyzer->SetParam(((TObjString*)paramNames.At(i))->GetString().Data(),((TObjString*)paramVals.At(i))->GetString().Data());
	}

	//printf("SPD DA: SPD phys STEP 2: Dead search for eq %d\n",deadAnalyzer->GetEqNr());

	// search for dead pixels:
	deadAnalyzer->ProcessDeadPixels();
	//Printf("Number of dead pixels =  %d",deadAnalyzer->ProcessDead());

	//printf("SPD DA: Dead pixels processed\n");
	if(deadAnalyzer!=NULL) delete deadAnalyzer;
	//printf("SPD DA: deadAnalyzer deleted\n\n");
	return 0;
}
//_______________________________________________________________________________________________________________________________________
void SendRefPhys(TString tarfiles){

	if  (tarfiles.Length() > 0) {   // make sure there are any files to send
		TString send_command = Form("cd %s; tar -cf ref_phys.tar %s",saveDirNoisyRef,tarfiles.Data());
		system(send_command.Data());
		TString fileName = Form("%s/ref_phys.tar",saveDirNoisyRef);
		TString id = "SPD_ref_phys";
		Int_t send_status = daqDA_FES_storeFile(fileName.Data(),id.Data());
		if (send_status!=0) {
			printf("SPD DA ERROR: Failed to export file %s , status %d\n",fileName.Data(),send_status);
			nrErrors++;
		}
		idsFXSfile << Form("%s\n",id.Data());
	}

}
//_______________________________________________________________________________________________________________________________________
void SendDead(TString tarfiles){

	if  (tarfiles.Length() > 0) {   // make sure there are any files to send
		TString send_command = Form("cd %s; tar -cf dead_phys.tar %s",saveDirDeadToFXS,tarfiles.Data());
		system(send_command.Data());
		TString fileName = Form("%s/dead_phys.tar",saveDirDeadToFXS);
		TString id = "SPD_phys_dead";
		Int_t send_status = daqDA_FES_storeFile(fileName.Data(),id.Data());
		if (send_status!=0) {
			printf("SPD DA ERROR: Failed to export file %s , status %d\n",fileName.Data(),send_status);
			nrErrors++;
		}
		idsFXSfile << Form("%s\n",id.Data());
	}
}
//_______________________________________________________________________________________________________________________________________
void SendNoisy(TString tarfiles){

	if  (tarfiles.Length() > 0) {   // make sure there are any files to send
		TString command = Form("cd %s; tar -cf noisy_phys.tar %s",saveDirNoisyToFXS,tarfiles.Data());
		system(command.Data());
		TString fileName = Form("%s/noisy_phys.tar",saveDirNoisyToFXS);
		TString id = "SPD_phys_noisy";
		Int_t send_status = daqDA_FES_storeFile(fileName.Data(),id.Data());
		if (send_status!=0) {
			printf("SPD DA ERROR: Failed to export file %s , status %d\n",fileName.Data(),send_status);
			nrErrors++;
		}
		idsFXSfile << Form("%s\n",id.Data());
	}
}


