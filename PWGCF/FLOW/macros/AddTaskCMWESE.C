#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#endif

AliAnalysisTaskCMWESE* AddTaskCMWESE(
      int                     debug=0, // debug level controls amount of output statements
      double              Harmonic=2,
      TString             trigger="kINT7",
      float                  chi2lo=0.1,
      float                  ptmin=0.2, // minimum pt for Q-vector components
      float                  ptmax=2.0, // maximum pt for Q-vector components
      int                     cbinlo=0, // lower centrality bin for histogram array
      int                     cbinhg=8, // higher centrality bin for histogram array
      double              etaGap=0.3,  
      bool                  v0calibOn=true,
      bool                  doNUE=true,
      bool                  doNUA=true,
      float                  centcut=7.5, // centrality restriction for V0M and TRK
      TString             period="LHC15o",
      TString	       uniqueID=""
      )	
{	
	// Creates a pid task and adds it to the analysis manager
	// Get the pointer to the existing analysis manager via the static
	//access method
	//=========================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskCMWESE.C", "No analysis manager to connect to.");
		return NULL;
	}

	// Check the analysis type using the event handlers connected to the
	// analysis manager The availability of MC handler can also be
	// checked here.
	// =========================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTaskCMWESE.C", "This task requires an input event handler");
		return NULL;
	}
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  	// --- instantiate analysis task
  	AliAnalysisTaskCMWESE *task = new AliAnalysisTaskCMWESE("TaskCMWESE", period, doNUE, doNUA, v0calibOn);
	task->SetDebug(debug);
	task-> SetHarmonic(Harmonic);
	task->SetTrigger(trigger);
  	task->SetChi2Low(chi2lo);
	task->SetPtMin(ptmin);
	task->SetPtMax(ptmax);
	task->SetCentBinLow(cbinlo);
	task->SetCentBinHigh(cbinhg);
	task->SetEtaGap(etaGap);
	task->SetV0CalibOn(v0calibOn);
	task->SetNUEOn(doNUE);
	task->SetNUAOn(doNUA);	
	task->SetCentCut(centcut);
	// task->SelectCollisionCandidates(AliVEvent::kINT7);
	mgr->AddTask(task);

	// Create ONLY the output containers for the data produced by the
	// task.  Get and connect other common input/output containers via
	// the manager as below
	//======================================================================
    	AliAnalysisDataContainer* cinput  = mgr->GetCommonInputContainer();
    	const char* outputFileName = mgr->GetCommonFileName();
  	AliAnalysisDataContainer* coutput = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), 
                                                           AliAnalysisManager::kOutputContainer,                                                          
                                                           Form("%s:%s", outputFileName, uniqueID.Data()));
   	mgr->ConnectInput (task, 0, cinput);
  	mgr->ConnectOutput(task, 1, coutput);

	Int_t inSlotCounter=1;
	TGrid::Connect("alien://");
        	TObjArray *AllContainers = mgr->GetContainers();

	if(task->GetNUEOn() || doNUE) {
                	if (period.EqualTo("LHC10h") || period.EqualTo("LHC11h")) {
                		TFile *inNUE;
			if(!AllContainers->FindObject("NUE")) {
				inNUE = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/Run1NUE.root");
				AliAnalysisDataContainer *cin_NUE = mgr->CreateContainer(Form("NUE"), TList::Class(), AliAnalysisManager::kInputContainer);             
				TList* wNUE_list = NULL;
				wNUE_list = dynamic_cast<TList*>(inNUE->Get("listNUE"));
				if (!wNUE_list) printf("Read TList wrong!\n");
            		    	cin_NUE->SetData(wNUE_list); 			
				mgr->ConnectInput(task,inSlotCounter,cin_NUE);
				inSlotCounter++;
			}
			else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUE"));
				inSlotCounter++;
				printf("NUE already loaded\n");
			}
		}	
		else if (period.EqualTo("LHC15o")) {
			TFile *inNUE;
			if(!AllContainers->FindObject("NUE")) {
				inNUE = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/efficiencyBothpol.root");
				AliAnalysisDataContainer *cin_NUE = mgr->CreateContainer(Form("NUE"), TList::Class(), AliAnalysisManager::kInputContainer);
				TList* wNUE_list = NULL;
				wNUE_list = dynamic_cast<TList*>(inNUE->Get("fMcEffiHij"));
				if (!wNUE_list) printf("Read TList wrong!\n");
            		    	cin_NUE->SetData(wNUE_list); 			
				mgr->ConnectInput(task,inSlotCounter,cin_NUE);
				inSlotCounter++;				
			}
			else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUE"));
				inSlotCounter++;
				printf("NUE already loaded\n");
			}
		} 
		else if (period.EqualTo("LHC18q") || period.EqualTo("LHC18r")) {
			TFile *inNUE;
			if(!AllContainers->FindObject("NUE")) {
				inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mhaque/calib2021/efficiencyBothpol18qnew.root");
				AliAnalysisDataContainer *cin_NUE = mgr->CreateContainer(Form("NUE"), TList::Class(), AliAnalysisManager::kInputContainer);
				TList* wNUE_list = NULL;
				wNUE_list = dynamic_cast<TList*>(inNUE->FindObjectAny("fMcEffiHij"));
				if (!wNUE_list) printf("Read TList wrong!\n");
            		    	cin_NUE->SetData(wNUE_list); 			
				mgr->ConnectInput(task,inSlotCounter,cin_NUE);
				inSlotCounter++;				
			}
			else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUE"));
				inSlotCounter++;
				printf("NUE already loaded\n");
			}
		} 
	}



	TString filenameNUA = "";
	if (period.EqualTo("LHC10h")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB1.root";
	if (uniqueID.EqualTo("10h_ChiHg3") && period.EqualTo("LHC10h")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB1ChiHg3.root";
	if (uniqueID.EqualTo("10h_Nhits60") && period.EqualTo("LHC10h")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB1Nhits60.root";
	if (uniqueID.EqualTo("10h_Nhits80") && period.EqualTo("LHC10h")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB1Nhits80.root";
	if (uniqueID.EqualTo("10h_FB768")   && period.EqualTo("LHC10h")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB768.root";
	if (uniqueID.EqualTo("10h_FB272")   && period.EqualTo("LHC10h")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB272.root";
	if (uniqueID.EqualTo("10h_FB768_ChiHg2")   && period.EqualTo("LHC10h")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB768_ChiHg2.root";
	if (uniqueID.EqualTo("10h_FB768_ChiHg3")   && period.EqualTo("LHC10h")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB768_ChiHg3.root";
	if (uniqueID.EqualTo("10h_FB768_Nhits80")   && period.EqualTo("LHC10h")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB768_Nhits80.root";
	if (uniqueID.EqualTo("10h_FB768_Nhits60")   && period.EqualTo("LHC10h")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB768_Nhits60.root";

	if (period.EqualTo("LHC15o")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/wgtPion_NUAFB768DeftwPUcut_LHC15op2_24Aug2021.root";
	if (uniqueID.EqualTo("15o_ChiHg3")  && period.EqualTo("LHC15o")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/LHC15o_pass2_NUA_ChiHg3.root";
	if (uniqueID.EqualTo("15o_Nhits60")  && period.EqualTo("LHC15o")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/LHC15o_pass2_NUA_Nhits60.root";
	if (uniqueID.EqualTo("15o_Nhits80")  && period.EqualTo("LHC15o")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/LHC15o_pass2_NUA_Nhits80.root";
	if (uniqueID.EqualTo("15o_FB96")      && period.EqualTo("LHC15o")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/LHC15o_pass2_NUA_FB96.root";	
	if (uniqueID.EqualTo("15o_NUA_fromEmil")  && period.EqualTo("LHC15o")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/LHC15o_pass2_NUA_wSyst.root";
	if (uniqueID.EqualTo("15o_ChiHg2")  && period.EqualTo("LHC15o")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/LHC15o_pass2_NUA_ChiHg2.root";

	if (period.EqualTo("LHC18q")) filenameNUA = "alien:///alice/cern.ch/user/m/mhaque/calib2021/WgtsNUAChargeAndPion_LHC18qPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root";
	if (uniqueID.EqualTo("18q_ChiHg2")     && period.EqualTo("LHC18q")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_ChiHg2.root";
	if (uniqueID.EqualTo("18q_ChiHg2d5") && period.EqualTo("LHC18q")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_ChiHg2d5.root";
	if (uniqueID.EqualTo("18q_Nhits60")     && period.EqualTo("LHC18q")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_Nhits60.root";
	if (uniqueID.EqualTo("18q_Nhits80")     && period.EqualTo("LHC18q")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_Nhits80.root";
	if (uniqueID.EqualTo("18q_Dcaz2")       && period.EqualTo("LHC18q")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_Dcaz2.root";
	if (uniqueID.EqualTo("18q_Dcaz2d5")   && period.EqualTo("LHC18q")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_Dcaz2d5.root";

	if (period.EqualTo("LHC18r")) filenameNUA = "alien:///alice/cern.ch/user/m/mhaque/calib2021/WgtsNUAChargeAndPion_LHC18rPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root";
	if (uniqueID.EqualTo("18r_ChiHg2")     && period.EqualTo("LHC18r")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_ChiHg2.root";
	if (uniqueID.EqualTo("18r_ChiHg2d5") && period.EqualTo("LHC18r")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_ChiHg2d5.root";
	if (uniqueID.EqualTo("18r_Nhits60")     && period.EqualTo("LHC18r")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_Nhits60.root";
	if (uniqueID.EqualTo("18r_Nhits80")     && period.EqualTo("LHC18r")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_Nhits80.root";
	if (uniqueID.EqualTo("18r_Dcaz2")       && period.EqualTo("LHC18r")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_Dcaz2.root";
	if (uniqueID.EqualTo("18r_Dcaz2d5")   && period.EqualTo("LHC18r")) filenameNUA = "alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_Dcaz2d5.root";

	if(task->GetNUAOn() ||doNUA) {
		if (period.EqualTo("LHC10h") ) { // NUA for 10h is too large to read, we separate them into 3 TList*s.
			TFile *inNUA;
			inNUA = TFile::Open(filenameNUA);
			if (!inNUA) return task;
			AliAnalysisDataContainer *cin_NUA;
			if (period.EqualTo("LHC10h")) cin_NUA = mgr->CreateContainer(Form("NUA"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("10h_ChiHg3") && period.EqualTo("LHC10h"))  cin_NUA = mgr->CreateContainer(Form("NUA_ChiHg3"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("10h_Nhits60") && period.EqualTo("LHC10h"))  cin_NUA = mgr->CreateContainer(Form("NUA_Nhits60"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("10h_Nhits80") && period.EqualTo("LHC10h"))  cin_NUA = mgr->CreateContainer(Form("NUA_Nhits80"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("10h_FB768")   && period.EqualTo("LHC10h"))  cin_NUA = mgr->CreateContainer(Form("NUA_FB768"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("10h_FB272")   && period.EqualTo("LHC10h"))  cin_NUA = mgr->CreateContainer(Form("NUA_FB272"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("10h_FB768_ChiHg2")   && period.EqualTo("LHC10h")) cin_NUA = mgr->CreateContainer(Form("NUA_FB768_ChiHg2"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("10h_FB768_ChiHg3")   && period.EqualTo("LHC10h")) cin_NUA = mgr->CreateContainer(Form("NUA_FB768_ChiHg3"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("10h_FB768_Nhits80")   && period.EqualTo("LHC10h")) cin_NUA = mgr->CreateContainer(Form("NUA_FB768_Nhits80"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("10h_FB768_Nhits60")   && period.EqualTo("LHC10h")) cin_NUA = mgr->CreateContainer(Form("NUA_FB768_Nhits60"), TList::Class(), AliAnalysisManager::kInputContainer);

	                	TList* wNUA_list = NULL;
			wNUA_list = dynamic_cast<TList*>(inNUA->Get("10hListNUA"));
		            cin_NUA->SetData(wNUA_list); 
			mgr->ConnectInput(task,inSlotCounter,cin_NUA);
			inSlotCounter++;
		} 
		else if (period.EqualTo("LHC15o")) {
			TFile *inNUA;
			inNUA = TFile::Open(filenameNUA);
			if (!inNUA) return task;
			AliAnalysisDataContainer *cin_NUA;
			if (period.EqualTo("LHC15o")) cin_NUA = mgr->CreateContainer(Form("NUA"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("15o_ChiHg3")  && period.EqualTo("LHC15o")) cin_NUA = mgr->CreateContainer(Form("NUA_ChiHg3"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("15o_Nhits60")  && period.EqualTo("LHC15o")) cin_NUA = mgr->CreateContainer(Form("NUA_Nhits60"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("15o_Nhits80")  && period.EqualTo("LHC15o")) cin_NUA = mgr->CreateContainer(Form("NUA_Nhits80"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("15o_FB96")      && period.EqualTo("LHC15o")) cin_NUA = mgr->CreateContainer(Form("NUA_FB96"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("15o_NUA_fromEmil")  && period.EqualTo("LHC15o")) cin_NUA = mgr->CreateContainer(Form("NUA_FromEmil"), TList::Class(), AliAnalysisManager::kInputContainer);
			if (uniqueID.EqualTo("15o_ChiHg2")  && period.EqualTo("LHC15o")) cin_NUA = mgr->CreateContainer(Form("NUA_ChiHg2"), TList::Class(), AliAnalysisManager::kInputContainer);

			TList* wNUA_list = NULL;
			wNUA_list = dynamic_cast<TList*>(inNUA->Get("15oListNUA"));
		            cin_NUA->SetData(wNUA_list); 
			mgr->ConnectInput(task,inSlotCounter,cin_NUA);
			inSlotCounter++;			
		}
		else if (period.EqualTo("LHC18q") || period.EqualTo("LHC18r")) {
			TFile *inNUA;
			inNUA = TFile::Open(filenameNUA);
			if (!inNUA) return task;

			if (	uniqueID.EqualTo("18q_ChiHg2") || uniqueID.EqualTo("18q_ChiHg2d5") || uniqueID.EqualTo("18q_Nhits60") || uniqueID.EqualTo("18q_Nhits80") || uniqueID.EqualTo("18q_Dcaz2") || uniqueID.EqualTo("18q_Dcaz2d5") 
				|| uniqueID.EqualTo("18r_ChiHg2") || uniqueID.EqualTo("18r_ChiHg2d5") || uniqueID.EqualTo("18r_Nhits60") || uniqueID.EqualTo("18r_Nhits80") || uniqueID.EqualTo("18r_Dcaz2") || uniqueID.EqualTo("18r_Dcaz2d5") 	
				// The ID of track-wise syst is unique!!			
			)  {
				AliAnalysisDataContainer *cin_NUA;
				if (uniqueID.EqualTo("18q_ChiHg2")     && period.EqualTo("LHC18q")) cin_NUA = mgr->CreateContainer(Form("NUA_ChiHg2"), TList::Class(), AliAnalysisManager::kInputContainer);
				if (uniqueID.EqualTo("18q_ChiHg2d5") && period.EqualTo("LHC18q")) cin_NUA = mgr->CreateContainer(Form("NUA_ChiHg2d5"), TList::Class(), AliAnalysisManager::kInputContainer);
				if (uniqueID.EqualTo("18q_Nhits60")     && period.EqualTo("LHC18q")) cin_NUA = mgr->CreateContainer(Form("NUA_Nhits60"), TList::Class(), AliAnalysisManager::kInputContainer);
				if (uniqueID.EqualTo("18q_Nhits80")     && period.EqualTo("LHC18q")) cin_NUA = mgr->CreateContainer(Form("NUA_Nhits80"), TList::Class(), AliAnalysisManager::kInputContainer);
				if (uniqueID.EqualTo("18q_Dcaz2")       && period.EqualTo("LHC18q")) cin_NUA = mgr->CreateContainer(Form("NUA_Dcaz2"), TList::Class(), AliAnalysisManager::kInputContainer);
				if (uniqueID.EqualTo("18q_Dcaz2d5")   && period.EqualTo("LHC18q")) cin_NUA = mgr->CreateContainer(Form("NUA_Dcaz2d5"), TList::Class(), AliAnalysisManager::kInputContainer);
				
				if (uniqueID.EqualTo("18r_ChiHg2")     && period.EqualTo("LHC18r")) cin_NUA = mgr->CreateContainer(Form("NUA_ChiHg2"), TList::Class(), AliAnalysisManager::kInputContainer); 
				if (uniqueID.EqualTo("18r_ChiHg2d5") && period.EqualTo("LHC18r")) cin_NUA = mgr->CreateContainer(Form("NUA_ChiHg2d5"), TList::Class(), AliAnalysisManager::kInputContainer); 
				if (uniqueID.EqualTo("18r_Nhits60")     && period.EqualTo("LHC18r")) cin_NUA = mgr->CreateContainer(Form("NUA_Nhits60"), TList::Class(), AliAnalysisManager::kInputContainer); 
				if (uniqueID.EqualTo("18r_Nhits80")     && period.EqualTo("LHC18r")) cin_NUA = mgr->CreateContainer(Form("NUA_Nhits80"), TList::Class(), AliAnalysisManager::kInputContainer); 
				if (uniqueID.EqualTo("18r_Dcaz2")       && period.EqualTo("LHC18r")) cin_NUA = mgr->CreateContainer(Form("NUA_Dcaz2"), TList::Class(), AliAnalysisManager::kInputContainer); 
				if (uniqueID.EqualTo("18r_Dcaz2d5")   && period.EqualTo("LHC18r")) cin_NUA = mgr->CreateContainer(Form("NUA_Dcaz2d5"), TList::Class(), AliAnalysisManager::kInputContainer); 
								
				TList* wNUA_list = NULL;
				wNUA_list = dynamic_cast<TList*>(inNUA->FindObjectAny("fNUA_ChPosChNeg"));
		            	cin_NUA->SetData(wNUA_list); 
				mgr->ConnectInput(task,inSlotCounter,cin_NUA);
				inSlotCounter++;
			}else
			{	if (!AllContainers->FindObject("NUA")){
					AliAnalysisDataContainer *cin_NUA;
					if (period.EqualTo("LHC18q") && !AllContainers->FindObject("NUA") )  cin_NUA = mgr->CreateContainer(Form("NUA"), TList::Class(), AliAnalysisManager::kInputContainer);
					if (period.EqualTo("LHC18r") && !AllContainers->FindObject("NUA") ) cin_NUA = mgr->CreateContainer(Form("NUA"), TList::Class(), AliAnalysisManager::kInputContainer);

					TList* wNUA_list = NULL;
					wNUA_list = dynamic_cast<TList*>(inNUA->FindObjectAny("fNUA_ChPosChNeg"));
		            		cin_NUA->SetData(wNUA_list); 
					mgr->ConnectInput(task,inSlotCounter,cin_NUA);
					inSlotCounter++;
				}else 
				{
					mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUA"));
					inSlotCounter++;
					printf("NUA(except for Trk-wise Syst) already loaded\n");
				}
			}

		}

	}; // Input NUA done




	if(task->GetV0CalibOn() || v0calibOn){
		if (period.EqualTo("LHC10h") ) {
			// GainEQ & Recenter
                		TFile *v0calib;
                		if(!AllContainers->FindObject("V0Calib")) {
				AliAnalysisDataContainer *cin_V0Calib = mgr->CreateContainer(Form("V0Calib"), TList::Class(), AliAnalysisManager::kInputContainer);
				v0calib = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hQnCalib.root");
                			TList* qncalib_list = NULL;
				qncalib_list = dynamic_cast<TList*>(v0calib->Get("10hlistqncalib"));
                		    	cin_V0Calib->SetData(qncalib_list); 
                		    	mgr->ConnectInput(task,inSlotCounter,cin_V0Calib);
                		    	inSlotCounter++;
                		}else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("V0Calib"));
				inSlotCounter++;
				printf("V0Calib already loaded\n");
			}
		}
		else if (period.EqualTo("LHC15o")){
                		TFile *qnSp;
                		if(!AllContainers->FindObject("qnSp")) {
				AliAnalysisDataContainer *cin_qnPercSp = mgr->CreateContainer(Form("qnSp"), TList::Class(), AliAnalysisManager::kInputContainer);                			
                			qnSp = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/calibSpq2V0C15oP2.root");
				// Ref V0 qn percentail data copied from alien:////alice/cern.ch/user/a/adobrin/cmeESE15oP2/calibSpq2V0C15oP2.root
                			TList* spperc_list = NULL;
				spperc_list = dynamic_cast<TList*>(qnSp->Get("15olistspPerc"));
                		    	cin_qnPercSp->SetData(spperc_list); 
                		    	mgr->ConnectInput(task,inSlotCounter,cin_qnPercSp);
                		    	inSlotCounter++;                		 
                		}else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("qnSp"));
				inSlotCounter++;
				printf("qnSp already loaded\n");
			}

                		TFile *wCent;
                		if(!AllContainers->FindObject("wCent")) {
				AliAnalysisDataContainer *cin_wCent = mgr->CreateContainer(Form("wCent"), TList::Class(), AliAnalysisManager::kInputContainer);                			
                			wCent = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/weightCent.root");
				// Ref V0 qn percentail data copied from alien:////alice/cern.ch/user/a/adobrin/cmeESE15oP2/calibSpq2V0C15oP2.root
                			TList* wCent_list = NULL;
				wCent_list = dynamic_cast<TList*>(wCent->Get("15owCent"));
                		    	cin_wCent->SetData(wCent_list); 
                		    	mgr->ConnectInput(task,inSlotCounter,cin_wCent);
                		    	inSlotCounter++;                		 
                		}else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("wCent"));
				inSlotCounter++;
				printf("wCent already loaded\n");
			}
			
		}
		else if (period.EqualTo("LHC18q")){
                		TFile *v0calib;
                		if(!AllContainers->FindObject("v0calib")) {
				AliAnalysisDataContainer *cin_v0calib = mgr->CreateContainer(Form("v0calib"), TList::Class(), AliAnalysisManager::kInputContainer);                			
                			v0calib = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/calibq2V0C18qP3.root");

                			TList* v0calib_list = NULL;
				v0calib_list = dynamic_cast<TList*>(v0calib->FindObjectAny("18qlistspPerc"));
                		    	cin_v0calib->SetData(v0calib_list); 
                		    	mgr->ConnectInput(task,inSlotCounter,cin_v0calib);
                		    	inSlotCounter++;                		 
                		}else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("v0calib"));
				inSlotCounter++;
				printf("v0calib already loaded\n");
			}
		}
		else if (period.EqualTo("LHC18r")){
                		TFile *v0calib;
                		if(!AllContainers->FindObject("v0calib")) {
				AliAnalysisDataContainer *cin_v0calib = mgr->CreateContainer(Form("v0calib"), TList::Class(), AliAnalysisManager::kInputContainer);                			
                			v0calib = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/calibq2V0C18rP3.root");
                			
                			TList* v0calib_list = NULL;
				v0calib_list = dynamic_cast<TList*>(v0calib->FindObjectAny("18rlistspPerc"));
                		    	cin_v0calib->SetData(v0calib_list); 
                		    	mgr->ConnectInput(task,inSlotCounter,cin_v0calib);
                		    	inSlotCounter++;                		 
                		}else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("v0calib"));
				inSlotCounter++;
				printf("v0calib already loaded\n");
			}
		}
	}
	


	// Return task pointer at the end
	return task;
}	

