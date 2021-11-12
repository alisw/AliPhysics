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
      int                     trigger=1,
      int                     filterBit=1, // AOD filter bit selection
      int                     nclscut=70, // ncls cut for all tracks 
      float                  chi2hg=4.0,
      float                  chi2lo=0.1,
      float                  dcacutz=3.2, // dcaz cut for all tracks
      float                  dcacutxy=2.4, // dcaxy cut for all tracks
      float                  ptmin=0.2, // minimum pt for Q-vector components
      float                  ptmax=5.0, // maximum pt for Q-vector components
      int                     cbinlo=0, // lower centrality bin for histogram array
      int                     cbinhg=8, // higher centrality bin for histogram array
      TString             period="LHC15o", // period
      TString             multComp="pileupByEDSTPC128", // multiplicity comparison
      double              etaGap=0.3,  
      bool                  v0calibOn=true,
      bool                  tpccalibOn=false,
      bool                  QAV0=true,
      bool                  QATPC=false,
      bool                  doNUE=true,
      bool                  doNUA=true,
      bool                  diffv2QA=true,
      bool                  calcv24ach=true,
      bool                  doTrivalCor=false,
      float                  centcut=7.5 // centrality restriction for V0M and TRK
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
	task->SetFilterBit(filterBit);
	task->SetNclsCut(nclscut);
  	task->SetChi2High(chi2hg);
  	task->SetChi2Low(chi2lo);
	task->SetDCAcutZ(dcacutz);
	task->SetDCAcutXY(dcacutxy);
	task->SetPtMin(ptmin);
	task->SetPtMax(ptmax);
	task->SetCentBinLow(cbinlo);
	task->SetCentBinHigh(cbinhg);
	task->SetPeriod(period);
	task->SetMultComp(multComp);
	task->SetEtaGap(etaGap);
	task->SetV0CalibOn(v0calibOn);
	task->SetTPCCalibOn(tpccalibOn);
	task->SetV0QAOn(QAV0);
	task->SetTPCQAOn(QATPC);
	task->SetNUEOn(doNUE);
	task->SetNUAOn(doNUA);
  	task->SetDiffV2QA(diffv2QA);
  	task->SetCalcV24VsAch(calcv24ach);
  	task->SetDoTrivalCor(doTrivalCor);
	task->SetCentCut(centcut);
	// task->SelectCollisionCandidates(AliVEvent::kINT7);
	mgr->AddTask(task);

	// Create ONLY the output containers for the data produced by the
	// task.  Get and connect other common input/output containers via
	// the manager as below
	//======================================================================
    	AliAnalysisDataContainer* cinput  = mgr->GetCommonInputContainer();
  	AliAnalysisDataContainer* coutput = mgr->CreateContainer("output", TList::Class(), 
                                                           AliAnalysisManager::kOutputContainer, 
                                                           mgr->GetCommonFileName());
   	mgr->ConnectInput (task, 0, cinput);
  	mgr->ConnectOutput(task, 1, coutput);


	Int_t inSlotCounter=1;
   	if (!gGrid) {
   	    TGrid::Connect("alien://");
   	}
        	TObjArray *AllContainers = mgr->GetContainers();

	if(task->GetNUEOn() || doNUE) {
		if(!AllContainers->FindObject("NUE")) {
			AliAnalysisDataContainer *cin_NUE = mgr->CreateContainer(Form("NUE"), TList::Class(), AliAnalysisManager::kInputContainer);
                		TFile *inNUE;
                		if (period.EqualTo("LHC10h") || period.EqualTo("LHC11h")) {
				inNUE = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/Run1NUE.root");
				TList* wNUE_list = NULL;
				wNUE_list = dynamic_cast<TList*>(inNUE->Get("listNUE"));
				if (!wNUE_list) printf("Read TList wrong!\n");
            		    	cin_NUE->SetData(wNUE_list); 			
				mgr->ConnectInput(task,inSlotCounter,cin_NUE);
				inSlotCounter++;
			}
			else if (period.EqualTo("LHC15o")) {
				inNUE = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/efficiencyBothpol.root");
				// Ref NUE data from alien:///alice/cern.ch/user/p/prottay/nuarootfiles_p5_one_two_two_FB768_15op2_withpileup/efficiencyBothpol.root
				TList* wNUE_list = NULL;
				wNUE_list = dynamic_cast<TList*>(inNUE->Get("fMcEffiHij"));
				if (!wNUE_list) printf("Read TList wrong!\n");
            		    	cin_NUE->SetData(wNUE_list); 			
				mgr->ConnectInput(task,inSlotCounter,cin_NUE);
				inSlotCounter++;				
			}
		} else {
			printf("Run2 NUE not been calculated!\n");
			mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUE"));
			inSlotCounter++;
			printf("NUE already loaded\n");
		}
	}




	if(task->GetNUAOn() ||doNUA) {
		TFile *inNUA;
		if (period.EqualTo("LHC10h") ) { // NUA for 10h is too large to read, we separate them into 3 TList*s.
			if (task->GetFilterBit() ==1) inNUA = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/10hNUAFB1.root");
			else if (task->GetFilterBit() ==128) inNUA = TFile::Open("alien:///alice/cern.ch/user/w/wenya/wNUA/weightNUA_LHC10hMAX_Cent_FB128.root");
	    			
			// Run 139510-138653
			if(!AllContainers->FindObject("NUA1")) {
				AliAnalysisDataContainer *cin_NUA1 = mgr->CreateContainer(Form("NUA1"), TList::Class(), AliAnalysisManager::kInputContainer);				
	                		TList* wNUA_list1 = NULL;
				wNUA_list1 = dynamic_cast<TList*>(inNUA->Get("listNUA_139510to138653"));
		                	cin_NUA1->SetData(wNUA_list1); 
				mgr->ConnectInput(task,inSlotCounter,cin_NUA1);
				delete wNUA_list1;
				inSlotCounter++;
				
			} else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUA1"));
				inSlotCounter++;
				printf("NUA1 already loaded\n");
			}

			// Run 138652-137693
			if(!AllContainers->FindObject("NUA2")) {
				AliAnalysisDataContainer *cin_NUA2 = mgr->CreateContainer(Form("NUA2"), TList::Class(), AliAnalysisManager::kInputContainer);		
	                		TList* wNUA_list2 = NULL;
				wNUA_list2 = dynamic_cast<TList*>(inNUA->Get("listNUA_138652to137693"));
		                	cin_NUA2->SetData(wNUA_list2); 
				mgr->ConnectInput(task,inSlotCounter,cin_NUA2);
				inSlotCounter++;
				delete wNUA_list2;
			} else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUA2"));
				inSlotCounter++;
				printf("NUA2 already loaded\n");
			}

			// Run 137692-137161
			if(!AllContainers->FindObject("NUA3")) {		
				AliAnalysisDataContainer *cin_NUA3 = mgr->CreateContainer(Form("NUA3"), TList::Class(), AliAnalysisManager::kInputContainer);
	                		TList* wNUA_list3 = NULL;
				wNUA_list3 = dynamic_cast<TList*>(inNUA->Get("listNUA_137692to137161"));
		                	cin_NUA3->SetData(wNUA_list3); 
				mgr->ConnectInput(task,inSlotCounter,cin_NUA3);
				inSlotCounter++;
				delete wNUA_list3;				
			} else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUA3"));
				inSlotCounter++;
				printf("NUA3 already loaded\n");
			}
		} else if (period.EqualTo("LHC11h")){
    		  		if(!AllContainers->FindObject("NUA")) {
				AliAnalysisDataContainer *cin_NUA = mgr->CreateContainer(Form("NUA"), TList::Class(), AliAnalysisManager::kInputContainer);				
    		    		if (task->GetFilterBit() ==1) inNUA = TFile::Open("alien:///alice/cern.ch/user/w/wenya/wNUA/weightNUA_LHC11hMAX_ReBined.root");
    		    		TList* wNUA_list = NULL;
    		    		wNUA_list = dynamic_cast<TList*>(inNUA->Get("listNUA"));
    		    		cin_NUA->SetData(wNUA_list);
    		    		mgr->ConnectInput(task,inSlotCounter,cin_NUA);
    		    		inSlotCounter++;
    		  	} else {
    		    		mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUA"));
    		    		inSlotCounter++;
    		    		printf("NUA already loaded\n");
    		  	}
    		} else if (period.EqualTo("LHC15o")) {
			if(!AllContainers->FindObject("NUA")) {
				if (task->GetFilterBit() ==768) inNUA = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/wgtPion_NUAFB768DeftwPUcut_LHC15op2_24Aug2021.root");
				// Ref NUA data from alien:///alice/cern.ch/user/p/prottay/nuarootfiles_p5_one_two_two_FB768_15op2_withpileup
				// /wgtPion_NUAFB768DeftwPUcut_LHC15op2_24Aug2021.root  (15o_pass2)
				TDirectoryFile* wNUA_directoryfile = (TDirectoryFile*)inNUA->Get("ZDCgains");
				AliAnalysisDataContainer *cin_NUA = mgr->CreateContainer(Form("NUA"), TList::Class(), AliAnalysisManager::kInputContainer);				
				TList* wNUA_list = NULL;
				wNUA_list = dynamic_cast<TList*>(wNUA_directoryfile->Get("fNUA_ChPosChNeg"));
		                	cin_NUA->SetData(wNUA_list); 
				mgr->ConnectInput(task,inSlotCounter,cin_NUA);
				inSlotCounter++;			
			} else {
				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUA"));
				inSlotCounter++;
				printf("NUA already loaded\n");
			}
		}

	}




	if(task->GetV0CalibOn() || v0calibOn){
		if (period.EqualTo("LHC10h") ) {
			// GainEQ & Recenter
                		TFile *v0calib;
                		if(!AllContainers->FindObject("V0Calib")) {
				AliAnalysisDataContainer *cin_V0Calib = mgr->CreateContainer(Form("V0Calib"), TList::Class(), AliAnalysisManager::kInputContainer);
				v0calib = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/10hQnCalib.root");
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
		}else if (period.EqualTo("LHC11h")) {
			TFile *qnSp;
                    		if(!AllContainers->FindObject("qnSp")) {
                    			AliAnalysisDataContainer *cin_qnPercSp = mgr->CreateContainer(Form("qnPercSp"), TList::Class(), AliAnalysisManager::kInputContainer);
                      			qnSp = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc11h/calibSpq2V0C11hP2.root");
                      			TList* spperc_list = NULL;
        				spperc_list = dynamic_cast<TList*>(qnSp->Get("11hlistspPerc"));
                          		cin_qnPercSp->SetData(spperc_list);
                          		mgr->ConnectInput(task,inSlotCounter,cin_qnPercSp);
                          		inSlotCounter++;
                    		}else {
        				mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("qnSp"));
        				inSlotCounter++;
        				printf("qnSp already loaded\n");
      			}
                	} else if (period.EqualTo("LHC15o")){
                		TFile *qnSp;
                		if(!AllContainers->FindObject("qnSp")) {
				AliAnalysisDataContainer *cin_qnPercSp = mgr->CreateContainer(Form("qnPercSp"), TList::Class(), AliAnalysisManager::kInputContainer);                			
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
		}
	}
	


	// Return task pointer at the end
	return task;
}	

