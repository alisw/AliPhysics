#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#endif

AliAnalysisTaskCMWESETrkSyst* AddTaskCMWESETrkSyst(
      int                     debug=0, // debug level controls amount of output statements
      double              Harmonic=2,
      TString             trigger="kINT7",
      int                     filterBit=1, // AOD filter bit selection
      int                     FltbitSyst=1, // AOD filter bit selection
      float                  ptmin=0.2, // minimum pt for Q-vector components
      float                  ptmax=5.0, // maximum pt for Q-vector components
      int                     cbinhg=8, // higher centrality bin for histogram array
      int                     cbinlo=0, // lower centrality bin for histogram array
      TString             period="LHC15o", // period
      bool                  v0calibOn=true,
      bool                  QAV0=true,
      bool                  doNUE=true,
      bool                  doNUA=true,
      bool                  MCSyst=true,
      bool                  FBSyst=true
      )	
{	
	// Creates a pid task and adds it to the analysis manager
	// Get the pointer to the existing analysis manager via the static
	//access method
	//=========================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskCMWESEsyst.C", "No analysis manager to connect to.");
		return NULL;
	}

	// Check the analysis type using the event handlers connected to the
	// analysis manager The availability of MC handler can also be
	// checked here.
	// =========================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTaskCMWESEsyst.C", "This task requires an input event handler");
		return NULL;
	}
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  	// --- instantiate analysis task
  	AliAnalysisTaskCMWESETrkSyst *task = new AliAnalysisTaskCMWESETrkSyst("TaskCMWESETrkSyst", period, doNUE, doNUA, v0calibOn, MCSyst, FBSyst);
	task->SetDebug(debug);
	task->SetHarmonic(Harmonic);
	task->SetTrigger(trigger);
	task->SetFilterBit(filterBit);
	task->SetFilterBitSyst(FltbitSyst);
	task->SetPtMin(ptmin);
	task->SetPtMax(ptmax);
	task->SetCentBinHigh(cbinhg);
	task->SetCentBinLow(cbinlo);
	task->SetPeriod(period);
	task->SetV0CalibOn(v0calibOn);
	task->SetV0QAOn(QAV0);
	task->SetNUEOn(doNUE);
	task->SetNUAOn(doNUA);
	task->SetMCSystOn(MCSyst);
	task->SetFBSystOn(FBSyst);
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
	
   	TGrid::Connect("alien://");

        	TObjArray *AllContainers = mgr->GetContainers();

	if(task->GetNUEOn()) {
		if(!AllContainers->FindObject("NUE") || !AllContainers->FindObject("NUESyst") || !AllContainers->FindObject("NUESystFB")) {
                		if (period.EqualTo("LHC10h") || period.EqualTo("LHC11h")) {
				AliAnalysisDataContainer *cin_NUE = mgr->CreateContainer(Form("NUE"), TList::Class(), AliAnalysisManager::kInputContainer);
                			TFile *inNUE;
				inNUE = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/Run1NUE.root");
				TList* wNUE_list = NULL;
				wNUE_list = dynamic_cast<TList*>(inNUE->Get("listNUE"));
				if (!wNUE_list) printf("Read TList wrong!\n");
            		    	cin_NUE->SetData(wNUE_list); 			
				mgr->ConnectInput(task,inSlotCounter,cin_NUE);
				inSlotCounter++;
				if (MCSyst){
					AliAnalysisDataContainer *cin_NUESyst = mgr->CreateContainer(Form("NUESyst"), TList::Class(), AliAnalysisManager::kInputContainer);
                				TFile *inNUESyst;
					inNUESyst = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/Run1NUE_AMPT.root");
					TList* wNUE_Systlist = NULL;
					wNUE_Systlist = dynamic_cast<TList*>(inNUESyst->Get("listNUEAMPT"));
					if (!wNUE_Systlist) printf("Read TList wrong!\n");
            		    		cin_NUESyst->SetData(wNUE_Systlist); 			
					mgr->ConnectInput(task,inSlotCounter,cin_NUESyst);
					inSlotCounter++;					
				}
				if (FBSyst){
					AliAnalysisDataContainer *cin_NUESystFB = mgr->CreateContainer(Form("NUESystFB"), TList::Class(), AliAnalysisManager::kInputContainer);
                				TFile *inNUESystFB;
					inNUESystFB = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/Run1NUEFB272.root");
					TList* wNUE_SystFBlist = NULL;
					wNUE_SystFBlist = dynamic_cast<TList*>(inNUESystFB->Get("listNUE"));
					if (!wNUE_SystFBlist) printf("Read TList wrong!\n");
            		    		cin_NUESystFB->SetData(wNUE_SystFBlist); 			
					mgr->ConnectInput(task,inSlotCounter,cin_NUESystFB);
					inSlotCounter++;						
				}

			}
			else if (period.EqualTo("LHC15o")) {
				AliAnalysisDataContainer *cin_NUE = mgr->CreateContainer(Form("NUE"), TList::Class(), AliAnalysisManager::kInputContainer);
				TFile *inNUE;
				inNUE = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/efficiencyBothpol.root");
				// Ref NUE data from alien:///alice/cern.ch/user/p/prottay/nuarootfiles_p5_one_two_two_FB768_15op2_withpileup/efficiencyBothpol.root
				TList* wNUE_list = NULL;
				wNUE_list = dynamic_cast<TList*>(inNUE->Get("fMcEffiHij"));
				if (!wNUE_list) printf("Read TList wrong!\n");
            		    	cin_NUE->SetData(wNUE_list); 			
				mgr->ConnectInput(task,inSlotCounter,cin_NUE);
				inSlotCounter++;
				if (MCSyst){
					AliAnalysisDataContainer *cin_NUESyst = mgr->CreateContainer(Form("NUESyst"), TList::Class(), AliAnalysisManager::kInputContainer);
                				TFile *inNUESyst;
					inNUESyst = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/efficiencyBothpolAMPT.root");
					TList* wNUE_Systlist = NULL;
					wNUE_Systlist = dynamic_cast<TList*>(inNUESyst->Get("fMcEffiHij"));
					if (!wNUE_Systlist) printf("Read TList wrong!\n");
            		    		cin_NUESyst->SetData(wNUE_Systlist); 			
					mgr->ConnectInput(task,inSlotCounter,cin_NUESyst);
					inSlotCounter++;					
				}
				if (FBSyst){
					AliAnalysisDataContainer *cin_NUESystFB = mgr->CreateContainer(Form("NUESystFB"), TList::Class(), AliAnalysisManager::kInputContainer);
                				TFile *inNUESystFB;
					inNUESystFB = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/Efficiency_LHC20j6a_wSyst.root");
					TList* wNUE_SystFBlist = NULL;
					wNUE_SystFBlist = dynamic_cast<TList*>(inNUESystFB->Get("EffAndFD"));
					if (!wNUE_SystFBlist) printf("Read TList wrong!\n");
            		    		cin_NUESystFB->SetData(wNUE_SystFBlist); 			
					mgr->ConnectInput(task,inSlotCounter,cin_NUESystFB);
					inSlotCounter++;						
				}				
			}
		} else {
			printf("Run2 NUE not been calculated!\n");
			mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUE"));
			inSlotCounter++;
			if (MCSyst) {mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUESyst")); inSlotCounter++;}
			if (FBSyst)  {mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUESystFB")); inSlotCounter++;}
			printf("NUE already loaded\n");
		}
	}





	if(task->GetNUAOn()) {
		if(!AllContainers->FindObject("NUA") || !AllContainers->FindObject("NUASystFB") ){
			if (period.EqualTo("LHC10h") ) { // NUA for 10h is too large to read, we separate them into 3 TList*s.
				AliAnalysisDataContainer *cin_NUA = mgr->CreateContainer(Form("NUA"), TList::Class(), AliAnalysisManager::kInputContainer);				
				TFile* inNUA = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB1.root");
	                		TList* wNUA_list = NULL;
				wNUA_list = dynamic_cast<TList*>(inNUA->Get("10hListNUAFB1"));
		                	cin_NUA->SetData(wNUA_list); 
				mgr->ConnectInput(task,inSlotCounter,cin_NUA);
				inSlotCounter++;
				if (FBSyst){
					AliAnalysisDataContainer *cin_NUASystFB = mgr->CreateContainer(Form("NUASystFB"), TList::Class(), AliAnalysisManager::kInputContainer);
                				TFile *inNUASystFB;
					inNUASystFB = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc10h/10hNUAFB272.root");
					TList* wNUA_SystFBlist = NULL;
					wNUA_SystFBlist = dynamic_cast<TList*>(inNUASystFB->Get("10hListNUAFB272"));
					if (!wNUA_SystFBlist) printf("Read TList wrong!\n");
            		    		cin_NUASystFB->SetData(wNUA_SystFBlist); 			
					mgr->ConnectInput(task,inSlotCounter,cin_NUASystFB);
					inSlotCounter++;						
				}				
			}else if (period.EqualTo("LHC11h")){
				AliAnalysisDataContainer *cin_NUA = mgr->CreateContainer(Form("NUA"), TList::Class(), AliAnalysisManager::kInputContainer);				
				TFile* inNUA = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc11h/11hNUAFB1.root");
    		    		TList* wNUA_list = NULL;
    		    		wNUA_list = dynamic_cast<TList*>(inNUA->Get("11hListNUAFB1"));
    		    		cin_NUA->SetData(wNUA_list);
    		    		mgr->ConnectInput(task,inSlotCounter,cin_NUA);
    		    		inSlotCounter++;
				if (FBSyst){
					AliAnalysisDataContainer *cin_NUASystFB = mgr->CreateContainer(Form("NUASystFB"), TList::Class(), AliAnalysisManager::kInputContainer);
                				TFile *inNUASystFB;
					inNUASystFB = TFile::Open("");
					TList* wNUA_SystFBlist = NULL;
					wNUA_SystFBlist = dynamic_cast<TList*>(inNUASystFB->Get(""));
					if (!wNUA_SystFBlist) printf("Read TList wrong!\n");
            		    		cin_NUASystFB->SetData(wNUA_SystFBlist); 			
					mgr->ConnectInput(task,inSlotCounter,cin_NUASystFB);
					inSlotCounter++;						
				}   		 
    			}else if (period.EqualTo("LHC15o")) {
				AliAnalysisDataContainer *cin_NUA = mgr->CreateContainer(Form("NUA"), TList::Class(), AliAnalysisManager::kInputContainer);				
    				TFile* inNUA = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/wgtPion_NUAFB768DeftwPUcut_LHC15op2_24Aug2021.root");
				TDirectoryFile* wNUA_directoryfile = (TDirectoryFile*)inNUA->Get("ZDCgains");
				TList* wNUA_list = NULL;
				wNUA_list = dynamic_cast<TList*>(wNUA_directoryfile->Get("fNUA_ChPosChNeg"));
		                	cin_NUA->SetData(wNUA_list); 
				mgr->ConnectInput(task,inSlotCounter,cin_NUA);
				inSlotCounter++;
				if (FBSyst){
					AliAnalysisDataContainer *cin_NUASystFB = mgr->CreateContainer(Form("NUASystFB"), TList::Class(), AliAnalysisManager::kInputContainer);
                				TFile *inNUASystFB;
					inNUASystFB = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/LHC15o_pass2_NUA_wSyst.root");
					TList* wNUA_SystFBlist = NULL;
					wNUA_SystFBlist = dynamic_cast<TList*>(inNUASystFB->Get("WeightList"));
					if (!wNUA_SystFBlist) printf("Read TList wrong!\n");
            		    		cin_NUASystFB->SetData(wNUA_SystFBlist); 			
					mgr->ConnectInput(task,inSlotCounter,cin_NUASystFB);
					inSlotCounter++;						
				}
			}
		} else {
			mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUA"));
			inSlotCounter++;
			if (FBSyst)  {mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUASystFB")); inSlotCounter++;}
			printf("NUE already loaded\n");
		}
	}





	if(task->GetV0CalibOn()){
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

