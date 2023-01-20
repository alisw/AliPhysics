#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliAnalysisTaskFlowPPTask.h"
#endif

AliAnalysisTaskFlowPPTask* AddFlowPPTask(
    Int_t		fFilterbit 		= 96,
		Double_t	fMinPt			= 0.2,
		Double_t	fMaxPt			= 3.0,
        Int_t           trigger                 = 0,
        Int_t           fSystFlag               = 0,
        TString         fPeriod                 = "LHC17",
        TString         fNtrksName              = "Mult",
		Bool_t		fNUE 				= true,
		Bool_t		fNUA				= true,
		Bool_t		UseCorrectedNTracks = true,
		TString		uniqueID        	= "Default"
    )
{
	TString name = "MyFlowPPTask";
     // The common parameters
	Double_t	fEtaCut 			= 0.8;
	Double_t	fVtxCut				= 10.0;
	Int_t		TPCclusters		        = 70;
	Double_t        chi2PerTPCcluster               = 10000;
	Int_t		fMinITSClus		        = 5;
	Double_t	fMaxChi2			= 2.5;
	Bool_t		fUseDCAzCut		        = false;
	Double_t	fDCAz				= 1.0;
	Bool_t		fUseDCAxyCut	                = false;
	Double_t	fDCAxy				= 0.2;
	Int_t		IsSample			= 10;
	Short_t		nCentFl				= 0;
	Bool_t		fLS				= false;
	Bool_t		fAddTPCPileupCuts = false;
	Double_t	fESDvsTPConlyLinearCut = 15000.;
	//Bool_t		fNUE 				= true;
	//Bool_t		fNUA				= true;
	


    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }

    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyResults";      
	//fileName += uniqueID; // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskFlowPPTask* task = new AliAnalysisTaskFlowPPTask(name.Data());   
    if(!task) return 0x0;
	task->SetDebugLevel(3);
	task->SetFilterbit(fFilterbit); // For systematics
	task->SetFilterbitDefault(fFilterbit);
	task->SetEtaCut(fEtaCut);
	task->SetVtxCut(fVtxCut); // For systematics
	task->SetVtxCutDefault(fVtxCut);
	task->SetMinPt(fMinPt);
	task->SetMaxPt(fMaxPt);
	task->SetTPCclusters(TPCclusters); // For systematics
	task->SetTPCclustersDefault(TPCclusters);
	task->SetChi2PerTPCcluster(chi2PerTPCcluster); // max. chi2 per TPC cluster
	task->SetMinITSClusters(fMinITSClus);
	task->SetMaxChi2(fMaxChi2);
	task->SetUseDCAzCut(fUseDCAzCut);
	task->SetDCAzCut(fDCAz); // For systematics
	task->SetDCAzCutDefault(fDCAz); 
	task->SetUseDCAxyCut(fUseDCAxyCut);
	task->SetDCAxyCut(fDCAxy); // For systematics
	task->SetDCAxyCutDefault(fDCAxy); 
	task->SetIsSample(IsSample);
	task->SetCentFlag(nCentFl);
	task->SetTrigger(trigger);
	task->SetLSFlag(fLS);
	task->SetNUEFlag(fNUE);
	task->SetNUA(fNUA);
	task->SetNtrksName(fNtrksName);
    task->SetSystFlag(fSystFlag);
    task->SetUseWeigthsRunByRun(false);
    task->SetUsePeriodWeigths(false);
    task->SetUseWeights3D(false); 
	task->SetPeriod(fPeriod);
	task->SetOnlineTrackCorrection(UseCorrectedNTracks);
	task->SetAddTPCPileupCuts(fAddTPCPileupCuts);
	task->SetESDvsTPConlyLinearCut(fESDvsTPConlyLinearCut);
	//task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    
    //task->SelectCollisionCandidates(AliVEvent::kINT7);
    /*
    task->SetDebugLevel(3);
	task->SetFilterbit(fFilterbit); // For systematics
	task->SetFilterbitDefault(fFilterbit);
	task->SetEtaCut(fEtaCut);
	task->SetVtxCut(fVtxCut); // For systematics
	task->SetVtxCutDefault(fVtxCut);
	task->SetMinPt(fMinPt);
	task->SetMaxPt(fMaxPt);
	task->SetTPCclusters(TPCclusters); // For systematics
	task->SetTPCclustersDefault(TPCclusters);
	task->SetChi2PerTPCcluster(chi2PerTPCcluster); // max. chi2 per TPC cluster
	task->SetMinITSClusters(fMinITSClus);
	task->SetMaxChi2(fMaxChi2);
	task->SetUseDCAzCut(fUseDCAzCut);
	task->SetDCAzCut(fDCAz); // For systematics
	task->SetDCAzCutDefault(fDCAz); 
	task->SetUseDCAxyCut(fUseDCAxyCut);
	task->SetDCAxyCut(fDCAxy); // For systematics
	task->SetDCAxyCutDefault(fDCAxy); 
	task->SetIsSample(IsSample);
	task->SetCentFlag(nCentFl);
	task->SetTrigger(trigger);
	task->SetLSFlag(fLS);
	task->SetNUEFlag(fNUE);
	task->SetNUA(fNUA);
	task->SetNtrksName(fNtrksName);
    task->SetSystFlag(fSystFlag);
    task->SetUseWeigthsRunByRun(false);
    task->SetUsePeriodWeigths(false);
    task->SetUseWeights3D(false); 
	task->SetPeriod(fPeriod);
    */
    
    // add your task to the manager
    mgr->AddTask(task);


    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid


	//Config Weights Input
	Int_t inSlotCounter=1;
	if(fNUA || fNUE)TGrid::Connect("alien:");

	TObjArray *AllContainers = mgr->GetContainers();

	//NUA
	if(fNUA) {
		//if(AllContainers->FindObject("NUA")){}(from AddTaskNonlinearFlowC)
		if(!AllContainers->FindObject("NUA")){
			AliAnalysisDataContainer *cin_NUA = mgr->CreateContainer(Form("NUA"), TFile::Class(), AliAnalysisManager::kInputContainer);
               
            TFile *inNUA=nullptr;

            if (fPeriod.EqualTo("LHC15o")) {
				// inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/LHC15o/RBRweights.root");
				inNUA = TFile::Open("alien:///alice/cern.ch/user/e/enielsen/Weights/NUA/WeightsPbPb15o.root");
				//inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zhlu/NUA_from_enielsen/WeightsPbPb15o.root");
				//inNUA = TFile::Open("/Users/lisck/Workdir/Alice_NonlinearFlow/Analysis/AnaFlow_LEGO/WeightsPbPb15o.root");
				//task->SetUseWeigthsRunByRun(true);
            } 
			else if (fPeriod.EqualTo("LHC15o_pass2")){
				inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/WeightsPbPb15o_pass2.root");
			}
			else if (fPeriod.EqualTo("LHC18qr_pass3")){
				inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/WeightsPbPb18qr_pass3.root");
			}
			else if (fPeriod.EqualTo("LHC17")) {
	            task->SetUsePeriodWeigths(true);
				inNUA = TFile::Open("alien:///alice/cern.ch/user/v/vvislavi/Weights/pp_NUA/Weights_pp17.root");
                //if (trigger == 0) {
                //    inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC17/weights_LHC17_MB_periods.root");
                //} else {
                //    inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC17/weights_LHC17_HM_periods.root");
                //}
            }
			// XeXe dataset
			else if (fPeriod.EqualTo("LHC17n")) {
	            task->SetUsePeriodWeigths(true);
				inNUA = TFile::Open("alien:///alice/cern.ch/user/e/enielsen/WeightsXeXe.root");

            } else if (fPeriod.EqualTo("LHC15i")) {
	            task->SetUsePeriodWeigths(true);
		    	if (trigger == 0) {
                    inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15i_MB.root");
                    } else {
                    inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15i_HM.root");
                    }
            } else if (fPeriod.EqualTo("LHC15l")) {
	            task->SetUsePeriodWeigths(true);
		    	if (trigger == 0) {
                    inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15l_MB.root");
                    } else {
                    inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15l_HM.root");
                    }
            } else if (fPeriod.EqualTo("LHC16")) {
	            task->SetUsePeriodWeigths(true);
				inNUA = TFile::Open("alien:///alice/cern.ch/user/v/vvislavi/Weights/pp_NUA/Weights_pp16.root");
		    	//if (trigger == 0) {
                        //inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC16/weights_LHC16_MB_periods.root");
                    //} else {
                     //   inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC16/weights_LHC16_periods.root");
                    //}
            } else if (fPeriod.EqualTo("LHC18")) {
	            task->SetUsePeriodWeigths(true);
				inNUA = TFile::Open("alien:///alice/cern.ch/user/v/vvislavi/Weights/pp_NUA/Weights_pp18.root");
		  		//if (trigger == 0) {
                    //    inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC18/weights_LHC18_MB_periods.root");
                    //} else {
                        //inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC18/weights_LHC18_periods.root");
					//	inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC18/weights_LHC18_allHM.root");
					//}
            } 

			if(!inNUA){
				printf("Can't Access NUA File\n");
				return 0;
			}
            TList* weight_list = NULL;
			weight_list = dynamic_cast<TList*>(inNUA->Get("WeightList"));
			//if (fSystFlag == 0 && !fPeriod.EqualTo("LHC15o")) {
			//	weight_list = dynamic_cast<TList*>(inNUA->Get("weights"));
			//} else {
			//	weight_list = dynamic_cast<TList*>(inNUA->Get("WeightList"));
			//}
			if(!weight_list) {
				printf("Could not get NUA!\n");
				return 0;
			}
			cin_NUA->SetData(weight_list);
			//Connect the weight input with task
			mgr->ConnectInput(task,inSlotCounter,cin_NUA);
			inSlotCounter++;
		}
		else{
			//Already have NUA container, avoid loading repeatedly
			mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUA"));
		    inSlotCounter++;
		    printf("NUA already loaded\n");
		}
		


	}
	//NUE
	if(fNUE) {

		if(!AllContainers->FindObject("NUE")){

			AliAnalysisDataContainer *cin_NUE = mgr->CreateContainer(Form("NUE"), TFile::Class(), AliAnalysisManager::kInputContainer);
			TFile *inNUE=nullptr;
			// Xe-Xe Dataset
			if (fPeriod.EqualTo("LHC17n")){
				inNUE =TFile::Open("alien:///alice/cern.ch/user/e/enielsen/LHC17nEfficiency_tmp.root");
			}
			// Pb-Pb Dataset
			else if(fPeriod.EqualTo("LHC15o")){
				inNUE =TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/LHC18e1_MBEff_FD_wSyst_v2.root");
				//inNUE = TFile::Open("/Users/lisck/Workdir/Alice_NonlinearFlow/Analysis/AnaFlow_LEGO/LHC18e1_MBEff_FD_wSyst_v2.root");
			}
			else if(fPeriod.EqualTo("LHC15o_pass2")){
				inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/Efficiency_LHC20j6a_wSyst.root");
			}
			else if(fPeriod.EqualTo("LHC18qr_pass3")){
				inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/Efficiency_LHC20e3a_wSyst.root");
			}
			// p-p Dataset
			else{
				inNUE =TFile::Open("alien:///alice/cern.ch/user/v/vvislavi/Aux/LHC17d20a1_WithModEff_Syst.root");
			}		
			TList* weight_listEff = NULL;
			weight_listEff = dynamic_cast<TList*>(inNUE->Get("EffAndFD"));
			//TFile *inNUE = (fFilterbit==96)?TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_HIR.root"): TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_HIR_FB768.root");
			if(!weight_listEff) {
				printf("Could not get efficiency!\n");
				return 0;
			}
			cin_NUE->SetData(weight_listEff);
			mgr->ConnectInput(task,inSlotCounter,cin_NUE);
			inSlotCounter++;
		}
		else{
			mgr->ConnectInput(task,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUE"));
		    inSlotCounter++;
		    printf("NUE already loaded\n");
		}

		
	}


    return task;
}
