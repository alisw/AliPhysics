#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliAnalysisTaskXDeptFlow.h"
#endif

AliAnalysisTaskXDeptFlow* AddTaskXDeptFlow(
		Int_t		fFilterbit 		= 96,
		Double_t	fMinPt			= 0.2,
		Double_t	fMaxPt			= 3.0,
                Int_t           trigger                 = 0,
                Int_t           fSystFlag               = 0,
                TString         fPeriod                 = "LHC15o",
                TString         fNtrksName              = "Mult",
	        Bool_t		fNUA			= true,
	        Bool_t		fNUE 			= true,
		TString		uniqueID        	= ""
		)
{
        // The common parameters
	Double_t	fEtaCut 			= 0.8;

	// Creates a pid task and adds it to the analysis manager
	// Get the pointer to the existing analysis manager via the static
	//access method
	//=========================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskXDeptFlow.C", "No analysis manager to connect to.");
		return NULL;
	}

	// Check the analysis type using the event handlers connected to the
	// analysis manager The availability of MC handler can also be
	// checked here.
	// =========================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTaskNonLinearFlow.C", "This task requires an input event handler");
		return NULL;
	}

	// Create the task and configure it
	//========================================================================
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

	AliAnalysisTaskXDeptFlow* taskFlowEp = new AliAnalysisTaskXDeptFlow("taskFlowEp", fNUA, fNUE, fPeriod);
	taskFlowEp->SetDebugLevel(3);
	taskFlowEp->SetEtaCut(fEtaCut);
	taskFlowEp->SetMinPt(fMinPt);
	taskFlowEp->SetMaxPt(fMaxPt);
	taskFlowEp->SetTrigger(trigger);
	taskFlowEp->SetNUEFlag(fNUE);
	taskFlowEp->SetNUA(fNUA);
	taskFlowEp->SetNtrksName(fNtrksName);
        taskFlowEp->SetSystFlag(fSystFlag);

	//....
	taskFlowEp->SetPeriod(fPeriod);
	mgr->AddTask(taskFlowEp);


	// Create ONLY the output containers for the data produced by the
	// task.  Get and connect other common input/output containers via
	// the manager as below
	//=======================================================================
	//TString fileName = AliAnalysisManager::GetCommonFileName();
	//fileName+=suffixName;
	AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *cout_hist = mgr->CreateContainer(Form("QA_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:%s", uniqueID.Data()));
	mgr->ConnectInput (taskFlowEp, 0, cinput);
	mgr->ConnectOutput(taskFlowEp, 1, cout_hist);
	AliAnalysisDataContainer *physics_hist = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:%s", uniqueID.Data()));
	mgr->ConnectOutput(taskFlowEp, 2, physics_hist);
        int outSlotCounter=2;
        for (int i = 0; i < 30; i++) {
            outSlotCounter++;
	    AliAnalysisDataContainer *physics_hist_i = mgr->CreateContainer(Form("output_%s_%d", uniqueID.Data(), i), TList::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:%s", uniqueID.Data()));
	    mgr->ConnectOutput(taskFlowEp, outSlotCounter, physics_hist_i);
        }
        
	Int_t inSlotCounter=1;
	TGrid::Connect("alien:");

        TObjArray *AllContainers = mgr->GetContainers();
	if(fNUA) {

                if(!AllContainers->FindObject("NUA")) {

		AliAnalysisDataContainer *cin_NUA = mgr->CreateContainer(Form("NUA"), TList::Class(), AliAnalysisManager::kInputContainer);
                TFile *inNUA;

		// PbPb periods
                if (fPeriod.EqualTo("LHC15o")) {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/WeightsPbPb15o.root");
		} else if (fPeriod.EqualTo("LHC15o_pass2")) {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/WeightsPbPb15o_pass2.root");
		} else if (fPeriod.EqualTo("LHC18qr_pass3")) {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/WeightsPbPb18qr_pass3.root");
                } else if (fPeriod.EqualTo("LHC15oKatarina")) {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/PhiWeight_Katarina.root");
                } 
		// XeXe Periods
		else if (fPeriod.EqualTo("LHC17n")) {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/WeightsXeXeLHC17n_wCLX.root");
		}
		// pPb Periods
		else if (fPeriod.EqualTo("LHC16qt")) {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/Weights_pPb16qt_pass1.root");
		}
	        // pp Periods
		else if (fPeriod.EqualTo("LHC15i")) { 
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15i_MB.root");
                    } else {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15i_HM.root");
                    }
                } else if (fPeriod.EqualTo("LHC15l")) {
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15l_MB.root");
                    } else {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15l_HM.root");
                    }
                } else if (fPeriod.EqualTo("LHC16Preview")) {
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC16/weights_LHC16_MB_periods.root");
                    } else {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/Weights_pp16.root");
                    }
                } else if (fPeriod.EqualTo("LHC17Preview")) {
                    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC17/weights_LHC17_MB_periods.root");
                    } else {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/Weights_pp17.root");
                    }
                } else if (fPeriod.EqualTo("LHC18Preview")) {
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC18/weights_LHC18_MB_periods.root");
                    } else {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/Weights_pp18.root");
                    }
                } else if (fPeriod.EqualTo("LHC16")) {
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC16/weights_LHC16_MB_periods.root");
                    } else {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/PeriodNUA_META_LHC16_AOD234.root");
                    }
                } else if (fPeriod.EqualTo("LHC17")) {
                    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC17/weights_LHC17_MB_periods.root");
                    } else {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/PeriodNUA_META_LHC17_pass1_AOD234.root");
                    }
                } else if (fPeriod.EqualTo("LHC18")) {
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC18/weights_LHC18_MB_periods.root");
                    } else {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/PeriodNUA_META_LHC18_pass1_AOD264.root");
                    }
		} else if (fPeriod.EqualTo("LHC16ZM")) {
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC16/weights_LHC16_MB_periods.root");
                    } else {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/weights_LHC16_periods.root");
                    }
                } else if (fPeriod.EqualTo("LHC17ZM")) {
                    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC17/weights_LHC17_HM_periods.root");
                    } else {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/Weights_pp17.root");
                    }
                } else if (fPeriod.EqualTo("LHC18ZM")) {
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC18/weights_LHC18_MB_periods.root");
                    } else {
			inNUA = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUA/weights_LHC18_periods.root");
                    }
		}
					
                TList* weight_list = NULL;
		if (fPeriod.EqualTo("LHC15oKatarina")) {
		    weight_list = dynamic_cast<TList*>(inNUA->Get("weightList"));
                    cin_NUA->SetData(weight_list); 
                } else if (fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18")) {
		    weight_list = dynamic_cast<TList*>(inNUA->Get("WeightList_Default"));
		    cin_NUA->SetData(weight_list);
		} else {
		    weight_list = dynamic_cast<TList*>(inNUA->Get("WeightList"));
		    cin_NUA->SetData(weight_list);
                }
		mgr->ConnectInput(taskFlowEp,inSlotCounter,cin_NUA);
		inSlotCounter++;
		} else {
		    mgr->ConnectInput(taskFlowEp,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUA"));
		    inSlotCounter++;
		    printf("NUA already loaded\n");
		}
	}

	if(fNUE) {

                if(!AllContainers->FindObject("NUE")) {
                TFile *inNUE = NULL;
		AliAnalysisDataContainer *cin_NUE_feeddown = NULL;

		AliAnalysisDataContainer *cin_NUE = mgr->CreateContainer(Form("NUE"), TList::Class(), AliAnalysisManager::kInputContainer);
		if (fPeriod.EqualTo("LHC16qt") ||
			fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18") ||
                        fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview")
                                ) {
		 cin_NUE_feeddown = mgr->CreateContainer(Form("NUE_feeddown"), TList::Class(), AliAnalysisManager::kInputContainer);
		}

		// PbPb Periods
                if (fPeriod.EqualTo("LHC15o")) {
			inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/LHC18e1_MBEff_FD_wSyst_v2.root");
                } else if (fPeriod.EqualTo("LHC15o_pass2")) {
			inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/Efficiency_LHC20j6a_wSyst.root");
                } else if (fPeriod.EqualTo("LHC18qr_pass3")) {
			inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/Efficiency_LHC20e3a_wSyst.root");
                } else if (fPeriod.EqualTo("LHC15oKatarina")) {
			inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/TrackingEfficiency_Katarina.root");
                } 
		// XeXe Periods
		else if (fPeriod.EqualTo("LHC17n")) {
			inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/LHC17nEfficiency.root");
		}
		// pPb Periods
		else if (fPeriod.EqualTo("LHC16qt")) {
			inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/pPb_DPMJet.root");
		}
		// pp Periods
		else { // pp, with preview or updated
	             if (fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18") || 
	                 fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview")
				     ) {
			inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/pp_HM.root");
		    } else { // pp, using Zuzana's previous weight
			inNUE = TFile::Open("alien:///alice/cern.ch/user/m/mzhao/Weights/NUE/LHC17d20a1_WithModEff_Syst.root");
		    }
		}

		if(!inNUE) {
			printf("Could not open efficiency file!\n");
			return 0;
		}

		// Connect the weight
                TList* weight_list = NULL;
                if (fPeriod.EqualTo("LHC15oKatarina")) {
		    weight_list = dynamic_cast<TList*>(inNUE->Get("weightList"));
		    cin_NUE->SetData(weight_list);
                } else if (fPeriod.EqualTo("LHC16qt")) {
		    weight_list = dynamic_cast<TList*>(inNUE->Get("EfficiencyMB"));
		    cin_NUE->SetData(weight_list);
		    TList* feeddown_list = dynamic_cast<TList*>(inNUE->Get("FeeddownMB"));
		    cin_NUE_feeddown->SetData(feeddown_list);
		} else {
	            if (fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18") || 
	                fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview")
			     ) {
		        weight_list = dynamic_cast<TList*>(inNUE->Get("EfficiencyMB"));
		        cin_NUE->SetData(weight_list);
		        TList* feeddown_list = dynamic_cast<TList*>(inNUE->Get("FeeddownMB"));
		        cin_NUE_feeddown->SetData(feeddown_list);
	            } else {
		        weight_list = dynamic_cast<TList*>(inNUE->Get("EffAndFD"));
		        cin_NUE->SetData(weight_list);
		    }
                }
		mgr->ConnectInput(taskFlowEp,inSlotCounter,cin_NUE);
		inSlotCounter++;
		if (fPeriod.EqualTo("LHC16qt")) {
		  mgr->ConnectInput(taskFlowEp,inSlotCounter,cin_NUE_feeddown);
		  inSlotCounter++;
		}
	        if (fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18") || 
	            fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview")) {
		  mgr->ConnectInput(taskFlowEp,inSlotCounter,cin_NUE_feeddown);
		  inSlotCounter++;
		}
            } else {
		    mgr->ConnectInput(taskFlowEp,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUE"));
		    inSlotCounter++;
		    if (fPeriod.EqualTo("LHC16qt")) {
		      mgr->ConnectInput(taskFlowEp,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUE_feeddown"));
		      inSlotCounter++;
		    }
	            if (fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18") || 
	              fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview")) {
		      mgr->ConnectInput(taskFlowEp,inSlotCounter,(AliAnalysisDataContainer*)AllContainers->FindObject("NUE_feeddown"));
		      inSlotCounter++;
		    }
		    printf("NUE already loaded\n");
	    }

	} 

	// Return task pointer at the end
	return taskFlowEp;
}
//
// EOF
//
