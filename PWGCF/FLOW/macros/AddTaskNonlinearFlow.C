#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliAnalysisTaskNonlinearFlow.h"
#endif

AliAnalysisTaskNonlinearFlow* AddTaskNonlinearFlow(
		Int_t		fFilterbit 		= 96,
		Double_t	fMinPt			= 0.2,
		Double_t	fMaxPt			= 3.0,
                Int_t           trigger                 = 0,
                Int_t           fSystFlag               = 0,
                TString         fPeriod                 = "LHC15o",
                TString         fNtrksName              = "Mult",
		TString		uniqueID        	= ""
		)
{
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
	Bool_t		fNUE 				= true;
	Bool_t		fNUA				= true;

	// Creates a pid task and adds it to the analysis manager
	// Get the pointer to the existing analysis manager via the static
	//access method
	//=========================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskNonlinearFlow.C", "No analysis manager to connect to.");
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

	AliAnalysisTaskNonlinearFlow* taskFlowEp = new AliAnalysisTaskNonlinearFlow("taskFlowEp");
	taskFlowEp->SetDebugLevel(3);
	taskFlowEp->SetFilterbit(fFilterbit); // For systematics
	taskFlowEp->SetFilterbitDefault(fFilterbit);
	taskFlowEp->SetEtaCut(fEtaCut);
	taskFlowEp->SetVtxCut(fVtxCut); // For systematics
	taskFlowEp->SetVtxCutDefault(fVtxCut);
	taskFlowEp->SetMinPt(fMinPt);
	taskFlowEp->SetMaxPt(fMaxPt);
	taskFlowEp->SetTPCclusters(TPCclusters); // For systematics
	taskFlowEp->SetTPCclustersDefault(TPCclusters);
	taskFlowEp->SetChi2PerTPCcluster(chi2PerTPCcluster); // max. chi2 per TPC cluster
	taskFlowEp->SetMinITSClusters(fMinITSClus);
	taskFlowEp->SetMaxChi2(fMaxChi2);
	taskFlowEp->SetUseDCAzCut(fUseDCAzCut);
	taskFlowEp->SetDCAzCut(fDCAz); // For systematics
	taskFlowEp->SetDCAzCutDefault(fDCAz); 
	taskFlowEp->SetUseDCAxyCut(fUseDCAxyCut);
	taskFlowEp->SetDCAxyCut(fDCAxy); // For systematics
	taskFlowEp->SetDCAxyCutDefault(fDCAxy); 
	taskFlowEp->SetIsSample(IsSample);
	taskFlowEp->SetCentFlag(nCentFl);
	taskFlowEp->SetTrigger(trigger);
	taskFlowEp->SetLSFlag(fLS);
	taskFlowEp->SetNUEFlag(fNUE);
	taskFlowEp->SetNUA(fNUA);
	taskFlowEp->SetNtrksName(fNtrksName);
        taskFlowEp->SetSystFlag(fSystFlag);

        taskFlowEp->SetUseWeigthsRunByRun(false);
        taskFlowEp->SetUsePeriodWeigths(false);
        taskFlowEp->SetUseWeights3D(false); 

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
	AliAnalysisDataContainer *cout_hist = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:%s", uniqueID.Data()));
	mgr->ConnectInput (taskFlowEp, 0, cinput);
	mgr->ConnectOutput(taskFlowEp, 1, cout_hist);
	Int_t inSlotCounter=1;
	if(fNUA || fNUE)
		TGrid::Connect("alien:");
	if(fNUA) {
		AliAnalysisDataContainer *cin_NUA = mgr->CreateContainer(Form("NUA%s", uniqueID.Data()), TFile::Class(), AliAnalysisManager::kInputContainer);
               
                TFile *inNUA;

                if (fPeriod.EqualTo("LHC15o")) {
			// inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/LHC15o/RBRweights.root");
			inNUA = TFile::Open("alien:///alice/cern.ch/user/e/enielsen/WeightsPbPb15o.root");
			taskFlowEp->SetUseWeigthsRunByRun(true);
                } else if (fPeriod.EqualTo("LHC17")) {
	            taskFlowEp->SetUsePeriodWeigths(true);
                    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC17/weights_LHC17_MB_periods.root");
                    } else {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC17/weights_LHC17_HM_periods.root");
                    }
                } else if (fPeriod.EqualTo("LHC15i")) {
	            taskFlowEp->SetUsePeriodWeigths(true);
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15i_MB.root");
                    } else {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15i_HM.root");
                    }
                } else if (fPeriod.EqualTo("LHC15l")) {
	            taskFlowEp->SetUsePeriodWeigths(true);
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15l_MB.root");
                    } else {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC15/weights_LHC15l_HM.root");
                    }
                } else if (fPeriod.EqualTo("LHC16")) {
	            taskFlowEp->SetUsePeriodWeigths(true);
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC16/weights_LHC16_MB_periods.root");
                    } else {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC16/weights_LHC16_periods.root");
                    }
                } else if (fPeriod.EqualTo("LHC18")) {
	            taskFlowEp->SetUsePeriodWeigths(true);
		    if (trigger == 0) {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC18/weights_LHC18_MB_periods.root");
                    } else {
                        inNUA = TFile::Open("alien:///alice/cern.ch/user/z/zumoravc/weights/pp_LHC18/weights_LHC18_allHM.root");
                    }
                } 
					
                TList* weight_list = NULL;
		if (fSystFlag == 0 && !fPeriod.EqualTo("LHC15o")) {
			weight_list = dynamic_cast<TList*>(inNUA->Get("weights"));
		} else {
			weight_list = dynamic_cast<TList*>(inNUA->Get("WeightList"));
		}
		cin_NUA->SetData(weight_list);
		mgr->ConnectInput(taskFlowEp,inSlotCounter,cin_NUA);
		inSlotCounter++;
	}
	if(fNUE) {
		AliAnalysisDataContainer *cin_NUE = mgr->CreateContainer(Form("NUE%s", uniqueID.Data()), TFile::Class(), AliAnalysisManager::kInputContainer);
		TFile *inNUE = (fFilterbit==96)?TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_HIR.root"): TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_HIR_FB768.root");
		if(!inNUE) {
			printf("Could not open efficiency file!\n");
			return 0;
		}
		cin_NUE->SetData(inNUE);
		mgr->ConnectInput(taskFlowEp,inSlotCounter,cin_NUE);
		inSlotCounter++;
	}

	// Return task pointer at the end
	return taskFlowEp;
}
//
// EOF
//
