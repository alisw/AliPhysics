///////////////////////////////////////////////////////////////////
//                                                               //            
// AddCaloHFEpp                                                  //
// Author: T. Suzuki Univ. of Tsukuba                            //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskCaloHFEpp* AddTaskCaloHFEpp(TString name = "name",
		TString dataname = "dataname",
		Bool_t flagEG1 = kTRUE,
		Bool_t flagEG2 = kFALSE,
		Bool_t flagDG1 = kFALSE,
		Bool_t flagDG2 = kFALSE,
		Bool_t SetFlagClsTypeEMC = kTRUE,
		Bool_t SetFlagClsTypeDCAL = kFALSE,
		Bool_t fMC = kFALSE,
		Double_t TrackEtaMin = -0.6,
		Double_t TrackEtaMax = 0.6,
		Int_t NTPCClust = 80,
		Int_t NITSClust = 3,
		Int_t NCrossedRow = 100,
		Double_t DCAxy = 2.4,
		Double_t DCAz = 3.2,
		Double_t NsigmaMin = -1.0,
		Double_t NsigmaMax = 3.0,
		Double_t M20Min = 0.015,
		Double_t M20Max = 0.35,
		Double_t EopMin = 0.85,
		Double_t EopMax = 1.3,
		Double_t coneR = 0.3,
		Double_t ptAsso = 0.2,
		Double_t mimcle = 0.2,
		TString  pte = "pte",
		Double_t MassMin = 0.15,
		Double_t nref = 12,
		Double_t nrefV0 = 258,
		TString estimatorFilename = "alien:///alice/cern.ch/user/s/ssakai/Multiplicity_pp13/estimator.root",
		Int_t minNtr = 0,
		Int_t maxNtr = 200,
                Int_t mtype = 0)
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function
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



    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":CaloHFEpp";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskCaloHFEpp* task = new AliAnalysisTaskCaloHFEpp(name.Data());   
    task -> SetEG1(flagEG1);
    task -> SetEG2(flagEG2);
    task -> SetDG1(flagDG1);
    task -> SetDG2(flagDG2);
    task -> SetfFlagClsTypeEMC(SetFlagClsTypeEMC);
    task -> SetfFlagClsTypeDCAL(SetFlagClsTypeDCAL);

    task -> SetTrackEta(TrackEtaMin,TrackEtaMax);
    task -> SetTrackClust(NTPCClust,NITSClust,NCrossedRow);
    task -> SetDCA(DCAxy,DCAz);
    task -> SetNsigma(NsigmaMin,NsigmaMax);
    task -> SetM20(M20Min,M20Max);
    task -> SetEop(EopMin,EopMax);
    task -> SetConeR(coneR);
    task -> SetptAsso(ptAsso);
    task -> SetMimClE(mimcle);
    task -> SetptCut(pte);
    task -> SetMassMin(MassMin);
    task -> SetNref(nref);
    task -> SetNrefV0(nrefV0);
    task -> SetMinNtr(minNtr);
    task -> SetMaxNtr(maxNtr);
    task -> SetEstimatorFile(estimatorFilename);
    task -> SetMultType(mtype);

    /*
    TFile* fEstimator=TFile::Open(estimatorFilename.Data());
    if(!fEstimator){
	    AliFatal("File with multiplicity estimator not found\n");
	    return;
    }

    // Get weight for N_{tracklet}
    TH1D* weightNtrkl = (TH1D*)fEstimator->Get("weightNtrkl")->Clone("weightNtrkl_clone");
    if(!weightNtrkl){
	    AliFatal("Multiplicity estimator for weight not found! Please check your estimator file");
	    return;
    }
    task->SetWeightNtrkl(weightNtrkl);
   */

    if(!task) return 0x0;

    // add your task to the manager
    mgr->AddTask(task);

    TString containerName = mgr->GetCommonFileName();
    containerName += ":PWGHF_hfeCalpp";
    TString SubcontainerName = Form("hfeCalpp");
    SubcontainerName += name;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName, TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput1); 

    return task;
}
