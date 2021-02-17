///////////////////////////////////////////////////////////////////
//                                                               //            
// AddCaloHFEpp                                                  //
// Author: T. Suzuki Univ. of Tsukuba                            //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskCaloHFEpp* AddTaskCaloHFEpp(TString name = "name",
		TString dataname = "dataname",
		Bool_t flagEG1,
		Bool_t flagEG2,
		Bool_t flagDG1,
		Bool_t flagDG2,
		Bool_t SetFlagClsTypeEMC,
		Bool_t SetFlagClsTypeDCAL,
		Bool_t fMC,
		Double_t TrackEtaMin,
		Double_t TrackEtaMax,
		Int_t NTPCClust,
		Int_t NITSClust,
		Int_t NCrossedRow,
		Double_t DCAxy,
		Double_t DCAz,
		Double_t NsigmaMin,
		Double_t NsigmaMax,
		Double_t M20Min,
		Double_t M20Max,
		Double_t EopMin,
		Double_t EopMax,
		Double_t coneR,
		Double_t ptAsso,
		Double_t mimcle,
		TString  pte = "pte",
		Double_t MassMin,
		Double_t nref,
		TString estimatorFilename,
		Int_t minNtr,
		Int_t maxNtr)
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
    task -> SetMinNtr(minNtr);
    task -> SetMaxNtr(maxNtr);

    TFile* fEstimator=TFile::Open(estimatorFilename.Data());
    if(!fEstimator){
	    AliFatal("File with multiplicity estimator not found\n");
	    return;
    }

    // MB get estimator file
    if(SetFlagClsTypeEMC && !flagEG1 && !flagEG2){
	    const Char_t* profilebasename="SPDTrklMB";
	    const Char_t* periodNames[11] = {"LHC16i", "LHC16j","LHC16k","LHC16o","LHC17h","LHC17i","LHC17k","LHC17l","LHC17m","LHC17o","LHC17r"};

	    TProfile* multEstimatorAvgMB[11];

	    for(Int_t ip=0; ip<11; ip++) {
		    cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
		    multEstimatorAvgMB[ip] = (TProfile*)(fEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
		    if(!multEstimatorAvgMB[ip]){
			    AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
			    return;
		    }
	    }
	    task->SetMultiProfileLHC16i(multEstimatorAvgMB[0]);
	    task->SetMultiProfileLHC16j(multEstimatorAvgMB[1]);
	    task->SetMultiProfileLHC16k(multEstimatorAvgMB[2]);
	    task->SetMultiProfileLHC16o(multEstimatorAvgMB[3]);
	    task->SetMultiProfileLHC17h(multEstimatorAvgEG1[4]);
	    task->SetMultiProfileLHC17i(multEstimatorAvgEG1[5]);
	    task->SetMultiProfileLHC17k(multEstimatorAvgEG1[6]);
	    task->SetMultiProfileLHC17l(multEstimatorAvgEG1[7]);
	    task->SetMultiProfileLHC17m(multEstimatorAvgEG1[8]);
	    task->SetMultiProfileLHC17o(multEstimatorAvgEG1[9]);
	    task->SetMultiProfileLHC17r(multEstimatorAvgEG1[10]);
    }

    // EG1 get estimator file
    if(SetFlagClsTypeEMC && flagEG1 && !flagEG2){
	    const Char_t* profilebasename="SPDTrklEG1";
	    const Char_t* periodNames[11] = {"LHC16i", "LHC16j","LHC16k","LHC16o","LHC17h","LHC17i","LHC17k","LHC17l","LHC17m","LHC17o","LHC17r"};

	    TProfile* multEstimatorAvgEG1[11];

	    for(Int_t ip=0; ip<11; ip++) {
		    cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
		    multEstimatorAvgEG1[ip] = (TProfile*)(fEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
		    if(!multEstimatorAvgEG1[ip]){
			    AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
			    return;
		    }
	    }
	    task->SetMultiProfileLHC16i(multEstimatorAvgEG1[0]);
	    task->SetMultiProfileLHC16j(multEstimatorAvgEG1[1]);
	    task->SetMultiProfileLHC16k(multEstimatorAvgEG1[2]);
	    task->SetMultiProfileLHC16o(multEstimatorAvgEG1[3]);
	    task->SetMultiProfileLHC17h(multEstimatorAvgEG1[4]);
	    task->SetMultiProfileLHC17i(multEstimatorAvgEG1[5]);
	    task->SetMultiProfileLHC17k(multEstimatorAvgEG1[6]);
	    task->SetMultiProfileLHC17l(multEstimatorAvgEG1[7]);
	    task->SetMultiProfileLHC17m(multEstimatorAvgEG1[8]);
	    task->SetMultiProfileLHC17o(multEstimatorAvgEG1[9]);
	    task->SetMultiProfileLHC17r(multEstimatorAvgEG1[10]);
    }

    // EG2 get estimator file
    if(SetFlagClsTypeEMC && !flagEG1 && flagEG2){
	    const Char_t* profilebasename="SPDTrklEG2";
	    const Char_t* periodNames[11] = {"LHC16i", "LHC16j","LHC16k","LHC16o"};

	    TProfile* multEstimatorAvgEG2[11];

	    for(Int_t ip=0; ip<11; ip++) {
		    cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
		    multEstimatorAvgEG2[ip] = (TProfile*)(fEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
		    if(!multEstimatorAvgEG2[ip]){
			    AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
			    return;
		    }
	    }
	    task->SetMultiProfileLHC16i(multEstimatorAvgEG2[0]);
	    task->SetMultiProfileLHC16j(multEstimatorAvgEG2[1]);
	    task->SetMultiProfileLHC16k(multEstimatorAvgEG2[2]);
	    task->SetMultiProfileLHC16o(multEstimatorAvgEG2[3]);
    }


    // MC get estimator file
    const Char_t* profilebasenameMC="SPDTrklMC";
    const Char_t* periodNamesMC[2] = {"LHC16k","LHC16l"};
    TProfile* multEstimatorAvgMC[2];

    for(Int_t ip=0; ip<2; ip++) {
	    cout<< " Trying to get "<<Form("%s_%s",profilebasenameMC,periodNamesMC[ip])<<endl;
	    multEstimatorAvgMC[ip] = (TProfile*)(fEstimator->Get(Form("%s_%s",profilebasenameMC,periodNamesMC[ip]))->Clone(Form("%s_%s_clone",profilebasenameMC,periodNamesMC[ip])));
	    if(!multEstimatorAvgMC[ip]){
		    AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNamesMC[ip]));
		    return;
	    }
    }
    task->SetMultiProfileMCLHC16k(multEstimatorAvgMC[0]);
    task->SetMultiProfileMCLHC16l(multEstimatorAvgMC[1]);


    // Get weight for N_{tracklet}
    TH1D* weightNtrkl = (TH1D*)fEstimator->Get("weightNtrkl")->Clone("weightNtrkl_clone");
    if(!weightNtrkl){
	    AliFatal("Multiplicity estimator for weight not found! Please check your estimator file");
	    return;
    }
    task->SetWeightNtrkl(weightNtrkl);


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
