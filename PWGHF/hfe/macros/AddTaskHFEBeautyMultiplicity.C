///////////////////////////////////////////////////////////////////
//                                                               //            
// AddTaskHFEBeautyMultiplicity                                  //
// Author: Shunya Chiba,                                         //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTask* AddTaskHFEBeautyMultiplicity(
	TString name = "name",
	Double_t EtaMin = -0.6,
	Double_t EtaMax = 0.6,
	Double_t NsigmaMin = -1.0,
	Double_t NsigmaMax = 3.0,
	Double_t HadNsigma = -3.5,
	Double_t M20Min = 0.015,
	Double_t M20Max = 0.3,
	Double_t EopMin = 0.8,
	Double_t EopMax = 1.2,
	Double_t DCAxy = 2.4,
	Double_t DCAz = 3.2,
	Double_t Diff = 0.05,
	Int_t NTPCClust = 100,
	Int_t NITSClust = 3,
	Int_t NCrossedRow = 100,
	Double_t TPCdEdx = 80.0,
	Double_t PhotInvMass = 0.15,
	Double_t PhotMinPt = 0.1,
	Double_t nref = 30.1,
	Double_t minNtrklet = 0,
	Double_t maxNtrklet = 9999,
	Bool_t   iGPMC = kFALSE,
	TString estimatorFilename = "alien:///alice/cern.ch/user/s/schiba/Mult_pPb16qt/estimatorAvg.root",
	//TString pTWeightFilename = "alien:///alice/cern.ch/user/s/schiba/Mult_pPb16qt/pTWeight.root"
	TString pTWeightFilename = "alien:///alice/cern.ch/user/s/schiba/Mult_pPb16qt/pTWeight_Fix.root"

	)
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
    fileName += ":MyTask";      // create a subfolder in the file
    

    //---- my task
    AliAnalysisTaskHFEBeautyMultiplicity* task = new AliAnalysisTaskHFEBeautyMultiplicity(name.Data());
    task->SetTrackEta(EtaMin, EtaMax);
    task->SetNsigma(NsigmaMin, NsigmaMax, HadNsigma);
    task->SetM20(M20Min, M20Max);
    task->SetEop(EopMin, EopMax);
    task->SetDCA(DCAxy, DCAz);
    task->SetTrackClust(NTPCClust, NITSClust, NCrossedRow, TPCdEdx);
    task->SetDiff(Diff);
    task->SetMass(PhotInvMass, PhotMinPt);
    task->SetNref(nref);
    task->SetNtrkletMin(minNtrklet);
    task->SetNtrkletMax(maxNtrklet);
    task->SetMCtype(iGPMC);



    //---- Get Estimator 
    TFile* fEstimator = TFile::Open(estimatorFilename.Data());
    if(!fEstimator){
		AliFatal("File with multiplicity estimator is not found!\n");
		return;
    }

    TProfile* multiEstimatorAvgMB[2];
    multiEstimatorAvgMB[0] = (TProfile*)(fEstimator->Get("mean_Trklet"))->Clone("multiEstimatorAvgMB_Data");	// Get Estimator file (Data)
    multiEstimatorAvgMB[1] = (TProfile*)(fEstimator->Get("mean_Trklet_MC"))->Clone("multiEstimatorAvgMB_MC");	// Get Estimator file (MC)

    task->SetMultiProfileLHC16qt(multiEstimatorAvgMB[0]);
    task->SetMultiProfileLHC16qt_MC(multiEstimatorAvgMB[1]);


    //---- Get weight for N tracklet
    TH1D* WeightNtrklet = (TH1D*)fEstimator->Get("weightNtrkl")->Clone("WeightNtrklet");
    if(!WeightNtrklet){
	    AliFatal("Multiplicity estimator for weight not found! Please check your estimator file.\n");
	    return;
    }
    task->SetWeightNtrklet(WeightNtrklet);


    //---- Get Pt Weight
    TFile* fpTWeight = TFile::Open(pTWeightFilename.Data());
    if(!fpTWeight){
		AliFatal("File with pT weight is not found!\n");
		return;
    }
    TGraphErrors* WeightPt_Dmeson = (TGraphErrors*)fpTWeight->Get("pTWeight_Dmeson")->Clone("WeightPt_Dmeson");
    TGraphErrors* WeightPt_Lc     = (TGraphErrors*)fpTWeight->Get("pTWeight_Lc")->Clone("WeightPt_Lc");
    TGraphErrors* WeightPt_Bmeson = (TGraphErrors*)fpTWeight->Get("pTWeight_Bmeson")->Clone("WeightPt_Bmeson");
    //TF1* 	  WeightPt_Bmeson = (TF1*)fpTWeight->Get("pTWeight_Bmeson")->Clone("WeightPt_Bmeson");
    //TF1* 	  WeightPt_Pi0 	  = (TF1*)fpTWeight->Get("pTWeight_Pi0")->Clone("WeightPt_Pi0");
    //TF1* 	  WeightPt_Eta    = (TF1*)fpTWeight->Get("pTWeight_Eta")->Clone("WeightPt_Eta");

    task->SetWeightDmeson(WeightPt_Dmeson);
    task->SetWeightLc(WeightPt_Lc);
    task->SetWeightBmeson(WeightPt_Bmeson);
    //task->SetWeightPi0(WeightPt_Pi0);
    //task->SetWeightEta(WeightPt_Eta);




    if(!task) return 0x0;

    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid

    return task;
}
