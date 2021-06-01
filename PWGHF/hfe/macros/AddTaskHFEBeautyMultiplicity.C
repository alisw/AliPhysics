///////////////////////////////////////////////////////////////////
//                                                               //            
// AddTaskHFEBeautyMultiplicity                                  //
// Author: Shunya Chiba,                                         //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTask* AddTaskHFEBeautyMultiplicity(
	TString name = "name",
	Double_t nref = 30.1,
	Double_t minNtrklet = 0,
	Double_t maxNtrklet = 9999,
	Bool_t   iGPMC = kFALSE,
	TString estimatorFilename = "alien:///alice/cern.ch/user/s/schiba/Mult_pPb16qt/estimatorAvg.root"

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
    task->SetNref(nref);
    task->SetNtrkletMin(minNtrklet);
    task->SetNtrkletMax(maxNtrklet);
    task->SetMCtype(iGPMC);




    //---- Get Estimator 
    TFile* fEstimator = TFile::Open(estimatorFilename.Data());
    if(!fEstimator){
		AliFatal("File with multiplicity estimator not found!\n");
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
