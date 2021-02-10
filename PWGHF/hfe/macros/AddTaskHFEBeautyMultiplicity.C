///////////////////////////////////////////////////////////////////
//                                                               //            
// AddTaskHFEBeautyMultiplicity                                  //
// Author: Shunya Chiba,                                         //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTask* AddTaskHFEBeautyMultiplicity(
	TString name = "name",
	TString estimatorFilename = "alien:///alice/cern.ch/user/s/schiba/Mult_pPb16qt/estimatorArg.root",
	Double_t nref = 30.1

	)
{


	TFile* fEstimator = TFile::Open(estimatorFilename.Data());
	if(!fEstimator){
		AliFatal("File with multiplicity estimator not found\n");
		return;
	}

	TProfile* multEstimatorAvgMB;
	multEstimatorAvgMB = (TProfile*)(fEstimator->Get("mean_Trklet"))->Clone("multEstimatorAvgMB");
	task->SetMultiProfileLHC16qt(multEstimatorAvgMB);

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
    // now we create an instance of your task

    AliAnalysisTaskHFEBeautyMultiplicity* task = new AliAnalysisTaskHFEBeautyMultiplicity(name.Data());
    task -> SetNref(nref);


    if(!task) return 0x0;
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
