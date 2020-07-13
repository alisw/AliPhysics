class AliAnalysisDataContainer;
class AliAnalysisTaskESEFlow;

AliAnalysisTaskESEFlow* AddESEFlowTask(AliAnalysisTaskESEFlow::ColSystem colSys, TString sWeightsFile = "", TString sVWeights = "",TString V0Calib = "", TString qSplines="", const char* suffix = "")
{
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
    TString name = Form("ESEFlow%s",suffix);

    Bool_t bUseOwnWeights = kFALSE;
    if(!sWeightsFile.IsNull()) { bUseOwnWeights = kTRUE; }

    Bool_t bUseEseSp       = kFALSE;
    if(!qSplines.IsNull()) { bUseEseSp = kTRUE; }

    Bool_t bUseV0Calibration = kFALSE;
    if(!V0Calib.IsNull()) { bUseV0Calibration = kTRUE; } 

    Bool_t bUseRBRWeights = kFALSE;
    if(!sVWeights.IsNull() || !sWeightsFile.IsNull()) { bUseRBRWeights = kTRUE; }
         


    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += Form(":%s", suffix);      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskESEFlow* task = new AliAnalysisTaskESEFlow(name.Data(), colSys, bUseV0Calibration);   
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("%s:EseTestFlow",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("%s:Observables",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("%s:Flow",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,4,mgr->CreateContainer(Form("%s:DiffFlow",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,5,mgr->CreateContainer(Form("%s:q_n",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,6,mgr->CreateContainer(Form("%s:DiffFlowESETPC",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,7,mgr->CreateContainer(Form("%s:FlowESETPC",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,8,mgr->CreateContainer(Form("%s:DiffFlowESEV0C",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,9,mgr->CreateContainer(Form("%s:FlowESEV0C",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,10,mgr->CreateContainer(Form("%s:DiffFlowESEV0A",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,11,mgr->CreateContainer(Form("%s:FlowESEV0A",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,12,mgr->CreateContainer(Form("%s:SP_Flow",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,13,mgr->CreateContainer(Form("%s:SP_FlowESE",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,14,mgr->CreateContainer(Form("%s:fQAEvents",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid


    task->SetFilterBit(96);
    task->SetVtxZCut(10.0);
    task->SetPhiBins(60);
    task->SetEtaBins(32);
    const int nPtBins = 28;
    Double_t PtEdges[nPtBins+1] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,7.0,8.0,9.0,10.0};
    task->SetPtBins(nPtBins,PtEdges);
    const Int_t nCentBins = 10;
    Double_t CentEdges[nCentBins+1] = {0, 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
    task->SetCentBin(nCentBins,CentEdges);
    task->SetReadMC(kFALSE); //activate monte carlo analysis
    task->SetAbsEta(0.8);
    task->SetRejectAddPileUp(kTRUE);
    task->SetFlowRFPsPt(0.2,5.0);
    task->SetFlowPOIsPt(0.0,10.0);
    task->SetRedFlowPt(0.2,20.0);
    task->SetTPCEse(kTRUE);
    task->SetV0CEse(kTRUE);
    task->SetV0AEse(kTRUE);
    task->SetSampling(kFALSE,1); //(kFALSE,1) for no sampling and only 1 sample size
    task->SetHasEtaGap(kTRUE);
    task->SetEtaGap(1.0);
    task->SetTPCEseqnBins(100,0.0,8.0);
    task->SetV0EseqnBins(100,0.0,15.0);
    task->SetChi2TPCFl(kFALSE, 4.0); // change to kTRUE for systematic, default in track cut is 4.0, so change to i.e. 3
    task->SetChi2ITSFl(kFALSE, 36.0); // change to kTRUE for systematic, default in track cut is 36.0, so change to i.e. 35
    task->SetQARejFiller(kFALSE);
    task->SetNUEWeights(kFALSE, 1);
    task->Set2018(kFALSE);
    task->SetBayesUnfoldingInput(kFALSE);

    const Int_t nEsePercentiles = 10;
    Double_t EseEdges[nEsePercentiles+1] = {0, 10., 20., 30., 40., 50., 60., 70., 80., 90.,100.};
    task->SetEventShapeBins(nEsePercentiles,EseEdges);

    if( colSys == AliAnalysisTaskESEFlow::ColSystem::kPbPb){
      task->SetCentralityEst("V0M"); // V0M
      task->SetChargedNumTPCclsMin(70);
    }
    else if ( colSys == AliAnalysisTaskESEFlow::ColSystem::kPPb){
      task->SetCentralityEst("V0A");
      task->SetChargedNumTPCclsMin(70);
    }
    else if ( colSys == AliAnalysisTaskESEFlow::ColSystem::kXeXe){
      task->SetCentralityEst("V0M");
      task->SetChargedNumTPCclsMin(70);
    }

    // RUN BY RUN own

    if(bUseRBRWeights){
      if(bUseOwnWeights) 
      {
        task->SetWeights(bUseOwnWeights);

        TObjArray* taskContainers = mgr->GetContainers();
        if(!taskContainers) { printf("E-AddTaskESEFlow: Task containers does not exists!\n"); return NULL; }

        // check if the input weights are already loaded (e.g. in different subwagon)
        AliAnalysisDataContainer* weights = (AliAnalysisDataContainer*) taskContainers->FindObject("inputWeights");
        if(!weights) 
        {  
          // if it does not exists create it

          // in case of non-local run, establish connection to ALiEn for loading the weights
          if(sWeightsFile.Contains("alien://")) { gGrid->Connect("alien://"); }

          TFile* weights_file = TFile::Open(sWeightsFile.Data(),"READ");
          if(!weights_file) { printf("E-AddTaskESEFlow: Input file with weights not found!\n"); return NULL; }

          TList* weights_list = static_cast<TList*>(weights_file->Get("weights"));
          if(!weights_list) { printf("E-AddTaskESEFlow: Input list with weights not found!\n"); weights_file->ls(); return NULL; }

          AliAnalysisDataContainer* cInputWeights = mgr->CreateContainer("inputWeights",TList::Class(), AliAnalysisManager::kInputContainer);
          cInputWeights->SetData(weights_list);
          mgr->ConnectInput(task,1,cInputWeights);
        }
        else 
        {
          // connect existing container
          mgr->ConnectInput(task,1,weights);
        }
      }
      if(!bUseOwnWeights)
      {
        TObjArray* taskContainersVy = mgr->GetContainers();
        if(!taskContainersVy) { printf("E-AddTaskESEFlow: Task containers does not exists!\n"); return NULL; }

        // check if the input weights are already loaded (e.g. in different subwagon)
        AliAnalysisDataContainer* weightsVy = (AliAnalysisDataContainer*) taskContainersVy->FindObject("inputWeights");
        if(!weightsVy) {  
          // if it does not exists create it
          // in case of non-local run, establish connection to ALiEn for loading the weights
          if(sVWeights.Contains("alien://")) { gGrid->Connect("alien://"); }

          TFile* weights_fileVy = TFile::Open(sVWeights.Data(),"READ");
          if(!weights_fileVy) { printf("E-AddTaskESEFlow: Input file with weights not found!\n"); return NULL; }

          TList* weights_listVy = static_cast<TList*>(weights_fileVy->Get("WeightList"));
          if(!weights_listVy) { printf("E-AddTaskESEFlow: Input list with weights not found!\n"); weights_fileVy->ls(); return NULL; }

          AliAnalysisDataContainer* cInputWeightsVy = mgr->CreateContainer("inputWeights",TList::Class(), AliAnalysisManager::kInputContainer);
          cInputWeightsVy->SetData(weights_listVy);
          mgr->ConnectInput(task,1,cInputWeightsVy);
        }
        else
        {
          // connect existing container
          mgr->ConnectInput(task,1,weightsVy);
        }
      }
    }
    else{
      task->SetMakeRBRweights(kTRUE);
    }
    

    if(bUseV0Calibration){
      TObjArray* taskContainersV0 = mgr->GetContainers();
      if(!taskContainersV0) { printf("E-AddTaskESEFlow: Task containers does not exists!\n"); return NULL; }

      AliAnalysisDataContainer* V0Cal = (AliAnalysisDataContainer*) taskContainersV0->FindObject("inputV0Cal");
      if(!V0Cal) {  
        if(V0Calib.Contains("alien://")) { gGrid->Connect("alien://"); }

        TFile* V0Cal_file = TFile::Open(V0Calib.Data(),"READ");
        if(!V0Cal_file) { printf("E-AddTaskESEFlow: Input file with V0 Calibration not found!\n"); return NULL; }

        TList* V0Cal_list = static_cast<TList*>(V0Cal_file->Get("Calibration"));
        if(!V0Cal_list) { printf("E-AddTaskESEFlow: Input list with V0 Calibration not found!\n"); V0Cal_file->ls(); return NULL; }

        AliAnalysisDataContainer* cInputV0Cal = mgr->CreateContainer("inputV0Cal",TList::Class(), AliAnalysisManager::kInputContainer);
        cInputV0Cal->SetData(V0Cal_list);
        mgr->ConnectInput(task,2,cInputV0Cal);
      }
      else
      {
        mgr->ConnectInput(task,2,V0Cal);
      }
    }

    if(bUseEseSp){
      TObjArray* taskContainersSp = mgr->GetContainers();
      if(!taskContainersSp) { printf("E-AddTaskESEFlow: Task containers does not exists!\n"); return NULL; }

      AliAnalysisDataContainer* qSelSp = (AliAnalysisDataContainer*) taskContainersSp->FindObject("inputqSp");
      if(!qSelSp) {  
        if(qSplines.Contains("alien://")) { gGrid->Connect("alien://"); }

        TFile* qSp_file = TFile::Open(qSplines.Data(),"READ");
        if(!qSp_file) { printf("E-AddTaskESEFlow: Input file with q-selection Splines not found!\n"); return NULL; }

        TList* qSp_list = static_cast<TList*>(qSp_file->Get("qSplines"));
        if(!qSp_list) { printf("E-AddTaskESEFlow: Input list with q Splines not found!\n"); qSp_file->ls(); return NULL; }

        AliAnalysisDataContainer* cInputqSp = mgr->CreateContainer("inputqSp",TList::Class(), AliAnalysisManager::kInputContainer);
        cInputqSp->SetData(qSp_list);
        mgr->ConnectInput(task,3,cInputqSp);
      }
      else
      {
        mgr->ConnectInput(task,3,qSelSp);
      }
    }
    else{ 
      task->SetMakeqSelectionRun(kTRUE); 
    }
    

  return task;
}
