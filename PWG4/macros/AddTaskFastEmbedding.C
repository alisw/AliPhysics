
AliAnalysisTaskFastEmbedding* AddTaskFastEmbedding(){

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if(!mgr){
	::Error("AddTaskCentralitySelection", "No analysis manager to connect ot.");
	return NULL;
    }
    if(!mgr->GetInputEventHandler()){
        ::Error("AddTaskCentralitySelection", "This task requires an input event handler.");
	return NULL;
    }


    AliAnalysisTaskFastEmbedding *task = new AliAnalysisTaskFastEmbedding("FastEmbedding");
    // ## set embedding mode ##
    // kAODFull=0, kAODJetTracks, kAODJet4Mom, kToySingle4Mom
    task->SetEmbedMode(AliAnalysisTaskFastEmbedding::kToyTracks);

    // ## set ranges for toy ##
    //SetToyTrackRanges(
    Double_t minPt = 20.;   Double_t maxPt = 200.;
    Double_t minEta = -0.5; Double_t maxEta = 0.5;
    Double_t minPhi = 0.;   Double_t maxPhi = 2*TMath::Pi();
    //fToyDistributionTrackPt: 0 = uniform distribution
    //                         else = exponential / power law (not implemented yet)
    //task->SetToyNumberOfTrackRange(4,4);
    //task->SetToyTrackRanges(0.15, 300., 5,-.9, .9, 0., 2*TMath::Pi());
    task->SetToyTrackRanges(minPt,maxPt,0.,minEta,maxEta,minPhi,maxPhi);
    task->SetToyFilterMap((1<<32)-1);

    // ## set event selection for events of the addition AOD ##
    // kEventsAll=0; kEventsJetPt
    task->SetEvtSelecMode(AliAnalysisTaskFastEmbedding::kEventsJetPt);

    // ## set jet pT range for event selection ##
    // SetEvtSelJetPtRange(Float_t minPt, Float_t maxPt)
    task->SetEvtSelJetPtRange(20.,-1.);

    mgr->AddTask(task);

    // ## create the output containers ##
    AliAnalysisDataContainer *coutputFastEmbedding = mgr->CreateContainer(
         "fastembedding", TList::Class(), AliAnalysisManager::kOutputContainer,
         Form("%s:PWG4_FastEmbedding", AliAnalysisManager::GetCommonFileName()));

    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
    mgr->ConnectOutput(task, 1, coutputFastEmbedding);


    return task;

}

AliAnalysisTaskFastEmbedding* AddTaskFastEmbedding(TObjArray* aodarray){

    AliAnalysisTaskFastEmbedding *task = AddTaskFastEmbedding();
    if(aodarray){
      task->SetArrayOfAODPaths(aodarray);
      task->SetEmbedMode(AliAnalysisTaskFastEmbedding::kAODFull);
    }

    return task;
}


AliAnalysisTaskFastEmbedding* AddTaskFastEmbedding(const char* filepath, Int_t mode = 0){

    AliAnalysisTaskFastEmbedding *task = AddTaskFastEmbedding();
    if(strlen(filepath)){

       if(mode==0){ // path to single AOD
          task->SetAODPath(filepath);
       }
       if(mode==1){ // path to text file with list of paths of multiple AODs
           Printf("Read aod paths from file %s", filepath);
           TObjArray* array = new TObjArray();
           TObjString* ostr = 0;
           TString line;
           ifstream in;
           in.open(filepath);
           while(in.good()){
              in >> line;
              if(line.Length() == 0) continue;
              Printf("found aod path %s", line.Data());
              ostr = new TObjString(line.Data());
              array->Add(ostr);
           }
           Printf("-> %d aod paths found", array->GetEntries());
           
           task->SetArrayOfAODPaths(array);
       }

       task->SetEmbedMode(AliAnalysisTaskFastEmbedding::kAODFull);
    }

    return task;
}

