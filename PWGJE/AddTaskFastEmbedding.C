
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
   Double_t minPt = 50.;   Double_t maxPt = 100.;
   Double_t minEta = -0.5; Double_t maxEta = 0.5;
   //Double_t minEta = -0.4; Double_t maxEta = 0.4;  // for LHC10h pass1
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
   task->SetEvtSelJetPtRange(50.,110.);
   //task->SetEvtSelJetEtaRange(-0.4, 0.4); // smaller eta window for LHC10h pass1
   task->SetEvtSelJetEtaRange(-0.5, 0.5);

   task->SetTrackFilterMap(272);

   // event selection
   task->SetOfflineTrgMask(AliVEvent::kMB);
   task->SetCentMin(0.);
   task->SetCentMax(80.);

   //task->SetVtxMin(-10.);
   //task->SetVtxMax(10.);

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

AliAnalysisTaskFastEmbedding* AddTaskFastEmbedding(const char* aodpath, const char* entriespath){

   AliAnalysisTaskFastEmbedding *task = AddTaskFastEmbedding(aodpath, 1);

   Printf("Read entries of aod files from %s", entriespath);
   TArrayI* array = new TArrayI();
   Int_t count = 0;
   Int_t iEntry = -1;
   Int_t iEntrySum = 0;
   Int_t iEntryMax = 0;
   TString line;
   ifstream in;
   in.open(entriespath);
   while(in.good()){
      in >> line;
      iEntry = line.Atoi();

      array->Set(count+1);
      array->AddAt(iEntry,count);
      count++;
      iEntrySum += iEntry;
      if(iEntry>iEntryMax) iEntryMax = iEntry;
   }

   task->SetArrayOfAODEntries(array);
   task->SetAODEntriesSum(iEntrySum);
   task->SetAODEntriesMax(iEntryMax);

   return task;
}

