
//_____________________________________________________________________________________________________
AliAnalysisTask* AddTask_iarsene_TrainTreeAnalysis(Bool_t isGrid=kFALSE, TString prod="LHC10h", Int_t reducedEventType=-1, Bool_t writeTree=kTRUE, TString tasks="dst", TString alienPathForMacros="alien:///alice/cern.ch/user/i/iarsene/analysisMacros") {
   //
   //  AddTask macro for the TreeMaker analysis task and eventual other dependent tasks
   //
   
   TString alienPath = alienPathForMacros;
   TString alirootPath("$ALICE_PHYSICS/PWGDQ/reducedTree/macros");
   
   TObjArray* tasksArray = tasks.Tokenize(";");
   if(tasksArray->GetEntries()==0) {
      printf("ERROR: In AddTask_iarsene_TrainTreeAnalysis(), no tasks are provided, nothing to add to the AliAnalysisManager! \n");
      return 0x0;
   }
   
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // check if the dst TreeMaker is included in the list of tasks and add it to the train if needed
   Bool_t makeTrees = kFALSE;
   AliAnalysisTask* treeMaker = 0x0;
   TString task = tasksArray->At(0)->GetName();   // NOTE: if the tree maker needs to be in the train, then this has to be the first task specified in the tasks string
   if(task.CompareTo("dst")==0) makeTrees = kTRUE;
   
   if(makeTrees) {
      TString addTaskFullPath = "";
      if(isGrid) {
         Int_t errCode = gSystem->Exec(Form("alien_cp %s/AddTask_iarsene_dst.C .", alienPath.Data()));
         if(errCode) {
            printf("ERROR: In AddTask_iarsene_TrainTreeAnalysis(), could not copy %s/AddTask_iarsene_dst.C from grid \n", alienPath.Data());
            return 0x0;
         }
         else addTaskFullPath = gSystem->pwd();
      }
      else {
         addTaskFullPath = alirootPath.Data();
      }
      addTaskFullPath += "/AddTask_iarsene_dst.C";
      if (!gROOT->GetListOfGlobalFunctions()->FindObject("AddTask_iarsene_dst")){
         gROOT->LoadMacro(addTaskFullPath.Data());
      }
      treeMaker = AddTask_iarsene_dst(reducedEventType, writeTree, prod);
      if(tasksArray->GetEntries()==1)     // no other tasks to be included, return
         return treeMaker;
   } // end if(makeTrees)
   
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // check if there are consumer tasks for the tree maker and add them to the train
   AliAnalysisTask* firstConsumerTask = 0x0;
   for(Int_t i=(makeTrees ? 1 : 0); i<tasksArray->GetEntries(); ++i) {
      TString task = tasksArray->At(i)->GetName();
      TString addTaskFullPath = "";
      if(isGrid) {
         Int_t errCode = gSystem->Exec(Form("alien_cp %s/AddTask_iarsene_%s.C .", alienPath.Data(), task.Data()));
         if(errCode) {
            printf("ERROR: In AddTask_iarsene_TrainTreeAnalysis(), could not copy %s/AddTask_iarsene_%s.C from grid \n", alienPath.Data(), task.Data());
            return 0x0;
         }
         else addTaskFullPath = gSystem->pwd();
      }
      else {
         addTaskFullPath = alirootPath.Data();
      }
      addTaskFullPath = Form("%s/AddTask_iarsene_%s.C", addTaskFullPath.Data(), task.Data());
      if (!gROOT->GetListOfGlobalFunctions()->FindObject(Form("AddTask_iarsene_%s", task.Data()))){
         gROOT->LoadMacro(addTaskFullPath.Data());
      }
      
      Int_t runMode = AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents;
      if(!makeTrees) runMode = AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree;
      if(!task.CompareTo("testTask")) {
         if(!firstConsumerTask) firstConsumerTask = AddTask_iarsene_testTask(kTRUE, runMode, prod);
         else AddTask_iarsene_testTask(kTRUE, runMode, prod);
      }
      if(!task.CompareTo("jpsi2ee")) {
         if(!firstConsumerTask) firstConsumerTask = AddTask_iarsene_jpsi2ee(kTRUE, runMode, prod);
         else AddTask_iarsene_jpsi2ee(kTRUE, runMode, prod);
      }
      // NOTE: here add the same for any other potential task
      /*
      if(!firstConsumerTask) {
        gROOT->ProcessLine(Form("firstConsumerTask = AddTask_iarsene_%s(kTRUE,%d,\"%s\")", task.Data(), runMode, prod.Data()));
      }
      else
        gROOT->ProcessLine(Form("AddTask_iarsene_%s(kTRUE,%d,\"%s\")", task.Data(), runMode, prod.Data()));  */
   }
   
   if(treeMaker) return treeMaker;
   if(firstConsumerTask) return firstConsumerTask;
   return 0x0;
}
