
//_____________________________________________________________________________________________________
AliAnalysisTask* AddTask_TrainTreeAnalysis(Bool_t isGrid=kFALSE, TString prod="LHC10h", Int_t reducedEventType=-1, Bool_t writeTree=kTRUE, TString tasks="dst", TString pathForMacros="$ALICE_PHYSICS/PWGDQ/reducedTree/macros") {
   //
   //  AddTask macro for the TreeMaker analysis task and eventual other dependent tasks
   //
   
   //TString alienPath = alienPathForMacros;
   //TString alirootPath("$ALICE_PHYSICS/PWGDQ/reducedTree/macros");
   
   TObjArray* tasksArray = tasks.Tokenize(";");
   if(tasksArray->GetEntries()==0) {
      printf("ERROR: In AddTask_TrainTreeAnalysis(), no tasks are provided as argument, nothing to add to the AliAnalysisManager! \n");
      return 0x0;
   }
   
   AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
   
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // check if the dst TreeMaker is included in the list of tasks and add it to the train if needed
   Bool_t makeTrees = kFALSE;
   AliAnalysisTask* treeMaker = 0x0;
   // NOTE: if the tree maker needs to be in the train, then this has to be the first task specified in the tasks string
   //         Also, the name specified in the "tasks" argument has to contain the keyword "dst", case insensitive
   TString treeTask = tasksArray->At(0)->GetName();   
   treeTask.ToLower();
   if(treeTask.Contains("dst")) makeTrees = kTRUE;
   
   if(makeTrees) {
      TString addTaskFullPath = "";
      if(isGrid) {
         //Int_t errCode = gSystem->Exec(Form("alien_cp %s/AddTask_%s.C .", alienPath.Data(), tasksArray->At(0)->GetName()));
         Int_t errCode = gSystem->Exec(Form("alien_cp %s/AddTask_%s.C .", pathForMacros.Data(), tasksArray->At(0)->GetName()));
         if(errCode) {
            //printf("ERROR: In AddTask_TrainTreeAnalysis(), could not copy %s/AddTask_%s.C from the grid directory %s \n", alienPath.Data(), tasksArray->At(0)->GetName(), alienPath.Data());
            printf("ERROR: In AddTask_TrainTreeAnalysis(), could not copy %s/AddTask_%s.C from the grid directory %s \n", pathForMacros.Data(), tasksArray->At(0)->GetName(), pathForMacros.Data());
            return 0x0;
         }
         else addTaskFullPath = gSystem->pwd();
      }
      else {
         //addTaskFullPath = alirootPath.Data();
         addTaskFullPath = pathForMacros.Data();
      }
      addTaskFullPath += "/AddTask_";
      addTaskFullPath += tasksArray->At(0)->GetName();
      addTaskFullPath += ".C";
      if (!gROOT->GetListOfGlobalFunctions()->FindObject(Form("AddTask_%s", tasksArray->At(0)->GetName()))){
#ifdef __CLING__
				// ROOT6 version of the Config macro. It cannot handle load and execute macro (compiler error) - need to call via gROOT->ProcessLine(...)
				std::stringstream addtaskload;
				addtaskload << ".L " << addTaskFullPath;
				std::string addtaskloadstring = addtaskload.str();
				std::cout << "Calling Load macro using command string " << addtaskloadstring << std::endl;
				gROOT->ProcessLine(addtaskloadstring.c_str());
#else
				// ROOT5 version, allows loading a macro
				gROOT->LoadMacro(addTaskFullPath.Data());
#endif
			}
      gROOT->ProcessLine(Form("AddTask_%s(%d,%d,\"%s\")", tasksArray->At(0)->GetName(), reducedEventType, (writeTree ? 1 : 0), prod.Data()));
      treeMaker = (AliAnalysisTask*)mgr->GetTasks()->FindObject("DSTTreeMaker");
      if(!treeMaker) {
         printf("ERROR: In AddTask_TrainTreeAnalysis(), something went wrong! The tree maker could not be added to the analysis manager!");
         return 0x0;
      }
      
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
         //Int_t errCode = gSystem->Exec(Form("alien_cp %s/AddTask_%s.C .", alienPath.Data(), task.Data()));
         Int_t errCode = gSystem->Exec(Form("alien_cp %s/AddTask_%s.C .", pathForMacros.Data(), task.Data()));
         if(errCode) {
            //printf("ERROR: In AddTask_TrainTreeAnalysis(), could not copy %s/AddTask_%s.C from the grid directory %s \n", alienPath.Data(), task.Data(), alienPath.Data());
            printf("ERROR: In AddTask_TrainTreeAnalysis(), could not copy %s/AddTask_%s.C from the grid directory %s \n", pathForMacros.Data(), task.Data(), pathForMacros.Data());
            return 0x0;
         }
         else addTaskFullPath = gSystem->pwd();
      }
      else {
         //addTaskFullPath = alirootPath.Data();
         addTaskFullPath = pathForMacros.Data();
      }
      addTaskFullPath = Form("%s/AddTask_%s.C", addTaskFullPath.Data(), task.Data());
      if (!gROOT->GetListOfGlobalFunctions()->FindObject(Form("AddTask_%s", task.Data()))){
#ifdef __CLING__
				// ROOT6 version of the Config macro. It cannot handle load and execute macro (compiler error) - need to call via gROOT->ProcessLine(...)
				std::stringstream addtaskload;
				addtaskload << ".L " << addTaskFullPath;
				std::string addtaskloadstring = addtaskload.str();
				std::cout << "Calling Load macro using command string " << addtaskloadstring << std::endl;
				gROOT->ProcessLine(addtaskloadstring.c_str());
#else
				// ROOT5 version, allows loading a macro
				gROOT->LoadMacro(addTaskFullPath.Data());
#endif
			}
      
      Int_t runMode = AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents;
      if(!makeTrees) runMode = AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree;
      
      gROOT->ProcessLine(Form("AddTask_%s(kTRUE,%d,\"%s\")", tasksArray->At(i)->GetName(), runMode, prod.Data()));
      if(!firstConsumerTask) firstConsumerTask = (AliAnalysisTask*)mgr->GetTasks()->At(mgr->GetTasks()->GetEntries()-1);
   }
   
   if(treeMaker) return treeMaker;
   if(firstConsumerTask) return firstConsumerTask;
   return 0x0;
}
