AliAnalysisTaskSE *AddTaskNuclMult(Bool_t isMC=kFALSE,
				     Bool_t isNTPCclsVaried=kFALSE, Bool_t isNsigmaTPCVaried=kFALSE,
				     Bool_t isDCAzMaxVaried=kFALSE){

  /*
    I can enable ONE of the last 3 arguments at a time
   */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  //for ESDs
  //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);

  AliPPVsMultUtils *fAliPPVsMultUtils = new AliPPVsMultUtils();

  Int_t NtaskT = GetNtask(isNTPCclsVaried, isNsigmaTPCVaried, isDCAzMaxVaried);
  const Int_t Ntask=NtaskT;
  
  Int_t nTPCclustersMin[Ntask];
  Float_t nsigmaTPCMax[Ntask];
  Float_t DCAxyMax[Ntask];
  Float_t DCAzMax[Ntask];
  
  GetCuts(Ntask, nTPCclustersMin, nsigmaTPCMax, DCAxyMax, DCAzMax,
	  isNTPCclsVaried, isNsigmaTPCVaried, isDCAzMaxVaried);
  
  AliAnalysisNuclMult *task[Ntask];
  for(Int_t i=0;i<Ntask;i++) {
    task[i] = new AliAnalysisNuclMult("AliAnalysisNuclMult");
    mgr->AddTask(task[i]);
  }
  
  for(Int_t i=0;i<Ntask;i++) {
    task[i]->SetIsMC(isMC);
    task[i]->SetPPVsMultUtilsObj(fAliPPVsMultUtils);
    task[i]->SetNTPCclustersMin(nTPCclustersMin[i]);
    task[i]->SetNsigmaTPCcut(nsigmaTPCMax[i]);
    task[i]->SetDCAxyMax(DCAxyMax[i]);
    task[i]->SetDCAzMax(DCAzMax[i]);
  }
  
  AliAnalysisDataContainer *cinput[Ntask];
  AliAnalysisDataContainer *cOutputL[Ntask];
  for(Int_t i=0;i<Ntask;i++) {
    Char_t nameIn[1000];
    Char_t nameOut[1000];
    
    if(!isNTPCclsVaried && !isNsigmaTPCVaried && !isDCAzMaxVaried) {
      snprintf(nameIn,1000,"cchain1_out_DCAxyMax=%.1f",DCAxyMax[i]);
      snprintf(nameOut,1000,"out_DCAxyMax=%.1f",DCAxyMax[i]);
    }
    else if(isNTPCclsVaried) {
      snprintf(nameIn,1000,"cchain1_out_DCAxyMax=%.1f_nTPCclsMin=%d",DCAxyMax[i],nTPCclustersMin[i]);
      snprintf(nameOut,1000,"out_DCAxyMax=%.1f_nTPCclsMin=%d",DCAxyMax[i],nTPCclustersMin[i]);
    }
    else if(isNsigmaTPCVaried) {
      snprintf(nameIn,1000,"cchain1_out_DCAxyMax=%.1f_nsigmaTPCMax=%.1f",DCAxyMax[i],nsigmaTPCMax[i]);
      snprintf(nameOut,1000,"out_DCAxyMax=%.1f_nsigmaTPCMax=%.1f",DCAxyMax[i],nsigmaTPCMax[i]);
    }
    else if(isDCAzMaxVaried) {
      snprintf(nameIn,1000,"cchain1_out_DCAxyMax=%.1f_DCAzMax=%.1f",DCAxyMax[i],DCAzMax[i]);
      snprintf(nameOut,1000,"out_DCAxyMax=%.1f_DCAzMax=%.1f",DCAxyMax[i],DCAzMax[i]);
    }
    
    cinput[i] = mgr->CreateContainer(nameIn,TChain::Class(),AliAnalysisManager::kInputContainer);
    mgr->ConnectInput(task[i],0,mgr->GetCommonInputContainer());
    
    cOutputL[i] = mgr->CreateContainer(nameOut,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(task[i],1,cOutputL[i]);
  }
  
  return task[0];
}
//_____________________________________________________________________________
Int_t GetNtask(Bool_t isNTPCclsVaried, Bool_t isNsigmaTPCVaried, Bool_t isDCAzMaxVaried) {
  
  if(!isNTPCclsVaried && !isNsigmaTPCVaried && !isDCAzMaxVaried) return 3;//3 dcaxy cuts 
  
  else if(isNTPCclsVaried) return 8;//4cuts x 2 (dcaxy cuts)

  else if(isNsigmaTPCVaried) return 4;//2cuts x 2 (dcaxy cuts)

  else if(isDCAzMaxVaried) return 8;//4cuts x 2 (dcaxy cuts)
  
}
//_____________________________________________________________________________
void GetCuts(const Int_t Ntask, Int_t *nTPCclustersMin, Float_t *nsigmaTPCMax, Float_t *DCAxyMax, Float_t *DCAzMax,
	     Bool_t isNTPCclsVaried, Bool_t isNsigmaTPCVaried, Bool_t isDCAzMaxVaried) {
  
  Float_t stdCuts[] = {70, 3, 0.1, 1.};//NTPCcls, NsigmaTPCMax, DCAxyMax, DCAxMax
 
  Int_t NTPCclustersMin[] = {60, 65, 75, 80};
  Float_t NsigmaTPCMax[] = {2.5, 3.5};
  Float_t dcaxyMax[] = {0.5, 1.};
  Float_t dcazMax[] = {0.5, 0.75, 1.5, 2};
  
  if(!isNTPCclsVaried && !isNsigmaTPCVaried && !isDCAzMaxVaried) {
    nTPCclustersMin[0]=(Int_t)stdCuts[0];
    nsigmaTPCMax[0]=stdCuts[1];
    DCAxyMax[0]=stdCuts[2];
    DCAzMax[0]=stdCuts[3];
    
    nTPCclustersMin[1]=(Int_t)stdCuts[0];
    nsigmaTPCMax[1]=stdCuts[1];
    DCAxyMax[1]=dcaxyMax[0];//!
    DCAzMax[1]=stdCuts[3];

    nTPCclustersMin[2]=(Int_t)stdCuts[0];
    nsigmaTPCMax[2]=stdCuts[1];
    DCAxyMax[2]=dcaxyMax[1];//!
    DCAzMax[2]=stdCuts[3];
  }

  else if(isNTPCclsVaried) {
    for(Int_t i=0;i<(Ntask/2);i++) {
      nTPCclustersMin[i]=NTPCclustersMin[i];//!
      nsigmaTPCMax[i]=stdCuts[1];
      DCAxyMax[i]=stdCuts[2];
      DCAzMax[i]=stdCuts[3];
    }
    for(Int_t i=(Ntask/2);i<Ntask;i++) {
      nTPCclustersMin[i]=nTPCclustersMin[i-Ntask/2];//!
      nsigmaTPCMax[i]=stdCuts[1];
      DCAxyMax[i]=dcaxyMax[0];//!
      DCAzMax[i]=stdCuts[3];
    }
  }

  else if(isNsigmaTPCVaried) {
    for(Int_t i=0;i<(Ntask/2);i++) {
      nTPCclustersMin[i]=(Int_t)stdCuts[0];
      nsigmaTPCMax[i]=NsigmaTPCMax[i];//!
      DCAxyMax[i]=stdCuts[2];
      DCAzMax[i]=stdCuts[3];
    }
    for(Int_t i=(Ntask/2);i<Ntask;i++) {
      nTPCclustersMin[i]=(Int_t)stdCuts[0];
      nsigmaTPCMax[i]=NsigmaTPCMax[i-Ntask/2];//!
      DCAxyMax[i]=dcaxyMax[0];//!
      DCAzMax[i]=stdCuts[3];
    }
  }

  else if(isDCAzMaxVaried) {
    for(Int_t i=0;i<(Ntask/2);i++) {
      nTPCclustersMin[i]=(Int_t)stdCuts[0];
      nsigmaTPCMax[i]=stdCuts[1];
      DCAxyMax[i]=stdCuts[2];
      DCAzMax[i]=dcazMax[i];//!
    }
    for(Int_t i=(Ntask/2);i<Ntask;i++) {
      nTPCclustersMin[i]=(Int_t)stdCuts[0];
      nsigmaTPCMax[i]=stdCuts[1];
      DCAxyMax[i]=dcaxyMax[0];//!
      DCAzMax[i]=dcazMax[i-Ntask/2];//!
    }
  }

  return;
}
