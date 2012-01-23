// $Id$

AliEmcalPhysicsSelectionTask* 
AddTaskEmcalPhysicsSelelection(Bool_t exFOnly, Bool_t rejectBG=kTRUE, Bool_t computeBG=kTRUE)
{
  // Add EMCAL physics selection task.

  //get the current analysis manager  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask1PhysSel", "No analysis manager found.");
    return 0;
  }
  isMC = (mgr->GetMCtruthEventHandler()) ? kTRUE:kFALSE; 

  AliEmcalPhysicsSelectionTask *pseltask = new AliEmcalPhysicsSelectionTask("PhysSel");
  //pseltask->SetDoWriteHistos(kFALSE);
  AliEmcalPhysicsSelection *physSel = pseltask->GetPhysicsSelection();
  physSel->SetExcludeFastOnly(exFOnly);
  if (rejectBG) 
    physSel->AddBackgroundIdentification(new AliBackgroundSelection());
  if (computeBG)
    physSel->SetComputeBG(computeBG);
  if (isMC)      
    physSel->SetAnalyzeMC();
  mgr->AddTask(pseltask);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(pseltask, 0, cinput);
  TString oname("EventStat.root");
  AliAnalysisDataContainer *co1 = 
    mgr->CreateContainer("PhysSel",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         oname);
  mgr->ConnectOutput(pseltask,1,co1);
  cout << " *** AliEmcalPhysicsTask configured *** " << endl;
  return task;
}   
