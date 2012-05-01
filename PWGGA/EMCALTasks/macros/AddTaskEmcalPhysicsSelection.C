// $Id$

AliEmcalPhysicsSelectionTask* AddTaskEmcalPhysicsSelelection(
  Bool_t exFOnly, 
  UInt_t computeBG = 0,
  Bool_t wHistos = kTRUE
)
{
  // Add EMCAL physics selection task.

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalPhysicsSelection", "No analysis manager found.");
    return 0;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEmcalPhysicsSelection", "This task requires an input event handler");
    return NULL;
  }

  Bool_t isMC = (mgr->GetMCtruthEventHandler()) ? kTRUE:kFALSE; 
  AliEmcalPhysicsSelectionTask *pseltask = new AliEmcalPhysicsSelectionTask("EmcalPSel");
  pseltask->SetDoWriteHistos(wHistos);

  AliEmcalPhysicsSelection *physSel = pseltask->GetPhysicsSelection();
  physSel->SetSkipFastOnly(exFOnly);
  if (computeBG)
    physSel->SetComputeBG(computeBG);
  if (isMC)      
    physSel->SetAnalyzeMC();
  mgr->AddTask(pseltask);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("cstatsout",
                TList::Class(),
                AliAnalysisManager::kOutputContainer,
                "EventStat_temp.root");
		
  mgr->ConnectInput(pseltask,  0, cinput);
  mgr->ConnectOutput(pseltask, 1, coutput);

  return pseltask;
}   
