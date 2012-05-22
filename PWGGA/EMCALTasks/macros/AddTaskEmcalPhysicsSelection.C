// $Id$

AliEmcalPhysicsSelectionTask* AddTaskEmcalPhysicsSelection(
  Bool_t exFOnly, 
  Bool_t wHistos   = kTRUE,
  UInt_t triggers  = 0,
  Double_t minE    = -1,
  Double_t minPt   = -1
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
  if (physSel) {
    physSel->SetSkipFastOnly(exFOnly);
    if (isMC)      
      physSel->SetAnalyzeMC();
    physSel->SetClusMinE(minE);
    physSel->SetTrackMinPt(minPt);
    physSel->SetTriggers(triggers);
  } else {
    ::Error("AddTaskEmcalPhysicsSelection", "No AliEmcalPhysicsSelection object found.");
  }

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
