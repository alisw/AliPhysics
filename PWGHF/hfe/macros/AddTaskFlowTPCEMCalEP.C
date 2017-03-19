AliAnalysisTask *AddTaskFlowTPCEMCalEP(Double_t SigmaITScut, Double_t SigmaTOFcut, Double_t SigmaTPCcut, Double_t AssPtCut, Int_t AssTPCnCut, Int_t ITSncut, Bool_t AssITSrefitCut, Int_t TPCnCut, Int_t period, Bool_t UseNewEP, Bool_t UseTender, TString ID="ContName", TString passV0, TString passTPC, Bool_t TimeCut,Bool_t WeightSyst,Bool_t SystTOFcut, Double_t CutM02, Double_t CutM20, Bool_t SScut, Bool_t EnablePileupRejVZEROTPCout)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskFlowTPCEMCalEP", "No analysis manager found.");
        return NULL;
    }
    
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskFlowTPCEMCalEP", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    //    if (type=="AOD"){
    //        ::Error("AddTaskFlowTPCEMCalEP", "The tasks exits because AODs are in input");
    //        return NULL;
    //    }
    Bool_t MCthere=kFALSE;
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
    if(!mcH){
        MCthere=kFALSE;
    }else{
        MCthere=kTRUE;
    }
    
    if (!UseNewEP){
        //Event plane task
        AliEPSelectionTask *eventplaneTask = new AliEPSelectionTask("EventplaneSelection");
        eventplaneTask->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kSemiCentral | AliVEvent::kCentral | AliVEvent::kEMCEGA | AliVEvent::kEMCEJE);
        
        eventplaneTask->SetTrackType("TPC");
        eventplaneTask->SetUsePtWeight();
        eventplaneTask->SetUsePhiWeight();
        eventplaneTask->SetSaveTrackContribution();
        
        mgr->AddTask(eventplaneTask);
        
        TString containerName0 = mgr->GetCommonFileName();
        containerName0 += ":PWGHF_hfeCalEventPlane";
        containerName0 += ID;
        
        TString name0 = "EPStat";
        name0 += ID;
        
        AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name0.Data(),TList::Class(), AliAnalysisManager::kOutputContainer,containerName0.Data());
        mgr->ConnectInput(eventplaneTask, 0, mgr->GetCommonInputContainer());
        mgr->ConnectOutput(eventplaneTask,1,coutput1);
    }
    
    Bool_t Is2015 = kTRUE;
    
    if(Is2015){
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/ConfigHFE_FLOW_TPCEMCal_EP.C");
        AliAnalysisTaskFlowTPCEMCalEP *taskMB = ConfigHFE_FLOW_TPCEMCal_EP(MCthere,SigmaITScut,SigmaTOFcut,SigmaTPCcut,AssPtCut,AssTPCnCut,ITSncut,AssITSrefitCut,TPCnCut,UseNewEP,UseTender,period,passV0,passTPC,TimeCut,WeightSyst,SystTOFcut,CutM02,CutM20,SScut,EnablePileupRejVZEROTPCout);
        
        mgr->AddTask(taskMB);
        
        taskMB->SelectCollisionCandidates(AliVEvent::kINT7);
        
        TString containerName1 = mgr->GetCommonFileName();
        containerName1 += ":PWGHF_hfeCalcorrINT7V2";
        containerName1 += ID;
        
        TString name1 = "histMB";
        name1 += ID;
        
        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
        mgr->ConnectInput(taskMB, 0, cinput);
        //        mgr->ConnectInput(taskMB, 1, corrTask);
        mgr->ConnectOutput(taskMB, 1, coutput1);
    }
    
    if(!Is2015){
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/ConfigHFE_FLOW_TPCEMCal_EP.C");
        
        AliAnalysisTaskFlowTPCEMCalEP *taskMB = ConfigHFE_FLOW_TPCEMCal_EP(MCthere,SigmaITScut,SigmaTOFcut,SigmaTPCcut,AssPtCut,AssTPCnCut,ITSncut,AssITSrefitCut,TPCnCut,UseNewEP,UseTender,period,passV0,passTPC,TimeCut,WeightSyst,SystTOFcut,CutM02,CutM20,SScut,EnablePileupRejVZEROTPCout);
        AliAnalysisTaskFlowTPCEMCalEP *taskcorrMB = ConfigHFE_FLOW_TPCEMCal_EP(MCthere,SigmaITScut,SigmaTOFcut,SigmaTPCcut,AssPtCut,AssTPCnCut,ITSncut,AssITSrefitCut,TPCnCut,UseNewEP,UseTender,period,passV0,passTPC,TimeCut,WeightSyst,SystTOFcut,CutM02,CutM20,SScut,EnablePileupRejVZEROTPCout);
        AliAnalysisTaskFlowTPCEMCalEP *taskTR = ConfigHFE_FLOW_TPCEMCal_EP(MCthere,SigmaITScut,SigmaTOFcut,SigmaTPCcut,AssPtCut,AssTPCnCut,ITSncut,AssITSrefitCut,TPCnCut,UseNewEP,UseTender,period,passV0,passTPC,TimeCut,WeightSyst,SystTOFcut,CutM02,CutM20,SScut,EnablePileupRejVZEROTPCout);
        
        mgr->AddTask(taskcorrMB);
        mgr->AddTask(taskMB);
        mgr->AddTask(taskTR);
        
        // Flattened semi central trigger
        
        taskcorrMB->SelectCollisionCandidates(AliVEvent::kAny);
        
        TString containerName1 = mgr->GetCommonFileName();
        containerName1 += ":PWGHF_hfeCalcorrSemiCentralV2";
        containerName1 += ID;
        
        TString name1 = "histcorrMB";
        name1 += ID;
        
        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
        mgr->ConnectInput(taskcorrMB, 0, cinput);
        mgr->ConnectOutput(taskcorrMB, 1, coutput1);
        
        // Central trigger
        taskMB->SelectCollisionCandidates(AliVEvent::kSemiCentral | AliVEvent::kCentral);
        
        TString containerName2 = mgr->GetCommonFileName();
        containerName2 += ":PWGHF_hfeCalCentralV2";
        containerName2 += ID;
        
        TString name2 = "histMB";
        name2 += ID;
        
        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name2.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
        mgr->ConnectInput(taskMB, 0, cinput);
        mgr->ConnectOutput(taskMB, 1, coutput1);
        
        //L1 gamma trigger
        taskTR->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        
        TString containerName3 = mgr->GetCommonFileName();
        containerName3 += ":PWGHF_hfeCalL1GammaV2";
        containerName3 += ID;
        
        TString name3 = "histTR";
        name3 += ID;
        
        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name3.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName3.Data());
        mgr->ConnectInput(taskTR, 0, cinput);
        mgr->ConnectOutput(taskTR, 1, coutput1);
        
    }
    
    if(MCthere){
        
        AliAnalysisTaskFlowTPCEMCalEP *taskMC = ConfigHFE_FLOW_TPCEMCal_EP(MCthere,SigmaITScut,SigmaTOFcut,SigmaTPCcut,AssPtCut,AssTPCnCut,ITSncut,AssITSrefitCut,TPCnCut,UseNewEP,UseTender,period,passV0,passTPC,TimeCut,WeightSyst,SystTOFcut,CutM02,CutM20,SScut,EnablePileupRejVZEROTPCout);
        mgr->AddTask(taskMC);
        
        taskMC->SelectCollisionCandidates(AliVEvent::kMB);
        
        TString containerName4 = mgr->GetCommonFileName();
        containerName4 += ":PWGHF_hfeCalMCV2";
        containerName4 += ID;
        
        TString name4 = "histMC";
        name4 += ID;
        
        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name4.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName4.Data());
        mgr->ConnectInput(taskMC, 0, cinput);
        mgr->ConnectOutput(taskMC, 1, coutput1);
    }
    
    
    return NULL;
}

