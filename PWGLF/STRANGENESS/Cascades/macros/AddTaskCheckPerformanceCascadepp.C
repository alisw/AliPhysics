//////////////////////////////////////////////////////////////////
//
//      Author: dcolella (domenico.colella@cern.ch)
//
//      Notes:
//       - collidingSystem: 0 --> "pp", 1 --> "pPb"
//       - SDD selection
//         This selection has been introduced for the pp@2.76TeV analysis
//         and allow to select the events wrt the status of the SDD during
//         the data acquisition. There are two different data re-productions:
//         one that use the infos from the SDD and another that do not use.
//         We have two parameters to select:
//           1)  ksddselection = kTRUE and kwithsdd = kTRUE   --> Sample withSDD_SDDon
//           2)  ksddselection = kTRUE and kwithsdd = kFALSE  --> Sample withSDD_SDDoff
//           3)  ksddselection = kFALSE and kwithsdd = kFALSE --> Sample withoutSDD
//           4)  ksddselection = kFALSE and kwithsdd = kTRUE  --> No selection applied
//
///////////////////////////////////////////////////////////////

AliAnalysisTaskCheckPerformanceCascadepp *AddTaskCheckPerformanceCascadepp( Int_t  collidingSystem                       = 0,
                                                                            AliVEvent::EOfflineTriggerTypes triggerclass = AliVEvent::kINT7,
                                                                            Int_t    minnTPCcls                          = 70,
                                                                            Float_t  minTPCcrossrawoverfindable          = 0.8,
                                                                            Float_t  vtxlimmax                           = 10.,
                                                                            Float_t  vtxlimmin                           = 0.,
                                                                            Bool_t   ksddselection                       = kFALSE,                                                                                  
                                                                            Bool_t   kwithsdd                            = kFALSE,
                                                                            Float_t  minptondaughtertracks               = 0.0,
                                                                            Float_t  etacutondaughtertracks              = 0.8,
                                                                            Bool_t   kacccut                             = kFALSE,
                                                                            TString  suffix                              = "" ) {
    
   //______________________________________________________________________________
   // Creates, configures and attaches to the train a cascades check task
   // Get the pointer to the existing analysis manager via the static access method
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckPerformanceCascade", "No analysis manager to connect to.");
      return NULL;
   }   


   //___________________________________________________________________________________
   // Check the analysis type using the event handlers connected to the analysis manager
   //===================================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckPerformanceCascade", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"


   //______________________________
   // Create and configure the task
   //==============================
   const Char_t *sddstatus = "";
   if      (collidingSystem == 0 && ksddselection && kwithsdd)   sddstatus = "_wSDDon";
   else if (collidingSystem == 0 && ksddselection && !kwithsdd)  sddstatus = "_wSDDoff";
   TString tasknameperf = Form("TaskCheckPerformanceCascade_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minptondaughtertracks,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   tasknameperf.Append(Form("%s",suffix.Data()));
   AliAnalysisTaskCheckPerformanceCascadepp *taskCheckPerfCascadepp = new AliAnalysisTaskCheckPerformanceCascadepp(tasknameperf);
     taskCheckPerfCascadepp->SetAnalysisType                 (type);                   // "ESD" or "AOD"
     taskCheckPerfCascadepp->SetCollidingSystem              (collidingSystem);        // choose the collision system to run on: "pp" and "pPb"
     taskCheckPerfCascadepp->SetSelectedTriggerClass         (triggerclass);           // trigger selection
     taskCheckPerfCascadepp->SetEventSelSDDstatus            (ksddselection);
     taskCheckPerfCascadepp->SetWithSDDOn                    (kwithsdd);
     taskCheckPerfCascadepp->SetQualityCutMinnTPCcls         (minnTPCcls);             // which value do you want apply for the minTPCcls cut?
     taskCheckPerfCascadepp->SetQualityCutClusterOverFindable(minTPCcrossrawoverfindable);
     taskCheckPerfCascadepp->SetVertexRange                  (vtxlimmin,vtxlimmax);
     taskCheckPerfCascadepp->SetMinptCutOnDaughterTracks     (minptondaughtertracks);  // which value do you want apply for cut on min pt daughter track?
     taskCheckPerfCascadepp->SetEtaCutOnDaughterTracks       (etacutondaughtertracks); // which value do you want apply for cut on eta daughter track?
     taskCheckPerfCascadepp->SetApplyAccCut                  (kacccut);                // choose if apply acceptance cut
     taskCheckPerfCascadepp->SetSuffix                       (suffix);

   mgr->AddTask(taskCheckPerfCascadepp);


   //______________________________________________________________________________
   // Create ONLY the output containers for the data produced by the task
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   // - Directory name
   TString outputFileNamePerf = Form("%s:PWGLFStrangeness.outputCheckPerformanceCascade", AliAnalysisManager::GetCommonFileName());
   // - Objects name
   TString outputnameperf0 = Form("clistCascPerf_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputnameperf1 = Form("cfcontPIDAsXiM_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputnameperf2 = Form("cfcontPIDAsXiP_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputnameperf3 = Form("cfcontPIDAsOmegaM_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputnameperf4 = Form("cfcontPIDAsOmegaP_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputnameperf5 = Form("cfcontAsCuts_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   outputnameperf0.Append(Form("%s",suffix.Data()));
   outputnameperf1.Append(Form("%s",suffix.Data()));
   outputnameperf2.Append(Form("%s",suffix.Data()));
   outputnameperf3.Append(Form("%s",suffix.Data()));
   outputnameperf4.Append(Form("%s",suffix.Data()));
   outputnameperf5.Append(Form("%s",suffix.Data()));
   // - Save objects into the train common file
   AliAnalysisDataContainer *coutputperf1 = mgr->CreateContainer(outputnameperf0, TList::Class(),          AliAnalysisManager::kOutputContainer, outputFileNamePerf);
   AliAnalysisDataContainer *coutputperf2 = mgr->CreateContainer(outputnameperf1, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileNamePerf);
   AliAnalysisDataContainer *coutputperf3 = mgr->CreateContainer(outputnameperf2, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileNamePerf);
   AliAnalysisDataContainer *coutputperf4 = mgr->CreateContainer(outputnameperf3, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileNamePerf);
   AliAnalysisDataContainer *coutputperf5 = mgr->CreateContainer(outputnameperf4, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileNamePerf);
   AliAnalysisDataContainer *coutputperf6 = mgr->CreateContainer(outputnameperf5, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileNamePerf);
   // - Connect output/input
   mgr->ConnectInput( taskCheckPerfCascadepp, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskCheckPerfCascadepp, 1, coutputperf1);
   mgr->ConnectOutput(taskCheckPerfCascadepp, 2, coutputperf2);
   mgr->ConnectOutput(taskCheckPerfCascadepp, 3, coutputperf3);
   mgr->ConnectOutput(taskCheckPerfCascadepp, 4, coutputperf4);
   mgr->ConnectOutput(taskCheckPerfCascadepp, 5, coutputperf5);
   mgr->ConnectOutput(taskCheckPerfCascadepp, 6, coutputperf6);
   
   return taskCheckPerfCascadepp;
}   
