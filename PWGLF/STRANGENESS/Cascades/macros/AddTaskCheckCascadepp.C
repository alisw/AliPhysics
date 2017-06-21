//////////////////////////////////////////////////////////////////
//    
//      Author: dcolella (domenico.colella@cern.ch)
//
//      Notes:
//       - collidingSystem: "pp", "pPb"
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

AliAnalysisTaskCheckCascadepp *AddTaskCheckCascadepp( Int_t    collidingSystem                     = 0,
                                                      AliVEvent::EOfflineTriggerTypes triggerclass = AliVEvent::kINT7, 
                                                      Int_t    minnTPCcls                          = 70,
                                                      Float_t  minTPCcrossrawoverfindable          = 0.8,
                                                      Float_t  vtxlimmax                           = 10.0,
                                                      Float_t  vtxlimmin                           = 0.0,
                                                      Bool_t   ksddselection                       = kFALSE,
                                                      Bool_t   kwithsdd                            = kFALSE,
                                                      Float_t  minptondaughtertracks               = 0.0,
                                                      Float_t  etacutondaughtertracks              = 0.8,
                                                      TString  suffix                              = "" ) {

   //______________________________________________________________________________
   // Creates, configures and attaches to the train a cascades check task
   // Get the pointer to the existing analysis manager via the static access method
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckCascadepp", "No analysis manager to connect to");
      return NULL;
   }   

   //___________________________________________________________________________________
   // Check the analysis type using the event handlers connected to the analysis manager
   //===================================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckCascadepp", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   //______________________________
   // Create and configure the task
   //==============================
   const Char_t *sddstatus = "";
   if      (collidingSystem == 0 && ksddselection && kwithsdd)   sddstatus = "_wSDDon";
   else if (collidingSystem == 0 && ksddselection && !kwithsdd)  sddstatus = "_wSDDoff";
   TString taskname = Form("TaskCheckCascadepp_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   taskname.Append(Form("%s",suffix.Data()));

   AliAnalysisTaskCheckCascadepp *taskcheckcascadepp = new AliAnalysisTaskCheckCascadepp(taskname);
     taskcheckcascadepp->SetAnalysisType                 (type);                   // "ESD" or "AOD"
     taskcheckcascadepp->SetCollidingSystem              (collidingSystem);        // choose the collision system to run on: "pp" and "pPb"
     taskcheckcascadepp->SetSelectedTriggerClass         (triggerclass);           // trigger selection
     taskcheckcascadepp->SetEventSelSDDstatus            (ksddselection);
     taskcheckcascadepp->SetWithSDDOn                    (kwithsdd);
     taskcheckcascadepp->SetQualityCutMinnTPCcls         (minnTPCcls);             // which value do you want apply for the minTPCcls cut?
     taskcheckcascadepp->SetQualityCutClusterOverFindable(minTPCcrossrawoverfindable); //
     taskcheckcascadepp->SetVertexRange                  (vtxlimmin,vtxlimmax);
     taskcheckcascadepp->SetMinptCutOnDaughterTracks     (minptondaughtertracks);  // which value do you want apply for cut on min pt daughter track?
     taskcheckcascadepp->SetEtaCutOnDaughterTracks       (etacutondaughtertracks); // which value do you want apply for cut on eta daughter track?
     taskcheckcascadepp->SetSuffix                       (suffix);

   mgr->AddTask(taskcheckcascadepp);

   //______________________________________________________________________________
   // Create ONLY the output containers for the data produced by the task
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   // Directory name
   TString outputFileName = Form("%s:PWGLFStrangeness.outputCheckCascade", AliAnalysisManager::GetCommonFileName());
   // Objects name
   TString outputname0 = Form("clistCasc_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputname1 = Form("cfcontPIDXiM_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputname2 = Form("cfcontPIDXiP_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputname3 = Form("cfcontPIDOmegaM_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputname4 = Form("cfcontPIDOmegaP_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputname5 = Form("cfcontCuts_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,minTPCcrossrawoverfindable,vtxlimmax,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   outputname0.Append(Form("%s",suffix.Data()));
   outputname1.Append(Form("%s",suffix.Data()));
   outputname2.Append(Form("%s",suffix.Data()));
   outputname3.Append(Form("%s",suffix.Data()));
   outputname4.Append(Form("%s",suffix.Data()));
   outputname5.Append(Form("%s",suffix.Data()));
   //Save objects into the train common file
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outputname0, TList::Class(),          AliAnalysisManager::kOutputContainer, outputFileName);
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(outputname1, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(outputname2, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(outputname3, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(outputname4, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   AliAnalysisDataContainer *coutput6 = mgr->CreateContainer(outputname5, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   mgr->ConnectInput( taskcheckcascadepp, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascadepp, 1, coutput1);
   mgr->ConnectOutput(taskcheckcascadepp, 2, coutput2);
   mgr->ConnectOutput(taskcheckcascadepp, 3, coutput3);
   mgr->ConnectOutput(taskcheckcascadepp, 4, coutput4);
   mgr->ConnectOutput(taskcheckcascadepp, 5, coutput5);
   mgr->ConnectOutput(taskcheckcascadepp, 6, coutput6);
   
   return taskcheckcascadepp;
}   
