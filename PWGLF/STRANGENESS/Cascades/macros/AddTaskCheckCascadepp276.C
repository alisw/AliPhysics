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

AliAnalysisTaskCheckCascadepp276 *AddTaskCheckCascadepp276( TString  collidingSystem        = "pp",
                                                            Int_t    minnTPCcls             = 70,
                                                            Float_t  vtxlim                 = 10.0,
                                                            Float_t  vtxlimmin              = 0.0,
                                                            Bool_t   ksddselection          = kFALSE,
                                                            Bool_t   kwithsdd               = kFALSE,
                                                            Bool_t   kextrasel              = kFALSE,
                                                            Bool_t   krelaunchvertexers     = kFALSE,
                                                            Float_t  minptondaughtertracks  = 0.0,
                                                            Float_t  etacutondaughtertracks = 0.8) {

   //______________________________________________________________________________
   // Creates, configures and attaches to the train a cascades check task
   // Get the pointer to the existing analysis manager via the static access method
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckCascadepp276", "No analysis manager to connect to");
      return NULL;
   }   


   //___________________________________________________________________________________
   // Check the analysis type using the event handlers connected to the analysis manager
   //===================================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckCascadepp276", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"


   //______________________________
   // Create and configure the task
   //==============================
   const Char_t *sddstatus = "";
   if      (collidingSystem == "pp" && ksddselection && kwithsdd)   sddstatus = "_wSDDon";
   else if (collidingSystem == "pp" && ksddselection && !kwithsdd)  sddstatus = "_wSDDoff";
   else if (collidingSystem == "pp" && !ksddselection && !kwithsdd) sddstatus = "_woSDD";
   TString taskname = Form("TaskCheckCascadepp276_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);

   AliAnalysisTaskCheckCascadepp276 *taskcheckcascadepp276 = new AliAnalysisTaskCheckCascadepp276(taskname);
     taskcheckcascadepp276->SetAnalysisType               (type);                   // "ESD" or "AOD"
     taskcheckcascadepp276->SetCollidingSystem            (collidingSystem);        // choose the collidiond system to run on: "pp" and "pPb"
     taskcheckcascadepp276->SetRelaunchV0CascVertexers    (krelaunchvertexers);     // choose if reconstruct the vertex of V0 in the cascades
     taskcheckcascadepp276->SetSDDSelection               (ksddselection);          // choose if apply SDD event selection
     taskcheckcascadepp276->SetWithSDDOn                  (kwithsdd);               // choose which SDD selection apply [if kTRUE select SDDon events]
     taskcheckcascadepp276->SetQualityCutZprimVtxPos      (kTRUE);                  // choose if apply Z vtx PV position event selection
     taskcheckcascadepp276->SetQualityCutNoTPConlyPrimVtx (kTRUE);                  // choose if apply no TPC only event selection
     taskcheckcascadepp276->SetQualityCutTPCrefit         (kTRUE);                  // choose if apply TPC refit on daughter tracks
     taskcheckcascadepp276->SetQualityCutnTPCcls          (kTRUE);                  // choose if apply n TPC cluster selection on daughter tracks
     taskcheckcascadepp276->SetQualityCutPileup           (kTRUE);                  // choose if apply no Pileup event selection
     taskcheckcascadepp276->SetQualityCutMinnTPCcls       (minnTPCcls);             // which value do you want apply for the minTPCcls cut?
     taskcheckcascadepp276->SetExtraSelections            (kextrasel);              // choose if apply the extra selection of cascade reco.
     taskcheckcascadepp276->SetVertexRange                (vtxlim);                 // which higher value do you want apply for vtx Z cut?
     taskcheckcascadepp276->SetVertexRangeMin             (vtxlimmin);              // which lower value do you want apply for vtx Z cut?
     taskcheckcascadepp276->SetMinptCutOnDaughterTracks   (minptondaughtertracks);  // which value do you want apply for cut on min pt daughter track?
     taskcheckcascadepp276->SetEtaCutOnDaughterTracks     (etacutondaughtertracks); // which value do you want apply for cut on eta daughter track?

   mgr->AddTask(taskcheckcascadepp276);


   //______________________________________________________________________________
   // Create ONLY the output containers for the data produced by the task
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   // Directory name
   TString outputFileName = Form("%s:PWGLFStrangeness.outputCheckCascade", AliAnalysisManager::GetCommonFileName());
   // Objects name
   TString outputname0 = Form("clistCasc_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputname1 = Form("cfcontPIDXiM_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputname2 = Form("cfcontPIDXiP_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputname3 = Form("cfcontPIDOmegaM_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputname4 = Form("cfcontPIDOmegaP_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   TString outputname5 = Form("cfcontCuts_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f%s",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks,sddstatus);
   //Save objects into the train common file
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outputname0, TList::Class(),          AliAnalysisManager::kOutputContainer, outputFileName);
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(outputname1, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(outputname2, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(outputname3, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(outputname4, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   AliAnalysisDataContainer *coutput6 = mgr->CreateContainer(outputname5, AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   mgr->ConnectInput( taskcheckcascadepp276, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascadepp276, 1, coutput1);
   mgr->ConnectOutput(taskcheckcascadepp276, 2, coutput2);
   mgr->ConnectOutput(taskcheckcascadepp276, 3, coutput3);
   mgr->ConnectOutput(taskcheckcascadepp276, 4, coutput4);
   mgr->ConnectOutput(taskcheckcascadepp276, 5, coutput5);
   mgr->ConnectOutput(taskcheckcascadepp276, 6, coutput6);
   
   return taskcheckcascadepp276;
}   
