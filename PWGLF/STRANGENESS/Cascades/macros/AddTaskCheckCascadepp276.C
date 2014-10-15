AliAnalysisTaskCheckCascadepp276 *AddTaskCheckCascadepp276( TString  collidingSystem        = "pp",
                                                            Int_t    minnTPCcls             = 70,
                                                            Float_t  vtxlim                 = 10.0,
                                                            Float_t  vtxlimmin              = 0.0,
                                                            Bool_t   fwithsdd               = kFALSE,
                                                            Bool_t   kextrasel              = kFALSE,
                                                            Bool_t   krelaunchvertexers     = kFALSE,
                                                            Bool_t   ksddonselection        = kTRUE,
                                                            Float_t  minptondaughtertracks  = 0.0,
                                                            Float_t  etacutondaughtertracks = 0.8) {

   // Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckCascadepp276", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckCascadepp276", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
   //==============================================================================
   TString taskname = Form("TaskCheckCascadepp276_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
     if(fwithsdd && collidingSystem == "pp") {
          if (ksddonselection)       taskname += "_wSDDon";
          else if (!ksddonselection) taskname += "_wSDDoff";
     } else if (!fwithsdd && collidingSystem == "pp") {
          taskname += "_woSDD";
     }
   AliAnalysisTaskCheckCascadepp276 *taskcheckcascadepp276 = new AliAnalysisTaskCheckCascadepp276(taskname);
     taskcheckcascadepp276->SetAnalysisType               (type);                   // "ESD" or "AOD"
     taskcheckcascadepp276->SetCollidingSystem            (collidingSystem);        // choose the collidiond system to run on: "pp" and "pPb"
     taskcheckcascadepp276->SetRelaunchV0CascVertexers    (krelaunchvertexers);     // choose if reconstruct the vertex of V0 in the cascades
     taskcheckcascadepp276->SetSDDSelection               (fwithsdd);               // choose if apply SDD event selection
     taskcheckcascadepp276->SetQualityCutZprimVtxPos      (kTRUE);                  // choose if apply Z vtx PV position event selection
     taskcheckcascadepp276->SetQualityCutNoTPConlyPrimVtx (kTRUE);                  // choose if apply no TPC only event selection
     taskcheckcascadepp276->SetQualityCutTPCrefit         (kTRUE);                  // choose if apply TPC refit on daughter tracks
     taskcheckcascadepp276->SetQualityCutnTPCcls          (kTRUE);                  // choose if apply n TPC cluster selection on daughter tracks
     taskcheckcascadepp276->SetQualityCutPileup           (kTRUE);                  // choose if apply no Pileup event selection
     taskcheckcascadepp276->SetWithSDDOn                  (ksddonselection);        // which SDD selection do you want apply? [if kTRUE select SDDon events]
     taskcheckcascadepp276->SetQualityCutMinnTPCcls       (minnTPCcls);             // which value do you want apply for the minTPCcls cut?
     taskcheckcascadepp276->SetExtraSelections            (kextrasel);              // choose if apply the extra selection of cascade reco.
     taskcheckcascadepp276->SetVertexRange                (vtxlim);                 // which higher value do you want apply for vtx Z cut?
     taskcheckcascadepp276->SetVertexRangeMin             (vtxlimmin);              // which lower value do you want apply for vtx Z cut?
     taskcheckcascadepp276->SetMinptCutOnDaughterTracks   (minptondaughtertracks);  // which value do you want apply for cut on min pt daughter track?
     taskcheckcascadepp276->SetEtaCutOnDaughterTracks     (etacutondaughtertracks); // which value do you want apply for cut on eta daughter track?

   mgr->AddTask(taskcheckcascadepp276);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   // Directory name
   TString outputFileName = Form("%s:PWGLFStrangeness.outputCheckCascadepp276", AliAnalysisManager::GetCommonFileName());
   // Objects name
   TString outputname0 = "clistCasc";
   TString outputname1 = "cfcontPIDAsXiM";
   TString outputname2 = "cfcontPIDAsXiP";
   TString outputname3 = "cfcontPIDAsOmegaM";
   TString outputname4 = "cfcontPIDAsOmegaP";
   TString outputname5 = "cfcontAsCuts";
      outputname0 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
      outputname1 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
      outputname2 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
      outputname3 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
      outputname4 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
      outputname5 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
      if(fwithsdd && collidingSystem == "pp") {
          if (ksddonselection) {
              outputname0 += "_wSDDon";
              outputname1 += "_wSDDon";
              outputname2 += "_wSDDon";
              outputname3 += "_wSDDon";
              outputname4 += "_wSDDon";
              outputname5 += "_wSDDon";
          } else if (!ksddonselection && collidingSystem == "pp") {
              outputname0 += "_wSDDoff";
              outputname1 += "_wSDDoff";
              outputname2 += "_wSDDoff";
              outputname3 += "_wSDDoff";
              outputname4 += "_wSDDoff";
              outputname5 += "_wSDDoff";
          }
      } else if (!fwithsdd && collidingSystem == "pp") {
          outputname0 += "_woSDD";
          outputname1 += "_woSDD";
          outputname2 += "_woSDD";
          outputname3 += "_woSDD";
          outputname4 += "_woSDD";
          outputname5 += "_woSDD";
      }

   //Save objects into the train common file
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outputname0,
			                                     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(outputname1,
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(outputname2,
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(outputname3,
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(outputname4,
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   AliAnalysisDataContainer *coutput6 = mgr->CreateContainer(outputname5,
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   mgr->ConnectInput( taskcheckcascadepp276, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascadepp276, 1, coutput1);
   mgr->ConnectOutput(taskcheckcascadepp276, 2, coutput2);
   mgr->ConnectOutput(taskcheckcascadepp276, 3, coutput3);
   mgr->ConnectOutput(taskcheckcascadepp276, 4, coutput4);
   mgr->ConnectOutput(taskcheckcascadepp276, 5, coutput5);
   mgr->ConnectOutput(taskcheckcascadepp276, 6, coutput6);
   
   return taskcheckcascadepp276;
}   
