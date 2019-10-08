AliAnalysisTaskStrangeCascadesDiscrete *AddTaskStrangeCascadesDiscrete(Bool_t lSaveEventTree = kTRUE,
                                                                                 Bool_t lSaveV0Tree = kFALSE,
                                                                                 Bool_t lSaveCascadeTree = kTRUE,
                                                                                 Bool_t lRunVertexers = kFALSE,
                                                                                 Bool_t lUseLightVertexer = kFALSE,
                                                                                 Bool_t lMomentaNeeded = kFALSE, //9 branches of p,pz,pt of daughters
                                                                                 Bool_t lExtraInfoCascadeKinematics = kTRUE,
                                                                                 Bool_t lDiscreteSymmetryInfo = kFALSE, //extra branch in cascade tree
                                                                                 Bool_t lMassesAndBackgroundFromScratch = kFALSE, //extra branch in cascade tree
                                                                                 Double_t lCascaderMaxChi2 = 33.,
                                                                                 Double_t lCascaderV0MinImpactParam = 0.05, //0.050 per def
                                                                                 Double_t lCascaderV0MassWindow = 0.01, //0.010 per def
                                                                                 Double_t lCascaderBachMinImpactParam = 0.03, //0.03 per def
                                                                                 Double_t lCascaderMaxDCAV0andBach = 2.0, //2.0 per def
                                                                                 Double_t lCascaderMinCosAngle = 0.95, //0.95 per def
                                                                                 Double_t lCascaderMinRadius = 0.4, //0.4 per def
                                                                                 Double_t lCascaderMaxRadius =100., //100. per def
                                                                                 TString lExtraOptions = "",
                                                                                 const TString lMasterJobSessionFlag = "",
                                                                                 TString lExtraOutputName = ""
                                                                                 )
{
    
    
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskStrangeCascadesDiscrete", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskStrangeCascadesDiscrete", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliAnalysisTaskStrangeCascadesDiscrete *taskAuxiliary =
    new AliAnalysisTaskStrangeCascadesDiscrete(lSaveEventTree,
                                                    lSaveV0Tree,
                                                    lSaveCascadeTree,
                                                    lRunVertexers,
                                                    lUseLightVertexer,
                                                    lMomentaNeeded,
                                                    lExtraInfoCascadeKinematics,
                                                    lDiscreteSymmetryInfo,
                                                    lMassesAndBackgroundFromScratch, //extra branch in cascade tree
                                                    lCascaderMaxChi2,
                                                    lCascaderV0MinImpactParam,
                                                    lCascaderV0MassWindow,
                                                    lCascaderBachMinImpactParam,
                                                    lCascaderMaxDCAV0andBach,
                                                    lCascaderMinCosAngle,
                                                    lCascaderMinRadius,
                                                    lCascaderMaxRadius,
                                                    "taskAuxiliary",
                                                    lExtraOptions
                                                    );  // const*charname =  Form("taskAuxiliary%s",lExtraOutputName.Data())
    
    mgr->AddTask(taskAuxiliary);
    
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    //   TString results = "AnalysisResults.root";
    if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    
    
    AliAnalysisDataContainer* coutput = mgr->CreateContainer("skorodum_cascadetreeAsevent", TTree::Class(),
                                                             AliAnalysisManager::kOutputContainer,outputFileName);
    
    mgr->ConnectInput(taskAuxiliary, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskAuxiliary, 1, coutput);
    return taskAuxiliary;
    
}


