AliAnalysisTaskStrangeCascadesDiscrete *AddTaskStrangeCascadesDiscrete(
                                                                       Bool_t lRunV0Vertexers = kTRUE,
                                                                       Bool_t lRunVertexers = kTRUE, //to rerun the cascade vertexers
                                                                       Bool_t lUseLightVertexer = kFALSE, //use light cascade vertexers
                                                                       Bool_t lUseOnTheFlyV0Cascading = kFALSE,
                                                                       Bool_t lguard_CheckTrackQuality = kTRUE,
                                                                       Bool_t lguard_CheckCascadeQuality = kTRUE,
                                                                       Bool_t lguard_CheckTPCPID = kTRUE,
                                                                       Double_t lV0MaxChi2 = 33.,
                                                                       Double_t lV0minDCAfirst = 0.01,
                                                                       Double_t lV0minDCAsecond = 0.01,
                                                                       Double_t lV0maxDCAdaughters = 2.,
                                                                       Double_t lV0minCosAngle = 0.6,
                                                                       Double_t lV0minRadius = 0.4,
                                                                       Double_t lV0maxRadius = 200.,
                                                                       
                                                                       Double_t lCascaderMaxChi2 = 33.,
                                                                       Double_t lCascaderV0MinImpactParam = 0.05, //0.050 per def
                                                                       Double_t lCascaderV0MassWindow = 0.008, //0.010 per def
                                                                       Double_t lCascaderBachMinImpactParam = 0.03, //0.03 per def
                                                                       Double_t lCascaderMaxDCAV0andBach = 2.0, //2.0 per def
                                                                       Double_t lCascaderMinCosAngle = 0.95, //0.95 per def
                                                                       Double_t lCascaderMinRadius = 0.4, //0.4 per def
                                                                       Double_t lCascaderMaxRadius =100., //100. per def
                                                                       Float_t sigmaRangeTPC = 3.,
                                                                       Float_t lOmegaCleanMassWindow = 0.1,
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
    new AliAnalysisTaskStrangeCascadesDiscrete(lRunV0Vertexers,
                                               lRunVertexers, lUseLightVertexer,
                                               lUseOnTheFlyV0Cascading,
                                               lguard_CheckTrackQuality,
                                               lguard_CheckCascadeQuality,
                                               lguard_CheckTPCPID,
                                               lV0MaxChi2,
                                               lV0minDCAfirst,
                                               lV0minDCAsecond,
                                               lV0maxDCAdaughters,
                                               lV0minCosAngle,
                                               lV0minRadius,
                                               lV0maxRadius,
                                               lCascaderMaxChi2,
                                               lCascaderV0MinImpactParam,
                                               lCascaderV0MassWindow,
                                               lCascaderBachMinImpactParam,
                                               lCascaderMaxDCAV0andBach,
                                               lCascaderMinCosAngle,
                                               lCascaderMinRadius,
                                               lCascaderMaxRadius, sigmaRangeTPC,lOmegaCleanMassWindow, "taskAuxiliary");  // const*charname =  Form("taskAuxiliary%s",lExtraOutputName.Data())
    
    mgr->AddTask(taskAuxiliary);
    
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    //   TString results = "AnalysisResults.root";
    if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    
    AliAnalysisDataContainer* coutput = mgr->CreateContainer("skorodum_cascadetreeAsevent_scratch", TTree::Class(),
                                                             AliAnalysisManager::kOutputContainer, outputFileName);
    
    mgr->ConnectInput(taskAuxiliary, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskAuxiliary, 1, coutput);
    
    return taskAuxiliary;
    
}


