

///////////////////////////////////////////////////////////////////
//                                                               //
// AddTaskFlowTPCEMCalQCSP2 macro                                 //
// Author: Andrea Dubla, Utrecht University, 2012                //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowTrackSimpleCuts;
class AliFlowEventCuts;
class AliFlowEventSimpleCuts;
class AliAnalysisDataContainer;
class AliHFEextraCuts;

AliAnalysisTaskHFEEfficiency*  AddTaskHFEEfficiency(
                                                    TString uniqueID = "",
                                                    Float_t centrMin,
                                                    Float_t centrMax,
                                                    Double_t etamin = -0.8,
                                                    Double_t etamax = 0.8,
                                                    Bool_t WeightsHF = kTRUE,
                                                    Bool_t Weights = kTRUE,
                                                    Bool_t CentralWeights = kTRUE,
                                                    Bool_t SemicentralWeights = kFALSE,
                                                    Bool_t UpWeights = kFALSE,
                                                    Bool_t DwWeights = kFALSE,
                                                    Bool_t SetSTACK = kFALSE,
                                                    Int_t minTPCCluster,
                                                    AliHFEextraCuts::ITSPixel_t pixel,
                                                    Double_t minTOFnSigma,
                                                    Double_t maxTOFnSigma,
                                                    Double_t minITSnsigmaLowpT,
                                                    Double_t maxITSnsigmaLowpT,
                                                    Double_t minITSnsigmaHighpT,
                                                    Double_t maxITSnsigmaHighpT,
                                                    Double_t minTPCnsigmaLowpT,
                                                    Double_t maxTPCnsigmaLowpT,
                                                    Double_t minTPCnsigmaHighpT,
                                                    Double_t maxTPCnsigmaHighpT,
                                                    Int_t TPCS,
                                                    Double_t MassCut = 0.07,
                                                    Int_t TPCClusterforAsso = 80,
                                                    Bool_t AssoITSref = kTRUE,
                                                    Double_t pTminAssoTrack = 0.0,
                                                    //const char *Cent = "V0M",
                                                    Bool_t shrinkSP = kTRUE,
                                                    Bool_t debug = kFALSE
                                                    )

{
    
    
    
    if(debug) cout << " === Adding Task ElectronEfficiency === " << endl;
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":ElectroID_";
    fileName += uniqueID;
    if(debug) cout << "    --> Reconstruction data container: " << fileName << endl;
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        if(debug) cout << " Fatal error: no analysis manager found! " << endl;
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        if(debug) cout << " Fatal error: no imput event handler found!" << endl;
        return 0x0;
    }
    
    //create a task
    AliAnalysisTaskHFEEfficiency *taskHFE = ConfigHFEEff(kTRUE, minTPCCluster, pixel);    //kTRUE if MC
    
    if(debug) cout << " === AliAnalysisTaskHFEEfficiency === " << taskHFE << endl;
    if(!taskHFE) {
        if(debug) cout << " --> Unexpected error occurred: NO TASK WAS CREATED! (could be a library problem!) " << endl;
        return 0x0;
    }
    
    
    // Set centrality percentiles and method V0M, FMD, TRK, TKL, CL0, CL1, V0MvsFMD, TKLvsV0M, ZEMvsZDC
    taskHFE->SetCentralityParameters(centrMin, centrMax, "V0M");
    taskHFE->SetTPCS(TPCS);
    taskHFE->SetIDCuts(minTOFnSigma, maxTOFnSigma, minITSnsigmaLowpT, maxITSnsigmaLowpT, minITSnsigmaHighpT, maxITSnsigmaHighpT, minTPCnsigmaLowpT, maxTPCnsigmaLowpT, minTPCnsigmaHighpT, maxTPCnsigmaHighpT);
    taskHFE->SetAssoTPCCluster(TPCClusterforAsso);
    taskHFE->SetAssopTmin(pTminAssoTrack);
    taskHFE->SetAssoITSRefit(AssoITSref);
    taskHFE->SetStackLoop(SetSTACK);
    taskHFE->SetMassCut(MassCut);
    taskHFE->SetWeights(Weights);
    taskHFE->SetWeightsHF(WeightsHF);
    taskHFE->SetCentralWeights(CentralWeights);
    taskHFE->SetSemicentralWeights(SemicentralWeights);
    
    taskHFE->SetTiltUpWeights(UpWeights);
    taskHFE->SetTiltDwWeights(DwWeights);
    taskHFE->SetEtaRange(etamin,etamax);//-0.8,0.8
    
    //set RP cuts for flow package analysis
    // TString foutputName = "PbPbPhotonicElecEfficiency010ITSTOFTPCWeights.root";
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("ccontainer0_%s",uniqueID.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
    
    mgr->ConnectInput(taskHFE,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHFE,1,coutput3);
    mgr->AddTask(taskHFE);
    
    
    return taskHFE;
    
}

//_____________________________________________________________________________



//_____________________________________________________________________________

AliAnalysisTaskHFEEfficiency* ConfigHFEEff(Bool_t useMC,Int_t minTPCCulster,AliHFEextraCuts::ITSPixel_t pixel){
    //
    // HFE standard task configuration
    //
    
    Bool_t kAnalyseTaggedTracks = kTRUE;
    
    AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsEffITSTOFTPC","HFE Standard Cuts");  //TODO....change the cuts values to PbPb
    //  hfecuts->CreateStandardCuts();
    hfecuts->SetMinNClustersTPC(minTPCCulster);
    hfecuts->SetMinNClustersITS(5);
    hfecuts->SetMinNTrackletsTRD(0);
    hfecuts->SetMinRatioTPCclusters(0.6);
    
    //   hfecuts->SetEtaRange(-0.9,0.9);
    //   hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    hfecuts->SetRequireITSPixel();
    hfecuts->SetCutITSpixel(pixel);//
    hfecuts->SetMaxChi2perClusterITS(-1);
    hfecuts->SetMaxChi2perClusterTPC(3.5);
    hfecuts->SetCheckITSLayerStatus(kFALSE); // shud be put back
    //  hfecuts->UnsetVertexRequirement();
    hfecuts->SetVertexRange(10.);
    hfecuts->SetRequireSigmaToVertex();
    //hfecuts->SetSigmaToVertex(10);
    hfecuts->SetTOFPIDStep(kFALSE);
    //  hfecuts->SetQAOn();
    hfecuts->SetPtRange(0, 30);
    
    AliAnalysisTaskHFEEfficiency *task = new AliAnalysisTaskHFEEfficiency("HFE_Eff");
    printf("task ------------------------ %p\n ", task);
    
    
    task->SetHFECuts(hfecuts);
    
    //   task->SetInvariantMassCut(0.05);
    task->SetRejectKinkMother(kTRUE);
    //  task->SetRemovePileUp(kTRUE);
    
    // Define PID
    AliHFEpid *pid = task->GetPID();
    if(useMC) pid->SetHasMCData(kTRUE);
    pid->AddDetector("ITS", 0);
    pid->AddDetector("TOF", 1);
    pid->AddDetector("TPC", 2);
    
    printf("*************************************\n");
    printf("Configuring standard Task:\n");
    //  task->PrintStatus();
    pid->PrintStatus();
    printf("*************************************\n");
    return task;
    
    
}

//_____________________________________________________________________________

