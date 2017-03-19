///////////////////////////////////////////////////////////////////
//                                                               //
// AddTaskFlowITSTPCTOFQCSP macro                                     //
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

AliAnalysisTaskFlowITSTPCTOFQCSP* AddTaskFlowITSTPCTOFQCSP(
                                                           TString uniqueID = "",
                                                           Float_t centrMin ,
                                                           Float_t centrMax ,
                                                           Bool_t CentFlatMine,
                                                           Double_t etamin = -0.8,
                                                           Double_t etamax = 0.8,
                                                           Double_t InvmassCut,
                                                           Int_t Trigger,
                                                           Bool_t multCorrcut,
                                                           Double_t pTCutmin,
                                                           Double_t pTCutmax,
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
                                                           Int_t minTPCCluster,
                                                           Int_t TPCS,
                                                           AliHFEextraCuts::ITSPixel_t pixel,
                                                           Int_t TPCClusterforAsso = 80,
                                                         //  Bool_t AssoITSref = kTRUE,
                                                           Double_t ptminassocut = 0.0,
                                                           Bool_t Weight = kFALSE,
                                                           Bool_t withmultetacorrection=kFALSE,
                                                           Double_t etaminpos = 0,
                                                           Double_t etaminneg = 0,
                                                     //      Bool_t PhiCut = kFALSE,
                                                           Bool_t PhotonicElectronDCA = kFALSE,
                                                           // Bool_t QaPidSparse = kFALSE,
                                                           const char *Cent = "V0M",
                                                           Bool_t QC = kTRUE, // use qc2 and qc4
                                                           Bool_t SP_TPC = kTRUE, //use tpc sp method
                                                           Bool_t VZERO_SP = kFALSE, // use vzero sp method
                                                           Int_t harmonic = 2,
                                                           Bool_t shrinkSP = kTRUE,
                                                           Bool_t debug = kFALSE,
                                                           Int_t RPFilterBit = 1
//Bool_t op_ang = kFALSE
//Int_t Vz = 10,
//Double_t op_angle_cut = 3.
)

{
    
    
    
    if(debug) cout << " === Adding Task ElectFlow === " << endl;
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
    AliAnalysisTaskFlowITSTPCTOFQCSP *taskHFE = ConfigHFEStandardCuts(kFALSE, minTPCCluster, pixel, withmultetacorrection);    //kTRUE if MC
    
    if(debug) cout << " === AliAnalysisTaskFlowITSTPCTOFQCSP === " << taskHFE << endl;
    if(!taskHFE) {
        if(debug) cout << " --> Unexpected error occurred: NO TASK WAS CREATED! (could be a library problem!) " << endl;
        return 0x0;
    }
    taskHFE->SetTrigger(Trigger);
    taskHFE->SetEPWeight(Weight);
    
    
    TString histoflatname = "alien:///alice/cern.ch/user/a/adubla/CentrDistrBins005.root";
    TString histoflatnameMine = "alien:///alice/cern.ch/user/a/adubla/CentFlat010_Smart.root";
    
    if(!CentFlatMine){
        if(Trigger==0 || Trigger==4){
            TFile *fFlat=TFile::Open(histoflatname.Data());
            TCanvas *c=fFlat->Get("cintegral");
            TH1F *hfl=(TH1F*)c->FindObject("hint");
            taskHFE->SetHistoForCentralityFlattening(hfl,centrMin,centrMax,0.,0);
        }
    }
    if(CentFlatMine){
        if(Trigger==0 || Trigger==4){
            TFile *fFlat=TFile::Open(histoflatnameMine.Data());
            TCanvas *c=fFlat->Get("CentFlat_mine");
            TH1F *hfl=(TH1F*)c->FindObject("cent");
            taskHFE->SetHistoForCentralityFlattening_Bis(hfl,centrMin,centrMax,0);
        }
    }
    taskHFE->SetCentralityMine(CentFlatMine);
    
    
    TString histoflatnameEP;
    if(centrMax == 10.) histoflatnameEP = "alien:///alice/cern.ch/user/a/adubla/EPVZero010_Smart.root";
    if(centrMax == 5.)  histoflatnameEP = "alien:///alice/cern.ch/user/a/adubla/EPVZero05_Smart.root";
    if(centrMin == 5.)  histoflatnameEP = "alien:///alice/cern.ch/user/a/adubla/EPVZero510_Smart.root";
    
    
    if(Weight){
        TFile *fFlatEP=TFile::Open(histoflatnameEP,"READ");
        TCanvas *cEP=fFlatEP->Get("c1_n7");
        TH1D *hEPfl=(TH1D*)cEP->FindObject("EPVz");
        taskHFE->SetHistoForEPFlattWeights(hEPfl);
    }
    
    
    // Set centrality percentiles and method V0M, FMD, TRK, TKL, CL0, CL1, V0MvsFMD, TKLvsV0M, ZEMvsZDC
    taskHFE->SetCentralityParameters(centrMin, centrMax, Cent);
    taskHFE->SetInvariantMassCut(InvmassCut);
    taskHFE->SetpTCuttrack(pTCutmin, pTCutmax);
    taskHFE->SetTPCS(TPCS);
    taskHFE->SetVz(10);
    taskHFE->SetIDCuts(minTOFnSigma, maxTOFnSigma, minITSnsigmaLowpT, maxITSnsigmaLowpT, minITSnsigmaHighpT, maxITSnsigmaHighpT, minTPCnsigmaLowpT, maxTPCnsigmaLowpT, minTPCnsigmaHighpT, maxTPCnsigmaHighpT);
    //  taskHFE->SetQAPIDSparse(QaPidSparse);
    taskHFE->SelectPhotonicElectronMethod(PhotonicElectronDCA);
    taskHFE->SetOpeningAngleflag(kFALSE);
    taskHFE->SetOpeningAngleCut(3);
    taskHFE->SetMultCorrelationCut(multCorrcut);
    taskHFE->SetPtMinAssoCut(ptminassocut);
    taskHFE->SetAssoTPCCluster(TPCClusterforAsso);
    taskHFE->SetAssoITSRefit(kTRUE);
    taskHFE->SetPhiCut(kFALSE);
    taskHFE->SetEtaMinPos(etaminpos); //0.2
    taskHFE->SetEtaMinNeg(etaminneg);//-0.2
    taskHFE->SetEtaRange(etamin,etamax);//-0.8,0.8
    
    
    
    //set RP cuts for flow package analysis
    cutsRP = new AliFlowTrackCuts(Form("RFPcuts%s",uniqueID));
    if(!cutsRP) {
        if(debug) cout << " Fatal error: no RP cuts found, could be a library problem! " << endl;
        return 0x0;
    }
    
    if(!VZERO_SP) {
        AliFlowTrackCuts::trackParameterType rptype = AliFlowTrackCuts::kGlobal;
        cutsRP->SetParamType(rptype);
        cutsRP->SetAODfilterBit(RPFilterBit);
        cutsRP->SetPtRange(0.2, 5.0);
        cutsRP->SetEtaRange(-0.8, 0.8);
        cutsRP->SetMinNClustersTPC(70);
        cutsRP->SetMinChi2PerClusterTPC(0.1);
        cutsRP->SetMaxChi2PerClusterTPC(4.0);
        cutsRP->SetRequireTPCRefit(kTRUE);
        cutsRP->SetMaxDCAToVertexXY(0.3);
        cutsRP->SetMaxDCAToVertexZ(0.3);
        cutsRP->SetAcceptKinkDaughters(kFALSE);
        cutsRP->SetMinimalTPCdedx(10.);
        if(debug) cout << "    --> kGlobal RP's " << cutsRP << endl;
    }
    if(VZERO_SP) { // use vzero sub analysis
        cutsRP = cutsRP->GetStandardVZEROOnlyTrackCuts(); // select vzero tracks
        SP_TPC = kFALSE; // disable other methods
        QC = kFALSE;
        if(debug) cout << "    --> VZERO RP's " << cutsRP << endl;
    }
    
    AliFlowTrackSimpleCuts *POIfilterLeft = new AliFlowTrackSimpleCuts();
    AliFlowTrackSimpleCuts *POIfilterRight = new AliFlowTrackSimpleCuts();
    if(SP_TPC){
        POIfilterLeft->SetEtaMin(-0.8);
        POIfilterLeft->SetEtaMax(0.0);
        POIfilterLeft->SetMassMin(263731); POIfilterLeft->SetMassMax(263733);
        
        POIfilterRight->SetEtaMin(0.0);
        POIfilterRight->SetEtaMax(0.8);
        POIfilterRight->SetMassMin(263731); POIfilterRight->SetMassMax(263733);
    }
    
    
    AliFlowTrackSimpleCuts *POIfilterVZERO = new AliFlowTrackSimpleCuts();
    if(VZERO_SP || QC){
        POIfilterVZERO->SetEtaMin(-0.8);
        POIfilterVZERO->SetEtaMax(0.8);
        POIfilterVZERO->SetMassMin(263731); POIfilterVZERO->SetMassMax(263733);
        
    }
    
    
    AliFlowTrackSimpleCuts *POIfilterLeftH = new AliFlowTrackSimpleCuts();
    AliFlowTrackSimpleCuts *POIfilterRightH = new AliFlowTrackSimpleCuts();
    if(SP_TPC){
        POIfilterLeftH->SetEtaMin(-0.8);
        POIfilterLeftH->SetEtaMax(0.0);
        POIfilterLeftH->SetMassMin(2636); POIfilterLeftH->SetMassMax(2638);
        
        POIfilterRightH->SetEtaMin(0.0);
        POIfilterRightH->SetEtaMax(0.8);
        POIfilterRightH->SetMassMin(2636); POIfilterRightH->SetMassMax(2638);
    }
    
    
    AliFlowTrackSimpleCuts *POIfilterQCH = new AliFlowTrackSimpleCuts();
    if(QC){
        POIfilterQCH->SetEtaMin(-0.8);
        POIfilterQCH->SetEtaMax(0.8);
        POIfilterQCH->SetMassMin(2636); POIfilterQCH->SetMassMax(2638);
        
    }
    
    
    
    taskHFE->SetRPCuts(cutsRP);
    
    
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("ccontainer0_%s",uniqueID.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
    
    mgr->ConnectInput(taskHFE,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHFE,1,coutput3);
    
    
    if(debug) cout << " === RECEIVED REQUEST FOR FLOW ANALYSIS === " << endl;
    AliAnalysisDataContainer *flowEvent = mgr->CreateContainer(Form("FlowContainer_%s",uniqueID.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(taskHFE, 2, flowEvent);
    if(debug) cout << "    --> Created IO containers " << flowEvent << endl;
    
    
    mgr->AddTask(taskHFE);
    
    if (QC) {  // add qc tasks
        TPCTOFnew::AddQCmethod(Form("QCTPCin_%s",uniqueID.Data()), harmonic, flowEvent,  debug ,uniqueID, -0.8, -0.0, 0.0, 0.8,false,POIfilterVZERO);
        if(debug) cout << "    --> Hanging QC task ...succes! "<< endl;
    }
    if (SP_TPC) {  // add sp subevent tasks
        TPCTOFnew::AddSPmethod(Form("SPTPCQa_in_%s", uniqueID.Data()), -0.8, -.0, .0, +0.8, "Qa", harmonic, flowEvent, false, shrinkSP, debug,uniqueID, false, POIfilterRight);
        if(debug) cout << "    --> Hanging SP Qa task ... succes!" << endl;
        TPCTOFnew::AddSPmethod(Form("SPTPCQb_in_%s", uniqueID.Data()), -0.8, -.0, .0, +0.8, "Qb", harmonic, flowEvent,  false, shrinkSP, debug,uniqueID, false, POIfilterLeft);
        if(debug) cout << "    --> Hanging SP Qb task ... succes!"<< endl;
    }
    if (VZERO_SP) {  // add sp subevent tasks
        TPCTOFnew::AddSPmethod(Form("SPVZEROQa_in_%s", uniqueID.Data()), -0.8, -.0, .0, +0.8, "Qa", harmonic, flowEvent, false, shrinkSP, debug,uniqueID, true, POIfilterVZERO);
        if(debug) cout << "    --> Hanging SP Qa task ... succes!" << endl;
        TPCTOFnew::AddSPmethod(Form("SPVZEROQb_in_%s", uniqueID.Data()), -0.8, -.0, .0, +0.8, "Qb", harmonic, flowEvent, false, shrinkSP, debug,uniqueID, true, POIfilterVZERO);
        if(debug) cout << "    --> Hanging SP Qb task ... succes!"<< endl;
    }
    
    
    return taskHFE;
    
}

//_____________________________________________________________________________
//_____________________________________________________________________________

AliAnalysisTaskFlowITSTPCTOFQCSP* ConfigHFEStandardCuts(Bool_t useMC,Int_t minTPCCulster,AliHFEextraCuts::ITSPixel_t pixel, Bool_t withmultetacorrection1){
    //
    // HFE standard task configuration
    //
    
    Bool_t kAnalyseTaggedTracks = kTRUE;
    
    AliHFEcuts *hfecuts = new AliHFEcuts("hfeCuts","HFE Standard Cuts");  //TODO....change the cuts values to PbPb
    //  hfecuts->CreateStandardCuts();
    hfecuts->SetMinNClustersTPC(minTPCCulster);
    hfecuts->SetMinNClustersITS(5);//5 for ITS pid.....usually i used 3.
    hfecuts->SetMinNTrackletsTRD(0);
    hfecuts->SetMinRatioTPCclusters(0.6);
    
    //   hfecuts->SetEtaRange(-0.9,0.9);
    //   hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    hfecuts->SetRequireITSPixel();
    hfecuts->SetCutITSpixel(pixel);//kAny
    hfecuts->SetMaxChi2perClusterITS(36);  //new from ALberto
    hfecuts->SetMaxChi2perClusterTPC(3.5);
    hfecuts->SetCheckITSLayerStatus(kFALSE); // shud be put back
    //  hfecuts->UnsetVertexRequirement();
    hfecuts->SetVertexRange(10.);
    hfecuts->SetRequireSigmaToVertex();
    //hfecuts->SetSigmaToVertex(10);
    hfecuts->SetTOFPIDStep(kFALSE);
    //  hfecuts->SetQAOn();
    hfecuts->SetPtRange(0, 5.);
    
    AliAnalysisTaskFlowITSTPCTOFQCSP *task = new AliAnalysisTaskFlowITSTPCTOFQCSP("HFE_Flow_TPCTOF");
    printf("*************************************************************");
    printf("task -------------------------------------------- %p\n ", task);
    printf("*************************************************************\n");
    
    
    task->SetHFECuts(hfecuts);
    
    //   task->SetInvariantMassCut(0.05);
    //  task->SetRejectKinkMother(kTRUE);
    //  task->SetRemovePileUp(kTRUE);
    
    // Define PID
    AliHFEpid *pid = task->GetPID();
    if(useMC) pid->SetHasMCData(kTRUE);
    pid->AddDetector("ITS", 0);
    pid->AddDetector("TOF", 1);
    pid->AddDetector("TPC", 2);
    
    
    if(withmultetacorrection1) {
        AliHFEpidTPC *tpcpid = pid->GetDetPID(AliHFEpid::kTPCpid);
        // Theo
        //  task->GetPIDQAManager()->SetFillMultiplicity();
        TF1 *etaCorrMean = GetEtaCorrection("LHC11h_etaCorrMean");
        TF1 *etaCorrWdth = GetEtaCorrection("LHC11h_etaCorrWidth");
        if(etaCorrMean && etaCorrWdth && withmultetacorrection1){
            tpcpid->SetEtaCorrections(etaCorrMean, etaCorrWdth);
            printf("TPC dE/dx Eta correction %p %p\n",etaCorrMean,etaCorrWdth);
        }
        TF1 *centCorrMean = GetCentralityCorrection("LHC11h_multCorrMean");
        TF1 *centCorrWdth = GetCentralityCorrection("LHC11h_multCorrWidth");
        if(centCorrMean && centCorrWdth && withmultetacorrection1){
            tpcpid->SetCentralityCorrections(centCorrMean, centCorrWdth);
            printf("TPC dE/dx multiplicity correction %p %p\n",centCorrMean,centCorrWdth);
        }
        task->SetMultCorrectionTheo(withmultetacorrection1);
        task->SetTPCPID(tpcpid);
    }
    
    
    printf("*************************************\n");
    printf("Configuring standard Task:\n");
    //  task->PrintStatus();
    pid->PrintStatus();
    printf("*************************************\n");
    return task;
    
    
}

//_____________________________________________________________________________

namespace TPCTOFnew{
    //_____________________________________________________________________________
    void AddSPmethod(char *name, double minEtaA, double maxEtaA, double minEtaB, double maxEtaB, char *Qvector, int harmonic, AliAnalysisDataContainer *flowEvent, bool bEP, bool shrink = false, bool debug, TString uniqueID,Bool_t VZERO_SP = kFALSE,  AliFlowTrackSimpleCuts* POIfilter)
    {
        // add sp task and invm filter tasks
        if(debug) (bEP) ? cout << " ****** Reveived request for EP task ****** " << endl : cout << " ******* Switching to SP task ******* " << endl;
        TString fileName = AliAnalysisManager::GetCommonFileName();
        (bEP) ? fileName+=":EP_tpctof" : fileName+=":SP_tpctof";
        //          if(etagap) {
        //            fileName+="_SUBEVENTS";
        //          if(debug) cout << "    --> Setting up subevent analysis <-- " << endl;
        //    }
        if(debug) cout << "    --> fileName " << fileName << endl;
        TString myFolder = fileName;
        if(debug) cout << "    --> myFolder " << myFolder << endl;
        TString myNameSP;
        (bEP) ? myNameSP = Form("%sEPv%d%s", name, harmonic, Qvector): myNameSP = Form("%sSPv%d%s", name, harmonic, Qvector);
        if(debug) cout << " myNameSP " << myNameSP << endl;
        AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
        AliAnalysisDataContainer *flowEventOut = mgr->CreateContainer(Form("Filter_%s",myNameSP.Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
        AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s", myNameSP.Data()), NULL, POIfilter);
        tskFilter->SetSubeventEtaRange(minEtaA, maxEtaA, minEtaB, maxEtaB);
        if(VZERO_SP) tskFilter->SetSubeventEtaRange(-10, 0, 0, 10);
        mgr->AddTask(tskFilter);
        mgr->ConnectInput(tskFilter, 0, flowEvent);
        mgr->ConnectOutput(tskFilter, 1, flowEventOut);
        AliAnalysisDataContainer *outSP = mgr->CreateContainer(myNameSP.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
        AliAnalysisTaskScalarProduct *tskSP = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s", myNameSP.Data()), kFALSE);
        tskSP->SetApplyCorrectionForNUA(kTRUE);
        tskSP->SetHarmonic(harmonic);
        tskSP->SetTotalQvector(Qvector);
        if (bEP) tskSP->SetBehaveAsEP();
        if (shrink) tskSP->SetBookOnlyBasicCCH(kTRUE);
        mgr->AddTask(tskSP);
        mgr->ConnectInput(tskSP, 0, flowEventOut);
        mgr->ConnectOutput(tskSP, 1, outSP);
    }
    //_____________________________________________________________________________
    void AddQCmethod(char *name, int harmonic, AliAnalysisDataContainer *flowEvent, Bool_t debug, TString uniqueID,double minEtaA, double maxEtaA, double minEtaB, double maxEtaB,Bool_t VZERO_SP = kFALSE,  AliFlowTrackSimpleCuts* POIfilter)
    {
        // add qc task and invm filter tasks
        if(debug) cout << " ****** Received request for QC v" << harmonic << " task " << name << ", IO ****** " << flowEvent << endl;
        TString fileName = AliAnalysisManager::GetCommonFileName();
        fileName+=":QC_tpctof";
        if(debug) cout << "    --> Common filename: " << fileName << endl;
        TString myFolder = Form("v%d", harmonic);
        if(debug) cout << "    --> myFolder: " << myFolder << endl;
        TString myName = Form("%s", name);
        if(debug) cout << "    --> myName: " << myName << endl;
        AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
        AliAnalysisDataContainer *flowEventOut = mgr->CreateContainer(Form("Filter_%s", myName.Data()), AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
        AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s", myName.Data()), NULL, POIfilter);
        tskFilter->SetSubeventEtaRange(minEtaA, maxEtaA, minEtaB, maxEtaB);
        //    if(VZERO_SP) tskFilter->SetSubeventEtaRange(-10, 0, 0, 10);
        mgr->AddTask(tskFilter);
        mgr->ConnectInput(tskFilter, 0, flowEvent);
        mgr->ConnectOutput(tskFilter, 1, flowEventOut);
        
        AliAnalysisDataContainer *outQC = mgr->CreateContainer(myName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
        AliAnalysisTaskQCumulants *tskQC = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_%s", myName.Data()), kFALSE);
        tskQC->SetApplyCorrectionForNUA(kTRUE);
        tskQC->SetHarmonic(harmonic);
        tskQC->SetBookOnlyBasicCCH(kTRUE);
        mgr->AddTask(tskQC);
        mgr->ConnectInput(tskQC, 0, flowEventOut);
        mgr->ConnectOutput(tskQC, 1, outQC);
    }
    //_____________________________________________________________________________
}
//_____________________________________________________________________________
TF1* GetCentralityCorrection(TString listname="LHC11h"){
    
    TString etaMap="$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/CentCorrMapsTPC.root";
    
    if (gSystem->AccessPathName(gSystem->ExpandPathName(etaMap.Data()))){
        Error("ConfigHFEpbpb","Eta map not found: %s",etaMap.Data());
        return 0;
    }
    
    TFile f(etaMap.Data());
    if (!f.IsOpen()) return 0;
    gROOT->cd();
    TList *keys=f.GetListOfKeys();
    
    for (Int_t i=0; i<keys->GetEntries(); ++i){
        TString kName=keys->At(i)->GetName();
        TPRegexp reg(kName);
        if (reg.MatchB(listname)){
            printf("Using Eta Correction Function: %s\n",kName.Data());
            return (TF1*)f.Get(kName.Data());
        }
    }
    return 0;
}
//_____________________________________________________________________________
TF1* GetEtaCorrection(TString listname="LHC11h"){
    
    TString etaMap="$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/EtaCorrMapsTPC.root";
    
    if (gSystem->AccessPathName(gSystem->ExpandPathName(etaMap.Data()))){
        Error("ConfigHFEpbpb","Eta map not found: %s",etaMap.Data());
        return 0;
    }
    
    TFile f(etaMap.Data());
    if (!f.IsOpen()) return 0;
    gROOT->cd();
    TList *keys=f.GetListOfKeys();
    
    for (Int_t i=0; i<keys->GetEntries(); ++i){
        TString kName=keys->At(i)->GetName();
        TPRegexp reg(kName);
        if (reg.MatchB(listname)){
            printf("Using Eta Correction Function: %s\n",kName.Data());
            return (TF1*)f.Get(kName.Data());
        }
    }
    return 0;
}
//_____________________________________________________________________________
