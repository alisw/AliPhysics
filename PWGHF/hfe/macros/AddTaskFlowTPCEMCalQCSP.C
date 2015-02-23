///////////////////////////////////////////////////////////////////
//                                                               //
// AddTaskFlowTPCEMCalQCSP macro                                 //
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

AliAnalysisTaskFlowTPCEMCalQCSP*  AddTaskFlowTPCEMCalQCSP(
                                                          TString uniqueID = "",
                                                          Float_t centrMin ,
                                                          Float_t centrMax ,
                                                          Bool_t CentFlatMine,
                                                          Double_t InvmassCut,
                                                          Double_t pTCut,
                                                          Int_t Trigger,
                                                          Bool_t multCorrcut,
                                                          Double_t minTPC,
                                                          Double_t maxTPC,
                                                          Double_t minEovP,
                                                          Double_t maxEovP,
                                                          Double_t minM20,
                                                          Double_t maxM20,
                                                          Double_t minM02,
                                                          Double_t maxM02,
                                                          Double_t Dispersion,
                                                          Int_t minTPCCluster,
                                                          AliHFEextraCuts::ITSPixel_t pixel,
                                                          Bool_t Weight = kFALSE,
                                                          Bool_t withmultetacorrection=kFALSE,
                                                          Double_t etaminpos = 0,
                                                          Double_t etaminneg = 0,
                                                          Bool_t NUA = kTRUE,
                                                          Bool_t PhotonicElectronDCA = kFALSE,
                                                          Int_t TPCClusterforAsso = 80,
                                                          Bool_t AssoITSref = kTRUE,
                                                          Double_t ptminassocut = 0.3,
                                                          Bool_t SideBandsFlow = kFALSE,
                                                          Bool_t Phi_minus_psi = kFALSE,
                                                          const char *Cent = "V0M",
                                                          Bool_t QC = kTRUE, // use qc2 and qc4
                                                          Bool_t SP_TPC = kTRUE, //use tpc sp method
                                                          Bool_t VZERO_SP = kFALSE, // use vzero sp method
                                                          Bool_t BaseH = kFALSE, // base histo
                                                          Int_t harmonic = 2,
                                                          Bool_t shrinkSP = kTRUE,
                                                          Bool_t debug = kFALSE,
                                                          Int_t RPFilterBit = 1,
                                                        //Bool_t purity = kTRUE,
                                                        //Bool_t op_ang = kFALSE,
                                                        //Double_t op_angle_cut=3.

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
    AliAnalysisTaskFlowTPCEMCalQCSP *taskHFE = ConfigHFEemcalMod(kFALSE, minTPCCluster, pixel, withmultetacorrection);    //kTRUE if MC
    
    if(debug) cout << " === AliAnalysisElectronFlow === " << taskHFE << endl;
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
    taskHFE->SetIDCuts(minTPC, maxTPC, minEovP, maxEovP, minM20, maxM20, minM02, maxM02, Dispersion);
    taskHFE->SetFlowSideBands(SideBandsFlow);
    taskHFE->Setphiminuspsi(Phi_minus_psi);
    taskHFE->SetPurity(kTRUE);
    taskHFE->SetpTCuttrack(pTCut);
    taskHFE->SelectPhotonicElectronMethod(PhotonicElectronDCA);
    taskHFE->SetOpeningAngleflag(kFALSE);
    taskHFE->SetOpeningAngleCut(3);
    taskHFE->SetAssoTPCCluster(TPCClusterforAsso);
    taskHFE->SetAssoITSRefit(AssoITSref);
    taskHFE->SetMultCorrelationCut(multCorrcut);
    taskHFE->SetPtMinAssoCut(ptminassocut);
    taskHFE->SetEtaMinPos(etaminpos); //0.2
    taskHFE->SetEtaMinNeg(etaminneg);//-0.2
    
    
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
        cutsRP->SetEtaRange(-0.7, 0.7);
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
    if(VZERO_SP || SP_TPC){
        POIfilterLeft->SetEtaMin(-0.7);
        POIfilterLeft->SetEtaMax(0.0);
        POIfilterLeft->SetMassMin(263731); POIfilterLeft->SetMassMax(263733);
        
        POIfilterRight->SetEtaMin(0.0);
        POIfilterRight->SetEtaMax(0.7);
        POIfilterRight->SetMassMin(263731); POIfilterRight->SetMassMax(263733);
    }
    
    
    AliFlowTrackSimpleCuts *POIfilterQC = new AliFlowTrackSimpleCuts();
    if(QC){
        POIfilterQC->SetEtaMin(-0.7);
        POIfilterQC->SetEtaMax(0.7);
        POIfilterQC->SetMassMin(263731); POIfilterQC->SetMassMax(263733);
        
    }
    
    if(SideBandsFlow){
        
        AliFlowTrackSimpleCuts *POIfilterLeftH = new AliFlowTrackSimpleCuts();
        AliFlowTrackSimpleCuts *POIfilterRightH = new AliFlowTrackSimpleCuts();
        if(SP_TPC){
            POIfilterLeftH->SetEtaMin(-0.7);
            POIfilterLeftH->SetEtaMax(0.0);
            POIfilterLeftH->SetMassMin(2636); POIfilterLeftH->SetMassMax(2638);
            
            POIfilterRightH->SetEtaMin(0.0);
            POIfilterRightH->SetEtaMax(0.7);
            POIfilterRightH->SetMassMin(2636); POIfilterRightH->SetMassMax(2638);
        }
        
        
        AliFlowTrackSimpleCuts *POIfilterQCH = new AliFlowTrackSimpleCuts();
        if(QC){
            POIfilterQCH->SetEtaMin(-0.7);
            POIfilterQCH->SetEtaMax(0.7);
            POIfilterQCH->SetMassMin(2636); POIfilterQCH->SetMassMax(2638);
            
        }
        
    }
    
    taskHFE->SetRPCuts(cutsRP);
    
    
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("ccontainer0_%s",uniqueID.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
    
    mgr->ConnectInput(taskHFE,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHFE,1,coutput3);
    
    
    if(debug) cout << " === RECEIVED REQUEST FOR FLOW ANALYSIS === " << endl;
    AliAnalysisDataContainer *flowEvent = mgr->CreateContainer(Form("FlowContainer_%s",uniqueID.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(taskHFE, 2, flowEvent);
    if(debug) cout << "    --> Created IO containers " << flowEvent << endl;
    
    if(SideBandsFlow){
        if(debug) cout << " === RECEIVED REQUEST FOR FLOW ANALYSIS === " << endl;
        AliAnalysisDataContainer *flowEventCont = mgr->CreateContainer(Form("FlowContainer_Cont_%s",uniqueID.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
        mgr->ConnectOutput(taskHFE, 3, flowEventCont);
        if(debug) cout << "    --> Created IO containers " << flowEventCont << endl;
    }
    
    
    mgr->AddTask(taskHFE);
    
    if (QC) {  // add qc tasks
        AddQCmethod(Form("QCTPCin_%s",uniqueID.Data()), harmonic, flowEvent,  debug ,uniqueID, -0.7, -0.0, 0.0, 0.7,false,POIfilterQC, NUA, shrinkSP, BaseH);
        if(debug) cout << "    --> Hanging QC task ...succes! "<< endl;
    }
    if (SP_TPC) {  // add sp subevent tasks
        AddSPmethod(Form("SPTPCQa_in_%s", uniqueID.Data()), -0.7, -.0, .0, +0.7, "Qa", harmonic, flowEvent, false, shrinkSP, debug,uniqueID, false, POIfilterRight, NUA,BaseH);
        if(debug) cout << "    --> Hanging SP Qa task ... succes!" << endl;
        AddSPmethod(Form("SPTPCQb_in_%s", uniqueID.Data()), -0.7, -.0, .0, +0.7, "Qb", harmonic, flowEvent,  false, shrinkSP, debug,uniqueID, false, POIfilterLeft, NUA,BaseH);
        if(debug) cout << "    --> Hanging SP Qb task ... succes!"<< endl;
    }
    if (VZERO_SP) {  // add sp subevent tasks
        AddSPmethod(Form("SPVZEROQa_in_%s", uniqueID.Data()), -0.7, -.0, .0, +0.7, "Qa", harmonic, flowEvent, false, shrinkSP, debug,uniqueID, true, POIfilterRight, NUA,BaseH);
        if(debug) cout << "    --> Hanging SP Qa task ... succes!" << endl;
        AddSPmethod(Form("SPVZEROQb_in_%s", uniqueID.Data()), -0.7, -.0, .0, +0.7, "Qb", harmonic, flowEvent,  false, shrinkSP, debug,uniqueID, true, POIfilterLeft, NUA,BaseH);
        if(debug) cout << "    --> Hanging SP Qb task ... succes!"<< endl;
    }
    
    //=========================================Flow event for elctronContamination==============================================================================================
    if(SideBandsFlow){
        if (QC) {  // add qc tasks
            AddQCmethod(Form("QCTPCCont_%s",uniqueID.Data()), harmonic, flowEventCont,  debug ,uniqueID, -0.7, -0.0, 0.0, 0.7,false,POIfilterQCH, NUA, shrinkSP,BaseH);
            if(debug) cout << "    --> Hanging QC task ...succes! "<< endl;
        }
        if (SP_TPC) {  // add sp subevent tasks
            AddSPmethod(Form("SPTPCQa_Cont_%s", uniqueID.Data()), -0.7, -.0, .0, +0.7, "Qa", harmonic, flowEventCont, false, shrinkSP, debug,uniqueID, false, POIfilterRightH, NUA,BaseH);
            if(debug) cout << "    --> Hanging SP Qa task ... succes!" << endl;
            AddSPmethod(Form("SPTPCQb_Cont_%s", uniqueID.Data()), -0.7, -.0, .0, +0.7, "Qb", harmonic, flowEventCont,  false, shrinkSP, debug,uniqueID, false, POIfilterLeftH, NUA,BaseH);
            if(debug) cout << "    --> Hanging SP Qb task ... succes!"<< endl;
        }
    }
    //==========================================================================================================================================================================
    
    
    return taskHFE;
    
}

//_____________________________________________________________________________


//_____________________________________________________________________________
void AddSPmethod(char *name, double minEtaA, double maxEtaA, double minEtaB, double maxEtaB, char *Qvector, int harmonic, AliAnalysisDataContainer *flowEvent, bool bEP, bool shrink = false, bool debug, TString uniqueID, Bool_t VZERO_SP = kFALSE,  AliFlowTrackSimpleCuts* POIfilter, Bool_t kNUA, Bool_t kBaseH)
{
    // add sp task and invm filter tasks
    if(debug) (bEP) ? cout << " ****** Reveived request for EP task ****** " << endl : cout << " ******* Switching to SP task ******* " << endl;
    TString fileName = AliAnalysisManager::GetCommonFileName();
    (bEP) ? fileName+=":EP" : fileName+=":SP";
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
    tskSP->SetApplyCorrectionForNUA(kNUA);
    tskSP->SetHarmonic(harmonic);
    tskSP->SetTotalQvector(Qvector);
    if (bEP) tskSP->SetBehaveAsEP();
    if (shrink)tskSP->SetBookOnlyBasicCCH(kBaseH);
    mgr->AddTask(tskSP);
    mgr->ConnectInput(tskSP, 0, flowEventOut);
    mgr->ConnectOutput(tskSP, 1, outSP);
}
//_____________________________________________________________________________
void AddQCmethod(char *name, int harmonic, AliAnalysisDataContainer *flowEvent, Bool_t debug, TString uniqueID,double minEtaA, double maxEtaA, double minEtaB, double maxEtaB,Bool_t VZERO_SP = kFALSE,  AliFlowTrackSimpleCuts* POIfilter, Bool_t kNUA, bool shrink = false, Bool_t kBaseH)
{
    // add qc task and invm filter tasks
    if(debug) cout << " ****** Received request for QC v" << harmonic << " task " << name << ", IO ****** " << flowEvent << endl;
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName+=":QC";
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
    tskQC->SetApplyCorrectionForNUA(kNUA);
    tskQC->SetHarmonic(harmonic);
    if (shrink)tskQC->SetBookOnlyBasicCCH(kBaseH);
    mgr->AddTask(tskQC);
    mgr->ConnectInput(tskQC, 0, flowEventOut);
    mgr->ConnectOutput(tskQC, 1, outQC);
}
//_____________________________________________________________________________


//_____________________________________________________________________________

AliAnalysisTaskFlowTPCEMCalQCSP* ConfigHFEemcalMod(Bool_t useMC,Int_t minTPCCulster,AliHFEextraCuts::ITSPixel_t pixel, Bool_t withmultetacorrection1){
    //
    // HFE standard task configuration
    //
    
    Bool_t kAnalyseTaggedTracks = kTRUE;
    
    AliHFEcuts *hfecuts = new AliHFEcuts("hfeCutsEMCAL","HFE Standard Cuts");  //TODO....change the cuts values to PbPb
    //  hfecuts->CreateStandardCuts();
    hfecuts->SetMinNClustersTPC(minTPCCulster);
    hfecuts->SetMinNClustersITS(3);
    hfecuts->SetMinNTrackletsTRD(0);
    hfecuts->SetMinRatioTPCclusters(0.6);
    
    //   hfecuts->SetEtaRange(-0.9,0.9);
    //   hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    hfecuts->SetRequireITSPixel();
    hfecuts->SetCutITSpixel(pixel);//kAny
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
    
    AliAnalysisTaskFlowTPCEMCalQCSP *task = new AliAnalysisTaskFlowTPCEMCalQCSP("HFE_Flow_TPCEMCal");
    printf("task ------------------------ %p\n ", task);
    
    
    task->SetHFECuts(hfecuts);
    
    //   task->SetInvariantMassCut(0.05);
    //  task->SetRejectKinkMother(kTRUE);
    //  task->SetRemovePileUp(kTRUE);
    
    // Define PID
    AliHFEpid *pid = task->GetPID();
    if(useMC) pid->SetHasMCData(kTRUE);
    pid->AddDetector("TPC", 0);
    pid->AddDetector("EMCAL", 1);
    
    
    if(withmultetacorrection1) {
        AliHFEpidTPC *tpcpid = pid->GetDetPID(AliHFEpid::kTPCpid);
        // Theo
        //  task->GetPIDQAManager()->SetFillMultiplicity();
        TF1 *etaCorrMean = GetEtaCorrectionEMCal("LHC11h_etaCorrMean");
        TF1 *etaCorrWdth = GetEtaCorrectionEMCal("LHC11h_etaCorrWidth");
        if(etaCorrMean && etaCorrWdth && withmultetacorrection1){
            tpcpid->SetEtaCorrections(etaCorrMean, etaCorrWdth);
            printf("TPC dE/dx Eta correction %p %p\n",etaCorrMean,etaCorrWdth);
        }
        TF1 *centCorrMean = GetCentralityCorrectionEMCal("LHC11h_multCorrMean");
        TF1 *centCorrWdth = GetCentralityCorrectionEMCal("LHC11h_multCorrWidth");
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
//_____________________________________________________________________________
TF1* GetCentralityCorrectionEMCal(TString listname="LHC11h"){
    
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
TF1* GetEtaCorrectionEMCal(TString listname="LHC11h"){
    
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



