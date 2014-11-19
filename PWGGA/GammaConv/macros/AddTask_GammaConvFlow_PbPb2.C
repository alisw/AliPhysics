class AliAnalysisDataContainer;
class AliFlowTrackCuts;
class AliFlowTrackSimpleCuts;
class AliFlowEventCuts;
class AliFlowEventSimpleCuts;
class AliAnalysisDataContainer;



//===================from here Flow pakage================
void AddSPmethod(char *name, char *Qvector, int harmonic, AliAnalysisDataContainer *flowEvent,  bool shrink = false, bool debug, TString uniqueID,/*Bool_t VZERO_SP = kFALSE,*/  AliFlowTrackSimpleCuts* POIfilter, int trainConfig)
{
    // add sp task and invm filter tasks
    if(debug) (bEP) ? cout << " ****** Reveived request for EP task ****** " << endl : cout << " ******* Switching to SP task ******* " << endl;
    TString fileName = Form("GammaConvFlow_%i.root",trainConfig);
    //    (bEP) ? fileName+=":EP_tpctof" : fileName+=":SP_tpctof";
    //          if(etagap) {
    //            fileName+="_SUBEVENTS";
    //          if(debug) cout << "    --> Setting up subevent analysis <-- " << endl;
    //    }
    if(debug) cout << "    --> fileName " << fileName << endl;
    TString myFolder = fileName;
    if(debug) cout << "    --> myFolder " << myFolder << endl;
    TString myNameSP;
    //(bEP) ? myNameSP = Form("%sEPv%d%s", name, harmonic, Qvector): myNameSP = Form("%sSPv%d%s", name, harmonic, Qvector);
    myNameSP = Form("%sSPv%d%s", name, harmonic, Qvector);
    
    if(debug) cout << " myNameSP " << myNameSP << endl;
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliAnalysisDataContainer *flowEventOut = mgr->CreateContainer(Form("Filter_%s",myNameSP.Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s", myNameSP.Data()), NULL, POIfilter);
    //  tskFilter->SetSubeventEtaRange(minEtaA, maxEtaA, minEtaB, maxEtaB);
    //if(VZERO_SP) tskFilter->SetSubeventEtaRange(-10, 0, 0, 10);
    tskFilter->SetSubeventEtaRange(-10, -1, 1, 10);
    mgr->AddTask(tskFilter);
    mgr->ConnectInput(tskFilter, 0, flowEvent);
    mgr->ConnectOutput(tskFilter, 1, flowEventOut);
    AliAnalysisDataContainer *outSP = mgr->CreateContainer(myNameSP.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
    AliAnalysisTaskScalarProduct *tskSP = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s", myNameSP.Data()), kFALSE);
    tskSP->SetApplyCorrectionForNUA(kTRUE);
    tskSP->SetHarmonic(harmonic);
    tskSP->SetTotalQvector(Qvector);
    // if (bEP) tskSP->SetBehaveAsEP();
    if (shrink) tskSP->SetBookOnlyBasicCCH(kTRUE);
    mgr->AddTask(tskSP);
    mgr->ConnectInput(tskSP, 0, flowEventOut);
    mgr->ConnectOutput(tskSP, 1, outSP);
}
//_____________________________________________________________________________



void AddTask_GammaConvFlow_PbPb2(
                               TString uniqueID = "",
                               Int_t harmonic=2,
                               Int_t trainConfig = 1,  //change different set of cuts
                               //Bool_t isMC   = kFALSE, //run MC
                               Int_t enableQAMesonTask = 0, //enable QA in AddTask_GammaConvFlow_PbPb2
                               Int_t enableQAPhotonTask = 0, // enable additional QA task
                               //TString fileNameInputForWeighting = "MCSpectraInput.root", // path to file for weigting input
                               Bool_t doWeighting = kFALSE,  //enable Weighting
                               TString cutnumberAODBranch = "1000000060084000001500000",
                               Bool_t shrinkSP = kTRUE,
                               Bool_t debug = kFALSE
                               ) {
    
    // ================= Load Librariers =================================
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libMinuit");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCDB.so");
    gSystem->Load("libSTEER.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libTENDER.so");
    gSystem->Load("libTENDERSupplies.so");
	gSystem->Load("libPWGflowBase.so");
	gSystem->Load("libPWGflowTasks.so");
    gSystem->Load("libPWGGAGammaConv.so");
    
    Int_t isHeavyIon = 1;
    
    // ================== GetAnalysisManager ===============================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
        return ;
    }
    
    // ================== GetInputEventHandler =============================
    AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
    
    //========= Add PID Reponse to ANALYSIS manager ====
    if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
        gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
        AddTaskPIDResponse(isMC);
    }
    
    //=========  Set Cutnumber for V0Reader ================================
    TString cutnumberPhoton = "000084001001500000000";
    TString cutnumberEvent = "1000000";
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
    if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
        AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
        
        fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
        fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
        fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
        
        if (!mgr) {
            Error("AddTask_V0ReaderV1", "No analysis manager found.");
            return;
        }
        
        AliConvEventCuts *fEventCuts=NULL;
        if(cutnumberEvent!=""){
            fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
            fEventCuts->SetPreSelectionCutFlag(kTRUE);
            if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
                fV0ReaderV1->SetEventCuts(fEventCuts);
                fEventCuts->SetFillCutHistograms("",kTRUE);
            }
        }
        
        // Set AnalysisCut Number
        AliConversionPhotonCuts *fCuts=NULL;
        if(cutnumberPhoton!=""){
            fCuts= new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
            fCuts->SetPreSelectionCutFlag(kTRUE);
            fCuts->SetIsHeavyIon(isHeavyIon);
            if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
                fV0ReaderV1->SetConversionCuts(fCuts);
                fCuts->SetFillCutHistograms("",kTRUE);
            }
        }
        
        if(inputHandler->IsA()==AliAODInputHandler::Class()){
            // AOD mode
            fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
        }
        fV0ReaderV1->Init();
        
        AliLog::SetGlobalLogLevel(AliLog::kInfo);
        
        //connect input V0Reader
        mgr->AddTask(fV0ReaderV1);
        mgr->ConnectInput(fV0ReaderV1,0,cinput);
        
    }
    
    //================================================
    //========= Add task to the ANALYSIS manager =====
    //================================================
    AliAnalysisTaskGammaConvFlow *task=NULL;
    task= new AliAnalysisTaskGammaConvFlow(Form("GammaConvV1_%i",trainConfig));
    task->SetIsHeavyIon(isHeavyIon);
    //task->SetIsMC(isMC);
    
    
//=====================from here call/book the flow package stuff==============
    //set RP cuts for flow package analysis
    cutsRP = new AliFlowTrackCuts(Form("RFPcuts%s",uniqueID));
    if(!cutsRP) {
        if(debug) cout << " Fatal error: no RP cuts found, could be a library problem! " << endl;
        return 0x0;
    }
    cutsRP = cutsRP->GetStandardVZEROOnlyTrackCuts(); // select vzero tracks
    cutsRP->SetVZEROgainEqualizationPerRing(kFALSE);
    cutsRP->SetApplyRecentering(kTRUE);
    if(debug) cout << "    --> VZERO RP's " << cutsRP << endl;
    
    task->SetRPCuts(cutsRP);
    
    AliFlowTrackSimpleCuts *POIfilterVZERO = new AliFlowTrackSimpleCuts();
    POIfilterVZERO->SetEtaMin(-0.8);
    POIfilterVZERO->SetEtaMax(0.8);
    //  POIfilterVZERO->SetMassMin(263731); POIfilterVZERO->SetMassMax(263733);
    
    
    if(debug) cout << " === RECEIVED REQUEST FOR FLOW ANALYSIS === " << endl;
    AliAnalysisDataContainer *flowEvent = mgr->CreateContainer(Form("FlowContainer_%s",uniqueID.Data()), AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(task, 2, flowEvent);
    if(debug) cout << "    --> Created IO containers " << flowEvent << endl;
    
    AddSPmethod(Form("SPVZEROQa_in_%s", uniqueID.Data()), "Qa", harmonic, flowEvent, shrinkSP, debug, uniqueID, POIfilterVZERO, trainConfig);
    if(debug) cout << "    --> Hanging SP Qa task ... succes!" << endl;
    AddSPmethod(Form("SPVZEROQb_in_%s", uniqueID.Data()), "Qb", harmonic, flowEvent, shrinkSP, debug, uniqueID, POIfilterVZERO, trainConfig);
    if(debug) cout << "    --> Hanging SP Qb task ... succes!"<< endl;
//================================END flow stuff===========================================
    
    
    
    
    
    
    // Cut Numbers to use in Analysis
    Int_t numberOfCuts = 1;
    
    TString *eventCutArray = new TString[numberOfCuts];
    TString *photonCutArray = new TString[numberOfCuts];
    TString *mesonCutArray = new TString[numberOfCuts];
    
    if (trainConfig == 1){
        eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522045009000";
    } else if (trainConfig == 2) {
        eventCutArray[ 0] = "6120001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522045009000";
    } else if (trainConfig == 3) {
        eventCutArray[ 0] = "5010001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522045009000";
    } else if (trainConfig == 4) {
        eventCutArray[ 0] = "5020001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522045009000";
    } else if (trainConfig == 5) {
        eventCutArray[ 0] = "5120001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522045009000";
    } else if (trainConfig == 6) {
        eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522045009000";
    } else if (trainConfig == 7) {
        eventCutArray[ 0] = "5460001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
    } else if (trainConfig == 8) {
        eventCutArray[ 0] = "5480001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
    } else if (trainConfig == 9) {
        eventCutArray[ 0] = "5450001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
    } else if (trainConfig == 10) {
        eventCutArray[ 0] = "5560001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
    } else if (trainConfig == 11) {
        eventCutArray[ 0] = "5680001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
    } else if (trainConfig == 12) {
        eventCutArray[ 0] = "5670001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
    } else if (trainConfig == 13) {
        eventCutArray[ 0] = "5780001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
    } else if (trainConfig == 14) {
        eventCutArray[ 0] = "4690001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
    } else if (trainConfig == 15) {
        eventCutArray[ 0] = "5890001"; photonCutArray[ 0] = "042092970023220000000"; mesonCutArray[ 0] = "01522065009000";
    } else  if (trainConfig == 16){
        eventCutArray[ 0] = "6010001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 17) {
        eventCutArray[ 0] = "6120001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 18) {
        eventCutArray[ 0] = "5010001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 19) {
        eventCutArray[ 0] = "5020001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 20) {
        eventCutArray[ 0] = "5120001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 21) {
        eventCutArray[ 0] = "5240001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 22) {
        eventCutArray[ 0] = "5460001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 23) {
        eventCutArray[ 0] = "5480001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 24) {
        eventCutArray[ 0] = "5450001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 25) {
        eventCutArray[ 0] = "5560001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 26) {
        eventCutArray[ 0] = "5680001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 27) {
        eventCutArray[ 0] = "5670001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 28) {
        eventCutArray[ 0] = "5780001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 29) {
        eventCutArray[ 0] = "4690001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 30) {
        eventCutArray[ 0] = "5890001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 31) {
        eventCutArray[ 0] = "5080001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 32) {
        eventCutArray[ 0] = "5250001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 33) {
        eventCutArray[ 0] = "5350001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 34) {
        eventCutArray[ 0] = "5450001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 35) {
        eventCutArray[ 0] = "5340001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else if (trainConfig == 36) {
        eventCutArray[ 0] = "5230001"; photonCutArray[ 0] = "002092970028250400000"; mesonCutArray[ 0] = "01525065000000";
    } else {
        Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
        return;
    }
    
    TList *EventCutList = new TList();
    TList *ConvCutList = new TList();
    TList *MesonCutList = new TList();
    
    TList *HeaderList = new TList();
    TObjString *Header1 = new TObjString("BOX");
    HeaderList->Add(Header1);
    //    TObjString *Header3 = new TObjString("eta_2");
    //    HeaderList->Add(Header3);
    
    EventCutList->SetOwner(kTRUE);
    AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
    ConvCutList->SetOwner(kTRUE);
    AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
    MesonCutList->SetOwner(kTRUE);
    AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
    
    for(Int_t i = 0; i<numberOfCuts; i++){
        analysisEventCuts[i] = new AliConvEventCuts();
        analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
        EventCutList->Add(analysisEventCuts[i]);
        analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
        
        analysisCuts[i] = new AliConversionPhotonCuts();
        analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
        ConvCutList->Add(analysisCuts[i]);
        analysisCuts[i]->SetFillCutHistograms("",kFALSE);
        
        analysisMesonCuts[i] = new AliConversionMesonCuts();
        analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
        MesonCutList->Add(analysisMesonCuts[i]);
        analysisMesonCuts[i]->SetFillCutHistograms("");
        
        analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    }
    
    task->SetEventCutList(numberOfCuts,EventCutList);
    task->SetConversionCutList(numberOfCuts,ConvCutList);
//    task->SetMesonCutList(numberOfCuts,MesonCutList);
//    task->SetMoveParticleAccordingToVertex(kTRUE);
    task->SetDoMesonAnalysis(kFALSE);
    task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
    task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
    
    //connect containers
    AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
                         AliAnalysisManager::kOutputContainer,Form("GammaConvFlow_%i.root",trainConfig));
    
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,cinput);
    mgr->ConnectOutput(task,1,coutput);
    
    return;
    
}



