#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C>
#endif
/***************************************************************************
//            Modified by Bong-Hwi Lim - 08/05/2019
//            Based on AddTaskRare_pp13
//
// Macro to configure the Xi1530 analysis task
//
// Input parameters:
// isMC: kTRUE if MC
// system: 0 for pp, 1 for pPb, 2 for PbPb
// EventCuts:      1*tigger(0 for INT7, 1 for HighMultV0)
//            +   10*Multibins(10 for V0M, 20 for RefMult08)
//            + 1000*EventMixing(1000 for no mxing, n*1000 for n mixing,
default:5)
// eg1. 10: INT7, V0M Mutliplicity, 5 mixing
// eg2. 20011: HighMultV0, V0M Multiplicity, 20 mixing
****************************************************************************/
enum ERsnCollType_t { kPP = 0, kPPb, kPbPb };

void AddMonitorOutput_P(TString s = "",
                        TObjArray* m = 0,
                        TString o = "",
                        AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_Pt(TString s = "",
                         TObjArray* m = 0,
                         TString o = "",
                         AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_Eta(TString s = "",
                          TObjArray* m = 0,
                          TString o = "",
                          AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_DCAxy(TString s = "",
                            TObjArray* m = 0,
                            TString o = "",
                            AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_DCAz(TString s = "",
                           TObjArray* m = 0,
                           TString o = "",
                           AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0NPt(TString s = "",
                            TObjArray* m = 0,
                            TString o = "",
                            AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0PPt(TString s = "",
                            TObjArray* m = 0,
                            TString o = "",
                            AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0Mass(TString s = "",
                             TObjArray* m = 0,
                             TString o = "",
                             AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0DCA(TString n = "",
                            TObjArray* m = 0,
                            TString o = "",
                            AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0Radius(TString s = "",
                               TObjArray* m = 0,
                               TString o = "",
                               AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0Lifetime(TString s = "",
                                 TObjArray* m = 0,
                                 TString o = "",
                                 AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0DaughterDCA(TString s = "",
                                    TObjArray* m = 0,
                                    TString o = "",
                                    AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0DCA2TPV(TString s = "",
                                TObjArray* m = 0,
                                TString o = "",
                                AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0CPA(TString s = "",
                            TObjArray* m = 0,
                            TString o = "",
                            AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0TPCpim(TString s = "",
                               TObjArray* m = 0,
                               TString o = "",
                               AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_V0TPCpip(TString s = "",
                               TObjArray* m = 0,
                               TString o = "",
                               AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_LambdaProtonPID(TObjArray* m = 0,
                                      TString o = "",
                                      AliRsnLoopDaughter* l = 0);
void AddMonitorOutput_LambdaAntiProtonPID(TObjArray* m = 0,
                                          TString o = "",
                                          AliRsnLoopDaughter* l = 0);

Bool_t Config_Xipi(  // From Anders's master macro.
    AliRsnMiniAnalysisTask* task,
    TString lname = "Xipi",
    Bool_t isMC = kFALSE,
    Int_t system = 0,
    Int_t EventCuts = 0,
    Int_t TrackCutsXi = 0,
    Int_t TrackCutsPi = 0);
AliRsnMiniAnalysisTask* AddTaskXi1530_RsnPack(TString lname,
                                                 Bool_t isMC,
                                                 Int_t system,
                                                 Int_t EventCuts = 0,
                                                 Int_t TrackCuts1 = 0,
                                                 Int_t TrackCuts2 = 0) {
    // ----- INITIALIZATION -----

    // retrieve analysis manager
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskXi1530_RsnPack",
                "No analysis manager to connect to.");
        return NULL;
    }

    // create the task and configure
    AliRsnMiniAnalysisTask* task = new AliRsnMiniAnalysisTask(lname, isMC);

    // trigger
    int trigger = EventCuts % 10;
    if (!trigger)
        task->UseESDTriggerMask(AliVEvent::kINT7);
    else if (trigger == 1)
        task->UseESDTriggerMask(AliVEvent::kHighMultV0);

    // multiplicity
    bool isPP = false;
    if (!system)
        isPP = true;
    int MultBins = (EventCuts / 10) % 10;
    if (system == 1 || system == 2)
        MultBins = 1;

    if (isPP) {
        if (MultBins == 1)
            task->UseMultiplicity("AliMultSelection_V0M");
        else if (MultBins == 2)
            task->UseMultiplicity("AliMultSelection_RefMult08");
        else
            task->UseMultiplicity("QUALITY");
    } else if (system == 1)
        task->UseMultiplicity("AliMultSelection_V0A");
    else if (system == 2)
        task->UseMultiplicity("AliMultSelection_V0M");
    else
        task->UseCentrality("V0M");

    // set event mixing options
    int nmix = 5;
    if ((EventCuts % 10000) / 1000 == 1)
        nmix = 0;
    float maxDiffVzMix = 1;
    float maxDiffMultMix = 5;
    task->UseContinuousMix();
    task->SetNMix(nmix);
    task->SetMaxDiffVz(maxDiffVzMix);
    task->SetMaxDiffMult(maxDiffMultMix);
    ::Info("AddTaskXi1530_RsnPack", "%s",
           Form("Event mixing configuration: \n events to mix = %i \n max "
                "diff. vtxZ = cm %5.3f \n max diff multi = %5.3f",
                nmix, maxDiffVzMix, maxDiffMultMix));

    // vertex cuts
    float vtxZcut = 10;
    Bool_t rejectPileUp = kTRUE;
    AliRsnCutPrimaryVertex* cutVertex = 0;
    if (!MultBins || fabs(vtxZcut - 10.) > 1.e-10) {
        cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vtxZcut, 0, kFALSE);
        if (!MultBins) {
            cutVertex->SetCheckZResolutionSPD();
            cutVertex->SetCheckDispersionSPD();
            cutVertex->SetCheckZDifferenceSPDTrack();
        }
        if (0)
            cutVertex->SetCheckGeneratedVertexZ();
    }

    // other event selection cuts
    AliRsnCutEventUtils* cutEventUtils = 0;
    if (1) {
        cutEventUtils =
            new AliRsnCutEventUtils("cutEventUtils", kTRUE, rejectPileUp);
        if (!MultBins) {
            cutEventUtils->SetCheckIncompleteDAQ();
            cutEventUtils->SetCheckSPDClusterVsTrackletBG();
        } else {
            cutEventUtils->SetRemovePileUppA2013(kFALSE);
            cutEventUtils->SetCheckAcceptedMultSelection();
        }
    }

    // set the check for pileup
    if (isPP && (!isMC) && cutVertex) {
        cutVertex->SetCheckPileUp(rejectPileUp);
        ::Info("AddTaskXi1530_RsnPack", "%s",
               Form(":::::::::::::::::: Pile-up rejection mode: %s",
                    (rejectPileUp) ? "ON" : "OFF"));
    }

    // define and fill cut set for event cuts
    AliRsnCutSet* eventCuts = 0;
    if (cutEventUtils || cutVertex) {
        eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);

        if (cutEventUtils && cutVertex) {
            eventCuts->AddCut(cutEventUtils);
            eventCuts->AddCut(cutVertex);
            eventCuts->SetCutScheme(
                Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));
        } else if (cutEventUtils && !cutVertex) {
            eventCuts->AddCut(cutEventUtils);
            eventCuts->SetCutScheme(Form("%s", cutEventUtils->GetName()));
        } else if (!cutEventUtils && cutVertex) {
            eventCuts->AddCut(cutVertex);
            eventCuts->SetCutScheme(Form("%s", cutVertex->GetName()));
        }

        task->SetEventCuts(eventCuts);
    }

    // ----- EVENT-ONLY COMPUTATIONS -----

    Double_t multbins[1000];
    int j, nmult = 0;
    if (!MultBins) {
        for (j = 0; j <= 401; j++) {
            multbins[nmult] = j - 0.5;
            nmult++;
        }
    } else if (!trigger) {
        for (j = 0; j <= 100; j++) {
            multbins[nmult] = j;
            nmult++;
        }
    } else {
        for (j = 0; j < 10; j++) {
            multbins[nmult] = 0.0001 * j;
            nmult++;
        }
        for (j = 1; j < 10; j++) {
            multbins[nmult] = 0.001 * j;
            nmult++;
        }
        for (j = 1; j < 10; j++) {
            multbins[nmult] = 0.01 * j;
            nmult++;
        }
        for (j = 1; j < 10; j++) {
            multbins[nmult] = 0.1 * j;
            nmult++;
        }
        for (j = 1; j <= 100; j++) {
            multbins[nmult] = j;
            nmult++;
        }
    }
    nmult--;

    // vertex
    Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
    AliRsnMiniOutput* outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
    outVtx->AddAxis(vtxID, 240, -12.0, 12.0);

    // multiplicity or centrality
    Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
    AliRsnMiniOutput* outMult =
        task->CreateOutput("eventMult", "HIST", "EVENT");
    outMult->AddAxis(multID, nmult + 1, multbins);

    TH1F* hEventsVsMulti = new TH1F("hAEventsVsMulti", "", nmult, multbins);
    task->SetEventQAHist(
        "EventsVsMulti",
        hEventsVsMulti);  // custom binning for fHAEventsVsMulti

    double ybins[1000];
    for (j = 0; j <= 240; j++)
        ybins[j] = -12 + 0.1 * j;

    TH2F* hvz = new TH2F("hVzVsCent", "", nmult, multbins, 240, ybins);
    task->SetEventQAHist(
        "vz", hvz);  // plugs this histogram into the fHAEventVz data member

    for (j = 0; j <= 401; j++)
        ybins[j] = j - 0.5;

    TH2F* hmc = new TH2F("MultiVsCent", "", nmult, multbins, 401, ybins);
    hmc->GetYaxis()->SetTitle("QUALITY");
    task->SetEventQAHist(
        "multicent",
        hmc);  // plugs this histogram into the fHAEventMultiCent data member

    // ----- CONFIGURE -----

    cerr << "configuring" << endl;
    Config_Xipi(task, lname, isMC, system, EventCuts, TrackCuts1, TrackCuts2);
    cerr << "done configuring" << endl;

    // ----- CONTAINERS -----

    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    Printf("AddTaskXi1530_RsnPack - Set OutputFileName : \n %s\n",
           outputFileName.Data());

    AliAnalysisDataContainer* output = mgr->CreateContainer(
        Form("RsnOut_%s", lname.Data()), TList::Class(),
        AliAnalysisManager::kOutputContainer, outputFileName);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, output);

    return task;
}

Bool_t Config_Xipi(AliRsnMiniAnalysisTask* task, TString lname,
                   Bool_t isMC, Int_t system, Int_t EventCuts,
                   Int_t TrackCutsXi, Int_t TrackCutsPi) {
    bool isPP = false;
    if (!system)
        isPP = true;
    int trigger = EventCuts % 10;
    int MultBins = (EventCuts / 10) % 10;
    if (system == 1 || system == 2)
        MultBins = 1;

    char suffix[1000];
    sprintf(suffix, "_%s", lname.Data());
    Bool_t enableMonitor = kTRUE;

    // set cuts for primary pion
    if (!(TrackCutsPi % 10000))
        TrackCutsPi += 3020;  // default settings
    Float_t nsigmaPiTPC = 0.1 * (TrackCutsPi % 100);
    Float_t nsigmaPiTOF = 0.1 * ((TrackCutsPi / 100) % 100);
    Int_t CutTypePi =
        (TrackCutsPi / 10000) %
        100;  // 0=TPC+TOFveto (default), 1=TPC only, 2=TOF only, 3 TPC+TOFcut
    Int_t MisidentifiedAsKaon =
        (TrackCutsPi / 1000000) %
        10;  // 0=pion assigned pion mass, 1=pion assigned kaon mass (for
             // Xi(1820)- analysis)

    AliRsnCutTrackQuality* trkQualityCut =
        new AliRsnCutTrackQuality("myQualityCut");
    trkQualityCut->SetDefaults2011(kTRUE, kTRUE);

    AliRsnCutSetDaughterParticle* cutSetQ = new AliRsnCutSetDaughterParticle(
        "cutQ", trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2010,
        AliPID::kPion, -1.);

    AliRsnCutSetDaughterParticle* cutSetPi = 0;
    if (!CutTypePi)
        cutSetPi = new AliRsnCutSetDaughterParticle(
            Form("cutPi%i_%2.1fsigma",
                 AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,
                 nsigmaPiTPC),
            trkQualityCut, AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,
            AliPID::kPion, nsigmaPiTPC, nsigmaPiTOF);
    else if (CutTypePi == 1)
        cutSetPi = new AliRsnCutSetDaughterParticle(
            Form("cutPi%i_%2.1fsigma",
                 AliRsnCutSetDaughterParticle::kFastTPCpidNsigma, nsigmaPiTPC),
            trkQualityCut, AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
            AliPID::kPion, nsigmaPiTPC, -1.);
    else if (CutTypePi == 2)
        cutSetPi = new AliRsnCutSetDaughterParticle(
            Form("cutPi%i_%2.1fsigma",
                 AliRsnCutSetDaughterParticle::kFastTOFpidNsigma, nsigmaPiTOF),
            trkQualityCut, AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,
            AliPID::kPion, -1., nsigmaPiTOF);
    else if (CutTypePi == 3)
        cutSetPi = new AliRsnCutSetDaughterParticle(
            Form("cutPi%i_%2.1fsigma",
                 AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011,
                 nsigmaPiTPC),
            trkQualityCut,
            AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011,
            AliPID::kPion, nsigmaPiTPC, -1.);
    if (!cutSetPi) {
        cerr << "Error in AddTaskXi1530_RsnPack::Config_Xipi(): missing cutSetPi"
             << endl;
        return kFALSE;
    }

    Int_t iCutQ = task->AddTrackCuts(cutSetQ);
    Int_t iCutPi = task->AddTrackCuts(cutSetPi);

    // selections for Xi daughters

    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("qualityTrackcut");
    esdTrackCuts->SetEtaRange(-0.8, 0.8);
    esdTrackCuts->SetRequireTPCRefit();
    esdTrackCuts->SetAcceptKinkDaughters(0);
    esdTrackCuts->SetMinNCrossedRowsTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetMinDCAToVertexXY(0.06);

    // selections for Xi
    Float_t XiPIDcut = 3.;
    Float_t V0dDCA = 1.6;
    Float_t XidDCA = 1.6;
    Float_t XiMinDCA = 0.07;
    Float_t Xi_massTol = 0.007;
    Float_t Xi_massTolVeto = 0.007;
    Float_t V0CosPoinAn = 0.97;
    Float_t XiCosPoinAn = 0.97;

    AliRsnCutCascade* cutXi = new AliRsnCutCascade(
        "cutXi", kXiMinus, AliPID::kProton, AliPID::kPion, AliPID::kPion);
    cutXi->SetPIDCutV0Proton(XiPIDcut);
    cutXi->SetPIDCutV0Pion(XiPIDcut);
    cutXi->SetPIDCutBachelor(XiPIDcut);
    if (TrackCutsXi == 1)
        cutXi->SetESDtrackCuts(esdTrackCuts);
    cutXi->SetV0MaxDaughtersDCA(V0dDCA);
    cutXi->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXi->SetV0MaxDCAVertex(1e5);  // not using
    cutXi->SetV0MinDCAVertex(XiMinDCA);
    cutXi->SetCascadeMaxDCAVertex(1e5);  // not using
    cutXi->SetCascadeMinDCAVertex(0);    // not using
    cutXi->SetV0LowRadius(0);            // not using
    cutXi->SetV0HighRadius(1e5);         // not using
    cutXi->SetCascadeLowRadius(0);       // not using
    cutXi->SetCascadeHighRadius(1e5);    // not using
    cutXi->SetMassTolerance(Xi_massTol);
    cutXi->SetMassToleranceVeto(
        Xi_massTolVeto);       // Rejection range for Competing Xi Rejection
    cutXi->SetSwitch(kFALSE);  // not using
    cutXi->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXi->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXi->SetMaxRapidity(2.);
    cutXi->SetMinTPCcluster(-1);

    AliRsnCutSet* cutSetXi = new AliRsnCutSet("setXi", AliRsnTarget::kDaughter);
    cutSetXi->AddCut(cutXi);
    cutSetXi->SetCutScheme(cutXi->GetName());
    Int_t icutXi = task->AddTrackCuts(cutSetXi);

    AliRsnCutCascade* cutXibar = new AliRsnCutCascade(
        "cutXibar", kXiPlusBar, AliPID::kProton, AliPID::kPion, AliPID::kPion);
    cutXibar->SetPIDCutV0Proton(XiPIDcut);
    cutXibar->SetPIDCutV0Pion(XiPIDcut);
    cutXibar->SetPIDCutBachelor(XiPIDcut);
    if (TrackCutsXi == 1)
        cutXibar->SetESDtrackCuts(esdTrackCuts);
    cutXibar->SetV0MaxDaughtersDCA(V0dDCA);
    cutXibar->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXibar->SetV0MaxDCAVertex(1e5);  // not using
    cutXibar->SetV0MinDCAVertex(XiMinDCA);
    cutXibar->SetCascadeMaxDCAVertex(1e5);  // not using
    cutXibar->SetCascadeMinDCAVertex(0);    // not using
    cutXibar->SetV0LowRadius(0);            // not using
    cutXibar->SetV0HighRadius(1e5);         // not using
    cutXibar->SetCascadeLowRadius(0);       // not using
    cutXibar->SetCascadeHighRadius(1e5);    // not using
    cutXibar->SetMassTolerance(Xi_massTol);
    cutXibar->SetMassToleranceVeto(
        Xi_massTolVeto);          // Rejection range for Competing Xi Rejection
    cutXibar->SetSwitch(kFALSE);  // not using
    cutXibar->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXibar->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXibar->SetMaxRapidity(2.);
    cutXibar->SetMinTPCcluster(-1);

    AliRsnCutSet* cutSetXibar =
        new AliRsnCutSet("setXibar", AliRsnTarget::kDaughter);
    cutSetXibar->AddCut(cutXibar);
    cutSetXibar->SetCutScheme(cutXibar->GetName());
    Int_t icutXibar = task->AddTrackCuts(cutSetXibar);

    // monitoring
    TString pname = "Xim";
    if (enableMonitor) {
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
#ifdef __CINT__
        gROOT->LoadMacro(
            "$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
#endif
        AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput());

        AddMonitorOutput_P(pname, cutSetXi->GetMonitorOutput());
        AddMonitorOutput_Pt(pname, cutSetXi->GetMonitorOutput());

        pname.Form("Xip");
        AddMonitorOutput_P(pname, cutSetXibar->GetMonitorOutput());
        AddMonitorOutput_Pt(pname, cutSetXibar->GetMonitorOutput());
    }

    // pair cuts
    AliRsnCutMiniPair* cutY =
        new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if (system != 1)
        cutY->SetRangeD(-0.5, 0.5);
    else
        cutY->SetRangeD(-0.465, 0.035);

    AliRsnCutMiniPair* cutV0 =
        new AliRsnCutMiniPair("cutV0", AliRsnCutMiniPair::kContainsV0Daughter);

    AliRsnCutSet* cutsPairSame =
        new AliRsnCutSet("pairCutsSame", AliRsnTarget::kMother);
    cutsPairSame->AddCut(cutY);
    cutsPairSame->AddCut(cutV0);
    cutsPairSame->SetCutScheme(
        TString::Format("%s&(!%s)", cutY->GetName(), cutV0->GetName()).Data());

    AliRsnCutSet* cutsPairMix =
        new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
    cutsPairMix->AddCut(cutY);
    cutsPairMix->SetCutScheme(cutY->GetName());

    // multiplicity binning
    Double_t multbins[200];
    int j, nmult = 0;
    if (!MultBins) {
        multbins[nmult] = 0.;
        nmult++;
        multbins[nmult] = 1.e6;
        nmult++;
    } else if (!trigger) {
        multbins[nmult] = 0.;
        nmult++;
        multbins[nmult] = 1.;
        nmult++;
        multbins[nmult] = 5.;
        nmult++;
        multbins[nmult] = 10.;
        nmult++;
        multbins[nmult] = 15.;
        nmult++;
        for (j = 2; j <= 10; j++) {
            multbins[nmult] = j * 10;
            nmult++;
        }
    } else {
        multbins[nmult] = 0.;
        nmult++;
        multbins[nmult] = 0.001;
        nmult++;
        multbins[nmult] = 0.005;
        nmult++;
        multbins[nmult] = 0.01;
        nmult++;
        multbins[nmult] = 0.05;
        nmult++;
        multbins[nmult] = 0.1;
        nmult++;
        multbins[nmult] = 1.;
        nmult++;
    }

    // -- Values
    // ------------------------------------------------------------------------------------
    /* invariant mass   */ Int_t imID =
        task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
    /* mother mass      */ Int_t mmID =
        task->CreateValue(AliRsnMiniValue::kInvMassMother, kFALSE);
    /* transv. momentum */ Int_t ptID =
        task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
    /* IM difference    */ Int_t diffID =
        task->CreateValue(AliRsnMiniValue::kInvMassDiff, kTRUE);
    /* centrality       */ Int_t centID =
        task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
    /* pseudorapidity   */ Int_t etaID =
        task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
    /* rapidity         */ Int_t yID =
        task->CreateValue(AliRsnMiniValue::kY, kFALSE);
    /* 1st daughter pt  */ Int_t fdpt =
        task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt =
        task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);

    // -- Create all needed outputs
    // -----------------------------------------------------------------

    int i, cutID1, ipdg, xID;
    TString name, comp;
    char charge1, charge2;
    double mass = 1.53178;
    AliRsnMiniOutput* out;

    for (i = 0; i < 16; i++) {
        if (i >= 8 && !isMC)
            continue;

        if (!i) {
            name.Form("XimPip");
            comp.Form("PAIR");
            charge1 = '-';
            charge2 = '+';
            cutID1 = icutXi;
            ipdg = 3324;
            xID = 0;
        } else if (i == 1) {
            name.Form("XipPim");
            comp.Form("PAIR");
            charge1 = '+';
            charge2 = '-';
            cutID1 = icutXibar;
            ipdg = -3324;
            xID = 0;
        } else if (i == 2) {
            name.Form("XimPim");
            comp.Form("PAIR");
            charge1 = '-';
            charge2 = '-';
            cutID1 = icutXi;
            ipdg = 3324;
            xID = 0;
        } else if (i == 3) {
            name.Form("XipPip");
            comp.Form("PAIR");
            charge1 = '+';
            charge2 = '+';
            cutID1 = icutXibar;
            ipdg = -3324;
            xID = 0;
        } else if (i == 4) {
            name.Form("XimPipMix");
            comp.Form("MIX");
            charge1 = '-';
            charge2 = '+';
            cutID1 = icutXi;
            ipdg = 3324;
            xID = 0;
        } else if (i == 5) {
            name.Form("XipPimMix");
            comp.Form("MIX");
            charge1 = '+';
            charge2 = '-';
            cutID1 = icutXibar;
            ipdg = -3324;
            xID = 0;
        } else if (i == 6) {
            name.Form("XimPimMix");
            comp.Form("MIX");
            charge1 = '-';
            charge2 = '-';
            cutID1 = icutXi;
            ipdg = 3324;
            xID = 0;
        } else if (i == 7) {
            name.Form("XipPipMix");
            comp.Form("MIX");
            charge1 = '+';
            charge2 = '+';
            cutID1 = icutXibar;
            ipdg = -3324;
            xID = 0;
        } else if (i == 8) {
            name.Form("Xistar0p_gen");
            comp.Form("MOTHER");
            charge1 = '-';
            charge2 = '+';
            cutID1 = icutXi;
            ipdg = 3324;
            xID = 0;
        } else if (i == 9) {
            name.Form("Xistar0a_gen");
            comp.Form("MOTHER");
            charge1 = '+';
            charge2 = '-';
            cutID1 = icutXibar;
            ipdg = -3324;
            xID = 0;
        } else if (i == 10) {
            name.Form("Xistar0p_true");
            comp.Form("TRUE");
            charge1 = '-';
            charge2 = '+';
            cutID1 = icutXi;
            ipdg = 3324;
            xID = 0;
        } else if (i == 11) {
            name.Form("Xistar0p_trueMM");
            comp.Form("TRUE");
            charge1 = '-';
            charge2 = '+';
            cutID1 = icutXi;
            ipdg = 3324;
            xID = 1;
        } else if (i == 12) {
            name.Form("Xistar0p_res");
            comp.Form("TRUE");
            charge1 = '-';
            charge2 = '+';
            cutID1 = icutXi;
            ipdg = 3324;
            xID = 2;
        } else if (i == 13) {
            name.Form("Xistar0a_true");
            comp.Form("TRUE");
            charge1 = '+';
            charge2 = '-';
            cutID1 = icutXibar;
            ipdg = -3324;
            xID = 0;
        } else if (i == 14) {
            name.Form("Xistar0a_trueMM");
            comp.Form("TRUE");
            charge1 = '+';
            charge2 = '-';
            cutID1 = icutXibar;
            ipdg = -3324;
            xID = 1;
        } else if (i == 15) {
            name.Form("Xistar0a_res");
            comp.Form("TRUE");
            charge1 = '+';
            charge2 = '-';
            cutID1 = icutXibar;
            ipdg = -3324;
            xID = 2;
        }

        out = task->CreateOutput(Form("Xipi_%s%s", name.Data(), suffix), "HIST",
                                 comp.Data());
        out->SetCutID(0, cutID1);
        out->SetCutID(1, iCutPi);
        out->SetDaughter(0, AliRsnDaughter::kXi);
        out->SetDaughter(1, AliRsnDaughter::kPion);
        out->SetCharge(0, charge1);
        out->SetCharge(1, charge2);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        if (comp.EqualTo("PAIR"))
            out->SetPairCuts(cutsPairSame);
        else
            out->SetPairCuts(cutsPairMix);

        // axis X: invmass or resolution
        if (!xID)
            out->AddAxis(imID, 1050, 1.45, 2.5);
        else if (xID == 1)
            out->AddAxis(mmID, 1050, 1.45, 2.5);
        else
            out->AddAxis(diffID, 200, -0.02, 0.02);

        // axis Y: transverse momentum
        out->AddAxis(ptID, 200, 0., 20.);

        // axis Z: centrality-multiplicity
        out->AddAxis(centID, nmult, multbins);
    }

    if (isMC) {  // phase-space histograms
        // Xi(1530)0
        out = task->CreateOutput(Form("Xistar0p_mother_ps%s", suffix), "HIST",
                                 "TRUE");
        out->SetDaughter(0, AliRsnDaughter::kXi);
        out->SetCutID(0, icutXi);
        out->SetCharge(0, '-');

        out->SetDaughter(1, AliRsnDaughter::kPion);
        out->SetCutID(1, iCutPi);
        out->SetCharge(1, '+');

        out->SetMotherPDG(3324);
        out->SetMotherMass(mass);
        out->SetPairCuts(
            cutsPairMix);  // just rapidity, no autocorrelation check
        out->AddAxis(fdpt, 100, 0., 10.);
        out->AddAxis(sdpt, 100, 0., 10.);
        out->AddAxis(ptID, 40, 0., 20.);

        // anti-Xi(1530)0
        out = task->CreateOutput(Form("Xistar0a_mother_ps%s", suffix), "HIST",
                                 "TRUE");
        out->SetDaughter(0, AliRsnDaughter::kXi);
        out->SetCutID(0, icutXibar);
        out->SetCharge(0, '+');

        out->SetDaughter(1, AliRsnDaughter::kPion);
        out->SetCutID(1, iCutPi);
        out->SetCharge(1, '-');

        out->SetMotherPDG(-3324);
        out->SetMotherMass(mass);
        out->SetPairCuts(cutsPairMix);
        out->AddAxis(fdpt, 100, 0., 10.);
        out->AddAxis(sdpt, 100, 0., 10.);
        out->AddAxis(ptID, 40, 0., 20.);
    }

    return kTRUE;
}

void AddMonitorOutput_P(TString name,
                        TObjArray* mon,
                        TString opt,
                        AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_mom", name.Data()), AliRsnValueDaughter::kP);
    a->SetBins(0., 10.0, 0.05);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_Pt(TString name,
                         TObjArray* mon,
                         TString opt,
                         AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(Form("%s_pt", name.Data()),
                                                     AliRsnValueDaughter::kPt);
    a->SetBins(0., 10.0, 0.05);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_Eta(TString name,
                          TObjArray* mon,
                          TString opt,
                          AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_eta", name.Data()), AliRsnValueDaughter::kEta);
    a->SetBins(-2., 2., 0.01);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_DCAxy(TString name,
                            TObjArray* mon,
                            TString opt,
                            AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_dcaxy", name.Data()), AliRsnValueDaughter::kDCAXY);
    a->SetBins(-0.5, 0.5, 0.001);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_DCAz(TString name,
                           TObjArray* mon,
                           TString opt,
                           AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_dcaz", name.Data()), AliRsnValueDaughter::kDCAZ);
    a->SetBins(-2.5, 2.5, 0.005);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0NPt(TString name,
                            TObjArray* mon,
                            TString opt,
                            AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_v0npt", name.Data()), AliRsnValueDaughter::kV0NPt);
    a->SetBins(0., 10.0, 0.05);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0PPt(TString name,
                            TObjArray* mon,
                            TString opt,
                            AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_v0ppt", name.Data()), AliRsnValueDaughter::kV0PPt);
    a->SetBins(0., 10.0, 0.05);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0Mass(TString name,
                             TObjArray* mon,
                             TString opt,
                             AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_v0mass", name.Data()), AliRsnValueDaughter::kV0Mass);
    name.ToLower();
    if (name.Contains("k0"))
        a->SetBins(0.4, 0.6, 0.001);
    else if (name.Contains("lambda"))
        a->SetBins(1.08, 1.16, 0.001);
    else
        a->SetBins(0., 3., 0.01);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0DCA(TString name,
                            TObjArray* mon,
                            TString opt,
                            AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_v0dca", name.Data()), AliRsnValueDaughter::kV0DCA);
    a->SetBins(0.0, 0.4, 0.001);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0Radius(TString name,
                               TObjArray* mon,
                               TString opt,
                               AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_v0radius", name.Data()), AliRsnValueDaughter::kV0Radius);
    a->SetBins(0.0, 200, 0.2);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0Lifetime(TString name,
                                 TObjArray* mon,
                                 TString opt,
                                 AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_v0lifetime", name.Data()), AliRsnValueDaughter::kV0Lifetime);
    a->SetBins(0.0, 200, 0.1);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0DaughterDCA(TString name,
                                    TObjArray* mon,
                                    TString opt,
                                    AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_v0ddca", name.Data()), AliRsnValueDaughter::kDaughterDCA);
    a->SetBins(0.0, 2, 0.001);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0DCA2TPV(TString name,
                                TObjArray* mon,
                                TString opt,
                                AliRsnLoopDaughter* loop) {
    // DCA of secondary tracks to primary vertex
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_v0dca2tpv", name.Data()), AliRsnValueDaughter::kV0DCAXY);
    a->SetBins(-10., 10., 0.01);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0CPA(TString name,
                            TObjArray* mon,
                            TString opt,
                            AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a = new AliRsnValueDaughter(
        Form("%s_v0cpa", name.Data()), AliRsnValueDaughter::kCosPointAng);
    a->SetBins(0.96, 1., 0.001);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0TPCpim(TString name,
                               TObjArray* mon,
                               TString opt,
                               AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a =
        new AliRsnValueDaughter(Form("%s_v0TPCpim", name.Data()),
                                AliRsnValueDaughter::kLambdaPionPIDCut);
    a->SetBins(0., 5., 0.01);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_V0TPCpip(TString name,
                               TObjArray* mon,
                               TString opt,
                               AliRsnLoopDaughter* loop) {
    AliRsnValueDaughter* a =
        new AliRsnValueDaughter(Form("%s_v0TPCpip", name.Data()),
                                AliRsnValueDaughter::kAntiLambdaAntiPionPIDCut);
    a->SetBins(-0., 5., 0.01);
    AliRsnListOutput* o = new AliRsnListOutput(Form("out_%s", a->GetName()),
                                               AliRsnListOutput::kHistoDefault);
    o->AddValue(a);
    if (mon)
        mon->Add(o);
    if (loop)
        loop->AddOutput(o);
}

void AddMonitorOutput_LambdaProtonPID(TObjArray* mon,
                                      TString opt,
                                      AliRsnLoopDaughter* lpPID) {
    // Lambda Cosine of the Pointing Angle
    AliRsnValueDaughter* axisLambdaProtonPID = new AliRsnValueDaughter(
        "lambda_protonPID", AliRsnValueDaughter::kLambdaProtonPIDCut);
    axisLambdaProtonPID->SetBins(0.0, 5, 0.01);

    // output: 2D histogram
    AliRsnListOutput* outMonitorLambdaProtonPID = new AliRsnListOutput(
        "Lambda_ProtonPID", AliRsnListOutput::kHistoDefault);
    outMonitorLambdaProtonPID->AddValue(axisLambdaProtonPID);

    // add outputs to loop
    if (mon)
        mon->Add(outMonitorLambdaProtonPID);
    if (lpPID)
        lpPID->AddOutput(outMonitorLambdaProtonPID);
}

void AddMonitorOutput_LambdaAntiProtonPID(TObjArray* mon,
                                          TString opt,
                                          AliRsnLoopDaughter* lapPID) {
    // Lambda Cosine of the Pointing Angle
    AliRsnValueDaughter* axisLambdaAntiProtonPID = new AliRsnValueDaughter(
        "lambda_antiprotonPID",
        AliRsnValueDaughter::kAntiLambdaAntiProtonPIDCut);
    axisLambdaAntiProtonPID->SetBins(0.0, 5, 0.01);

    // output: 2D histogram
    AliRsnListOutput* outMonitorLambdaAntiProtonPID = new AliRsnListOutput(
        "Lambda_AntiProtonPID", AliRsnListOutput::kHistoDefault);
    outMonitorLambdaAntiProtonPID->AddValue(axisLambdaAntiProtonPID);

    // add outputs to loop
    if (mon)
        mon->Add(outMonitorLambdaAntiProtonPID);
    if (lapPID)
        lapPID->AddOutput(outMonitorLambdaAntiProtonPID);
}
