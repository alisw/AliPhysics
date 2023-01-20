#ifndef __CINT__
#include <vector>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNanoLKr.h"
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0Cuts.h"
#endif

AliAnalysisTaskSE *AddTaskFemtoNanoLKr(
    bool isMC = true,  // 2 MC
    int fFilterBit =
        128,  // 3 type of tracks we select from raw data. TPC only track with 128. With 96 (for sanity check)
    TString triggerData = "kINT7",  // 4 minimum bias  (for default is this) or sample with high multiplicity KHM
    int CutKaon = 0,                /// 7 decide which cut use (Ramona s or Oton s)
    const char *sTcut = "0",        /// 6 for the sphericity cuts. If it s 1, it does the sphericity cuts
    bool DoAncestors = true,        // for common or uncommon ancestors
    const char *cutVariation =
        "0") {  /// 8 to set subwagon, and vary random combination of different cuts (systematic variation)

    TString suffix = TString::Format("%s", cutVariation);  // convert into a string the information about the
                                                           // systematics
    TString sTsuffix = TString::Format("%s", sTcut);  /// convert into a string the information about the sphericity
                                                      /// cuts

    ///
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        printf("No analysis manager found.\n");
        return nullptr;
    }
    if (!mgr->GetInputEventHandler()) {
        printf("This task requires an input event handler!\n");
        return nullptr;
    }

    //========= Init subtasks and start analysis ============================
    // Event Cuts: we select the standard cuts, already used in run2
    AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
    evtCuts->CleanUpMult(false, false, false, true);  /// to select multiplicity estimator

    /// sphericity cuts
    if (sTsuffix == "1") {
        evtCuts->SetSphericityCuts(0.7, 1.0);
    }

    // Lambda Cuts
    AliFemtoDreamv0Cuts *v0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
    AliFemtoDreamTrackCuts *Posv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);  // PileUpRej, false
    AliFemtoDreamTrackCuts *Negv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
    // by default in AliFemtoDreamv0Cuts, the parameter SetCutInvMass is 0.004. It should stay as it is, but we want to
    // do an analysis with all the mass spectrum, so just for this analysis, we ll change this parameter. Afterwards we
    // will put again as it is
    // v0Cuts->SetCutInvMass(3);
    v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
    v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
    v0Cuts->SetPDGCodePosDaug(2212);  // Proton
    v0Cuts->SetPDGCodeNegDaug(211);   // Pion
    v0Cuts->SetPDGCodev0(3122);       // Lambda

    /// antilambda cuts
    AliFemtoDreamv0Cuts *Antiv0Cuts = AliFemtoDreamv0Cuts::LambdaCuts(isMC, true, true);
    AliFemtoDreamTrackCuts *PosAntiv0Daug = AliFemtoDreamTrackCuts::DecayPionCuts(isMC, true, false);
    PosAntiv0Daug->SetCutCharge(1);
    AliFemtoDreamTrackCuts *NegAntiv0Daug = AliFemtoDreamTrackCuts::DecayProtonCuts(isMC, true, false);
    NegAntiv0Daug->SetCutCharge(-1);
    // Antiv0Cuts->SetCutInvMass(3);
    Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
    Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
    Antiv0Cuts->SetPDGCodePosDaug(211);   // Pion
    Antiv0Cuts->SetPDGCodeNegDaug(2212);  // Proton
    Antiv0Cuts->SetPDGCodev0(-3122);      // Anti Lambda (plus or minus?)

    // Track Cuts. PrimKaonCuts is a function already defined
    AliFemtoDreamTrackCuts *TrackPosKaonCuts =
        AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, true);  // DCAplots,CombSigma,ContribSplitting
    TrackPosKaonCuts->SetFilterBit(fFilterBit);
    TrackPosKaonCuts->SetCutCharge(1);  /// positive particle
    TrackPosKaonCuts->SetPtRange(0.15, 1.4);
        
    if (CutKaon == 0) {
        TrackPosKaonCuts->SetPIDkd();  // Oton
    } else if (CutKaon == 1) {
        TrackPosKaonCuts->SetFilterBit(128);
        TrackPosKaonCuts->SetCutCharge(1);
        TrackPosKaonCuts->SetPtRange(0.15, 1.4);
        TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
        TrackPosKaonCuts->SetNClsTPC(70);
        TrackPosKaonCuts->SetDCAVtxZ(1.0); //large DCa, used by Ramona: we usually use 0.1
        TrackPosKaonCuts->SetDCAVtxXY(1.0);
        TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, 3); //only if we ignore the proton/kaon rejection

    }
    else if(CutKaon == 2){
        TrackPosKaonCuts->SetFilterBit(128);
        TrackPosKaonCuts->SetCutCharge(1);
        TrackPosKaonCuts->SetPtRange(0.15, 1.4);
        TrackPosKaonCuts->SetPtExclusion(0.3,0.4);
        TrackPosKaonCuts->SetNClsTPC(70);
        TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, 3); //only if we ignore the proton/kaon rejection
    }

    AliFemtoDreamTrackCuts *TrackNegKaonCuts = AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, true);
    TrackNegKaonCuts->SetFilterBit(fFilterBit);
    TrackNegKaonCuts->SetCutCharge(-1);  /// negative particle
    TrackNegKaonCuts->SetPtRange(0.15, 1.4);
    
    if (CutKaon == 0) {
        TrackNegKaonCuts->SetPIDkd();  // Oton
    } else if (CutKaon == 1) {
        //correct Ramona cuts 
        //TrackNegKaonCuts->SetPIDkd(true, true);  // Ramona
        TrackNegKaonCuts->SetFilterBit(128);
        TrackNegKaonCuts->SetCutCharge(-1);
        TrackNegKaonCuts->SetPtRange(0.15, 1.4);
        TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
        TrackNegKaonCuts->SetNClsTPC(70);
        TrackNegKaonCuts->SetDCAVtxZ(1.0); //large DCa, used by Ramona: we usually use 0.1
        TrackNegKaonCuts->SetDCAVtxXY(1.0);
        TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, 3); //only if we ignore the proton/kaon rejection

    }
    else if(CutKaon == 2){
        TrackNegKaonCuts->SetFilterBit(128);
        TrackNegKaonCuts->SetCutCharge(-1);
        TrackNegKaonCuts->SetPtRange(0.15, 1.4);
        TrackNegKaonCuts->SetPtExclusion(0.3,0.4);
        TrackNegKaonCuts->SetNClsTPC(70);
        TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, 3); //only if we ignore the proton/kaon rejection
    }

    /* if (suffix != "0") {
        evtCuts->SetMinimalBooking(true);
        TrackPosKaonCuts->SetMinimalBooking(false);
        TrackNegKaonCuts->SetMinimalBooking(false);
        v0Cuts->SetMinimalBooking(true);
        Antiv0Cuts->SetMinimalBooking(true);
    } */

    /// specify the particle we re mixing
    // Femto Collection
    std::vector<int> PDGParticles;
    PDGParticles.push_back(321);  // kaon
    PDGParticles.push_back(321);
    PDGParticles.push_back(3122);  // Lambda
    PDGParticles.push_back(-3122);

    std::vector<float> ZVtxBins;  /// distance from vertex
    ZVtxBins.push_back(-10);
    ZVtxBins.push_back(-8);
    ZVtxBins.push_back(-6);
    ZVtxBins.push_back(-4);
    ZVtxBins.push_back(-2);
    ZVtxBins.push_back(0);
    ZVtxBins.push_back(2);
    ZVtxBins.push_back(4);
    ZVtxBins.push_back(6);
    ZVtxBins.push_back(8);
    ZVtxBins.push_back(10);

    /// standard binning of multiplicity
    std::vector<int> MultBins;
    MultBins.push_back(0);
    MultBins.push_back(4);
    MultBins.push_back(8);
    MultBins.push_back(12);
    MultBins.push_back(16);
    MultBins.push_back(20);
    MultBins.push_back(24);
    MultBins.push_back(28);
    MultBins.push_back(32);
    MultBins.push_back(36);
    MultBins.push_back(40);
    MultBins.push_back(44);
    MultBins.push_back(48);
    MultBins.push_back(52);
    MultBins.push_back(56);
    MultBins.push_back(60);
    MultBins.push_back(64);
    MultBins.push_back(68);
    MultBins.push_back(72);
    MultBins.push_back(76);
    MultBins.push_back(80);
    MultBins.push_back(84);
    MultBins.push_back(88);
    MultBins.push_back(92);
    MultBins.push_back(96);
    MultBins.push_back(100);

    /// config->SetZBins(ZVtxBins);
    /// only the ones with same z and multiplicity will be compared

    std::vector<int> NBins;   /// number of bins
    std::vector<float> kMin;  /// minimum value of k*
    std::vector<float> kMax;  /// maximum value of k*
    std::vector<int> pairQA;
    std::vector<bool> closeRejection;

    /// we can have 10 different pairs:
    /// 0 k+ k+
    /// 1 k+ k-
    /// 2 k+ lambda;
    /// 3 k+ antilambda;
    /// 4 k- k-
    /// 5 k- lambda;
    /// 6 k- antilambda;
    /// 7 lambda lambda;
    /// 8 lambda antilambda;
    /// 9 antilambda antilambda

    const int nPairs = 10;
    for (int i = 0; i < nPairs; ++i) {
        pairQA.push_back(0);  /// I initialize all the vectors for the close rejection to zero
        closeRejection.push_back(false);
        NBins.push_back(1500);  /// why 1500
        kMin.push_back(0.);
        kMax.push_back(6.);
    }
    /* if (suffix != "0") {  /// what is Systematic
        pairQA[2] = 12;
        pairQA[3] = 12;
        pairQA[5] = 12;
        pairQA[6] = 12;
        pairQA[7] = 22;  /// the lambdas are just for a cross check
        pairQA[8] = 22;
        pairQA[9] = 22;
 */
    //} else {
        pairQA[0] = 11;
        pairQA[1] = 11;
        pairQA[2] = 12;
        pairQA[3] = 12;
        pairQA[4] = 11;
        pairQA[5] = 12;
        pairQA[6] = 12;
        pairQA[7] = 22;
        pairQA[8] = 22;
        pairQA[9] = 22;

        closeRejection[0] = true;  ///pos kaon - pos kaon
        closeRejection[2] = true;  ///pos kaon- lambda
        closeRejection[3] = true;  ///pos kaon- antilambda
        closeRejection[4] = true;  ///neg kaon-neg kaon
        closeRejection[5] = true;  ///neg kaon-lambda
        closeRejection[6] = true;  ///neg kaon-antilambda
        
    //}

    AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto");

    config->SetExtendedQAPairs(pairQA);  //(commented for systematics)
    config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);
    config->SetMultBinning(true);
    config->SetZBins(ZVtxBins);     /// only the ones with same z and multiplicity will be compared
    config->SetMultBins(MultBins);  // the event will be put in the correspondent multiplicity bin
    config->SetPDGCodes(PDGParticles);
    config->SetNBinsHist(NBins);
    config->SetPhiEtaBinnign(true);
    config->SetDeltaEtaMax(0.012);
    config->SetDeltaPhiMax(0.012);
    config->SetClosePairRejection(closeRejection);
    config->SetMinKRel(kMin);
    config->SetMaxKRel(kMax);
    config->SetUseEventMixing(true);
    config->SetMixingDepth(30);  /// how many events i want to mix. 10 is usually okay
    //config->SetMinimalBookingME(suffix != "0");
    
    /* if (suffix == "0") {
       config->SetPtQA(true);
       config->SetMassQA(true);
       config->SetkTBinning(true);
       config->SetMultBinning(true);
       config->SetmTBinning(true);
    }
 */
    if (isMC) {
        /// create a histogram, which says how e.g. k true (montecarlo s), vs the one we reconstruct
        config->SetMomentumResolution(true);
    } else {
        std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
    }

    // for common ancestors. It is tue only when we run MC
    if (isMC && DoAncestors) {
        config->SetAncestors(true);
        config->GetDoAncestorsPlots();
    }

    /// Setting all cuts for systematics

    // Variation cuts
    const float KaonPtlow = 0.1;
    const float KaonPtup = 0.2;
    const float KaonEtaLow = 0.75;
    const float KaonEtaUp = 0.85;
    const float KaonNClsLow = 70;
    const float KaonNClsUp = 90;
    const float KaonPtMax = 4.0;

    AliPID::EParticleType aliPIDParticle;
    aliPIDParticle = AliPID::kKaon;
    std::map<std::string, float> kaonPIDTight;
    std::map<std::string, float> kaonPIDLoose;

    kaonPIDTight = {
        {"COMB", 2.7},
        {"TPC", 2.7},
        {"EXCLUSION", 3.3},
    };  // for SetPIDkd() when using oton's K selection
    kaonPIDLoose = {
        {"COMB", 3.3},
        {"TPC", 3.3},
        {"EXCLUSION", 2.7},
    };

    if (suffix != "0") {
        if (suffix == "1") {  /// depending on the value of "suffix", we have a different random combination of kaons and lambdas
            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

        } else if (suffix == "2") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

        } else if (suffix == "3") {
            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

        } else if (suffix == "4") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

        } else if (suffix == "5") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

        } else if (suffix == "6") {
            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

        } else if (suffix == "7") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

        } else if (suffix == "8") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

        } else if (suffix == "9") {
            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

        } else if (suffix == "10") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

        } else if (suffix == "11") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

        } else if (suffix == "12") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

        } else if (suffix == "13") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

        } else if (suffix == "14") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

        } else if (suffix == "15") {
            TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

        } else if (suffix == "16") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

        } else if (suffix == "17") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

        } else if (suffix == "18") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

        } else if (suffix == "19") {
            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

        } else if (suffix == "20") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

        } else if (suffix == "21") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

        } else if (suffix == "22") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

        } else if (suffix == "23") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

        } else if (suffix == "24") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

            // XI

        } else if (suffix == "25") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-0.83, 0.83);
            TrackNegKaonCuts->SetEtaRange(-0.83, 0.83);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

        } else if (suffix == "26") {
            TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

        } else if (suffix == "27") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

        } else if (suffix == "28") {
            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);

        } else if (suffix == "29") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
        } else if (suffix == "30") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
        } else if (suffix == "31") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
        } else if (suffix == "32") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
        } else if (suffix == "33") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
        } else if (suffix == "34") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

        } else if (suffix == "35") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
        } else if (suffix == "36") {
            TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
        } else if (suffix == "37") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);

        } else if (suffix == "38") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);
        } else if (suffix == "39") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDTight["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
        } else if (suffix == "40") {
            TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

            TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
            TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
        } else if (suffix == "41") {
            TrackPosKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtlow, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetEtaRange(-0.77, 0.77);
            Negv0Daug->SetEtaRange(-0.77, 0.77);
            PosAntiv0Daug->SetEtaRange(-0.77, 0.77);
            NegAntiv0Daug->SetEtaRange(-0.77, 0.77);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
        } else if (suffix == "42") {
            TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDTight["EXCLUSION"]);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            v0Cuts->SetCutDCADaugTov0Vtx(1.2);
            Antiv0Cuts->SetCutDCADaugTov0Vtx(1.2);
        } else if (suffix == "43") {
            TrackPosKaonCuts->SetPtRange(KaonPtup, KaonPtMax);
            TrackNegKaonCuts->SetPtRange(KaonPtup, KaonPtMax);

            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDLoose["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);
        } else if (suffix == "44") {
            TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
            TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);

            TrackPosKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);
            TrackNegKaonCuts->SetPIDkd(true, false, kaonPIDLoose["COMB"], kaonPIDTight["TPC"],
                                       kaonPIDLoose["EXCLUSION"]);

            v0Cuts->SetCutCPA(0.995);
            Antiv0Cuts->SetCutCPA(0.995);

            Posv0Daug->SetPID(AliPID::kProton, 999.9, 4);
            Negv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            PosAntiv0Daug->SetPID(AliPID::kPion, 999.9, 4);
            NegAntiv0Daug->SetPID(AliPID::kProton, 999.9, 4);

            Posv0Daug->SetNClsTPC(80);
            Negv0Daug->SetNClsTPC(80);
            PosAntiv0Daug->SetNClsTPC(80);
            NegAntiv0Daug->SetNClsTPC(80);

            Posv0Daug->SetEtaRange(-0.83, 0.83);
            Negv0Daug->SetEtaRange(-0.83, 0.83);
            PosAntiv0Daug->SetEtaRange(-0.83, 0.83);
            NegAntiv0Daug->SetEtaRange(-0.83, 0.83);

            v0Cuts->SetCutDCADaugToPrimVtx(0.06);
            Antiv0Cuts->SetCutDCADaugToPrimVtx(0.06);
        }
        else if(suffix == "103"){
            TrackPosKaonCuts->SetPtRange(0.15, 0.5);
            TrackNegKaonCuts->SetPtRange(0.15, 0.5);

        }
        else if(suffix == "104"){
            TrackPosKaonCuts->SetPtRange(0.5, 1.0);
            TrackNegKaonCuts->SetPtRange(0.5, 1.0);

        }       
        else if(suffix == "105"){
            TrackPosKaonCuts->SetPtRange(1.0, 1.4);
            TrackNegKaonCuts->SetPtRange(1.0, 1.4);

        }
           
    }
    AliAnalysisTaskNanoLKr *task = new AliAnalysisTaskNanoLKr("AliAnalysisTaskNanoLKr", isMC);
    /// very important
    if (triggerData == "kINT7") {
        task->SetTrigger(AliVEvent::kINT7);
        task->SelectCollisionCandidates(AliVEvent::kINT7);
        std::cout << "Added kINT7 Trigger \n";
    } else if (triggerData == "kHM") {
        task->SetTrigger(AliVEvent::kHighMultV0);
        task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
        std::cout << "Added kHighMultV0 Trigger \n";
    } else {
        std::cout << "============================================================="
                     "========"
                  << std::endl;
        std::cout << "============================================================="
                     "========"
                  << std::endl;
        std::cout << "Centrality Estimator not set, fix it else your Results will "
                     "be empty!"
                  << std::endl;
        std::cout << "============================================================="
                     "========"
                  << std::endl;
        std::cout << "============================================================="
                     "========"
                  << std::endl;
    }
    /// all the settings go into the taskTrackPosKaonCuts
    task->SetEventCuts(evtCuts);
    task->SetPosKaonCuts(TrackPosKaonCuts);
    task->SetNegKaonCuts(TrackNegKaonCuts);
    task->Setv0Cuts(v0Cuts);
    task->SetAntiv0Cuts(Antiv0Cuts);
    task->SetCorrelationConfig(config);
    mgr->AddTask(task);
    TString file = AliAnalysisManager::GetCommonFileName();

    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

    mgr->ConnectInput(task, 0, cinput);

    TString addon = "";
    if (triggerData == "kINT7") {
        addon += "MB";
    } else if (triggerData == "kHM") {
        addon += "HM";
    }

    mgr->ConnectInput(task, 0, cinput);

    TString QAName = Form("%sQA%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
    AliAnalysisDataContainer *coutputQA = mgr->CreateContainer(
        QAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", file.Data(), QAName.Data()));
    mgr->ConnectOutput(task, 1, coutputQA);

    TString EvtCutsName = Form("%sEvtCuts%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
    AliAnalysisDataContainer *coutputEvtCuts =
        mgr->CreateContainer(EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), EvtCutsName.Data()));
    mgr->ConnectOutput(task, 2, coutputEvtCuts);

    TString TrackCutsName = Form("%sPosKaonCuts%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
    AliAnalysisDataContainer *couputTrkCuts =
        mgr->CreateContainer(TrackCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), TrackCutsName.Data()));
    mgr->ConnectOutput(task, 3, couputTrkCuts);

    TString AntiTrackCutsName = Form("%sNegKaonCuts%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
    AliAnalysisDataContainer *couputAntiTrkCuts =
        mgr->CreateContainer(AntiTrackCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
                             Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
    mgr->ConnectOutput(task, 4, couputAntiTrkCuts);

    AliAnalysisDataContainer *coutputv0Cuts;
    TString v0CutsName = Form("%sLambdaCuts%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
    coutputv0Cuts = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        v0CutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), v0CutsName.Data()));
    mgr->ConnectOutput(task, 5, coutputv0Cuts);

    AliAnalysisDataContainer *coutputAntiv0Cuts;
    TString Antiv0CutsName = Form("%sAntiLambdaCuts%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
    coutputAntiv0Cuts = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        Antiv0CutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
    mgr->ConnectOutput(task, 6, coutputAntiv0Cuts);

    AliAnalysisDataContainer *coutputResults;
    TString ResultsName = Form("%sResults%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
    coutputResults = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        ResultsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), ResultsName.Data()));
    mgr->ConnectOutput(task, 7, coutputResults);

    AliAnalysisDataContainer *coutputResultsQA;
    TString ResultsQAName = Form("%sResultsQA%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
    coutputResultsQA = mgr->CreateContainer(
        //@suppress("Invalid arguments") it works ffs
        ResultsQAName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s:%s", file.Data(), ResultsQAName.Data()));
    mgr->ConnectOutput(task, 8, coutputResultsQA);

    if (isMC) {
        AliAnalysisDataContainer *coutputTrkCutsMC;
        TString TrkCutsMCName = Form("%sTrkCutsMC%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
        coutputTrkCutsMC = mgr->CreateContainer(
            //@suppress("Invalid arguments") it works ffs
            TrkCutsMCName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
            Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
        mgr->ConnectOutput(task, 9, coutputTrkCutsMC);

        AliAnalysisDataContainer *coutputAntiTrkCutsMC;
        TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
        coutputAntiTrkCutsMC = mgr->CreateContainer(
            //@suppress("Invalid arguments") it works ffs
            AntiTrkCutsMCName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
            Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
        mgr->ConnectOutput(task, 10, coutputAntiTrkCutsMC);

        AliAnalysisDataContainer *coutputv0CutsMC;
        TString v0CutsMCName = Form("%sv0CutsMC%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
        coutputv0CutsMC = mgr->CreateContainer(
            //@suppress("Invalid arguments") it works ffs
            v0CutsMCName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
            Form("%s:%s", file.Data(), v0CutsMCName.Data()));
        mgr->ConnectOutput(task, 11, coutputv0CutsMC);

        AliAnalysisDataContainer *coutputAntiv0CutsMC;
        TString Antiv0CutsMCName = Form("%sAntiv0CutsMC%s%s", addon.Data(), suffix.Data(), sTsuffix.Data());
        coutputAntiv0CutsMC = mgr->CreateContainer(
            //@suppress("Invalid arguments") it works ffs
            Antiv0CutsMCName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
            Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
        mgr->ConnectOutput(task, 12, coutputAntiv0CutsMC);
    }
    return task;
}
