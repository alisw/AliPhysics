void InitHistograms(AliDielectron *die, Int_t cutDefinition);
TVectorD *GetVector(Int_t var);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts (AliDielectron *die, Int_t cutDefinition);
void SetTPCSigmaEleCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim);
void SetTOFSigmaEleCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim);
const AliDielectronEventCuts *GetEventCutsMinBias();
const AliDielectronEventCuts *GetEventCutsHighMult();

//
TString names=("pt200_TPCTOFcombITSshared");
TObjArray *arrNames=names.Tokenize(",");
const Int_t nDie=arrNames->GetEntriesFast();

enum {kPhiV = 0, kPt3D, kEta3D, kPhi3D, kMee, kMeeLinear, kMeeForDCAee, kPhiVRebinned, kP2D, kPtee3D, kMee3D, kPairDCA};

Bool_t kPairing  = 1;
Bool_t kMixing   = 1; // kPairing has a higher priority
Bool_t kPairCuts = 0;
Bool_t kTPCCorr  = 0;
Bool_t kTOFCorr  = 0;

Bool_t hasMC;

AliDielectron* Config_ivorobye_pp(Int_t cutDefinition, Bool_t kMinBias = kFALSE)
{
    //
    // Setup the instance of AliDielectron
    //
    
    //=== additional event cuts if needed ============================
    
//    AliDielectronVarCuts *evV0Cut = new AliDielectronVarCuts("evV0Cut","evV0Cut");
//    evV0Cut->AddCut(AliDielectronVarManager::kMultV0,450.,4000.,kFALSE);
//    
//    if (kMinBias == 0){
//        die->GetEventFilter().AddCuts(evV0Cut);
//    }
    
    //=== create the framework object =================================
    TString name=arrNames->At(cutDefinition)->GetName();
    
    AliDielectron *die = new AliDielectron(Form("%s",name.Data()),
                                           Form("Track cuts: %s",name.Data()));
    die->SetHasMC(hasMC);
    
    //=== Switch off pairing for single track analysis ===============
    if (!kPairing) die->SetNoPairing();
    
    //=== cut setup ===================================================
    //--- track cuts
    SetupTrackCuts(die, cutDefinition);
    //--- pair cuts
    if (kPairCuts) SetupPairCuts(die, cutDefinition);
    
    //--- TPCSigmaEle corrections vs P and Eta
    if (kTPCCorr) SetTPCSigmaEleCorrection(die, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
    
    //--- TOFSigmaEle corrections vs P and Eta
    if (kTOFCorr) SetTOFSigmaEleCorrection(die, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
    
    //=== add event mixing handler ====================================
    
    if (kMixing && kPairing){
        
        AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
        mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
        mix->SetMixType(AliDielectronMixingHandler::kAll);
        mix->SetDepth(20);
        mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
        die->SetMixingHandler(mix);
        
    }

    InitHistograms(die, cutDefinition);
    //  InitCF(die,cutDefinition);
    
    //=== Kalman Filter ===============================================
    //--- calculate pairing information not from the KF package,
    //    but directly from the 4-vectors (TLorentzVector)
    //    KF information might be unreliable for some cases, also for event mixing
    die->SetUseKF(kFALSE);
    
    return die;
}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
    //
    // Setup the track cuts
    //
    
    AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
    die->GetTrackFilter().AddCuts(cuts);
    
    //=== ESD tracks cuts =============================================
    
    AliESDtrackCuts* trackCuts=new AliESDtrackCuts();
    
    // pT and eta
    trackCuts->SetPtRange(0.2, 1e30);
    trackCuts->SetEtaRange(-0.8, 0.8);
    
    //TPC
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMinNCrossedRowsTPC(100);
    trackCuts->SetMinNClustersTPC(80);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    trackCuts->SetMaxChi2PerClusterTPC(4.0);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    
    //ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    trackCuts->SetMinNClustersITS(3);
    trackCuts->SetMaxChi2PerClusterITS(4.5);
    
    //primary selection
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetDCAToVertex2D(kFALSE);
    trackCuts->SetRequireSigmaToVertex(kFALSE);
    trackCuts->SetMaxDCAToVertexZ(3.0);
    trackCuts->SetMaxDCAToVertexXY(1.0);
    
    cuts->AddCut(trackCuts);
    
    // ITS shared cluster cut
    
    AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
    varCuts->AddCut(AliDielectronVarManager::kNclsSITS,  -0.1,   0.1);
    cuts->AddCut(varCuts);
    
    
    //=== PID cuts ====================================================
    AliDielectronPID *pidTPCTOFreq = new AliDielectronPID("pidTPCTOFreq","pidTPCTOFreq");
    //TPC
    pidTPCTOFreq->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pidTPCTOFreq->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    // TOF required
    pidTPCTOFreq->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.4, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kP);
    
    AliDielectronPID *pidTPCHadRejTOFif = new AliDielectronPID("pidTPCHadRejTOFif","pidTPCHadRejTOFif");
    //TPC
    pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kKaon,       -4. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kProton,       -4. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    // TOF if available
    pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.4, 1e30,    kFALSE,  AliDielectronPID::kIfAvailable, AliDielectronVarManager::kP);
    
    // combine 2 cut sets with OR option
    AliDielectronCutGroup* combinedPIDcuts = new AliDielectronCutGroup("combinedPIDcuts","combinedPIDcuts",AliDielectronCutGroup::kCompOR);
    combinedPIDcuts->AddCut(pidTPCTOFreq);
    combinedPIDcuts->AddCut(pidTPCHadRejTOFif);
    
    die->GetTrackFilter().AddCuts(combinedPIDcuts);
    
}


//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
    
    AliDielectronVarCuts *PhiV = new AliDielectronVarCuts("PhiV","PhiV");
    PhiV->AddCut(AliDielectronVarManager::kM, 0.0 , 0.05);
    PhiV->AddCut(AliDielectronVarManager::kPhivPair, 2.0 , 3.2 );
    die->GetPairPreFilter().AddCuts(PhiV);
    
    AliDielectronVarCuts *OpAngle = new AliDielectronVarCuts("OpAngle","OpAngle");
    OpAngle->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.050);
    die->GetPairPreFilter().AddCuts(OpAngle);
    
}


//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
    
    //=== Setup histogram Manager =====================================
    AliDielectronHistos *histos= new AliDielectronHistos(die->GetName(),
                                                         die->GetTitle());
    
    //=== Initialise histogram classes ================================
    //    Histograms for classes of Reserved Words are added for all
    //    sub classes
    histos->SetReservedWords("Track;Pair");
    
    //Event class
    histos->AddClass("Event");
    
    //--- Track sub-classes
    //    same event +-
    for (Int_t i=0; i<2; ++i){
        histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
    }
    
    if (kPairing){
        
        //--- Pair sub-classes
        //    same event ++, +-, --
        for (Int_t i=0; i<3; ++i){
            histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
        }
    }
    
    //ME and track rot
    if (kMixing && kPairing) {
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
    }
    
    //=== add histograms to event class ===============================
    
    histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0.,1.,AliDielectronVarManager::kNevents);
    histos->UserHistogram("Event","VtxZ","Vertex Z;Vertex Z [cm];N of events",150,-15.,15.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","NVtxContrib","Number of Vertex Contributor;N of Vertex Contributors;N of events",200,-0.5,199.5,AliDielectronVarManager::kNVtxContrib);
    
    // --- V0 info ---------------------------------------------
    histos->UserHistogram("Event","MultV0","Multiplicity V0;V0M amplitude",4000,-0.5,3999.5,AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","EqMultV0","Equalized Multiplicity V0;Equalized V0M amplitude",4000,-0.5,3999.5,AliDielectronVarManager::kEqMultV0);
    histos->UserHistogram("Event","ChMultV0","Charged Multiplicity V0;Charged V0M amplitude",1000,-0.5,999.5,AliDielectronVarManager::kVZEROchMult);
    histos->UserHistogram("Event","CentralityV0M","Centrality V0;V0M percentile",300,-50,250,AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Event","CentralityV0Mzoomed","Centrality V0 zoomed;V0M percentile",200,0,2,AliDielectronVarManager::kCentralityNew);
    
    //--- N tracks ---------------------------------------------
    histos->UserHistogram("Event","Ntracks","Number of tracks;N of tracks;N of events", 100, -0.5, 99.5,AliDielectronVarManager::kTracks);
    histos->UserHistogram("Event","NtracksVsVtxZ","N tracks vs VtxZ;Vertex Z [cm];N of tracks",150,-15,15,50,-0.5,49.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kTracks);
    histos->UserHistogram("Event","NaccTracks","AliDielectronHelper::GetNacc;N of tracks;N of events", 300, -0.5, 299.5,AliDielectronVarManager::kNacc);
    histos->UserHistogram("Event","NaccItsPureEsd05","ITS SA tracks(AliESDtrackCuts::GetReferenceMultiplicity) in |#eta| < 0.5;N of SPD tracklets in |#eta| < 0.5;N of events", 300, -0.5, 299.5,AliDielectronVarManager::kNaccItsPureEsd05);
    histos->UserHistogram("Event","NaccItsPureEsd10","ITS SA tracks(AliESDtrackCuts::GetReferenceMultiplicity) in |#eta| < 1.0;N of SPD tracklets in |#eta| < 1.0;N of events", 300, -0.5, 299.5,AliDielectronVarManager::kNaccItsPureEsd10);
    histos->UserHistogram("Event","NaccItsTpcEsd05","ITS-TPC tracks(AliESDtrackCuts::GetReferenceMultiplicity) in |#eta| < 0.5;Reference multiplicity in |#eta| < 0.5;N of events", 300, -0.5, 299.5,AliDielectronVarManager::kNaccItsTpcEsd05);
    histos->UserHistogram("Event","NaccItsTpcEsd10","ITS-TPC tracks(AliESDtrackCuts::GetReferenceMultiplicity) in |#eta| < 1.0;Reference multiplicity in |#eta| < 1.0;N of events", 300, -0.5, 299.5,AliDielectronVarManager::kNaccItsTpcEsd10);

    //--- Mult estimators vs vertex Z ---------------------------------------------
    histos->UserHistogram("Event","NaccItsPureEsd05VsVtxZ","NaccItsPureEsd05 vs VtxZ;Vertex Z [cm];N of SPD tracklets in |#eta| < 0.5",110,-11,11,250,-0.5,249.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsPureEsd05);
    histos->UserHistogram("Event","NaccItsPureEsd10VsVtxZ","NaccItsPureEsd10 vs VtxZ;Vertex Z [cm];N of SPD tracklets in |#eta| < 1.0",110,-11,11,250,-0.5,249.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsPureEsd10);
    histos->UserHistogram("Event","NaccItsTpcEsd05VsVtxZ","NaccItsTpcEsd05 vs VtxZ;Vertex Z [cm];Reference multiplicity in |#eta| < 0.5",110,-11,11,250,-0.5,249.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsTpcEsd05);
    histos->UserHistogram("Event","NaccItsTpcEsd10VsVtxZ","NaccItsTpcEsd10 vs VtxZ;Vertex Z [cm];Reference multiplicity in |#eta| < 1.0",110,-11,11,250,-0.5,249.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kNaccItsTpcEsd10);
    
    // --- correlation between V0 mult and ref mult ---------------------------------------------
    histos->UserHistogram("Event","MultV0vsRefMult05","Multiplicity V0 vs Ref.Mult in |#eta|<0.5;V0M amplitude;Reference multiplicity in |#eta| < 0.5",2500,-0.5,2499.5, 200, -0.5, 199.5, AliDielectronVarManager::kMultV0, AliDielectronVarManager::kNaccItsTpcEsd05);
    histos->UserHistogram("Event","MultV0vsSPDMult05","Multiplicity V0 vs SPD Mult in |#eta|<0.5;V0M amplitude;N of SPD tracklets in |#eta| < 0.5",2500,-0.5,2499.5, 200, -0.5, 199.5, AliDielectronVarManager::kMultV0, AliDielectronVarManager::kNaccItsPureEsd05);
    
    //=== add histograms to Track classes =============================
    //    name and title (2nd and 3rd parameter) may be omitted
    //    in that case it is take directly from the var manager
    histos->UserHistogram("Track","Pt","Pt;#it{p}_{T} (GeV/#it{c});N of tracks",500,0,50.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","P","P;#it{p} (GeV/#it{c});N of tracks",500,0.,50.,AliDielectronVarManager::kPIn);
    histos->UserHistogram("Track","P_PIn",";#it{p} (GeV/#it{c});#it{p}^{TPC}_{inner wall} (GeV/#it{c})", GetVector(kP2D), GetVector(kP2D), AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
     histos->UserHistogram("Track","Mass","Mass;ESD mass (GeV/#it{c}^{2});N of tracks",200,0.,2.,AliDielectronVarManager::kM);
    
    // --- Eta and Phi ---------------------------------------------
    histos->UserHistogram("Track","Eta","Eta;#eta;#tracks", 100,-1,1,AliDielectronVarManager::kEta);
    histos->UserHistogram("Track","Phi","Phi;#varphi;#tracks", 72,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    histos->UserHistogram("Track","Eta_Phi","Eta Phi Map;#eta;#varphi", 200,-1,1,360,0,TMath::TwoPi(),AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    histos->UserHistogram("Track","Pt_Phi","Pt Phi Map;#it{p}_{T} (GeV/#it{c});#varphi", 200,0,20,360,0,TMath::TwoPi(),AliDielectronVarManager::kPt,AliDielectronVarManager::kPhi);
    histos->UserHistogram("Track","Pt_eta_phi","Pt Eta Phi;#it{p}_{T} (GeV/#it{c});#eta;#varphi", GetVector(kPt3D),GetVector(kEta3D),GetVector(kPhi3D), AliDielectronVarManager::kPt, AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    
    //--- DCA ---------------------------------------------
    histos->UserHistogram("Track","dXY","dXY;DCA_{xy} (cm);N of tracks", 200,-2.,2.,AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","dZ","dZ;DCA_{z} (cm);N of tracks",  400,-4.,4.,AliDielectronVarManager::kImpactParZ);
    histos->UserHistogram("Track","dXY_dZ","dXY dZ Map;DCA_{xy} (cm);DCA_{z} (cm)", 150,-1.5.,1.5.,400,-4.,4.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
    histos->UserHistogram("Track","phi_dXY","dXY vs phi Map;#varphi;DCA_{xy} (cm)", 144,0.,TMath::TwoPi(), 150,-1.5,1.5, AliDielectronVarManager::kPhi,AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","phi_dZ","dZ vs phi Map;#varphi;DCA_{z} (cm)", 144,0.,TMath::TwoPi(), 400,-4.,4., AliDielectronVarManager::kPhi,AliDielectronVarManager::kImpactParZ);
    
    //--- track checks (TPC) ---------------------------------------------
    
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable;Ratio TPC crossed rows over findable;N of tracks",150,0.,1.5,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRows","Number of Crossed Rows TPC;N of TPC crossed rows;N of tracks",200,-0.5,199.5,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;N of TPC clusters;N of tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","TPCnFCls","Number of Findable Clusters TPC;N of TPC findable clusters;N of tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPC);
    histos->UserHistogram("Track","TPCchi2","TPC Chi2 value;TPC #chi^{2}/N of clusters;N of tracks",50,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
    histos->UserHistogram("Track","NclsSTPC","Number of shared clusters assigned in the TPC;N of TPC shared clusters;N of tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsSTPC);
    histos->UserHistogram("Track","NclsSFracTPC","Fraction of shared clusters assigned in the TPC;Fraction of TPC shared clusters;N of tracks",20,0,1.,AliDielectronVarManager::kNclsSFracTPC);
    histos->UserHistogram("Track","TPCsignalN","Number of PID Clusters TPC;N of TPC PID clusters;N of #tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
    histos->UserHistogram("Track","TPCsignalNfrac","Fraction TPCSignalN/TPCncls;Fraction of TPC PID clusters;N of tracks",60,0.,1.2,AliDielectronVarManager::kTPCsignalNfrac);
    
    //--- track checks (ITS) ---------------------------------------------
    histos->UserHistogram("Track","ITSnCls","Number of Clusters ITS;N of ITS clusters;N of tracks",10,-0.5,9.5,AliDielectronVarManager::kNclsITS);
    histos->UserHistogram("Track","ITSchi2","ITS Chi2 value;ITS #chi^{2}/N of ITS clusters;N of tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
    histos->UserHistogram("Track","ITSSharedClusters","N of ITS shared clusters;N of shared clusters ITS;N of tracks",20,-10.,10.,AliDielectronVarManager::kNclsSITS);
    
    //--- PID information ---------------------------------------------
    //-----------------------------------------------------------------
    //--- ITS ---------------------------------------------------------
//    histos->UserHistogram("Track","ITS_dEdx_P","ITS dEdx;#it{p} (GeV/#it{c});ITS signal (arb units)",
//                          1000,0.,50.,700,0.,700.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
    histos->UserHistogram("Track","ITSnSigmaEle_P","ITS number of sigmas Electrons;#it{p} (GeV/#it{c});ITS number of sigmas Electrons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
//    histos->UserHistogram("Track","ITSnSigmaPio_P","ITS number of sigmas Pions;#it{p} (GeV/#it{c});ITS number of sigmas Pions",
//                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
//    histos->UserHistogram("Track","ITSnSigmaKao_P","ITS number of sigmas Kaons;#it{p} (GeV/#it{c});ITS number of sigmas Kaons",
//                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
//    histos->UserHistogram("Track","ITSnSigmaPro_P","ITS number of sigmas Protons;#it{p} (GeV/#it{c});ITS number of sigmas Protons",
//                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
//    histos->UserHistogram("Track","ITSnSigmaEle_Eta_P_lin","ITS number of sigmas Electrons vs Eta and P;#eta;n#sigma_{ele}^{ITS};#it{p} (GeV/#it{c})",
//                          200, -1., 1., 200, -10., 10., 100, 0., 5.,
//                          AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kP);
    
    //--- TPC ----------------------------------------------------------
//    histos->UserHistogram("Track","TPC_dEdx_P","TPC dEdx;#it{p} (GeV/#it{c});TPC signal (arb units)",
//                          1000,0.,50.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;#it{p} (GeV/#it{c});TPC number of sigmas Electrons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
//    histos->UserHistogram("Track","TPCnSigmaEleRaw_P","raw TPC number of sigmas Electrons;#it{p} (GeV/#it{c});raw TPC number of sigmas Electrons",
//                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEleRaw);
    histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;#it{p} (GeV/#it{c});TPC number of sigmas Pions",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;#it{p} (GeV/#it{c});TPC number of sigmas Kaons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    histos->UserHistogram("Track","TPCnSigmaPro_P","TPC number of sigmas Protons;#it{p} (GeV/#it{c});TPC number of sigmas Protons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_P_lin",";Eta;n#sigma_{ele}^{TPC};#it{p}(GeV/#it{c})",
                          200, -1., 1., 200, -10., 10., 100, 0., 5.,
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kP);
//    histos->UserHistogram("Track","TPCnSigmaEleRaw_Eta_P_lin",";Eta;n#sigma_{ele}^{TPC};#it{p}(GeV/#it{c})",
//                          200, -1., 1., 200, -10., 10., 100, 0., 5.,
//                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw,AliDielectronVarManager::kP);

    //--- TOF ----------------------------------------------------------
//    histos->UserHistogram("Track","TOFbeta","TOF beta;#it{p} (GeV/#it{c});TOF #beta",
//                          1000,0.,50.,300,0.,3.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
//    histos->UserHistogram("Track","TOFbeta_Pt","TOF beta;#it{p}_{T} (GeV/#it{c});TOF #beta",
//                          1000,0.,50.,300,0.,3.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTOFbeta);
    histos->UserHistogram("Track","TOFnSigmaEle_P","TOF number of sigmas Electrons;#it{p} (GeV/#it{c});TOF number of sigmas Electrons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
//    histos->UserHistogram("Track","TOFnSigmaEleRaw_P","raw TOF number of sigmas Electrons;#it{p} (GeV/#it{c});raw TOF number of sigmas Electrons",
//                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEleRaw);
//    histos->UserHistogram("Track","TOFnSigmaEle_Pt","TOF number of sigmas Electrons;#it{p} (GeV/#it{c});TOF number of sigmas Electrons",
//                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTOFnSigmaEle);
    histos->UserHistogram("Track","TOFnSigmaPio_P","TOF number of sigmas Pions;#it{p} (GeV/#it{c});TOF number of sigmas Pions",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);
    histos->UserHistogram("Track","TOFnSigmaKao_P","TOF number of sigmas Kaons;#it{p} (GeV/#it{c});TOF number of sigmas Kaons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao);
    histos->UserHistogram("Track","TOFnSigmaPro_P","TOF number of sigmas Protons;#it{p} (GeV/#it{c});TOF number of sigmas Protons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro);
    histos->UserHistogram("Track","TOFnSigmaEle_Eta_P_lin","TOF number of sigmas Electrons vs Eta and P;Eta;n#sigma_{ele}^{TOF};#it{p} (GeV/#it{c})",
                          200, -1., 1., 200, -10., 10., 100, 0., 5.,
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle,AliDielectronVarManager::kP);
//    histos->UserHistogram("Track","TOFnSigmaEleRaw_Eta_P_lin","raw TOF number of sigmas Electrons vs Eta and P;Eta;n#sigma_{ele}^{TOF};#it{p} (GeV/#it{c})",
//                          200, -1., 1., 200, -10., 10., 100, 0., 5.,
//                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEleRaw,AliDielectronVarManager::kP);
    
    //=== add histograms to Pair classes ==============================
    
    if (kPairing){
        
        histos->UserHistogram("Pair","InvMass","Inv.Mass;#it{m}_{ee} (GeV/#it{c}^{2});N of pairs", 400,0.,4.,AliDielectronVarManager::kM);
        histos->UserHistogram("Pair","InvMass_low","Inv.Mass;#it{m}_{ee} (GeV/#it{c}^{2});N of pairs",500,0,0.5,AliDielectronVarManager::kM);
        histos->UserHistogram("Pair","PairPt","PairPt;Pair #it{p}_{T} (GeV/#it{c});N of pairs", 160,0.,8.,AliDielectronVarManager::kPt);
        histos->UserHistogram("Pair","Rapidity","Rapidity;#eta;N of pairs", 100,-1.,1.,AliDielectronVarManager::kY);
        histos->UserHistogram("Pair","OpeningAngle","Opening angle;Opening angle;n of pairs", 180,0.,TMath::Pi(),AliDielectronVarManager::kOpeningAngle);
        histos->UserHistogram("Pair","Phi_Pair","Phi;#varphi;N of pairs", 360, 0. , TMath::TwoPi(), AliDielectronVarManager::kPhi);
        histos->UserHistogram("Pair","PhiV_Pair","PhiV_Pair;#phi_{V};N of pairs", 180,0.,TMath::Pi(),AliDielectronVarManager::kPhivPair);
        histos->UserHistogram("Pair","kDeltaEta","kDeltaEta;#delta_{#eta};N of pairs",160,0.,1.6,AliDielectronVarManager::kDeltaEta);
        histos->UserHistogram("Pair","kDeltaPhi","kDeltaPhi;#delta_{#varphi};N of pairs",360,0.,TMath::TwoPi(),AliDielectronVarManager::kDeltaPhi);

//        // --- 2d ----------------------------------------------------------
//        histos->UserHistogram("Pair","InvMass_OpAngle","InvMass_OpAngle;#it{m}_{ee} (GeV/#it{c}^{2});Opening angle", 200, 0. , 4. , 180, 0. , TMath::Pi() ,AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
//        histos->UserHistogram("Pair","OpAngle_PhiV","OpAngle_PhiV;Opening angle;#phi_{V}",180, 0. , TMath::Pi(), 180 , 0. , TMath::Pi() ,AliDielectronVarManager::kOpeningAngle,AliDielectronVarManager::kPhivPair);
        
        // --- mass vs phiV ----------------------------------------------------------
        // flexible vs mass and phiV
        histos->UserHistogram("Pair","InvMass_PhivPair","InvMass:PhivPair;#it{m}_{ee} (GeV/#it{c}^{2});#phi_{V}",
                              400, 0., 4., 30 , 0. , TMath::Pi(),
                              AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
        
//        histos->UserHistogram("Pair","InvMassLow_PhivPair","InvMass:PhivPair;#it{m}_{ee} (GeV/#it{c}^{2});#phi_{V}",
//                              400, 0., 0.4, 30 , 0. , TMath::Pi(),
//                              AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
//
        
        // dark photons for Taku
        histos->UserHistogram("Pair","InvMassUltraLow_PhivPair","InvMass:PhivPair;#it{m}_{ee} (GeV/#it{c}^{2});#phi_{V}",
                              1500, 0., 0.15, 30 , 0. , TMath::Pi(),
                              AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
        
        // --- mass vs pT vs phiV ----------------------------------------------------------
        // mass 10 MeV, pT 1 GeV
        histos->UserHistogram("Pair","InvMass_PtRebinned_PhivPair","InvMass:Pt:PhivPair;#it{m}_{ee} (GeV/#it{c}^{2});Pair #it{p}_{T} (GeV/#it{c});#phi_{V};",GetVector(kMeeLinear),GetVector(kPtee3D),GetVector(kPhiV),AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
        
        // mass very coarse, fine pT
        histos->UserHistogram("Pair","InvMassRebinned_Pt_PhivPair","InvMass:Pt:PhivPair;#it{m}_{ee} (GeV/#it{c}^{2});Pair #it{p}_{T} (GeV/#it{c});#phi_{V};",GetVector(kMee3D),GetVector(kPt3D),GetVector(kPhiV),AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
        
        // --- mass vs pair DCA vs phiV ----------------------------------------------------------
        // mass 10 MeV, pT 1 GeV
        histos->UserHistogram("Pair","InvMass_PairDCAxy_PhivPair","InvMass:PairDCAxy:PhivPair;#it{m}_{ee} (GeV/#it{c}^{2});Pair DCA_{xy};#phi_{V};",GetVector(kMeeForDCAee),GetVector(kPairDCA),GetVector(kPhiV),AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPhivPair);
        
//        histos->UserHistogram("Pair","InvMass_PairDCAz_PhivPair","InvMass:PairDCAz:PhivPair;#it{m}_{ee} (GeV/#it{c}^{2});Pair DCA_{z};#phi_{V};",GetVector(kMeeForDCAee),GetVector(kPairDCA),GetVector(kPhiV),AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigZ, AliDielectronVarManager::kPhivPair);
    }
    
    //=== add histograms relevant for mixing ============================
    AliDielectronMixingHandler *mix=die->GetMixingHandler();
    if (kMixing && kPairing){
        const Int_t nMix=mix->GetNumberOfBins();
        histos->UserHistogram("Event","","", nMix,0.,(Double_t)nMix,AliDielectronVarManager::kMixingBin);
        histos->UserHistogram("Pair","","", nMix, 0, nMix, AliDielectronVarManager::kMixingBin);
    }
    
    die->SetHistogramManager(histos);
    
}

//______________________________________________________________________________________
TVectorD *GetVector(Int_t var)
{
    
    switch (var)
    {
        case kPhiV:   return AliDielectronHelper::MakeLinBinning(30, 0., TMath::Pi());
        case kPt3D:   return AliDielectronHelper::MakeArbitraryBinning(" 0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.30, 2.50, 3.00, 3.50, 4.00, 5.0, 6.0, 7.0, 10.0, 20.0");
        case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);
        case kPhi3D:  return AliDielectronHelper::MakeLinBinning( 60,0,6.2832); // same as in EffTask
        case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.32, 0.38, 0.54, 0.66, 0.72, 0.76, 0.78, 0.86, 0.98, 1.02, 1.2, 1.4, 1.64, 1.96, 2.34, 2.84, 2.95, 3.05, 3.1, 3.15, 3.3, 3.5, 3.75, 4.0");
        case kMeeLinear:      return AliDielectronHelper::MakeLinBinning(400, 0., 4.);
        case kMeeForDCAee:    return AliDielectronHelper::MakeArbitraryBinning("0., 0.08, 0.14, 0.2, 1.1, 2.7, 2.8, 3.2, 5.0");
        case kPhiVRebinned:   return AliDielectronHelper::MakeArbitraryBinning("0.0, 2.0, 3.2");
        case kP2D:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 9.00, 10.00");
        case kPtee3D:   return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0, 100.0");
        case kMee3D:   return AliDielectronHelper::MakeArbitraryBinning("0., 0.15, 0.75, 1.25, 2.75, 3.5, 4.0");
        case kPairDCA: return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 4.0, 5.0, 7.0, 10.0");
            
        default: cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << endl;
            break;
    }
    
}

//______________________________________________________________________________________
void SetTPCSigmaEleCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim) {
    //
    // eta correction for the centroid and width of electron sigmas in the TPC
    //
    
    //
    printf("starting SetTPCSigmaEleCorrection()\n");
    printf(" correction Xdim = %s\n", AliDielectronVarManager::GetValueName(corrXdim));
    printf(" correction Ydim = %s\n", AliDielectronVarManager::GetValueName(corrYdim));
    
    if (corrYdim!=AliDielectronVarManager::kEta) { printf(" correction only available for Ydim = eta!\n"); return; }
    
    if (corrXdim!=AliDielectronVarManager::kP)
    {
        printf(" no correction available for Xdim = %s!\n", AliDielectronVarManager::GetValueName(corrXdim));
        printf(" no correction applied!\n");
        return;
    }
    
    // die->SetCentroidCorrFunction(histMean2D, corrXdim, corrYdim);
    // die->SetWidthCorrFunction(histWidth2D, corrXdim, corrYdim);
    printf("no default TPC PID correction! Corrections were not applied!\n");
    
}

//______________________________________________________________________________________
void SetTOFSigmaEleCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim) {
    //
    // eta correction for the centroid and width of electron sigmas in the TOF
    //
    
    //
    printf("starting SetTOFSigmaEleCorrection()\n");
    printf(" correction Xdim = %s\n", AliDielectronVarManager::GetValueName(corrXdim));
    printf(" correction Ydim = %s\n", AliDielectronVarManager::GetValueName(corrYdim));
    
    if (corrYdim!=AliDielectronVarManager::kEta) { printf(" correction only available for Ydim = eta!\n"); return; }
    
    if (corrXdim!=AliDielectronVarManager::kP)
    {
        printf(" no correction available for Xdim = %s!\n", AliDielectronVarManager::GetValueName(corrXdim));
        printf(" no correction applied!\n");
        return;
    }
    
    // die->SetCentroidCorrFunctionTOF(histMean2DTOF, corrXdim, corrYdim);
    // die->SetWidthCorrFunctionTOF(histWidth2DTOF, corrXdim, corrYdim);
    printf("no default TOF PID correction! Corrections were not applied!\n");
    
}

//______________________________________________________________________________________
const AliDielectronEventCuts *GetEventCutsMinBias(){
    
    AliDielectronEventCuts *eventCutsMB=new AliDielectronEventCuts("eventCutsMB","Vertex Track && |vtxZ|<10 && ncontrib>0");
    eventCutsMB->SetVertexType(AliDielectronEventCuts::kVtxAny);
    eventCutsMB->SetRequireVertex();
    eventCutsMB->SetMinVtxContributors(1);
    eventCutsMB->SetVertexZ(-10.,10.);
    eventCutsMB->SetMaxSPDVertexResolution(0.25);
    eventCutsMB->SetMaxSPDVertexDispersion(0.03);
    eventCutsMB->SetMaxVertexDisplacement(0.5);
//    eventCutsMB->SetCentralityRange(0., 1., kTRUE); // optional centrality selection for pp run2
    return eventCutsMB;
}

//______________________________________________________________________________________
const AliDielectronEventCuts *GetEventCutsHighMult(){
    
    AliDielectronEventCuts *eventCutsHM=new AliDielectronEventCuts("eventCutsHM","Vertex Track && |vtxZ|<10 && ncontrib>0");
    eventCutsHM->SetVertexType(AliDielectronEventCuts::kVtxAny);
    eventCutsHM->SetRequireVertex();
    eventCutsHM->SetMinVtxContributors(1);
    eventCutsHM->SetVertexZ(-10.,10.);
    eventCutsHM->SetMaxSPDVertexResolution(0.25);
    eventCutsHM->SetMaxSPDVertexDispersion(0.03);
    eventCutsHM->SetMaxVertexDisplacement(0.5);
    eventCutsHM->SetCentralityRange(0., 0.05, kTRUE); // optional centrality selection for pp run2
    return eventCutsHM;
}
