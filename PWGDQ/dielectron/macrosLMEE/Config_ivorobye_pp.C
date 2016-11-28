void InitHistograms(AliDielectron *die, Int_t cutDefinition);
TVectorD *GetVector(Int_t var);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts (AliDielectron *die, Int_t cutDefinition);
void SetTPCSigmaEleCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim);
const AliDielectronEventCuts *GetEventCutsMinBias();
const AliDielectronEventCuts *GetEventCutsHighMult();

//
TString names=("pt200_TPCTOFreq");
TObjArray *arrNames=names.Tokenize(",");
const Int_t nDie=arrNames->GetEntriesFast();

enum {kPhiV = 0, kPt3D, kEta3D, kPhi3D, kMee, kPhiVRebinned, kP2D};

Bool_t kPairing  = 1;
Bool_t kMixing   = 1; // kPairing has a higher priority
Bool_t kPairCuts = 0;
Bool_t kTPCCorr  = 0;

Bool_t hasMC;

AliDielectron* Config_ivorobye_pp(Int_t cutDefinition)
{
    //
    // Setup the instance of AliDielectron
    //
    
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
    
    // TPC
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMinNCrossedRowsTPC(100);
    trackCuts->SetMinNClustersTPC(80);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    trackCuts->SetMaxChi2PerClusterTPC(4.);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // trackCuts->SetMaxChi2TPCConstrainedGlobal(36.);
    
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    trackCuts->SetMinNClustersITS(3);
    trackCuts->SetMaxChi2PerClusterITS(4.5);
    
    //
    // primary selection
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetDCAToVertex2D(kFALSE);
    trackCuts->SetRequireSigmaToVertex(kFALSE);
    trackCuts->SetMaxDCAToVertexZ(3.0);
    trackCuts->SetMaxDCAToVertexXY(1.0);
    cuts->AddCut(trackCuts);
    
    
    //=== PID cuts ====================================================
    AliDielectronPID *pid = new AliDielectronPID();
    
    //            Detector                 Particle            nSigma      pt Range    exclusion?   require?
    //                                                        min, max     min, max
    
    //TPC
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
    
    // TOF
    pid->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.4, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kP);
    
    cuts->AddCut(pid);
    
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
    histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",150,-15.,15.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","NVtxContrib","Number of Vertex Contributor;NVtx;Number events",200,-0.5,199.5,AliDielectronVarManager::kNVtxContrib);
    
    // --- V0 info ---------------------------------------------
    histos->UserHistogram("Event","MultV0","Multiplicity V0;V0 multiplicity",4000,-0.5,3999.5,AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","EqMultV0","Equalized Multiplicity V0;V0 multiplicity",4000,-0.5,3999.5,AliDielectronVarManager::kEqMultV0);
    histos->UserHistogram("Event","ChMultV0","Charged Multiplicity V0;V0 multiplicity",1000,-0.5,999.5,AliDielectronVarManager::kVZEROchMult);
    histos->UserHistogram("Event","CentralityV0M","Centrality V0;V0 centrality",300,-50,250,AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Event","CentralityV0Mzoomed","Centrality V0 zoomed;V0 centrality",200,0,2,AliDielectronVarManager::kCentralityNew);
    
    //--- N tracks ---------------------------------------------
    histos->UserHistogram("Event","Ntracks","Number of tracks;Ntracks;#events", 100, -0.5, 99.5,AliDielectronVarManager::kTracks);
    histos->UserHistogram("Event","NtracksVsVtxZ","N tracks vs VtxZ;VtxZ;N tracks",150,-15,15,50,-0.5,49.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kTracks);
    histos->UserHistogram("Event","NaccTracks","AliDielectronHelper::GetNacc;Ntracks;#events", 300, -0.5, 299.5,AliDielectronVarManager::kNacc);
    histos->UserHistogram("Event","NaccItsPureEsd05","ITS SA tracks(AliESDtrackCuts::GetReferenceMultiplicity) in eta 0.5;NItsTracks;#events", 300, -0.5, 299.5,AliDielectronVarManager::kNaccItsPureEsd05);
    histos->UserHistogram("Event","NaccItsPureEsd10","ITS SA tracks(AliESDtrackCuts::GetReferenceMultiplicity) in eta 1.0;NItsTracks;#events", 300, -0.5, 299.5,AliDielectronVarManager::kNaccItsPureEsd10);
    
    //=== add histograms to Track classes =============================
    //    name and title (2nd and 3rd parameter) may be omitted
    //    in that case it is take directly from the var manager
    histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",500,0,50.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","P","P;P [GeV];#tracks",500,0.,50.,AliDielectronVarManager::kPIn);
    histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)", GetVector(kP2D), GetVector(kP2D), AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
    histos->UserHistogram("Track","Mass","Mass;mass [GeV];#tracks",200,0.,2.,AliDielectronVarManager::kM);
    
    // --- Eta and Phi ---------------------------------------------
    histos->UserHistogram("Track","Eta","Eta;Eta;#tracks", 100,-1,1,AliDielectronVarManager::kEta);
    histos->UserHistogram("Track","Phi","Phi;Phi;#tracks", 72,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    histos->UserHistogram("Track","Eta_Phi","Eta Phi Map;Eta;Phi", 200,-1,1,360,0,TMath::TwoPi(),AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    histos->UserHistogram("Track","Pt_Phi","Pt Phi Map;Pt;Phi", 200,0,20,360,0,TMath::TwoPi(),AliDielectronVarManager::kPt,AliDielectronVarManager::kPhi);
    histos->UserHistogram("Track","Pt_eta_phi","Pt Eta Phi;Pt;Eta;Phi", GetVector(kPt3D),GetVector(kEta3D),GetVector(kPhi3D), AliDielectronVarManager::kPt, AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    
    //--- DCA ---------------------------------------------
    histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks", 200,-2.,2.,AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",  400,-4.,4.,AliDielectronVarManager::kImpactParZ);
    histos->UserHistogram("Track","dXY_dZ","dXY dZ Map;dXY;dZ", 150,-1.5.,1.5.,400,-4.,4.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
    histos->UserHistogram("Track","phi_dXY","dXY vs phi Map;phi;dXY", 144,0.,TMath::TwoPi(), 150,-1.5,1.5, AliDielectronVarManager::kPhi,AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","phi_dZ","dZ vs phi Map;phi;dZ", 144,0.,TMath::TwoPi(), 400,-4.,4., AliDielectronVarManager::kPhi,AliDielectronVarManager::kImpactParZ);
    
    //--- track checks (TPC) ---------------------------------------------
    
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable;TPC crossed rows over findable;#tracks",150,0.,1.5,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRows","Number of Crossed Rows TPC;TPC crossed rows;#tracks",200,-0.5,199.5,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","TPCnFCls","Number of Findable Clusters TPC;TPC number of findable clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPC);
    histos->UserHistogram("Track","TPCchi2","TPC Chi2 value;TPC chi2/Cl;#tracks",50,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
    histos->UserHistogram("Track","NclsSTPC","Number of shared clusters assigned in the TPC;TPC number of shared clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsSTPC);
    histos->UserHistogram("Track","NclsSFracTPC","Fraction of shared clusters assigned in the TPC;TPC fraction of shared clusters;#tracks",20,0,1.,AliDielectronVarManager::kNclsSFracTPC);
    histos->UserHistogram("Track","TPCsignalN","Number of PID Clusters TPC;TPC number PID clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
    histos->UserHistogram("Track","TPCsignalNfrac","Fraction TPCSignalN/TPCncls;TPCSignalN/TPCNcls;#tracks",60,0.,1.2,AliDielectronVarManager::kTPCsignalNfrac);
    
    //--- track checks (ITS) ---------------------------------------------
    histos->UserHistogram("Track","ITSnCls","Number of Clusters ITS;ITS number clusters;#tracks",10,-0.5,9.5,AliDielectronVarManager::kNclsITS);
    histos->UserHistogram("Track","ITSchi2","ITS Chi2 value;ITS chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
    
    //--- PID information ---------------------------------------------
    //--- ITS ---------------------------------------------------------
    histos->UserHistogram("Track","ITS_dEdx_P","ITS dEdx;P [GeV];ITS signal (arb units)",
                          1000,0.,50.,700,0.,700.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
    histos->UserHistogram("Track","ITSnSigmaEle_P","ITS number of sigmas Electrons;P [GeV];ITS number of sigmas Electrons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
    histos->UserHistogram("Track","ITSnSigmaPio_P","ITS number of sigmas Pions;P [GeV];ITS number of sigmas Pions",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
    histos->UserHistogram("Track","ITSnSigmaKao_P","ITS number of sigmas Kaons;P [GeV];ITS number of sigmas Kaons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
    histos->UserHistogram("Track","ITSnSigmaPro_P","ITS number of sigmas Protons;P [GeV];ITS number of sigmas Protons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
    histos->UserHistogram("Track","ITSnSigmaEle_Eta_P_lin","ITS number of sigmas Electrons vs Eta and P;Eta;n#sigma_{ele}^{ITS};p (GeV/c)",
                          200, -1., 1., 200, -10., 10., 100, 0., 5.,
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kP);
    
    
    //--- TPC ----------------------------------------------------------
    histos->UserHistogram("Track","TPC_dEdx_P","TPC dEdx;P [GeV];TPC signal (arb units)",
                          1000,0.,50.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    histos->UserHistogram("Track","TPCnSigmaPro_P","TPC number of sigmas Protons;P [GeV];TPC number of sigmas Protons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_P_lin",";Eta;n#sigma_{ele}^{TPC};p (GeV/c)",
                          200, -1., 1., 200, -10., 10., 100, 0., 5.,
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kP);

    //--- TOF ----------------------------------------------------------
    histos->UserHistogram("Track","TOFbeta","TOF beta;P [GeV];TOF beta",
                          1000,0.,50.,300,0.,3.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
    histos->UserHistogram("Track","TOFbeta_Pt","TOF beta;Pt [GeV];TOF beta",
                          1000,0.,50.,300,0.,3.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTOFbeta);
    histos->UserHistogram("Track","TOFnSigmaEle_P","TOF number of sigmas Electrons;P [GeV];TOF number of sigmas Electrons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
    histos->UserHistogram("Track","TOFnSigmaEle_Pt","TOF number of sigmas Electrons;Pt [GeV];TOF number of sigmas Electrons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTOFnSigmaEle);
    histos->UserHistogram("Track","TOFnSigmaPio_P","TOF number of sigmas Pions;P [GeV];TOF number of sigmas Pions",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);
    histos->UserHistogram("Track","TOFnSigmaKao_P","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao);
    histos->UserHistogram("Track","TOFnSigmaPro_P","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons",
                          1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro);
    histos->UserHistogram("Track","TOFnSigmaEle_Eta_P_lin","TOF number of sigmas Electrons vs Eta and P;Eta;n#sigma_{ele}^{TOF};p (GeV/c)",
                          200, -1., 1., 200, -10., 10., 100, 0., 5.,
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle,AliDielectronVarManager::kP);
    
    //=== add histograms to Pair classes ==============================
    
    if (kPairing){
        
        histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs", 250,0.,5.,AliDielectronVarManager::kM);
        histos->UserHistogram("Pair","InvMass_low","Inv.Mass;Inv. Mass [GeV];#pairs",250,0,0.5,AliDielectronVarManager::kM);
        histos->UserHistogram("Pair","PairPt","PairPt;Pair Pt [GeV];#pairs", 160,0.,8.,AliDielectronVarManager::kPt);
        histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs", 100,-1.,1.,AliDielectronVarManager::kY);
        histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle;#pairs", 180,0.,TMath::Pi(),AliDielectronVarManager::kOpeningAngle);
        histos->UserHistogram("Pair","Phi_Pair","Phi;counts;Phi", 360, 0. , TMath::TwoPi(), AliDielectronVarManager::kPhi);
        histos->UserHistogram("Pair","PhiV_Pair","PhiV_Pair;PhiV;#pairs", 180,0.,TMath::Pi(),AliDielectronVarManager::kPhivPair);
        histos->UserHistogram("Pair","kDeltaEta","kDeltaEta;kDeltaEta;#pairs",160,0.,1.6,AliDielectronVarManager::kDeltaEta);
        histos->UserHistogram("Pair","kDeltaPhi","kDeltaPhi;kDeltaPhi;#pairs",360,0.,TMath::TwoPi(),AliDielectronVarManager::kDeltaPhi);
        
        // --- 2d ----------------------------------------------------------
        histos->UserHistogram("Pair","InvMass_Pt","InvMass_Pt;InvMass;Pt",250, 0. , 5., 100 , 0., 5. , AliDielectronVarManager::kM , AliDielectronVarManager::kPt );
        histos->UserHistogram("Pair","OpAngle_InvMass","OpAngle_InvMass;Opening angle;Invariant Mass",180, 0. , TMath::Pi(), 250 , 0. , 5. ,AliDielectronVarManager::kOpeningAngle,AliDielectronVarManager::kM);
        histos->UserHistogram("Pair","OpAngle_PhiV","OpAngle_PhiV;Opening angle;PhiV",180, 0. , TMath::Pi(), 180 , 0. , TMath::Pi() ,AliDielectronVarManager::kOpeningAngle,AliDielectronVarManager::kPhivPair);
        
        // --- mass vs phiV ----------------------------------------------------------
        histos->UserHistogram("Pair","InvMass_PhivPair","InvMass:PhivPair;Inv. Mass [GeV];PhiV",
                              250, 0.,5., 32, 0., 3.2,
                              AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
        histos->UserHistogram("Pair","InvMassRebinned_PhivPair","InvMass:PhivPair;Inv. Mass [GeV];PhiV",GetVector(kMee),GetVector(kPhiV),AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
        histos->UserHistogram("Pair","InvMassRebinned_PhivPairRebinned","InvMass:PhivPair;Inv. Mass [GeV];PhiV",GetVector(kMee),GetVector(kPhiVRebinned),AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
        
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
        case kPhiV:   return AliDielectronHelper::MakeLinBinning(32, 0., 3.2);
        case kPt3D:   return AliDielectronHelper::MakeArbitraryBinning(" 0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.30, 2.50, 3.00, 3.50, 4.00, 5.0, 6.0, 7.0, 10.0, 20.0");
        case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);
        case kPhi3D:  return AliDielectronHelper::MakeLinBinning( 60,0,6.2832); // same as in EffTask
        case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0., 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.65, 0.688, 0.725, 0.75, 0.775, 0.8, 0.85, 0.95, 0.975, 1.0, 1.025, 1.05, 1.125, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 2.85, 2.95, 3.05, 3.1, 3.15, 3.3, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0");
        case kPhiVRebinned:   return AliDielectronHelper::MakeArbitraryBinning("0.0, 2.0, 3.2");
        case kP2D:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 9.00, 10.00");
        default: cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << endl;
            break;
    }
    
}

//______________________________________________________________________________________
void SetTPCSigmaEleCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim) {
    //
    // eta correction for the centroid and width of electron sigmas in the TPC
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
    printf("no default TPC PID eta correction! Corrections were not applied!\n");
    
}

//______________________________________________________________________________________
const AliDielectronEventCuts *GetEventCutsMinBias(){
    
    AliDielectronEventCuts *eventCutsMB=new AliDielectronEventCuts("eventCutsMB","Vertex Track && |vtxZ|<10 && ncontrib>0");
    eventCutsMB->SetRequireVertex();
    eventCutsMB->SetMinVtxContributors(1);
    eventCutsMB->SetVertexZ(-10.,10.);
    // eventCutsMB->SetCentralityRange(0., 1., kTRUE); // optional centrality selection for pp run2
    return eventCutsMB;
}

//______________________________________________________________________________________
const AliDielectronEventCuts *GetEventCutsHighMult(){
    
    AliDielectronEventCuts *eventCutsHM=new AliDielectronEventCuts("eventCutsHM","Vertex Track && |vtxZ|<10 && ncontrib>0");
    eventCutsHM->SetRequireVertex();
    eventCutsHM->SetMinVtxContributors(1);
    eventCutsHM->SetVertexZ(-10.,10.);
    // eventCutsHM->SetCentralityRange(0., 0.1, kTRUE); // optional centrality selection for pp run2
    return eventCutsHM;
}


