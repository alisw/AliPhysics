void InitHistograms(AliDielectron *die, Int_t cutDefinition);
TVectorD *GetVector(Int_t var);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition, Double_t ptMin, Bool_t reqTPCnTOF);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);
void SetTPCSigmaEleCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim);
void SetTOFSigmaEleCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim);
const AliDielectronEventCuts *GetEventCutsMinBias();

//TString names = ("pt200_TPCTOFcombITSshared");
TString names = ("pt200_TPCTOFcombITSshared,TPCtight,TPCloose,TOFtight,TOFloose,track1,track2,track3,track4,track5,track6,track7,track8,track9,track10,track11,track12,track13,track14,track15,track16,track17,track18,track19,track20");

TObjArray *arrNames = names.Tokenize(",");
const Int_t nDie = arrNames->GetEntriesFast();

Bool_t hasMC;

enum {kPhiV = 0, kPt3D, kEta3D, kPhi3D, kMee, kMeeLinear, kMeeLinear2, kMeeForDCAee, kPhiVRebinned, kP2D, kPtee3D, kMee3D, kPairDCAsig, kPairDCA, kDeltaPhiLin, kDeltaEtaLin};

Bool_t kPairing  = 1;
Bool_t kMixing   = 1; // kPairing has a higher priority
Bool_t kPairCuts = 0;
Bool_t kTPCCorr  = 1;
Bool_t kTOFCorr  = 1;
Bool_t copyCorr  = 1; // kTRUE to download the TPC & TOF correction maps

AliDielectron* Config_hdegenhardt_pp(Int_t cutDefinition, Bool_t kMinBias = kFALSE, char *period = "16d", Double_t ptMin = 0.2, Bool_t sysUnc = kFALSE, Bool_t reqTPCnTOF = kFALSE)
{
    //
    // Setup the instance of AliDielectron
    //

    //=== create the framework object =================================
    TString name = arrNames->At(cutDefinition)->GetName();
    printf("\n------------------------------------------------------\n");
    printf("Config %d: %s\n",cutDefinition,name.Data());

    AliDielectron *die = new AliDielectron(Form("%s",name.Data()),
                                           Form("Track cuts: %s",name.Data()));
    die->SetHasMC(hasMC);

    //=== Switch off pairing for single track analysis ===============
    if (!kPairing) die->SetNoPairing();

    //=== cut setup ===================================================
    //--- track cuts
    SetupTrackCuts(die, cutDefinition, ptMin, reqTPCnTOF);
    //--- pair cuts
    if (kPairCuts) SetupPairCuts(die, cutDefinition);

    //--- TPCSigmaEle corrections vs P and Eta
    if (kTPCCorr) SetNSigmaEleCorrection(die, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, period, kTRUE); //last kTRUE means TPC

    //--- TOFSigmaEle corrections vs P and Eta
    if (kTOFCorr) SetNSigmaEleCorrection(die, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, period, kFALSE); //last kFALSE means TOF

    //=== add event mixing handler ====================================

    if (kMixing && kPairing){

        AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
        mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
        mix->SetMixType(AliDielectronMixingHandler::kAll);
        mix->SetDepth(20);
        mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
        die->SetMixingHandler(mix);

    }

    InitHistograms(die, cutDefinition, sysUnc);
    //  InitCF(die,cutDefinition);

    //=== Kalman Filter ===============================================
    //--- calculate pairing information not from the KF package,
    //    but directly from the 4-vectors (TLorentzVector)
    //    KF information might be unreliable for some cases, also for event mixing
    die->SetUseKF(kFALSE);

    return die;
}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition, Double_t ptMin, Bool_t reqTPCnTOF)
{
    //
    // Setup the track cuts
    //

    printf("Track cuts: %s\n",trackCutsVar[cutDefinition]);
    Double_t impactParXY = 3.;
    if (trackCutsVar[cutDefinition][0] == '0') impactParXY = 5.;
    else if (trackCutsVar[cutDefinition][0] == '2') impactParXY = 2.;
    Double_t impactParZ = 1.;
    if (trackCutsVar[cutDefinition][1] == '0') impactParZ = 2.;
    else if (trackCutsVar[cutDefinition][1] == '2') impactParZ = 0.7;
    Double_t nClusITS = 3.;
    if (trackCutsVar[cutDefinition][2] == '0') nClusITS = 2.;
    else if (trackCutsVar[cutDefinition][2] == '2') nClusITS = 4.;
    Double_t chi2ITS = 4.5;
    if (trackCutsVar[cutDefinition][3] == '0') chi2ITS = 6.;
    else if (trackCutsVar[cutDefinition][3] == '2') chi2ITS = 3.5;
    Double_t nClusTPC = 80;
    Double_t nFclsTPCr = 100;
    if (trackCutsVar[cutDefinition][4] == '0'){
      nClusTPC = 60;
      nFclsTPCr = 80;
    }
    else if (trackCutsVar[cutDefinition][4] == '2'){
      nClusTPC = 100;
      nFclsTPCr = 120;
    }
    Double_t nFclsTPCfCross = 0.8;
    if (trackCutsVar[cutDefinition][5] == '0') nFclsTPCfCross = 0.6;
    else if (trackCutsVar[cutDefinition][5] == '2') nFclsTPCfCross = 0.9;
    Double_t nClsSFracTPC = 0.4;
    if (trackCutsVar[cutDefinition][6] == '0') nClsSFracTPC = 0.2;
    else if (trackCutsVar[cutDefinition][6] == '2') nClsSFracTPC = 1.0;
    Double_t chi2TPC = 4.;
    if (trackCutsVar[cutDefinition][7] == '0') chi2TPC = 6.;
    else if (trackCutsVar[cutDefinition][7] == '2') chi2TPC = 3.;

    AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
    die->GetTrackFilter().AddCuts(cuts);
    // for AOD analysis, one can cut on the DielectronVarManager variables directly
    AliDielectronVarCuts *etaRange08 = new AliDielectronVarCuts("etaRange08","etaRange08");
    etaRange08->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
    AliDielectronVarCuts *ptRangeMinTo100 = new AliDielectronVarCuts("ptRangeMinTo100","ptRangeMinTo100");
    ptRangeMinTo100->AddCut(AliDielectronVarManager::kPt, ptMin, 100.0);

    AliDielectronVarCuts* trackCutsAOD = new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");

    trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,   nClusTPC, 160.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,nFclsTPCr, 160.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross, nFclsTPCfCross,   1.5);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSFracTPC,   0.0,   nClsSFracTPC);
    trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0,   chi2TPC);

    trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,   nClusITS, 100.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl, 0.0, chi2ITS);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS, -0.1, 0.1);

    trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1*impactParXY, impactParXY);
    trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,   -1*impactParZ, impactParZ);

    cuts->AddCut(etaRange08);
    cuts->AddCut(ptRangeMinTo100);
    cuts->AddCut(trackCutsAOD);

    Double_t TPCv = 0.;
    Double_t TOFv = 0.;
    if (cutDefinition == 1) TPCv = -0.5;	//TPC tight
    if (cutDefinition == 2) TPCv = 0.5;		//TPC loose
    if (cutDefinition == 3) TOFv = -0.5;	//TOF tight
    if (cutDefinition == 4) TOFv = 0.5;		//TOF loose

    // first SPD layer
    AliDielectronTrackCuts *trackCutsSPDfirst = new AliDielectronTrackCuts("trackCutsSPDfirst","trackCutsSPDfirst");
         if (trackCutsVar[cutDefinition][8] == '0') trackCutsSPDfirst->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    else if (trackCutsVar[cutDefinition][8] == '1') trackCutsSPDfirst->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    else if (trackCutsVar[cutDefinition][8] == '2') trackCutsSPDfirst->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    cuts->AddCut(trackCutsSPDfirst);

    // ITS shared cluster cut
    AliDielectronVarCuts *NclsSITSCuts   = new AliDielectronVarCuts("NclsSITSCuts","NclsSITSCuts");
    NclsSITSCuts->AddCut(AliDielectronVarManager::kNclsSITS,  -0.1,   0.1);
    cuts->AddCut(NclsSITSCuts);

    //=== PID cuts ====================================================
    AliDielectronPID *pidTPCTOFreq = new AliDielectronPID("pidTPCTOFreq","pidTPCTOFreq");
    //TPC
    pidTPCTOFreq->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,    -3.-TPCv,  3.+TPCv,   0.0, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
    pidTPCTOFreq->AddCut(AliDielectronPID::kTPC, AliPID::kPion,           -100.,  4.-TPCv,   0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
    // TOF required
    pidTPCTOFreq->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,    -3.-TOFv,  3.+TOFv,   0.4, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kP);

    AliDielectronPID *pidTPCHadRejTOFif = new AliDielectronPID("pidTPCHadRejTOFif","pidTPCHadRejTOFif");
    //TPC
    pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3.-TPCv,  3.+TPCv,   0.0, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
    pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kPion,     -100.     ,  4.-TPCv,   0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
    pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kKaon,       -4.+TPCv,  4.-TPCv,   0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
    pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kProton,     -4.+TPCv,  4.-TPCv,   0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
    // TOF if available
    pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,   -3.-TOFv,  3.+TOFv,   0.4, 1e30, kFALSE,  AliDielectronPID::kIfAvailable, AliDielectronVarManager::kP);

    // combine 2 cut sets with OR option
    AliDielectronCutGroup* combinedPIDcuts = new AliDielectronCutGroup("combinedPIDcuts","combinedPIDcuts",AliDielectronCutGroup::kCompOR);
    combinedPIDcuts->AddCut(pidTPCTOFreq);
    combinedPIDcuts->AddCut(pidTPCHadRejTOFif);

	if (reqTPCnTOF) die->GetTrackFilter().AddCuts(pidTPCTOFreq);
    else die->GetTrackFilter().AddCuts(combinedPIDcuts);
    //die->GetTrackFilter().AddCuts(pidTPCHadRejTOFif);
    
    printf("Track cuts: %.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%c\n",impactParXY,impactParZ,nClusITS,chi2ITS,nClusTPC,nFclsTPCr,nFclsTPCfCross,nClsSFracTPC,chi2TPC,trackCutsVar[cutDefinition][8]);
}


//______________________________________________________________________________________
// not used atm
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
void InitHistograms(AliDielectron *die, Int_t cutDefinition, Bool_t sysUnc)
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
	if (!sysUnc){
		histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0.,1.,AliDielectronVarManager::kNevents);
		histos->UserHistogram("Event","VtxZ","Vertex Z;Vertex Z [cm];N of events",300,-15.,15.,AliDielectronVarManager::kZvPrim);

		// --- V0 info ---------------------------------------------
		histos->UserHistogram("Event","MultV0","Multiplicity V0;V0M amplitude",4000,-0.5,3999.5,AliDielectronVarManager::kMultV0);

		//--- N tracks ---------------------------------------------
		histos->UserHistogram("Event","Ntracks","Number of tracks;N of tracks;N of events", 100, -0.5, 99.5,AliDielectronVarManager::kTracks);

		//=== add histograms to Track classes =============================
		//    name and title (2nd and 3rd parameter) may be omitted
		//    in that case it is take directly from the var manager
		histos->UserHistogram("Track","Pt","Pt;#it{p}_{T} (GeV/#it{c});N of tracks",500,0,50.,AliDielectronVarManager::kPt);

		// --- Eta and Phi ---------------------------------------------
		histos->UserHistogram("Track","Eta","Eta;#eta;#tracks", 100,-1,1,AliDielectronVarManager::kEta);
		histos->UserHistogram("Track","Phi","Phi;#varphi;#tracks", 72,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
		histos->UserHistogram("Track","Eta_Phi","Eta Phi Map;#eta;#varphi", 200,-1,1,360,0,TMath::TwoPi(),AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
	}
    //--- DCA ---------------------------------------------
    histos->UserHistogram("Track","dXY","dXY;DCA_{xy} (cm);N of tracks", 200,-2.,2.,AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","dZ","dZ;DCA_{z} (cm);N of tracks",  400,-4.,4.,AliDielectronVarManager::kImpactParZ);
    histos->UserHistogram("Track","dXY_dZ","dXY dZ Map;DCA_{xy} (cm);DCA_{z} (cm)", 150,-1.5.,1.5.,400,-4.,4.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
	if (!sysUnc){
		histos->UserHistogram("Track","Pt_dcaXYres0","Pt dXYres Map;#it{p}_{T} (GeV/#it{c}); DCA_{xy}^{res} (cm)",150,0.,15.,1000,0.,0.4,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParXYres);
		histos->UserHistogram("Track","Pt_dcaXYres1","Pt dXYres Map;#it{p}_{T} (GeV/#it{c}); DCA_{xy}^{res} (cm)",150,0.,15.,1000,0.,0.04,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParXYres);
		histos->UserHistogram("Track","Pt_dcaXYres2","Pt dXYres Map;#it{p}_{T} (GeV/#it{c}); DCA_{xy}^{res} (cm)",150,0.,15.,1000,0.,0.004,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParXYres);

		histos->UserHistogram("Track","Pt_dXY_phi","#phi vs DCA and Pt;#phi;DCA_{xy}^{e} (cm);#it{p}_{T} (GeV/#it{c}))",72, 0., TMath::TwoPi(), 800, -1., 1., 100, 0., 10., AliDielectronVarManager::kPhi, AliDielectronVarManager::kImpactParXY, AliDielectronVarManager::kPt);
		histos->UserHistogram("Track","dXY_phi","DCA vs #phi;#phi;DCA_{xy}^{e} (cm)",72, 0., TMath::TwoPi(), 800, -1., 1., AliDielectronVarManager::kPhi, AliDielectronVarManager::kImpactParXY);
	}

    //--- track checks (TPC) ---------------------------------------------
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable;Ratio TPC crossed rows over findable;N of tracks",150,0.,1.5,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRows","Number of Crossed Rows TPC;N of TPC crossed rows;N of tracks",200,-0.5,199.5,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;N of TPC clusters;N of tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","TPCchi2","TPC Chi2 value;TPC #chi^{2}/N of clusters;N of tracks",50,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
    histos->UserHistogram("Track","NclsSFracTPC","Fraction of shared clusters assigned in the TPC;Fraction of TPC shared clusters;N of tracks",20,0,1.,AliDielectronVarManager::kNclsSFracTPC);

    //--- track checks (ITS) ---------------------------------------------
    histos->UserHistogram("Track","ITSnCls","Number of Clusters ITS;N of ITS clusters;N of tracks",10,-0.5,9.5,AliDielectronVarManager::kNclsITS);
    histos->UserHistogram("Track","ITSchi2","ITS Chi2 value;ITS #chi^{2}/N of ITS clusters;N of tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
    histos->UserHistogram("Track","ITSSharedClusters","N of ITS shared clusters;N of shared clusters ITS;N of tracks",20,-10.,10.,AliDielectronVarManager::kNclsSITS);
	if (!sysUnc){
		//--- PID information ---------------------------------------------
		//-----------------------------------------------------------------
		//--- ITS ---------------------------------------------------------
		histos->UserHistogram("Track","ITS_dEdx_P","ITS dEdx;#it{p} (GeV/#it{c});ITS signal (arb units)",
							  1000,0.,50.,200,0.,200.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
		histos->UserHistogram("Track","ITSnSigmaEle_P","ITS number of sigmas Electrons;#it{p} (GeV/#it{c});ITS number of sigmas Electrons",
							  1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);

		//--- TPC ----------------------------------------------------------
		histos->UserHistogram("Track","TPC_dEdx_P","TPC dEdx;#it{p} (GeV/#it{c});TPC signal (arb units)",
							  1000,0.,50.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
		histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;#it{p} (GeV/#it{c});TPC number of sigmas Electrons",
							  1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
		histos->UserHistogram("Track","TPCnSigmaEle_Eta_P_lin","TPC number of sigmas Electrons vs Eta and P;Eta;n#sigma_{ele}^{TPC};#it{p}(GeV/#it{c})",
							  200, -1., 1., 200, -10., 10., 100, 0., 5., AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kP);
		histos->UserHistogram("Track","TPCnSigmaEleRaw_Eta_P_lin","TPC number of sigmas Electrons vs Eta and P;Eta;n#sigma_{ele}^{TPC};#it{p}(GeV/#it{c})",
							  200, -1., 1., 200, -10., 10., 100, 0., 5., AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEleRaw, AliDielectronVarManager::kP);

		//--- TOF ----------------------------------------------------------
		histos->UserHistogram("Track","TOFbeta","TOF beta;#it{p} (GeV/#it{c});TOF #beta",
							  1000,0.,50.,300,0.,3.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
		histos->UserHistogram("Track","TOFnSigmaEle_P","TOF number of sigmas Electrons;#it{p} (GeV/#it{c});TOF number of sigmas Electrons",
							  1000,0.,50.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
		if (cutDefinition == 0)
			histos->UserHistogram("Track","TOFnSigmaEle_Eta_P_lin","TOF number of sigmas Electrons vs Eta and P;Eta;n#sigma_{ele}^{TOF};#it{p}(GeV/#it{c})",
							  200, -1., 1., 200, -10., 10., 100, 0., 5., AliDielectronVarManager::kEta, AliDielectronVarManager::kTOFnSigmaEle, AliDielectronVarManager::kP);
		if (cutDefinition == 0)
			histos->UserHistogram("Track","TOFnSigmaEleRaw_Eta_P_lin","TOF number of sigmas Electrons vs Eta and P;Eta;n#sigma_{ele}^{TOF};#it{p}(GeV/#it{c})",
							  200, -1., 1., 200, -10., 10., 100, 0., 5., AliDielectronVarManager::kEta, AliDielectronVarManager::kTOFnSigmaEleRaw, AliDielectronVarManager::kP);
	}
    //=== add histograms to Pair classes ==============================

    if (kPairing){
		if (!sysUnc){
				histos->UserHistogram("Pair","InvMass","Inv.Mass;#it{m}_{ee} (GeV/#it{c}^{2});N of pairs", 400,0.,4.,AliDielectronVarManager::kM);
				histos->UserHistogram("Pair","PairPt","PairPt;Pair #it{p}_{T} (GeV/#it{c});N of pairs", 160,0.,8.,AliDielectronVarManager::kPt);
				histos->UserHistogram("Pair","OpeningAngle","Opening angle;Opening angle;n of pairs", 180,0.,TMath::Pi(),AliDielectronVarManager::kOpeningAngle);
				histos->UserHistogram("Pair","PhiV_Pair","PhiV_Pair;#phi_{V};N of pairs", 180,0.,TMath::Pi(),AliDielectronVarManager::kPhivPair);

				// --- mass vs pT vs phiV ----------------------------------------------------------
				// mass 10 MeV, pT 1 GeV
				histos->UserHistogram("Pair","InvMass_PtRebinned_PhivPair","InvMass:Pt:PhivPair;#it{m}_{ee} (GeV/#it{c}^{2});Pair #it{p}_{T} (GeV/#it{c});#phi_{V}",GetVector(kMeeLinear),GetVector(kPtee3D),GetVector(kPhiV),AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
		}
        // --- mass vs pT vs deltaPhi ----------------------------------------------------------
        // mass 10 MeV, pT 1 GeV
		histos->UserHistogram("Pair","InvMass_PtRebinned_deltaPhi","InvMass:Pt:DeltaPhi;#it{m}_{ee} (GeV/#it{c}^{2});Pair #it{p}_{T} (GeV/#it{c});#delta#varphi",GetVector(kMeeLinear2),GetVector(kPtee3D),GetVector(kDeltaPhiLin),AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kDeltaPhi);
		histos->UserHistogram("Pair","InvMass_PtRebinned_dca","InvMass:Pt:DCA;#it{m}_{ee} (GeV/#it{c}^{2});Pair #it{p}_{T} (GeV/#it{c});DCA_{ee} (#sigma_{xy})",GetVector(kMeeLinear2),GetVector(kPtee3D),GetVector(kPairDCAsig),AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPairDCAsigXY);
		if (!sysUnc){
				// --- mass vs pT vs deltaPhi vs deltaEta ----------------------------------------------------------
				// 4D THnSparse
				const int nDim = 4;
				TObjArray *limits = new TObjArray(nDim);
				limits->AddFirst(GetVector(kMeeLinear));
				limits->AddLast(GetVector(kPtee3D));
				limits->AddLast(GetVector(kDeltaPhiLin));
				limits->AddLast(GetVector(kDeltaEtaLin));
				UInt_t vars[nDim] = {AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kDeltaPhi, AliDielectronVarManager::kDeltaEta};
				histos->UserSparse("Pair", nDim, limits, vars);

				// --- mass vs pT vs deltaPhi vs DCAee ----------------------------------------------------------
				// 4D THnSparse
				const int nDim2 = 4;
				TObjArray *limits2 = new TObjArray(nDim2);
				limits2->AddFirst(GetVector(kMeeForDCAee));
				limits2->AddLast(GetVector(kPtee3D));
				limits2->AddLast(GetVector(kDeltaPhiLin));
				limits2->AddLast(GetVector(kPairDCAsig));
				UInt_t vars2[nDim2] = {AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kDeltaPhi, AliDielectronVarManager::kPairDCAsigXY};
				histos->UserSparse("Pair", nDim2, limits2, vars2);
		}
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
        case kPhiV:         return AliDielectronHelper::MakeLinBinning(30, 0., TMath::Pi());
        case kPt3D:         return AliDielectronHelper::MakeArbitraryBinning(" 0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.30, 2.50, 3.00, 3.50, 4.00, 5.0, 6.0, 7.0, 10.0, 20.0");
        case kEta3D:        return AliDielectronHelper::MakeLinBinning( 50,-1,1);
        case kPhi3D:        return AliDielectronHelper::MakeLinBinning( 60,0,6.2832); // same as in EffTask
        case kMee:          return AliDielectronHelper::MakeArbitraryBinning("0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.32, 0.38, 0.54, 0.66, 0.72, 0.76, 0.78, 0.86, 0.98, 1.02, 1.2, 1.4, 1.64, 1.96, 2.34, 2.84, 2.95, 3.05, 3.1, 3.15, 3.3, 3.5, 3.75, 4.0");
        case kMeeLinear:    return AliDielectronHelper::MakeLinBinning(200, 0., 4.);
        case kMeeLinear2:    return AliDielectronHelper::MakeLinBinning(40, 0., 4.);
        case kMeeForDCAee:  return AliDielectronHelper::MakeArbitraryBinning("0., 0.08, 0.14, 0.2, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.2, 5.0");
        case kPhiVRebinned: return AliDielectronHelper::MakeArbitraryBinning("0.0, 2.0, 3.2");
        case kP2D:          return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 9.00, 10.00");
        case kPtee3D:       return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0, 100.0");
        case kMee3D:        return AliDielectronHelper::MakeArbitraryBinning("0., 0.15, 0.75, 1.25, 2.75, 3.5, 4.0");
        case kPairDCAsig:   return AliDielectronHelper::MakeLinBinning(200, 0., 20.);//MakeArbitraryBinning("0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0");
        case kPairDCA:      return AliDielectronHelper::MakeLinBinning(100, 0., 2.);
        case kDeltaPhiLin:  return AliDielectronHelper::MakeLinBinning(36, -TMath::Pi(), TMath::Pi());
        case kDeltaEtaLin:  return AliDielectronHelper::MakeLinBinning(36, 0., 1.6.);

        default: cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << endl;
            break;
    }

}

//______________________________________________________________________________________
void SetNSigmaEleCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, char *period, Bool_t isTPC) {
    //
    // eta correction for the centroid and width of electron sigmas in the TPC
    //
    Bool_t test = 0;
    //
    printf("-> Starting SetNSigmaEleCorrection() for %s\n", isTPC?"TPC":"TOF");
    printf("   Correction Xdim = %s\n", AliDielectronVarManager::GetValueName(corrXdim));
    printf("   Correction Ydim = %s\n", AliDielectronVarManager::GetValueName(corrYdim));

    if (corrYdim!=AliDielectronVarManager::kEta) { printf("#  Correction only available for Ydim = eta!\n"); return; }

    if (corrXdim!=AliDielectronVarManager::kP){
        printf("#  No correction available for Xdim = %s!\n", AliDielectronVarManager::GetValueName(corrXdim));
        printf("#  No correction applied!\n");
        return;
    }

    TH2F* histMean2D;
    TH2F* histWidth2D;

    if (!GetCorrectionsHisto(&histMean2D, &histWidth2D, period, isTPC)){
        printf("#  File with corrections could not be accessed!\n");
        printf("#  No correction applied for %s!\n",isTPC?"TPC":"TOF");
        return;
    }

    if (isTPC){
      die->SetCentroidCorrFunction(histMean2D, corrXdim, corrYdim);
      die->SetWidthCorrFunction(histWidth2D, corrXdim, corrYdim);
    }
    else{
      die->SetCentroidCorrFunctionTOF(histMean2D, corrXdim, corrYdim);
      die->SetWidthCorrFunctionTOF(histWidth2D, corrXdim, corrYdim);
    }
    
    printf(">> %s PID eta correction loaded!\n", isTPC?"TPC":"TOF");
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

Bool_t GetCorrectionsHisto(TH2F **mean, TH2F **width, char *period, Bool_t isTPC){
	char rootFile[100];
    sprintf(rootFile,"alien://alice/cern.ch/user/h/hfranzde/TPCnTOFcor/calMaps1%c.root",period[1]);
    
    TFile *f;
    char fileName[100];
    sprintf(fileName,"calMaps1%c.root",period[1]);
    f = TFile::Open(fileName);
  	
    if (copyCorr && !f){
		TGrid::Connect("alien://");
		gSystem->Exec(Form("alien_cp %s .",rootFile));
		f = TFile::Open(fileName);
		if (!f) return kFALSE;
	}
	
  	char *p = GetCorPeriodMap(period);

  	char hNameM[100]; //Mean histo name
  	char hNameW[100]; //Width histo name
  	sprintf(hNameM,"%sm%s",isTPC?"TPC":"TOF",p);
  	sprintf(hNameW,"%sw%s",isTPC?"TPC":"TOF",p);
  	
  	printf("Loading Recalibration Map of %s for %s!\n",p,period);
    TH2F *m;
    TH2F *w;
  	f->GetObject(hNameM,m);
  	f->GetObject(hNameW,w);

    for (Int_t i = 0; i <= m->GetNbinsX()+1; i++){
          for (Int_t k = 0; k <= m->GetNbinsY()+1; k++){
              if ( (i == 0) || (k == 0) || (i > m->GetNbinsX()) || (k > m->GetNbinsY())) { // under/overflows
                  m->SetBinContent(i, k, 0.0 );
                  w->SetBinContent(i, k, 1.0 );
              }
          }
    }

    *mean = m->Clone();
    *width = w->Clone();
    
    return kTRUE;
}

char *GetCorPeriodMap(char *p){ //If low stat period, use the maps for all periods
		 if (p[1] == '6') return "16ALL";
	else if (p[1] == '7') return "17ALL";
	else if (p == "18spl")   return "18ALLsplines";
	else if (p == "18noSpl") return "18ALLnoSplines";
	else if (p == "18b") return "18b";
	else if (p == "18d") return "18d";
	else if (p == "18e") return "18e";
	else if (p == "18f") return "18f";
	else if (p == "18h") return "18h";
	else if (p == "18j") return "18j";
	else if (p == "18l") return "18l";
	else if (p == "18g") return "18g";
	else if (p == "18i") return "18i";
	else if (p == "18k") return "18k";
	else if (p == "18m") return "18m";
	else if (p == "18n") return "18n";
	else if (p == "18o") return "18o";
	else if (p == "18p") return "18p";
	else return p;
}

char *trackCutsVar[25];
trackCutsVar[0] = "1111111111"; //Default
trackCutsVar[1] = "1111111111"; //TPC loose
trackCutsVar[2] = "1111111111"; //TPC tight
trackCutsVar[3] = "1111111111"; //TOF loose
trackCutsVar[4] = "1111111111"; //TOF tight

trackCutsVar[5] = "0120002021"; //track01
trackCutsVar[6] = "1212202100";
trackCutsVar[7] = "2121200002";
trackCutsVar[8] = "2012011220";
trackCutsVar[9] = "0112021112"; //track05
trackCutsVar[10] = "2022210010";
trackCutsVar[11] = "0202221220";
trackCutsVar[12] = "2111001020";
trackCutsVar[13] = "0200221122";
trackCutsVar[14] = "2200010100"; //track10
trackCutsVar[15] = "1020110110";
trackCutsVar[16] = "0001120212";
trackCutsVar[17] = "1222122110";
trackCutsVar[18] = "2201101100";
trackCutsVar[19] = "2222212002"; //track15
trackCutsVar[20] = "2212102120";
trackCutsVar[21] = "1212100100";
trackCutsVar[22] = "2101212101";
trackCutsVar[23] = "2002120100";
trackCutsVar[24] = "2201210011"; //track20
