class LMEECutLib {

  public:
    LMEECutLib():
      fCutName("")
  { 
  }

    LMEECutLib(TString cutname):
      fCutName("")
  {
    fCutName = cutname;
    Print();
  }
    virtual ~LMEECutLib(){}
    void Print(){ printf("cut name = %s\n",fCutName.Data()); }
    void SetCutName(TString cutname) {fCutName = cutname;}

    static TString GetGeneratorMCSignalName(){
      const TString generators = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6;Pythia CC_0;Pythia BB_0;Pythia B_0;";
      return generators;
    }

    static TString GetResolutionFileName(){
      TString filename = "";
      return filename;
    }

    //define cut configuration
    static AliDielectronEventCuts *SetupEventCuts(Float_t CenMin, Float_t CenMax, Bool_t isRun2, TString estimator){
      AliDielectronEventCuts *eventCuts = new AliDielectronEventCuts("eventCuts","Vertex Any && |vtxZ|<10 && ncontrib>0");
      eventCuts->SetRequireVertex(kTRUE);
      eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
      eventCuts->SetVertexZ(-10.,10.);
      eventCuts->SetMinVtxContributors(1);
      eventCuts->SetTimeRangeCut(kTRUE);
      if(-1 < CenMin && -1 < CenMax){
        eventCuts->SetCentralityEstimator(estimator);
        eventCuts->SetCentralityRange(CenMin,CenMax,isRun2);
      }
      return eventCuts;
    }

    AliESDtrackCuts *SetupESDtrackCuts(){//only for ESD
      AliESDtrackCuts *esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);//bit4
      esdTrackCuts->SetMaxDCAToVertexXY(2.4);
      esdTrackCuts->SetMaxDCAToVertexZ(3.2);
      esdTrackCuts->SetDCAToVertex2D(kTRUE);
      return esdTrackCuts;
      //further cuts are defined in SetupTrackCuts
    }

    AliAnalysisCuts *SetupPhiVPreFilter(){
      AliDielectronVarCuts *phiVcut = new AliDielectronVarCuts("phiVcut","phiVcut");
      phiVcut->AddCut(AliDielectronVarManager::kPhivPair, 2.0,3.2);
      phiVcut->AddCut(AliDielectronVarManager::kM       , 0.0,0.1);
      return phiVcut;
    }

    AliAnalysisCuts *SetupTrackCuts(Float_t PtMin, Float_t PtMax, Float_t EtaMin, Float_t EtaMax){
      AliDielectronCutGroup *trCG = new AliDielectronCutGroup("TrackCutsGroup","TrackCutsGroup",AliDielectronCutGroup::kCompAND);

      if(fCutName.Contains("PIDCalib",TString::kIgnoreCase) || fCutName.Contains("GammaConv",TString::kIgnoreCase)){
        AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
        gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly
        // gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
        gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
        gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);
        gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,  0.25, kFALSE);
        gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
        gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,  0.05, kFALSE);
        gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,  0.05, kFALSE);
        gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,  0.05, kFALSE);
        gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35, 0.35, kFALSE); // should increase purity...
        gammaV0Cuts->SetExcludeTracks(kFALSE);

        /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
        AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");

        varCuts->AddCut(AliDielectronVarManager::kPt,0.15,1e+10);
        varCuts->AddCut(AliDielectronVarManager::kEta,-0.9,+0.9);

        varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,1.0);
        varCuts->AddCut(AliDielectronVarManager::kImpactParZ, -3.0,3.0);

        varCuts->AddCut(AliDielectronVarManager::kNclsITS,  3.0,6.1);
        //varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,0.0,5.0);

        //varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,     80.0, 160.0);//crossed rows
        varCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,    2.);
        varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0,   4.0);

        //ITS shared cluster cut
        varCuts->AddCut(AliDielectronVarManager::kNclsSFracITS,0.2,1.1);

        trCG->AddCut(gammaV0Cuts);
        trCG->AddCut(varCuts);
      }
      else if(fCutName.Contains("Resolution",TString::kIgnoreCase)){
        AliDielectronTrackCuts *trCuts = new AliDielectronTrackCuts("TrackCuts","TrackCuts");
        trCuts->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
        //trCuts->SetRequireITSRefit(kTRUE);
        //trCuts->SetRequireTPCRefit(kTRUE);
        trCuts->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); // 1<<4
        trCG->AddCut(trCuts);

        AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");

        //varCuts->AddCut(AliDielectronVarManager::kPt, PtMin ,PtMax );
        //varCuts->AddCut(AliDielectronVarManager::kEta,EtaMin,EtaMax);

        varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,+1.0);
        varCuts->AddCut(AliDielectronVarManager::kImpactParZ, -3.0,+3.0);

        varCuts->AddCut(AliDielectronVarManager::kNclsITS,  3.0,6.1);
        varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,0.0,15.0);

        //varCuts->AddCut(AliDielectronVarManager::kNclsTPC,        80.0, 160.0);//should not be used in 2018 PbPb analyses
        varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,     80.0, 160.0);//crossed rows
        varCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,    2.);
        varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0,   4.0);

        //ITS shared cluster cut
        varCuts->AddCut(AliDielectronVarManager::kNclsSITS,-0.1,+4.1);//accept upto Nsc = 4 on ITS

        trCG->AddCut(varCuts);
      }
      else if(fCutName.Contains("Loose",TString::kIgnoreCase)){
				trCG = (AliDielectronCutGroup*)SetupPreTrackCutsForProbe();
			}
			else{//Default
        AliDielectronTrackCuts *trCuts = new AliDielectronTrackCuts("TrackCuts","TrackCuts");
        trCuts->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
        //trCuts->SetRequireITSRefit(kTRUE);
        //trCuts->SetRequireTPCRefit(kTRUE);
        trCuts->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); // 1<<4
        trCG->AddCut(trCuts);

        AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");

        varCuts->AddCut(AliDielectronVarManager::kPt, PtMin ,PtMax );
        varCuts->AddCut(AliDielectronVarManager::kEta,EtaMin,EtaMax);

        varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,+1.0);
        varCuts->AddCut(AliDielectronVarManager::kImpactParZ, -3.0,+3.0);

        varCuts->AddCut(AliDielectronVarManager::kNclsITS,  3.9,6.1);
        varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,0.0,5.0);

        //varCuts->AddCut(AliDielectronVarManager::kNclsTPC,        80.0, 160.0);//should not be used in 2018 PbPb analyses
        varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,     80.0, 160.0);//crossed rows
        varCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,    2.);
        varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0,   4.0);

        //ITS shared cluster cut
        //varCuts->AddCut(AliDielectronVarManager::kNclsSITS,-0.1,+0.1);//accept only Nsc = 0 on ITS

				AliDielectronCutGroup *ITSscCG = new AliDielectronCutGroup("ITSscCutsGroup","ITSscCutsGroup",AliDielectronCutGroup::kCompOR);
        AliDielectronVarCuts *varCuts_ITSsc0 = new AliDielectronVarCuts("VarCuts_ITSsc0","VarCuts_ITSsc0");
        varCuts_ITSsc0->AddCut(AliDielectronVarManager::kNclsSITS,-0.1,+0.1);//accept only Nsc = 0 on ITS

        AliDielectronVarCuts *varCuts_ITSsc1 = new AliDielectronVarCuts("VarCuts_ITSsc1","VarCuts_ITSsc1");
        varCuts_ITSsc1->AddCut(AliDielectronVarManager::kNclsSITS,0.9,1.1);//accept Nsc = 1, but not exactly on first SPD
        varCuts_ITSsc1->AddCut(AliDielectronVarManager::kClsS1ITS,-0.1,+0.1);//accept Nsc = 1, but not exactly on first SPD//be carefull! value is double

				ITSscCG->AddCut(varCuts_ITSsc0);
				if(fCutName.Contains("Nsc01",TString::kIgnoreCase)) ITSscCG->AddCut(varCuts_ITSsc1);
				trCG->AddCut(ITSscCG);

        trCG->AddCut(varCuts);
      }
      return trCG;

    }

    AliAnalysisCuts *SetupPIDCuts(){

      if(fCutName.Contains("DefaultPID",TString::kIgnoreCase)){
        AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts", "cuts", AliDielectronCutGroup::kCompAND);
        AliDielectronPID* cutsPID   = new AliDielectronPID("PID", "PID");
        cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -2. ,3. ,0.0, 1e+10, kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
        cutsPID->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 1e+10, kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
        cutsPID->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 1e+10, kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
        cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 1e+10, kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
        cuts->AddCut(cutsPID);
        return cuts;
      }
      else if(fCutName.Contains("ITSTPChadrejORTOFrec",TString::kIgnoreCase)){
        AliDielectronCutGroup* hadrej = new AliDielectronCutGroup("hadrej","hadrej", AliDielectronCutGroup::kCompOR);
        AliDielectronPID* cutsTPC = new AliDielectronPID("cutsTCP", "cutsTCP");
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 3.5, 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kKaon    ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kProton  ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kITS, AliPID::kElectron,   -3., 1., 0., 100., kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        AliDielectronPID* recoverTOF = new AliDielectronPID("recoverTOF", "recoverTOF");
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 3.5, 0., 100., kTRUE , AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kITS, AliPID::kElectron,   -3., 1., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        hadrej->AddCut(cutsTPC);
        hadrej->AddCut(recoverTOF);
        return hadrej;
      }
      else if(fCutName.Contains("TPChadrejORTOFrec",TString::kIgnoreCase)){
        AliDielectronCutGroup* hadrej = new AliDielectronCutGroup("hadrej","hadrej", AliDielectronCutGroup::kCompOR);
        AliDielectronPID* cutsTPC = new AliDielectronPID("cutsTCP", "cutsTCP");
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 3.5, 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kKaon    ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kProton  ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        AliDielectronPID* recoverTOF = new AliDielectronPID("recoverTOF", "recoverTOF");
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 3.5, 0., 100., kTRUE , AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        hadrej->AddCut(cutsTPC);
        hadrej->AddCut(recoverTOF);
        return hadrej;
      }
      else if(fCutName.Contains("ITSTPChadrej",TString::kIgnoreCase)){
        AliDielectronPID* cutsTPC = new AliDielectronPID("TPCHadRej", "TPCHadRej");
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 3.5, 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kKaon    ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kProton  ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kITS, AliPID::kElectron,   -3., 1., 0., 100., kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        return cutsTPC;
      }
      else if(fCutName.Contains("TPChadrej",TString::kIgnoreCase)){
        AliDielectronPID* cutsTPC = new AliDielectronPID("TPCHadRej", "TPCHadRej");
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 3.5, 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kKaon    ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kProton  ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        return cutsTPC;
      }
      else if(fCutName.Contains("ITSTOFrecover",TString::kIgnoreCase)){
        AliDielectronPID* recoverTOF = new AliDielectronPID("recoverTOF", "recoverTOF");
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 3.5, 0., 100., kTRUE , AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kITS, AliPID::kElectron,   -3., 1., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        return recoverTOF;
      }
      else if(fCutName.Contains("TOFrecover",TString::kIgnoreCase)){
        AliDielectronPID* recoverTOF = new AliDielectronPID("recoverTOF", "recoverTOF");
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 3.5, 0., 100., kTRUE , AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        return recoverTOF;
      }
      else if(fCutName.Contains("TightTPCTOF",TString::kIgnoreCase)){
        AliDielectronPID* cutsPID = new AliDielectronPID("PID", "PID");
        cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -2. ,2. ,0.0, 1e+10, kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
        cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 1e+10, kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
        cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        return cutsPID;
      }
      else if(fCutName.Contains("TightTPC",TString::kIgnoreCase)){
        AliDielectronPID* cutsPID = new AliDielectronPID("PID", "PID");
        cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -2. ,2. ,0.0, 1e+10, kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
        cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 1e+10, kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
        return cutsPID;
      }
      else if(fCutName.Contains("PIDCalib",TString::kIgnoreCase) || fCutName.Contains("noPID",TString::kIgnoreCase)){
        AliDielectronPID* cutsPID   = new AliDielectronPID("PID", "PID");
        return cutsPID;
      }
      else {
        AliDielectronPID* cutsPID   = new AliDielectronPID("PID", "PID");
        return cutsPID;
      }
    }

		AliAnalysisCuts *SetupPIDCutsForTag(){
			//AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts", "cuts", AliDielectronCutGroup::kCompAND);
			AliDielectronPID* cutsPIDtag   = new AliDielectronPID("PIDtag", "PIDtag");
			cutsPIDtag->AddCut(AliDielectronPID::kTPC,AliPID::kPion,      -100, 3.5,0.0, 1e+10, kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
			cutsPIDtag->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -3. ,3. ,0.0, 1e+10, kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
			cutsPIDtag->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 1e+10, kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
			return cutsPIDtag;
		}

		AliAnalysisCuts *SetupTrackCutsForTag(){
			AliDielectronCutGroup *trCG = new AliDielectronCutGroup("TrackCutsGroup","TrackCutsGroup",AliDielectronCutGroup::kCompAND);
			AliDielectronTrackCuts *trCuts = new AliDielectronTrackCuts("TrackCuts","TrackCuts");
			trCuts->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
			trCuts->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); // 1<<4
			trCG->AddCut(trCuts);

			AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");

			varCuts->AddCut(AliDielectronVarManager::kPt, 0.2,10.0);
			varCuts->AddCut(AliDielectronVarManager::kEta,-0.8,+0.8);

			//varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,+1.0);
			//varCuts->AddCut(AliDielectronVarManager::kImpactParZ, -3.0,+3.0);
			varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-0.5,+0.5);
			varCuts->AddCut(AliDielectronVarManager::kImpactParZ, -1.0,+1.0);

			varCuts->AddCut(AliDielectronVarManager::kNclsITS,  4.0,6.1);
			varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,0.0,5.0);

			varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,     80.0, 160.0);//crossed rows
			varCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,    2.);
			varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0,   4.0);
			trCG->AddCut(varCuts);

			AliDielectronCutGroup *ITSscCG = new AliDielectronCutGroup("ITSscCutsGroup","ITSscCutsGroup",AliDielectronCutGroup::kCompOR);

			AliDielectronVarCuts *varCuts_ITSsc0 = new AliDielectronVarCuts("VarCuts_ITSsc0","VarCuts_ITSsc0");
			varCuts_ITSsc0->AddCut(AliDielectronVarManager::kNclsSITS,-0.1,+0.1);//accept only Nsc = 0 on ITS

			//AliDielectronVarCuts *varCuts_ITSsc1 = new AliDielectronVarCuts("VarCuts_ITSsc1","VarCuts_ITSsc1");
			//varCuts_ITSsc1->AddCut(AliDielectronVarManager::kNclsSITS, 0.9, 1.1);//accept Nsc = 1, but not exactly on first SPD
			//varCuts_ITSsc1->AddCut(AliDielectronVarManager::kClsS1ITS,-0.1,+0.1);//accept Nsc = 1, but not exactly on first SPD//be carefull! value is double

			ITSscCG->AddCut(varCuts_ITSsc0);
			//ITSscCG->AddCut(varCuts_ITSsc1);
			trCG->AddCut(ITSscCG);
			return trCG;
		}

		AliAnalysisCuts *SetupPrePIDCutsForProbe(){
			AliDielectronPID* precutsPID = new AliDielectronPID("prePID", "prePID");
			precutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -3. ,3. ,0.0, 1e+10, kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
			return precutsPID;
		}
		AliAnalysisCuts *SetupPreTrackCutsForProbe(){
			//do not use AOD filter bit or ESDtrackCuts
			AliDielectronCutGroup *trCG = new AliDielectronCutGroup("TrackCutsGroup","TrackCutsGroup",AliDielectronCutGroup::kCompAND);
			AliDielectronTrackCuts *trCuts = new AliDielectronTrackCuts("TrackCuts","TrackCuts");
			//trCuts->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); // 1<<4
			trCuts->SetRequireITSRefit(kTRUE);
			trCuts->SetRequireTPCRefit(kTRUE);
			trCG->AddCut(trCuts);

			AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");

			varCuts->AddCut(AliDielectronVarManager::kPt, 0.2,10.0);
			varCuts->AddCut(AliDielectronVarManager::kEta,-0.8,+0.8);

			varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,+1.0);
			varCuts->AddCut(AliDielectronVarManager::kImpactParZ, -3.0,+3.0);

			varCuts->AddCut(AliDielectronVarManager::kNclsITS, -0.1,6.1);
			varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,0.0,36.0);

			varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,     70.0, 160.0);//crossed rows
			varCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8, 1e+10);
			varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0,   4.0);
			varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,    -0.1,  +0.1);//should be kFALSE for no kink
			trCG->AddCut(varCuts);

			//trCG->AddCut(SetupESDtrackCuts());
			return trCG;
		}

    AliAnalysisCuts *SetupPIDCutsForPassingProbe(){
			return SetupPIDCuts();
		}

    AliAnalysisCuts *SetupTrackCutsForPassingProbe(){
			return SetupTrackCuts(0.2,10.0,-0.8,+0.8);
		}

  protected:
    TString fCutName;

    ClassDef(LMEECutLib,1);
};

