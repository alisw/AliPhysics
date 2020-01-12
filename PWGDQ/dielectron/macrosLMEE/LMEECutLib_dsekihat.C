class LMEECutLib {

  public:
    LMEECutLib():
      fCutName(""),
      fIsAOD(kTRUE)
  { 
  }

    LMEECutLib(TString cutname, Bool_t isAOD):
      fCutName(""),
      fIsAOD(kTRUE)
  {
    fCutName = cutname;
    fIsAOD = isAOD;
    Print();
  }
    virtual ~LMEECutLib(){}
    void Print(){ printf("cut name = %s , isAOD = %d\n",fCutName.Data(),fIsAOD); }
    void SetCutName(TString cutname) {fCutName = cutname;}
    void SetAOD(Bool_t isAOD) {fIsAOD = isAOD;}

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
      phiVcut->AddCut(AliDielectronVarManager::kPhivPair, 2.5, 3.2);
      phiVcut->AddCut(AliDielectronVarManager::kM       , 0.0,0.05);
      return phiVcut;
    }

    AliAnalysisCuts *SetupTrackCuts(Float_t PtMin, Float_t PtMax, Float_t EtaMin, Float_t EtaMax){
      AliDielectronCutGroup *trCG = new AliDielectronCutGroup("TrackCutsGroup","TrackCutsGroup",AliDielectronCutGroup::kCompAND);

      AliDielectronTrackCuts *trCuts = new AliDielectronTrackCuts("TrackCuts","TrackCuts");
      trCuts->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
      //trCuts->SetRequireITSRefit(kTRUE);
      //trCuts->SetRequireTPCRefit(kTRUE);
      trCuts->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); // 1<<4

      //if(!fCutName.Contains("PIDCalib",TString::kIgnoreCase)){
      //  if(!fIsAOD){//for ESD
      //    trCG->AddCut(SetupESDtrackCuts());
      //  }
      //  else{//for AOD
      //    trCuts->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); // 1<<4
      //  }
      //}
      trCG->AddCut(trCuts);

      AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");

      varCuts->AddCut(AliDielectronVarManager::kPt, PtMin ,PtMax );
      varCuts->AddCut(AliDielectronVarManager::kEta,EtaMin,EtaMax);

      varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,1.0);
      varCuts->AddCut(AliDielectronVarManager::kImpactParZ, -3.0,3.0);

      varCuts->AddCut(AliDielectronVarManager::kNclsITS,  3.0,6.1);
      varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,0.0,5.0);

      //varCuts->AddCut(AliDielectronVarManager::kNclsTPC,        80.0, 160.0);//should not be used in 2018 PbPb analyses
      varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,     80.0, 160.0);//crossed rows
      varCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,    2.);
      varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0,   4.0);

      //ITS shared cluster cut
      varCuts->AddCut(AliDielectronVarManager::kNclsSITS,-0.1,+0.1);
      trCG->AddCut(varCuts);

      return trCG;
    }

    AliAnalysisCuts *SetupPIDCuts(){
      //AliDielectronPID *pid = new AliDielectronPID();
      //if(fCutName.Contains("noPID",TString::kIgnoreCase) || fCutName.Contains("postPIDCalib",TString::kIgnoreCase) ) return pid;
      //pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -2. ,3. ,0.0, 1e+10, kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
      //pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 1e+10, kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
      //pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 1e+10, kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
      //pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 1e+10, kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
      //return pid;

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
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 4., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kKaon    ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kProton  ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kITS, AliPID::kElectron,   -3., 1., 0., 100., kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        AliDielectronPID* recoverTOF = new AliDielectronPID("recoverTOF", "recoverTOF");
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 4., 0., 100., kTRUE , AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
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
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 4., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kKaon    ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kProton  ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        AliDielectronPID* recoverTOF = new AliDielectronPID("recoverTOF", "recoverTOF");
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 4., 0., 100., kTRUE , AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        hadrej->AddCut(cutsTPC);
        hadrej->AddCut(recoverTOF);
        return hadrej;
      }
      else if(fCutName.Contains("ITSTPChadrej",TString::kIgnoreCase)){
        AliDielectronPID* cutsTPC = new AliDielectronPID("TPCHadRej", "TPCHadRej");
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 4., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kKaon    ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kProton  ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kITS, AliPID::kElectron,   -3., 1., 0., 100., kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        return cutsTPC;
      }
      else if(fCutName.Contains("TPChadrej",TString::kIgnoreCase)){
        AliDielectronPID* cutsTPC = new AliDielectronPID("TPCHadRej", "TPCHadRej");
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 4., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kKaon    ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kProton  ,   -3., 3., 0., 100., kTRUE ,AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        return cutsTPC;
      }
      else if(fCutName.Contains("ITSTOFrecover",TString::kIgnoreCase)){
        AliDielectronPID* recoverTOF = new AliDielectronPID("recoverTOF", "recoverTOF");
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 4., 0., 100., kTRUE , AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kITS, AliPID::kElectron,   -3., 1., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        return recoverTOF;
      }
      else if(fCutName.Contains("TOFrecover",TString::kIgnoreCase)){
        AliDielectronPID* recoverTOF = new AliDielectronPID("recoverTOF", "recoverTOF");
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kPion    , -100., 4., 0., 100., kTRUE , AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        recoverTOF->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,   -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire,AliDielectronVarManager::kPt);
        return recoverTOF;
      }
      else if(fCutName.Contains("noPID",TString::kIgnoreCase)){
        AliDielectronPID* cutsPID   = new AliDielectronPID("PID", "PID");
        return cutsPID;
      }
      else {
        AliDielectronPID* cutsPID   = new AliDielectronPID("PID", "PID");
        return cutsPID;
      }
    }

  protected:
    TString fCutName;
    Bool_t fIsAOD;

    ClassDef(LMEECutLib,1);
};

