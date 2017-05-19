class LMEECutLib {
  
    public:

    enum LMMECutSet{
        kAllSpecies,
        kElectrons,
        kHighMult,
        kMidMult,
        kLowMult,
        kLowPt,
        kMidLowPt,
        kMidPt,
        kMidHighPt,
        kHighPt
    };
    
    
    LMEECutLib(Bool_t wSDD): wSDD(wSDD){}
    
    //Getters
    AliDielectronEventCuts*     GetEventCuts(Int_t cutSet);
    AliAnalysisCuts*            GetCentralityCuts(Int_t centSel);
    AliDielectronTrackRotator*  GetTrackRotator(Int_t cutSet);
    AliDielectronMixingHandler* GetMixingHandler(Int_t cutSet);
    
    AliAnalysisCuts* GetPairCutsAna(Int_t cutSet); //Bool_t togglePC=kFALSE
    AliAnalysisCuts* GetPairCutsPre(Int_t cutSet);
    
    AliAnalysisCuts* GetPIDCutsAna(Int_t cutSet);
    AliAnalysisCuts* GetPIDCutsPre(Int_t cutSet);
    
    AliAnalysisCuts* GetTrackCutsAna(Int_t cutSet);
    AliAnalysisCuts* GetTrackCutsPre(Int_t cutSet);

  
    private:
        Bool_t wSDD;
};


// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
// the selection is hardcoded in the AddTask
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Int_t cutSet) {

    AliDielectronEventCuts* eventCuts = 0x0;

    switch(cutSet){
        case kAllSpecies:
        case kElectrons:
        case kHighMult:
        case kMidMult:
        case kLowMult:
        case kLowPt:
        case kMidLowPt:
        case kMidPt:
        case kMidHighPt:
        case kHighPt:
            eventCuts = new AliDielectronEventCuts("eventCuts_acapon","Vertex Track && |vtxZ|<10 && ncontrib>0");
            eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
            eventCuts->SetRequireVertex();
            eventCuts->SetMinVtxContributors(1);
            eventCuts->SetVertexZ(-10.,10.);
            //eventCuts->SetCentralityRange(0,80,kTRUE); //Use Run2 centrality definitions
            break;
 
    default: cout << "No Event Cut defined" << endl;
    }
    return eventCuts;
}


//Centrality selection done in Event selection
AliAnalysisCuts* LMEECutLib::GetCentralityCuts(Int_t centSel) {
    AliAnalysisCuts* centCuts = 0x0;
    switch(centSel){
        case kAllSpecies:
        case kElectrons:
        case kLowPt:
        case kMidLowPt:
        case kMidPt:
        case kMidHighPt:
        case kHighPt:
            break;
        case kHighMult:
            AliDielectronVarCuts* centCut1 = new AliDielectronVarCuts("centCutsHigh","MultiplicitypPbLHC16qHigh");
            centCut1->AddCut(AliDielectronVarManager::kCentralityNew,0.,20.);
            centCuts = centCut1;
            break;
         case kMidMult:
            AliDielectronVarCuts* centCut2 = new AliDielectronVarCuts("centCutsMid","MultiplicitypPbLHC16qMid");
            centCut2->AddCut(AliDielectronVarManager::kCentralityNew,20.,80.);
            centCuts = centCut2;
            break;
        case kLowMult:
            AliDielectronVarCuts* centCut3 = new AliDielectronVarCuts("centCutsLow","MultiplicitypPbLHC16qLow");
            centCut3->AddCut(AliDielectronVarManager::kCentralityNew,80.,100.);
            centCuts = centCut3;
            break;
        default: cout << "No Centrality selected" << endl;
    }
    return centCuts;
}


AliDielectronMixingHandler* LMEECutLib::GetMixingHandler(Int_t cutSet) {
    AliDielectronMixingHandler* mixingHandler = 0x0;
    switch (cutSet) {
        case kAllSpecies: 
        case kElectrons:
        case kHighMult:
        case kMidMult:
        case kLowMult:
        case kLowPt:
        case kMidLowPt:
        case kMidPt:
        case kMidHighPt:
        case kHighPt:
            mixingHandler = new AliDielectronMixingHandler;
            mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
            mixingHandler->AddVariable(AliDielectronVarManager::kNacc,"0,500");
            // for using TPC event plane, uncorrected. (also, the old phi range was wrong, now same effective binning.)
            // mixingHandler->AddVariable(AliDielectronVarManager::kTPCrpH2uc, 6, TMath::Pi()/-2., TMath::Pi()/2.);
            mixingHandler->SetDepth(10);
            mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
            break;
        //[...]
        default: cout << "No Mixer defined" << endl;
    }
    return mixingHandler;
}



//Pair Cuts for Analysis step - take care of logic - inverted compared to other PairCuts!!
// cuts = SELECTION!!!
AliAnalysisCuts* LMEECutLib::GetPairCutsAna(Int_t cutSet)  {
    
    std::cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    AliAnalysisCuts* cuts = 0x0;
    switch(cutSet){
        case kAllSpecies: 
        case kElectrons: 
        case kHighMult:
        case kMidMult:
        case kLowMult:
            AliDielectronVarCuts* pairCuts = new AliDielectronVarCuts("pairCuts", "pairCuts");
            pairCuts->AddCut(AliDielectronVarManager::kM, 0.05, 5.0);
            pairCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.0, 2.5);
            pairCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.05, 3.2); 
            cuts = pairCuts;
            break;
        case kLowPt:
            AliDielectronVarCuts* lowPtCuts = new AliDielectronVarCuts("lowPtCuts ","lowPtCuts ");
            lowPtCuts->AddCut(AliDielectronVarManager::kPt, 0.05, 0.35);
            lowPtCuts->AddCut(AliDielectronVarManager::kM, 0.05, 5.0);
            lowPtCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.0, 2.5);
            lowPtCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.05, 3.2); 
            cuts = lowPtCuts;
            break;
        case kMidLowPt:
            AliDielectronVarCuts* midLowPtCuts = new AliDielectronVarCuts("midLowPtCuts ","midLowPtCuts ");
            midLowPtCuts->AddCut(AliDielectronVarManager::kPt, 0.35, 0.55);
            midLowPtCuts->AddCut(AliDielectronVarManager::kM, 0.05, 5.0);
            midLowPtCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.0, 2.5);
            midLowPtCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.05, 3.2); 
            cuts = midLowPtCuts;
            break;
        case kMidPt:
            AliDielectronVarCuts* midPtCuts = new AliDielectronVarCuts("midPtCuts ","midPtCuts ");
            midPtCuts->AddCut(AliDielectronVarManager::kPt, 0.55, 0.75);
            midPtCuts->AddCut(AliDielectronVarManager::kM, 0.05, 5.0);
            midPtCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.0, 2.5);
            midPtCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.05, 3.2); 
            cuts = midPtCuts;
            break;
         case kMidHighPt:
            AliDielectronVarCuts* midHighPtCuts = new AliDielectronVarCuts("midHighPtCuts ","midHighPtCuts ");
            midHighPtCuts->AddCut(AliDielectronVarManager::kPt, 0.75, 1.1);
            midHighPtCuts->AddCut(AliDielectronVarManager::kM, 0.05, 5.0);
            midHighPtCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.0, 2.5);
            midHighPtCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.05, 3.2); 
            cuts = midHighPtCuts;
            break;
        case kHighPt:
            AliDielectronVarCuts* highPtCuts = new AliDielectronVarCuts("highPtCuts ","highPtCuts ");
            highPtCuts->AddCut(AliDielectronVarManager::kPt, 1.1, 10.0);
            highPtCuts->AddCut(AliDielectronVarManager::kM, 0.05, 5.0);
            highPtCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.0, 2.5);
            highPtCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.05, 3.2); 
            cuts = highPtCuts;
            break;
        default: 
            std::cout << "No Pair Cuts defined " << std::endl;
    }
    return cuts;
}


//Pair Cuts for PREFILTER step
// cuts = REJECTION!!!
AliAnalysisCuts* LMEECutLib::GetPairCutsPre(Int_t cutSet)  {  
    cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
    AliAnalysisCuts* pairCutsPre = 0x0;
    switch(cutSet){
        case kAllSpecies: 
        case kElectrons: 
        case kLowPt:
        case kMidLowPt:
        case kMidPt:
        case kMidHighPt:
        case kHighPt:
            AliDielectronVarCuts* cutsPre = new AliDielectronVarCuts("cutsPre", "cutsPre");
            cutsPre->AddCut(AliDielectronVarManager::kM, 0.0, 0.05);
            cutsPre->AddCut(AliDielectronVarManager::kPhivPair, 2.5, 3.2);
            cutsPre->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.05);
            pairCutsPre = cutsPre;
            break;
    
        default: cout << "No Prefilter Pair Cuts defined " << endl;
    } 
    return pairCutsPre;
}



AliAnalysisCuts* LMEECutLib::GetPIDCutsAna(Int_t cutSet) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pidCuts=0x0;
  
  //-----------------------------------------------
  // Define different PID Cuts, that are used later
  //-----------------------------------------------
  // PID cuts depend on TPC_inner_p, if not specified
  // PID cut ranges correspond to global momentum P
  // check it again!!!
  //-----------------------------------------------
  
  //
  //
  //For electron PID
  //wSDD cuts
  AliDielectronPID* PID_wSDD_looseCuts = new AliDielectronPID("PID_wSDD_looseCuts","PID_wSDD_looseCuts");
  PID_wSDD_looseCuts->AddCut(AliDielectronPID::kITS, AliPID::kElectron, -3.0, 1.0, 0.2, 100., kFALSE);
  PID_wSDD_looseCuts->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -1.5, 4.0, 0.2, 100., kFALSE);
  PID_wSDD_looseCuts->AddCut(AliDielectronPID::kTPC, AliPID::kPion, -100., 3.5, 0.2, 100., kTRUE);
  PID_wSDD_looseCuts->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -3.0, 3.0, 0.2, 100., kFALSE, AliDielectronPID::kIfAvailable);

  //CENT_woSDD and pass1_FAST, loose cuts
  AliDielectronPID* PID_woSDD_FAST_looseCuts = new AliDielectronPID("PID_woSDD_FAST_looseCuts","PID_woSDD_FAST_looseCuts");
  PID_woSDD_FAST_looseCuts->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -3.0, 3.0, 0.2, 100., kFALSE);
  PID_woSDD_FAST_looseCuts->AddCut(AliDielectronPID::kTPC, AliPID::kPion, -100., 4.0, 0.2, 100., kTRUE);
  PID_woSDD_FAST_looseCuts->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -3.0, 3.0, 0.4, 100., kFALSE, AliDielectronPID::kRequire);


  ///////////////////////
  // LOOSE PID ITS+TPC+TOFif
  AliDielectronPID* pidTPCITS_TOFif_LOOSE = new AliDielectronPID("pidTPCITS_TOFif_LOOSE","pidTPCITS_TOFif_LOOSE");
  pidTPCITS_TOFif_LOOSE->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3. ,3. , 0. ,100., kFALSE);
  pidTPCITS_TOFif_LOOSE->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-3. ,3. , 0. ,100., kFALSE);
  pidTPCITS_TOFif_LOOSE->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  //

  // tighter PID ITS+TPC+TOFif
  // ITS only up to momentum where proton contamination is seen in TPC signal
  AliDielectronPID *pidTPCITS_TOFif56 = new AliDielectronPID("pidTPCITS_TOFif56","pidTPCITS_TOFif56");
  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 2.5, 0. ,100., kFALSE);
  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 0.5, 0. ,  2., kFALSE);
  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  
  // PID for V0 task
  AliDielectronPID *pid_V0select_1 = new AliDielectronPID("pid_V0select_1","pid_V0select_1");
  pid_V0select_1->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
  pid_V0select_1->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -1.5, 1.5, 0. ,100., kFALSE);
  
  //-----------------------------------------------
  // Now see what Config actually loads and assemble final cuts
  //-----------------------------------------------

    switch(cutSet){
        case kElectrons:
        case kHighMult:   
        case kMidMult:
        case kLowMult:
        case kLowPt:
        case kMidLowPt:
        case kMidPt:
        case kMidHighPt:
        case kHighPt:
            AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts", AliDielectronCutGroup::kCompAND);
            if(wSDD){
              cuts->AddCut(PID_wSDD_looseCuts);
            }else{
              cuts->AddCut(PID_woSDD_FAST_looseCuts);
            }
            pidCuts = cuts; 
            break;
        case kAllSpecies:;
            break;
        default: 
            std::cout << "No Analysis PID Cut defined " << std::endl;
        return 0x0;
    }
    return pidCuts;
}

AliAnalysisCuts* LMEECutLib::GetKineCutsAna(Int_t cutSet){

    cout << "--------------  Get Kinematic Cuts ---------------" << endl;
    AliDielectronVarCuts* kineCuts = static_cast<AliDielectronVarCuts*>(GetKineCutsPre(cutSet));
    kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
    if(!kineCuts){
        std::cout << "Kinemtaic cuts could not be setup!" << std::endl;
        return 0x0;
    }
    switch(cutSet){
        case kAllSpecies:
        case kElectrons:
        case kLowPt:
        case kMidLowPt:
        case kMidPt:
        case kMidHighPt:
        case kHighPt:     
            kineCuts->AddCut(AliDielectronVarManager::kPt, 0.2, 10.);
            break;
        
        default:
            kineCuts->AddCut(AliDielectronVarManager::kPt, 0.2, 10.);
  }

  return kineCuts;
  
}

AliAnalysisCuts* LMEECutLib::GetKineCutsPre(Int_t cutSet){

  cout << "--------------  Get Kinematic Cuts Prefilter ---------------" << endl;
  AliDielectronVarCuts* kineCuts = new AliDielectronVarCuts("kineCuts","kineCuts");
  
  switch(cutSet){
    case kAllSpecies:
    case kElectrons:
      kineCuts->AddCut(AliDielectronVarManager::kPt, 0.2, 10.);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      break;
    default:
      cout << "No kinematic cuts used for prefilter" << endl;
  }
  return kineCuts;
}

//Make/Tighten track Cuts that are *NOT* already
//done in the AOD production
//**IMPORTANT**: For AODs, select FilterBit
//the method is ignored for ESDs

AliAnalysisCuts* LMEECutLib::GetTrackCutsAna(Int_t cutSet) {
    cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
    AliDielectronCutGroup* trackCuts=0x0;
    switch(cutSet){
    //----------
    // these MAIN settings just load the main track selection directly below:
    //----------
        case kAllSpecies: 
        case kElectrons: 
        case kHighMult:
        case kMidMult:
        case kLowMult:
        case kLowPt:
        case kMidLowPt:
        case kMidPt:
        case kMidHighPt:
        case kHighPt:
            AliDielectronVarCuts* trackCutsAOD = new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY,  - 1.0, 1.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,   - 3.0, 3.0);
            if(wSDD){
                trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
            }else{
                trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      2.0, 100.0);
            }
            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0, 36.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,       70.0, 200.);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     60.0, 200.);
            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0, 6.0);
            //trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCrFrac,  0.3, 10.); //Number of found/findable
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.3, 1.1); //Crossed rows over findable 
            //Lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
            AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
            //trackCutsDiel->SetAODFilterBit(0<<0); // (=0) filterbit 0! //GetStandardITSTPCTrackCuts2010(kFALSE); loose DCA, 2D cut
            trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
            trackCutsDiel->SetRequireITSRefit(kTRUE);
            trackCutsDiel->SetRequireTPCRefit(kTRUE);
            
            AliDielectronCutGroup* cgPIDCuts = new AliDielectronCutGroup("cgPIDCuts","cgPIDCuts",AliDielectronCutGroup::kCompAND);
            cgPIDCuts->AddCut(trackCutsDiel);
            cgPIDCuts->AddCut(trackCutsAOD);
            trackCuts = cgPIDCuts;
            break;
        default: cout << "No Analysis Track Cut defined " << endl;
    }
    return trackCuts;
} 



//Relaxed PID cuts for additional rejection step, do not use blindly
AliAnalysisCuts* LMEECutLib::GetPIDCutsPre(Int_t cutSet) {
    std::cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    AliAnalysisCuts* pidCuts=0x0;

    switch (cutSet) {
    case kAllSpecies: 
        AliDielectronCutGroup* InitialFilterCG = new AliDielectronCutGroup("IntitialFilterCG","InitialFilterCG", AliDielectronCutGroup::kCompOR);
        //InitialFilterCG->AddCut(GetTrackCutsAna(cutSet));
        InitialFilterCG->AddCut(GetPIDCutsAna(cutSet));
        pidCuts = InitialFilterCG;
        break;
    case kElectrons: 
        AliDielectronCutGroup* InitialFilterCGelecs = new AliDielectronCutGroup("IntitialFilterCGelecs","InitialFilterCGelecs", AliDielectronCutGroup::kCompOR);
        //InitialFilterCGelecs->AddCut(GetTrackCutsAna(cutSet));
        InitialFilterCGelecs->AddCut(GetPIDCutsAna(cutSet));
        pidCuts = InitialFilterCGelecs;
        break;
      
    default: 
        std::cout << "No Prefilter PID Cut defined " << std::endl;
    }
    return pidCuts;
}


//Possibly different cut sets for Prefilter step
//Not used at the moment
AliAnalysisCuts* LMEECutLib::GetTrackCutsPre(Int_t cutSet) {
    std::cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    AliDielectronCutGroup* trackCuts=0x0;
    switch (cutSet) {
    case kAllSpecies: 
    case kElectrons: 
        trackCuts = LMEECutLib::GetTrackCutsAna(cutSet);
    break;
        /*
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1); 

        cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
        cgTrackCutsPre->AddCut(trackCutsDiel);
        cgTrackCutsPre->AddCut(trackCutsAOD);
        cgTrackCutsPre->AddCut(GetKineCutsAna(cutSet));
        trackCuts = cgTrackCutsPre;*/

    default: 
        std::cout << "No Prefilter Track Cut defined " << std::endl;
    }
    return trackCuts;
}



