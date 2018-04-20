class LMEECutLib {

    public:

    enum LMMECutSet{
        kAllSpecies,
        kElectrons,
        kHighMult,
        kMidMult,
        kLowMult,
        kTTreeCuts
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
				case kTTreeCuts:
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
        case kTTreeCuts
            break;
        case kHighMult:
            AliDielectronVarCuts* centCut1 = new AliDielectronVarCuts("centCutsHigh","MultiplicitypPbLHC16qHigh");
            centCut1->AddCut(AliDielectronVarManager::kCentralityNew,0.,20.);
            centCuts = centCut1;
            break;
         case kMidMult:
            AliDielectronVarCuts* centCut2 = new AliDielectronVarCuts("centCutsMid","MultiplicitypPbLHC16qMid");
            centCut2->AddCut(AliDielectronVarManager::kCentralityNew,20.,60.);
            centCuts = centCut2;
            break;
        case kLowMult:
            AliDielectronVarCuts* centCut3 = new AliDielectronVarCuts("centCutsLow","MultiplicitypPbLHC16qLow");
            centCut3->AddCut(AliDielectronVarManager::kCentralityNew,60.,100.);
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
        case kTTreeCuts
            mixingHandler = new AliDielectronMixingHandler;
            mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
            //mixingHandler->AddVariable(AliDielectronVarManager::kNacc,"0,500");
            mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100");
            // for using TPC event plane, uncorrected. (also, the old phi range was wrong, now same effective binning.)
            // mixingHandler->AddVariable(AliDielectronVarManager::kTPCrpH2uc, 6, TMath::Pi()/-2., TMath::Pi()/2.);
            mixingHandler->SetDepth(20);
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

     //AND cut group to select low mass pairs with large opening angle
    AliDielectronCutGroup* convRejCut = new AliDielectronCutGroup("convRejCut", "convRejCut", AliDielectronCutGroup::kCompAND);
    AliDielectronVarCuts* convMassCut = new AliDielectronVarCuts("convMassCut", "convMassCut");
    AliDielectronVarCuts* convPhiVCut = new AliDielectronVarCuts("convPhiVCut", "convPhiVCut");
    convMassCut->AddCut(AliDielectronVarManager::kM, 0.00, 0.1);
    convPhiVCut->AddCut(AliDielectronVarManager::kPhivPair, 0., 2.);
    convRejCut->AddCut(convMassCut);
    convRejCut->AddCut(convPhiVCut);

    //Mass cut to include any pairs with mass greater than 0.1 GeV
    AliDielectronVarCuts* pairMassCut = new AliDielectronVarCuts("pairMassCut", "pairMassCut");
    pairMassCut->AddCut(AliDielectronVarManager::kM, 0.1, 5.0);

    //OR cut group to accept for above cuts
    AliDielectronCutGroup* allCuts    = new AliDielectronCutGroup("allCuts", "allCuts", AliDielectronCutGroup::kCompOR);

    allCuts->AddCut(convRejCut);
    allCuts->AddCut(pairMassCut);

    //AND cut group for selecting pair momenta
    AliDielectronCutGroup* finalCuts  = new AliDielectronCutGroup("finalCuts", "finalCuts", AliDielectronCutGroup::kCompAND);


    switch(cutSet){
        case kAllSpecies:
        case kElectrons:
        case kHighMult:
        case kMidMult:
        case kLowMult:
				case kTTreeCuts:
            cuts = allCuts;
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
        case kTTreeCuts:
            AliDielectronVarCuts* pairCutsInvM =new AliDielectronVarCuts("pairCutsInvM","pairCutsInvM");
            pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.03); // in upgrade: 0.01
            AliDielectronVarCuts* pairCutsOpAng =new AliDielectronVarCuts("pairCutsOpAng","pairCutsOpAng");
            pairCutsOpAng->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.06); // in upgrade: 0.05

            AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
            pairCutsCG->AddCut(pairCutsInvM);
            pairCutsCG->AddCut(pairCutsOpAng);
            pairCuts = pairCutsCG;
            break;

        default: cout << "No Prefilter Pair Cuts defined " << endl;
    }
    return pairCutsPre;
}



AliAnalysisCuts* LMEECutLib::GetPIDCutsAna(Int_t cutSet) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pidCuts = 0x0;

  //-----------------------------------------------
  // Define different PID Cuts, that are used later
  //-----------------------------------------------
  // PID cuts depend on TPC_inner_p, if not specified
  // PID cut ranges correspond to global momentum P
  // check it again!!!
  //-----------------------------------------------
	
  //Loose cuts used during QA phase
  //wSDD cuts
  AliDielectronPID* PID_wSDD_looseCuts = new AliDielectronPID("PID_wSDD_looseCuts","PID_wSDD_looseCuts");
  PID_wSDD_looseCuts->AddCut(AliDielectronPID::kITS, AliPID::kElectron, -3.0, 1.0, 0.2, 100., kFALSE);
  PID_wSDD_looseCuts->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -1.5, 4.0, 0.2, 100., kFALSE);
  PID_wSDD_looseCuts->AddCut(AliDielectronPID::kTPC, AliPID::kPion, -100., 3.5, 0.2, 100., kTRUE);
  PID_wSDD_looseCuts->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -3.0, 3.0, 0.2, 100., kFALSE, AliDielectronPID::kIfAvailable);
  //CENT_woSDD and pass1_FAST
  AliDielectronPID* PID_woSDD_FAST_looseCuts = new AliDielectronPID("PID_woSDD_FAST_looseCuts","PID_woSDD_FAST_looseCuts");
  PID_woSDD_FAST_looseCuts->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -3.0, 3.0, 0.2, 100., kFALSE);
  PID_woSDD_FAST_looseCuts->AddCut(AliDielectronPID::kTPC, AliPID::kPion, -100., 4.0, 0.2, 100., kTRUE);
  PID_woSDD_FAST_looseCuts->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -3.0, 3.0, 0.4, 100., kFALSE, AliDielectronPID::kRequire);

  //PID cuts used during TTree creating
	//Momentum range relaxed as it cuts on P not Pt. Kinematic cuts applied
	//separately.
  AliDielectronPID* PID_TTreeCuts = new AliDielectronPID("PID_TTreeCuts", "PID_TTreeCuts");
  PID_TTreeCuts->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4. , 0.1, 100., kFALSE);
  
  //-----------------------------------------------
  // Now see what Config actually loads and assemble final cuts
  //-----------------------------------------------

	AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts", AliDielectronCutGroup::kCompAND);
	switch(cutSet){
		case kElectrons:
		case kHighMult:
		case kMidMult:
		case kLowMult:
			if(wSDD){
				cuts->AddCut(PID_wSDD_looseCuts);
			}else{
				cuts->AddCut(PID_woSDD_FAST_looseCuts);
			}
			pidCuts = cuts;
			break;
		case kAllSpecies:
			break;
		case kTTreeCuts:
			cuts->AddCut(PID_TTreeCuts);
			pidCuts = cuts;
			break;
		default:
			std::cout << "No Analysis PID Cut defined " << std::endl;
			return 0x0;
    }
    return pidCuts;
}

AliAnalysisCuts* LMEECutLib::GetKineCutsAna(Int_t cutSet){

	std::cout << "--------------  Get Kinematic Cuts ---------------" << std::endl;
	AliDielectronVarCuts* kineCutsAna = new AliDielectronVarCuts("kineCutsAna", "kineCutsAna");
	if(!kineCutsAna){
			std::cout << "Kinemtaic cuts could not be setup!" << std::endl;
			return 0x0;
	}

	kineCutsAna->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
	kineCutsAna->AddCut(AliDielectronVarManager::kPt, 0.2, 10.);

  return kineCutsAna;

}

AliAnalysisCuts* LMEECutLib::GetKineCutsPre(Int_t cutSet){

  std::cout << "--------------  Get Kinematic Cuts Prefilter ---------------" << std::endl;
  AliDielectronVarCuts* kineCuts = new AliDielectronVarCuts("kineCuts","kineCuts");

  switch(cutSet){
    case kAllSpecies:
    case kElectrons:
		case kTTreeCuts:
      kineCuts->AddCut(AliDielectronVarManager::kPt, 0.2, 10.);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      break;
    default:
      std::cout << "No kinematic cuts used for prefilter" << std::endl;
  }
  return kineCuts;
}

//Make/Tighten track Cuts that are *NOT* already
//done in the AOD production
//**IMPORTANT**: For AODs, select FilterBit
//the method is ignored for ESDs

AliAnalysisCuts* LMEECutLib::GetTrackCutsAna(Int_t cutSet){

	std::cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << std::endl;
	AliDielectronCutGroup* trackCuts = 0x0;
	switch(cutSet){
    //----------
    // these MAIN settings just load the main track selection directly below:
    //----------
		case kAllSpecies:
		case kElectrons:
		case kHighMult:
		case kMidMult:
		case kLowMult:
			AliDielectronVarCuts* trackCutsAOD = new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
			trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY,  - 1.0, 1.0);
			trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,   - 3.0, 3.0);
			if(wSDD){
					trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,  3.0, 100.0);
			}else{
					trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,  2.0, 100.0);
			}
			//trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0, 0.1);
			trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0, 4.5);
			trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      60.0, 200.); //clusters
			trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    70.0, 200.); //findable
			//trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0, 6.0);
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
		case kTTreeCuts:
			//Cuts appear in same order as written in SimpleTreeMaker
			AliDielectronVarCuts* trackCutsTTree = new AliDielectronVarCuts("trackCutsTTree","trackCutsTTree");
			//TPC cuts	
			trackCutsTTree->AddCut(AliDielectronVarManager::kNclsTPC,      70.0, 200.); //clusters
			trackCutsTTree->AddCut(AliDielectronVarManager::kNFclsTPCr,      60.0, 200.); //Crossed (60)
			trackCutsTTree->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.3, 1.1); //Crossed rows over findable
			//DCA
			trackCutsTTree->AddCut(AliDielectronVarManager::kImpactParXY,  - 3.0, 3.0);
			trackCutsTTree->AddCut(AliDielectronVarManager::kImpactParZ,   - 4.0, 4.0);
			//ITS cuts
			if(wSDD){
				trackCutsTTree->AddCut(AliDielectronVarManager::kNclsITS,  3.0, 100.0); // < 3
			}else{
				trackCutsTTree->AddCut(AliDielectronVarManager::kNclsITS,  1.0, 100.0); // < 1
			}
			trackCutsTTree->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0, 36.);
			
			AliDielectronTrackCuts *ttreeCutsDiel = new AliDielectronTrackCuts("ttreeCutsDiel","ttreeCutsDiel");
			//Select filterbit 4
			ttreeCutsDiel->SetAODFilterBit(16);//or 1<<4
			//Refits	
			ttreeCutsDiel->SetRequireITSRefit(kTRUE);
			ttreeCutsDiel->SetRequireTPCRefit(kTRUE);
	
			AliDielectronCutGroup* allTTreeCuts = new AliDielectronCutGroup("allTTreeCuts","allTTreeCuts",AliDielectronCutGroup::kCompAND);
			allTTreeCuts->AddCut(ttreeCutsDiel);
			allTTreeCuts->AddCut(trackCutsTTree);
			trackCuts = allTTreeCuts;
			break;
		default:
			std::cout << "No Analysis Track Cut defined" << std::endl;
		}
		return trackCuts;
}



//Relaxed PID cuts for additional rejection step, do not use blindly
AliAnalysisCuts* LMEECutLib::GetPIDCutsPre(Int_t cutSet){
    std::cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    AliAnalysisCuts* pidCuts = 0x0;

    switch(cutSet){
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
 -
 -
//Possibly different cut sets for Prefilter step
//Not used at the moment
AliAnalysisCuts* LMEECutLib::GetTrackCutsPre(Int_t cutSet){
    std::cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    AliDielectronCutGroup* trackCuts = 0x0;
    switch(cutSet){
        case kAllSpecies:
        case kElectrons:
        trackCuts = LMEECutLib::GetTrackCutsAna(cutSet);
        break;
    default:
        std::cout << "No Prefilter Track Cut defined " << std::endl;
    }
return trackCuts;
}



