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

	AliDielectronCutGroup* GetPairCuts(Int_t cutSet);
	AliDielectronCutGroup* GetPIDCuts(Int_t PIDcuts);
	AliDielectronCutGroup* GetTrackCuts(Int_t cutSet, Int_t PIDcuts);

	static TH3D SetEtaCorrectionTPC( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise);

	private:
			Bool_t wSDD;
};


static TH3D LMEECutLib::SetEtaCorrectionTPC( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise, int sel) {
	

}

// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
// the selection is hardcoded in the AddTask
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Int_t cutSet) {

    ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetEventCuts() >>>>>>>>>>>>>>>>>>>>>> ");
    AliDielectronEventCuts* eventCuts = new AliDielectronEventCuts("eventCuts_acapon","Vertex Track && |vtxZ|<10 && ncontrib>0");

    switch(cutSet){
        case kAllSpecies:
        case kElectrons:
        case kHighMult:
        case kMidMult:
        case kLowMult:
				case kTTreeCuts:
            eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
            eventCuts->SetRequireVertex();
            eventCuts->SetMinVtxContributors(1);
            eventCuts->SetVertexZ(-10.,10.);
            break;
				default: cout << "No Event Cut defined" << endl;
    }
    return eventCuts;
}


//Centrality selection done in Event selection
AliDielectronCutGroup* LMEECutLib::GetCentralityCuts(Int_t centSel) {
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
				case kTTreeCuts:
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
AliDielectronCutGroup* LMEECutLib::GetPairCuts(Int_t cutSet)  {

    ::Info("LMEECutLibg_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetPairCuts() >>>>>>>>>>>>>>>>>>>>>> ");
		//Final OR cut group to incorporate the following cuts (below)
    AliDielectronCutGroup* allCuts    = new AliDielectronCutGroup("allCuts", "allCuts", AliDielectronCutGroup::kCompOR);

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

    allCuts->AddCut(convRejCut);
    allCuts->AddCut(pairMassCut);

    return allCuts;
}


AliDielectronCutGroup* LMEECutLib::GetPIDCuts(Int_t PIDcuts) {
  
	::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetPIDCuts() >>>>>>>>>>>>>>>>>>>>>> ");

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
	switch(PIDcuts){
		case kElectrons:
		case kHighMult:
		case kMidMult:
		case kLowMult:
			if(wSDD){
				cuts->AddCut(PID_wSDD_looseCuts);
			}else{
				cuts->AddCut(PID_woSDD_FAST_looseCuts);
			}
			break;
		case kAllSpecies:
			break;
		case kTTreeCuts:
			cuts->AddCut(PID_TTreeCuts);
			break;
		default:
			std::cout << "No Analysis PID Cut defined " << std::endl;
			return 0x0;
    }
    return cuts;
}

//Make/Tighten track Cuts that are *NOT* already
//done in the AOD production
//**IMPORTANT**: For AODs, select FilterBit
//the method is ignored for ESDs

AliDielectronCutGroup* LMEECutLib::GetTrackCuts(Int_t cutSet, Int_t PIDcuts){

	::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetTrackCuts() >>>>>>>>>>>>>>>>>>>>>> ");
	AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("trackCuts", "trackCuts", AliDielectronCutGroup::kCompAND);

	AliDielectronVarCuts* varCutsFilter     = new AliDielectronVarCuts("varCutsFilter", "varCutsFilter");
	AliDielectronTrackCuts* trackCutsFilter = new AliDielectronTrackCuts("trackCutsFilter", "trackCutsFilter");

	switch(cutSet){
    //----------
    // these MAIN settings just load the main track selection directly below:
    //----------
		case kAllSpecies:
		case kElectrons:
		case kHighMult:
		case kMidMult:
		case kLowMult:
			varCutsFilter->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
			varCutsFilter->AddCut(AliDielectronVarManager::kPt, 0.2, 10.);
			varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,  - 1.0, 1.0);
			varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,   - 3.0, 3.0);
			if(wSDD){
					varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  3.0, 100.0);
			}else{
					varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  2.0, 100.0);
			}
			//varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0, 0.1);
			varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0, 4.5);
			varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,      60.0, 200.); //clusters
			varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,    70.0, 200.); //findable
			//varCutsFilter->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0, 6.0);
			//varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCrFrac,  0.3, 10.); //Number of found/findable
			varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.3, 1.1); //Crossed rows over findable

			trackCutsFilter->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
			trackCutsFilter->SetRequireITSRefit(kTRUE);
			trackCutsFilter->SetRequireTPCRefit(kTRUE);

			trackCuts->AddCut(trackCutsFilter);
			trackCuts->AddCut(varCutsFilter);
			
			trackCuts->AddCut(GetPIDCuts(PIDcuts));
			break;
		case kTTreeCuts:
			varCutsFilter->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
			varCutsFilter->AddCut(AliDielectronVarManager::kPt, 0.2, 10.);
			//TPC cuts
			//Clusters
			varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,      70.0, 200.); 
			//Crossed rows
			varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      60.0, 200.); 
			//Crossed rows over findable
			varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.3, 1.1); 
			//DCA
			varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,  - 3.0, 3.0);
			varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,   - 4.0, 4.0);
			//ITS cuts
			if(wSDD){
				varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  3.0, 100.0); // < 3
			}else{
				varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  1.0, 100.0); // < 1
			}
			varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0, 36.);
			
			//Select filterbit 4
			trackCutsFilter->SetAODFilterBit(16);//or 1<<4
			//Refits	
			trackCutsFilter->SetRequireITSRefit(kTRUE);
			trackCutsFilter->SetRequireTPCRefit(kTRUE);

			trackCuts->AddCut(trackCutsFilter);
			trackCuts->AddCut(varCutsFilter);

			trackCuts->AddCut(GetPIDCuts(PIDcuts));
			break;
		default:
			std::cout << "No Analysis Track Cut defined" << std::endl;
		}
		return trackCuts;
}

