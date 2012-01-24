
AliAnalysisTaskGammaConvDalitz* AddTaskGammaConvDalitz(
                                                    AliV0Reader *v0Reader,
                                                    Bool_t computeBkg,
                                                    Bool_t standalone = kFALSE
                                                  ){


    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddGammaConvDalitz", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddGammaConvDalitz", "This task requires an input event handler");
        return 0;
    }


	//------------------------- Dalitz task ---------------------------------
	AliGammaConversionHistograms* histogramsD = new AliGammaConversionHistograms();
	AddHistogramsDalitz( histogramsD );
	
	AliAnalysisTaskGammaConvDalitz* dalitzTask = new AliAnalysisTaskGammaConvDalitz("GammaConversionTaskDalitz");
	AliV0Reader* v0ReaderD = v0Reader;
	
	dalitzTask->AdoptHistograms( histogramsD );
	dalitzTask->SetComputeBackground( computeBkg );
	dalitzTask->SetUseBayesPID( kGCUseBayesPID );
	dalitzTask->SetNSigmasElecTPC( kGCNSigmaBelowElecTPCbethe, kGCNSigmaAboveElecTPCbethe);
	dalitzTask->SetNSigmasPionTPC( kGCNSigmaAbovePionTPCbethe );
	dalitzTask->SetNSigmasKaonTPC( kGCNSigmaAboveKaonTPCbethe );
	dalitzTask->SetNSigmasProtonTPC( kGCNSigmaAboveProtonTPCbethe );

	dalitzTask->SetUseESDtrackIndexCut( kGCUseTrackIndexCut );
	dalitzTask->SetUsePsiPairCut( kGCUsePsiPairCut );
	dalitzTask->SetPsiPairCut( kGCPsiPairCut, kGCDeltaPhiCutMin, kGCDeltaPhiCutMax, kGCReadMagFieldSign);
	dalitzTask->SetUseMassCut( kGCUseMassCut );
	dalitzTask->SetUseGammaCut( kGCUseGammaCut );
	dalitzTask->SetUseAliKF( kGCUseAliKF );
	
	dalitzTask->SetTrackSelectionCriteria(kGCTrkSelectionCriteria);
	
	//  track selection
	// -----------------------------
	
	AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("track cuts for pi0 Dalitz decay");
	
	esdTrackCuts->SetEtaRange( -0.9, 0.9 );
	esdTrackCuts->SetAcceptKinkDaughters(kFALSE);

	esdTrackCuts->SetMinNClustersITS(2);
	esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
	esdTrackCuts->SetRequireITSRefit(kTRUE);
	
	//esdTrackCuts->SetRequireSigmaToVertex(kTRUE);
	//esdTrackCuts->SetMaxNsigmaToVertex(3);
	
	esdTrackCuts->SetMaxDCAToVertexXY(1.);
	esdTrackCuts->SetMaxDCAToVertexZ(3.);
	esdTrackCuts->SetDCAToVertex2D(kTRUE);
	
	// ITSsa
	AliESDtrackCuts* esdITSsaTrackCuts = new AliESDtrackCuts(*esdTrackCuts);
	esdITSsaTrackCuts->SetRequireITSStandAlone(kTRUE);
	
	dalitzTask->AdoptITSsaTrackCuts( esdITSsaTrackCuts );
	
	esdTrackCuts->SetMinNClustersTPC(80);
	esdTrackCuts->SetMaxChi2PerClusterTPC(4.);
	esdTrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
	esdTrackCuts->SetRequireTPCRefit(kTRUE);
	
	dalitzTask->AdoptESDtrackCuts( esdTrackCuts );
	
	// --------------------- electrons pid cuts -----------------------------
        esdPidCuts = new AliESDpidCuts("Electrons", "Electron PID cuts");
        //   esdPidCuts->SetTPCnSigmaCut(AliPID::kElectron, 3.);
        esdPidCuts->SetTPCnSigmaCut(AliPID::kElectron, -3., 3.);
        //esdPidCuts->SetTPCnSigmaCut(AliPID::kPion, -1000., 0.);

    dalitzTask->AdoptESDpidCuts( esdPidCuts );

	if( standalone )
	{
		v0ReaderD = new AliV0Reader( *v0Reader ); // Copy constructor
		v0ReaderD->SetHistograms( histogramsD );
		v0ReaderD->SetDoMCTruth( kGCdoMCTruth );    // NOTE: not copied in Copy constructor!!!
		v0ReaderD->SetESDtrackCuts( esdTrackCuts );
		dalitzTask->SetRunStandalone( kTRUE );
	}
	
	dalitzTask->SetV0Reader( v0ReaderD );
	dalitzTask->SetDoMC(kGCdoMCTruth);
	
	dalitzTask->SetDebugLevel(AliLog::kFatal);
	
	// Select triggers events (needs AliPhysicsSelectionTask in the train)
	 dalitzTask->SelectCollisionCandidates();
	
	// Add task to the manager 
	mgr->AddTask( dalitzTask );

    // Connect I/O to the task
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGGA_GammaConversion";
    outputfile += Form("_%s",kGCAnalysisCutSelectionId.Data());
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("histogramsAliGammaConversionDalitz", TList::Class(),AliAnalysisManager::kOutputContainer, outputfile);
    mgr->ConnectInput ( dalitzTask, 0, mgr->GetCommonInputContainer() );
    mgr->ConnectOutput( dalitzTask, 1, coutput3 );

    return dalitzTask;
}

//----------------------------------------------------------------------------------

// Table_MC
Int_t kGCnElementsPi0Table  = 5;
const char* kGCPi0Table[] = {
    "Events",                            // 0
    "\\pi^{0}->\\gamma\\gamma",          // 1
    "\\pi^{0}->e+e-\\gamma",             // 2
    "\\pi^{0}->e+e-\\gamma acceptance",  // 3
    "\\pi^{0}->e+e-\\gamma converted"    // 4
};

//----------------------------------------------------------------------------------
// Table_PID
Int_t kGCnElementsPidTable  = 6;
const char* kGCPidTable[] = { "e", "\\mu", "\\pi", "K", "p", "sum" };

//----------------------------------------------------------------------------------

// Table_ElectronDalitz (ESD cuts)
Int_t kGCnElementsElectronDalitzTable = 8;
const char * kGCElectronDalitzTable[] = {
  "Events",                  // 0
  "e+ candidates",           // 1
  "e- candidates",           // 2
  "e+ Dalitz",               // 3
  "e- Dalitz",               // 4
  "dalitz pairs",            // 5
  "\\gamma \\pi^{0}",        // 6
  "\\pi^{0} Dalitz"          // 7
};

//----------------------------------------------------------------------------------
// Table_ElectronDalitz (ESD cuts)
Int_t kGCnElementsCutsTable  = 6;
const char * kGCCutsTable[] = {
  "Events",                   //0
  "No cut",                   //1
  "Index cut",                //2
  "PsiPair cut",              //3
  "Mass    cut",              //4
  "Angle   cut"               //5
};

//----------------------------------------------------------------------------------
void AddHistogramsDalitz( AliGammaConversionHistograms *histograms )
{
    // Generation
    //
  if(kGCdoMCTruth)
  {
    // pi0
    histograms->AddHistogram("MC_Pi0Dalitz_P", "\\pi^{0} Dalitz decay", kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt/5., "p (GeV/c)", "");
    histograms->AddHistogram("MC_Pi0Dalitz_Pt", "\\pi^{0} Dalitz decay", kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt/5., "p_{t} (GeV/c)", "");
    histograms->AddHistogram("MC_Pi0Dalitz_Eta", "\\pi^{0} Dalitz decay", kGCnXBinsEta, kGCfirstXBinEta, kGClastXBinEta, "\\eta" , "");
    histograms->AddHistogram("MC_Pi0Dalitz_Pt_vs_Y", "\\pi^{0} Dalitz decay", kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt/5., "y","p_{t} (GeV/c)");
    histograms->AddHistogram("MC_Acceptance_Pi0Dalitz_Pt_vs_Y", "\\pi^{0} Acceptance", kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt/5., "y","p_{t} (GeV/c)");
    histograms->AddHistogram("MC_Acceptance_GC_Pi0Dalitz_Pt_vs_Y", "\\pi^{0} Acceptance with converted \\gamma", kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, kGCnXBinsPt, kGCfirstXBinPt, kGClastXBinPt/5., "y","p_{t} (GeV/c)");
   
    // e+ from pi0
    histograms->AddHistogram("MC_EposDalitz_Pt", "e+ from Dalitz pair",              kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_EposDalitz_Eta", "e+ from Dalitz pair",              kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
    
    // e- from pi0
    histograms->AddHistogram("MC_EnegDalitz_Pt", "e- from Dalitz pair",              kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_EnegDalitz_Eta", "e- from Dalitz pair",              kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");

    // gamma from pi0
    histograms->AddHistogram("MC_GammaPi0Dalitz_Pt", "\\gamma from \\pi^{0} Dalitz decay",            kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_GammaPi0Dalitz_Eta", "\\gamma from \\pi^{0} Dalitz decay",            kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");

    // angles
    // e+e-
    histograms->AddHistogram("MC_EposEnegDalitz_Angle", "e+e- Dalitz pair", kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,"Angle(e+,e-) (rad)","");
    // e+e- gamma
    histograms->AddHistogram("MC_EposEnegDalitz_GammaPi0_Angle",  "Angle between the plane e+e- and \\gamma from Pi0 Dalitz decay", kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,"\\theta (rad)","");
    histograms->AddHistogram("MC_EposEnegDalitz_GammaPi0_Angle_vs_P","Angle between the plane e+e- and \\gamma from Pi0 Dalitz decay",kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"\\theta (rad)","p(\\pi^{0}) (GeV/c)");
    histograms->AddHistogram("MC_EposEnegDalitz_GammaPi0_Angle_vs_Pt","Angle between the plane e+e- and \\gamma from Pi0 Dalitz decay",kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"\\theta (rad)","p_{t}(\\pi^{0}) (GeV/c)");

    // Pi0 within acceptance
    histograms->AddHistogram("MC_Acceptance_Pi0Dalitz_Pt", "\\pi^{0} Dalitz decay within acceptance",             kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_Acceptance_Pi0Dalitz_Eta", "\\pi^{0} Dalitz decay within acceptance",             kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
    histograms->AddHistogram("MC_Acceptance_Pi0Dalitz_Radius", "\\pi^{0} Dalitz decay within acceptance",             kGCnXBinsR,kGCfirstXBinR,kGClastXBinR,"Radius(XY) cm","");
    
    histograms->AddHistogram("MC_Acceptance_GammaPi0Dalitz_Pt", "\\gamma from \\pi^{0} Dalitz,  \\pi^{0} within acceptance ",  kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_Acceptance_GammaPi0Dalitz_Eta", "\\gamma from \\pi^{0} Dalitz,  \\pi^{0} within acceptance ",  kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
    
    // e+
    histograms->AddHistogram("MC_Acceptance_EposPi0Dalitz_Pt",  "e+ from \\pi^{0} Dalitz,  \\pi^{0} within acceptance", kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_Acceptance_EposPi0Dalitz_Eta", "e+ from \\pi^{0} Dalitz,  \\pi^{0} within acceptance", kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
    
    // e-
    histograms->AddHistogram("MC_Acceptance_EnegPi0Dalitz_Pt",  "e- from \\pi^{0} Dalitz,  \\pi^{0} within acceptance",              kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_Acceptance_EnegPi0Dalitz_Eta", "e- from \\pi^{0} Dalitz,  \\pi^{0} within acceptance",              kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
    
    
    histograms->AddHistogram("MC_Acceptance_DalitzPair_OpeningAngle", "Dalitz pair within acceptance", kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,"Angle(e+,e-) (rad)","");
    histograms->AddHistogram("MC_GC_GammaPi0Dalitz_All_Z_vs_R",    "\\gamma from \\pi^{0} Dalitz decay", kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "Z (cm)"," R (cm)");

    // e+ within acceptance
    histograms->AddHistogram("MC_Acceptance_EposDalitz_Pt",  "e+ from Dalitz pair within acceptance",              kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_Acceptance_EposDalitz_Eta", "e+ from Dalitz pair within acceptance",              kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
    
    // e-
    histograms->AddHistogram("MC_Acceptance_EnegDalitz_Pt",  "e- from Dalitz pair within acceptance",              kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_Acceptance_EnegDalitz_Eta", "e- from Dalitz pair within acceptance",              kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
    // e+e-
    histograms->AddHistogram("MC_Acceptance_DalitzPair_EposPt_vs_EnegPt","Dalitz pair within acceptance",     kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t}(e+) (GeV/c)","p_{t}(e-) (GeV/c)");

    // within acceptance with gamma converted
    histograms->AddHistogram("MC_Acceptance_GC_Pi0Dalitz_Pt", "\\pi^{0} Dalitz decay within acceptance with \\gamma converted",  kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_Acceptance_GC_Pi0Dalitz_Eta", "\\pi^{0} Dalitz decay within acceptance with \\gamma converted", kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
    
    histograms->AddHistogram("MC_Acceptance_GC_EposPi0Dalitz_Pt", "e+ from Dalitz pair,  \\pi^{0}  within acceptance with \\gamma converted", kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_Acceptance_GC_EposPi0Dalitz_Eta", "e+ from Dalitz pair,  \\pi^{0}  within acceptance with \\gamma converted", kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");

    histograms->AddHistogram("MC_Acceptance_GC_EnegPi0Dalitz_Pt", "e- from Dalitz pair,  \\pi^{0}  within acceptance with \\gamma converted", kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_Acceptance_GC_EnegPi0Dalitz_Eta", "e- from Dalitz pair,  \\pi^{0}  within acceptance with \\gamma converted", kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
    
    histograms->AddHistogram("MC_Acceptance_GC_GammaPi0Dalitz_Pt", "\\gamma from \\pi^{0} Dalitz,  \\pi^{0}  within acceptance with \\gamma converted", kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_Acceptance_GC_GammaPi0Dalitz_Eta", "\\gamma from \\pi^{0} Dalitz,  \\pi^{0}  within acceptance with \\gamma converted", kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
    histograms->AddHistogram("MC_Acceptance_GC_GammaPi0Dalitz_Z_vs_R",  "\\gamma from \\pi^{0} Dalitz decay within acceptance GC", kGCnXBinsZR, kGCfirstXBinZR, kGClastXBinZR, kGCnYBinsZR, kGCfirstYBinZR, kGClastYBinZR, "Z (cm)","R (cm)");
  }
    //  Reconstruction
    // ---------------------------------------------------------------------------

    // reconstructed epos
    histograms->AddHistogram("ESD_Debug_Table","Table Debug",10,-5,5,"Quality","");
    histograms->AddHistogram("ESD_EposCandidates_Pt", "Reconstructed e+ candidates",   kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("ESD_EposCandidates_Eta", "Reconstructed e+ candidates",       kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");

    // reconstructed eneg
    histograms->AddHistogram("ESD_EnegCandidates_Pt", "Reconstructed e- candidates",   kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("ESD_EnegCandidates_Eta", "Reconstructed e- candidates",       kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");

    // Reconstructed dalitz pair
    histograms->AddHistogram("ESD_DalitzPairCandidates_Pt", "Reconstructed Dalitz pair candidates",   kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("ESD_DalitzPairCandidates_Mass", "Reconstructed Dalitz pair candidates", kGCnXBinsDalitzMass, kGCfirstXBinDalitzMass, kGClastXBinDalitzMass, "M(e+,e-) (GeV/c^{2})","");
    
    // gammas
    histograms->AddHistogram("ESD_GammaCandidates_Pt", "Reconstructed \\gamma candidates",   kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("ESD_GammaCandidates_Eta", "Reconstructed \\gamma candidates",       kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");

    histograms->AddHistogram("ESD_EposEneg_GammaCandidates_Angle", "Angle between e+, e- candidates and \\gamma candidates", kGCnXBinsEPosNegAngle,kGCfirstXBinEPosNegAngle,kGClastXBinEPosNegAngle,"Angle(e+e-, \\gamma) (rad)","");

    // Reconstructed pi0
    //
    histograms->AddHistogram("ESD_Pi0_P", "\\pi^{0} Dalitz decay", kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p (GeV/c)","");
    histograms->AddHistogram("ESD_Pi0_Pt", "\\pi^{0} Dalitz decay",               kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("ESD_Pi0_Eta", "\\pi^{0} Dalitz decay",               kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
     histograms->AddHistogram("ESD_Pi0_Y",        "\\pi^{0} Dalitz decay",               kGCnXBinsRapid,kGCfirstXBinRapid,kGClastXBinRapid,"y","");
    histograms->AddHistogram("ESD_Pi0_Phi",          "\\pi^{0} Dalitz decay",               kGCnXBinsPhi,kGCfirstXBinPhi,kGClastXBinPhi," \\phi (rad)","");
    histograms->AddHistogram("ESD_Pi0_Mass",      "\\pi^{0} Dalitz decay",                kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"M(e+,e-,\\gamma) (GeV/c^{2})","");
    histograms->AddHistogram("ESD_Pi0_Pt_vs_Y","\\pi^{0} Dalitz decay", kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"y","p_{t} (GeV/c)");
    histograms->AddHistogram("ESD_Pi0_Mass_vs_Pt","\\pi^{0} Dalitz decay",               kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5., kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"p_{t} (GeV/c)","M(e+,e-,\\gamma) (GeV/c^{2})");
    histograms->AddHistogram("ESD_Pi0_Mass_vs_Eta","\\pi^{0} Dalitz decay",               kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta, kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"\\eta","M(e+,e-,\\gamma) (GeV/c^{2})");
    histograms->AddHistogram("ESD_Pi0_Mass_vs_Y","\\pi^{0} Dalitz decay", kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid,  kGCnXBinsPi0DalitzMass,  kGCfirstXBinPi0DalitzMass,  kGClastXBinPi0DalitzMass, "y","M(e+,e-,\\gamma) (GeV/c^{2})");

    // Background
    histograms->AddHistogram("ESD_BKG_EposEneg", "Like sign background",              kGCnXBinsDalitzMass, kGCfirstXBinDalitzMass, kGClastXBinDalitzMass,"M(e+,e+) & M(e-,e-) (GeV/c^{2})","");
    histograms->AddHistogram("ESD_BKG_PrevGamma", "e+e- with \\gamma used in the signal",              kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"M(e+,e-,\\gamma) (GeV/c^{2})","");
    histograms->AddHistogram("ESD_BKG_PrevGamma_vs_Pt", "Backgroud: \\gamma used in the signal",  kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,  kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass, "p_{t} (GeV/c)", "M(e+,e-,\\gamma) (GeV/c^{2})");

    histograms->AddHistogram("ESD_BKG_BGHandler", "e+e- with \\gamma from a pool of events",              kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"M(e+,e-,\\gamma) (GeV/c^{2})","");
    histograms->AddHistogram("ESD_BKG_BGHandler_vs_Pt", "e+e- with \\gamma from a pool of events",  kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,      kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"p_{t} (GeV/c)","M(e+,e-,\\gamma) (GeV/c^{2})");

    histograms->AddHistogram("ESD_BKG_Electron", "e+ with e-,\\gamma from a pool of events",              kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"M(e+,e-,\\gamma) (GeV/c^{2})","");
    histograms->AddHistogram("ESD_BKG_Electron_vs_Pt", "e+ with e-,\\gamma from a pool of events", kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,       kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"p_{t} (GeV/c)","M(e+,e-,\\gamma) (GeV/c^{2})");


    // GEANT
    // ---------------------
  if(kGCdoMCTruth)
  {
    histograms->AddHistogram("MC_ESD_EposDalitz_Pt", "Reconstructed e+ from Dalitz pair (MC)",   kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_ESD_EposDalitz_Eta", "Reconstructed e+ from Dalitz pair (MC)",       kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");

    // reconstructed eneg
    histograms->AddHistogram("MC_ESD_EnegDalitz_Pt", "Reconstructed e- from Dalitz pair (MC)",   kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_ESD_EnegDalitz_Eta", "Reconstructed e- from Dalitz pair (MC)",       kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");

    // Reconstructed dalitz pair
    histograms->AddHistogram("MC_ESD_DalitzPair_Pt", "Reconstructed Dalitz pair (MC)",   kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_ESD_DalitzPair_Mass", "Reconstructed Dalitz pair (MC)",              kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"M(e+,e-) (GeV/c^{2})","");

    // reconstructed gammas
    histograms->AddHistogram("MC_ESD_Gamma_Pt", "Reconstructed \\gamma (MC)",   kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_ESD_Gamma_Eta", "Reconstructed \\gamma (MC)",       kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");

    // gamma from Dalitz
    histograms->AddHistogram("MC_ESD_GammaPi0Dalitz_Pt", "Reconstructed \\gamma from \\pi^{0} Dalitz (MC)",   kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_ESD_GammaPi0Dalitz_Eta", "Reconstructed \\gamma from \\pi^{0} Dalitz (MC)",       kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");

    // Reconstructed pi0 (MC)
    histograms->AddHistogram("MC_ESD_Pi0_P", "\\pi^{0} Dalitz decay (MC)",               kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p (GeV/c)","");
    histograms->AddHistogram("MC_ESD_Pi0_Pt", "\\pi^{0} Dalitz decay (MC)",               kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"p_{t} (GeV/c)","");
    histograms->AddHistogram("MC_ESD_Pi0_Eta", "\\pi^{0} Dalitz decay (MC)",               kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta,"\\eta","");
     histograms->AddHistogram("MC_ESD_Pi0_Y",        "\\pi^{0} Dalitz decay (MC)",               kGCnXBinsRapid,kGCfirstXBinRapid,kGClastXBinRapid,"y","");
    histograms->AddHistogram("MC_ESD_Pi0_Phi",          "\\pi^{0} Dalitz decay (MC)",               kGCnXBinsPhi,kGCfirstXBinPhi,kGClastXBinPhi," \\phi (rad)","");
    histograms->AddHistogram("MC_ESD_Pi0_Mass",      "\\pi^{0} Dalitz decay (MC)",                kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"M(e+,e-,\\gamma) (GeV/c^{2})","");
    
    histograms->AddHistogram("MC_ESD_Pi0_Pt_vs_Y","\\pi^{0} Dalitz decay (MC)", kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid, kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5.,"y","p_{t} (GeV/c)");
    histograms->AddHistogram("MC_ESD_Pi0_Mass_vs_Pt","\\pi^{0} Dalitz decay (MC)",               kGCnXBinsPt,kGCfirstXBinPt,kGClastXBinPt/5., kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"p_{t} (GeV/c)","M(e+,e-,\\gamma) (GeV/c^{2})");
    histograms->AddHistogram("MC_ESD_Pi0_Mass_vs_Eta","\\pi^{0} Dalitz decay (MC)",               kGCnXBinsEta,kGCfirstXBinEta,kGClastXBinEta, kGCnXBinsPi0DalitzMass, kGCfirstXBinPi0DalitzMass, kGClastXBinPi0DalitzMass,"\\eta","M(e+,e-,\\gamma) (GeV/c^{2})");
    histograms->AddHistogram("MC_ESD_Pi0_Mass_vs_Y","\\pi^{0} Dalitz decay (MC)", kGCnXBinsRapid, kGCfirstXBinRapid, kGClastXBinRapid,  kGCnXBinsPi0DalitzMass,  kGCfirstXBinPi0DalitzMass,  kGClastXBinPi0DalitzMass, "y","M(e+,e-,\\gamma) (GeV/c^{2})");
  }
  
  //if(kGCTrkSelectionCriteria != 1)
    // ITS standalone
    histograms->AddHistogram("ESD_ITSsa_dEdx_vs_P",          "ITS standalone tracks",    kGCnXBinsP*5, kGCfirstXBinP, kGClastXBinP/5.,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"p (GeV/c)","ITS Signal (a.u.)");
    histograms->AddHistogram("ESD_ITSsa_PidCut_dEdx_vs_P", "ITS standalone tracks after PID",     kGCnXBinsP*5, kGCfirstXBinP, kGClastXBinP/5.,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"p (GeV/c)","ITS Signal (a.u.)");
  
  if(kGCdoMCTruth)
  {
    histograms->AddHistogram("MC_ESD_ITSsa_PidCut_dEdx_vs_P", "ITS standalone tracks after PID (MC)",     kGCnXBinsP*5, kGCfirstXBinP, kGClastXBinP/5.,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"p (GeV/c)","ITS Signal (a.u.)");
    histograms->AddHistogram("MC_ESD_ITSsa_Electron_dEdx_vs_P",       "ITS standalone electrons",    kGCnXBinsP*5, kGCfirstXBinP, kGClastXBinP/5.,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"p (GeV/c)","ITS Signal (a.u.)");
  }
    // TPC dEdx
    histograms->AddHistogram("ESD_TPC_dEdx_vs_P" ,  "Global tracks",               kGCnXBinsP*5, kGCfirstXBinP, kGClastXBinP/5.,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"p (GeV/c)","TPC Signal (a.u.)");
    histograms->AddHistogram("ESD_TPC_PidCut_dEdx_vs_P",   "Global tracks after PID",     kGCnXBinsP*5, kGCfirstXBinP, kGClastXBinP/5.,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"p (GeV/c)","TPC Signal (a.u.)");
    
  if(kGCdoMCTruth)
  {
    histograms->AddHistogram("MC_ESD_TPC_Electron_dEdx_vs_P" ,        "Electron global tracks",     kGCnXBinsP*5, kGCfirstXBinP, kGClastXBinP/5.,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"p (GeV/c)","TPC Signal (a.u.)");
    histograms->AddHistogram("MC_ESD_TPC_PidCut_dEdx_vs_P",   "Global tracks after PID (MC)",     kGCnXBinsP*5, kGCfirstXBinP, kGClastXBinP/5.,kGCnYBinsdEdx, kGCfirstYBindEdx, kGClastYBindEdx,"p (GeV/c)","TPC Signal (a.u.)");
  }

// psi pair
     histograms->AddHistogram("ESD_EposEneg_PsiPair_vs_DPhi", "e+e- candidates", 2*kGCnXBinsOpeningAngle, -1*kGClastXBinOpeningAngle, kGClastXBinOpeningAngle, 2*kGCnXBinsOpeningAngle, -1*kGClastXBinOpeningAngle, kGClastXBinOpeningAngle, "\\Delta\\phi (rad)","\\psi_{pair}");
  
  if(kGCdoMCTruth)
  {
     histograms->AddHistogram("MC_ESD_DalitzPair_PsiPair_vs_DPhi", "Dalitz pairs (MC)", 2*kGCnXBinsOpeningAngle, -1*kGClastXBinOpeningAngle, kGClastXBinOpeningAngle, 2*kGCnXBinsOpeningAngle, -1*kGClastXBinOpeningAngle, kGClastXBinOpeningAngle, "\\Delta\\phi (rad)","\\psi_{pair}");

     histograms->AddHistogram("MC_ESD_EposEnegGamma_PsiPair_vs_DPhi", "e+e- from \\gamma conversion (MC)", 2*kGCnXBinsOpeningAngle, -1*kGClastXBinOpeningAngle, kGClastXBinOpeningAngle, 2*kGCnXBinsOpeningAngle, -1*kGClastXBinOpeningAngle, kGClastXBinOpeningAngle, "\\Delta\\phi (rad)","\\psi_{pair}");
  }

    // Tables
  if(kGCdoMCTruth)
  {
    histograms->AddTable("Table_Generation","Generation",kGCnElementsPi0Table,kGCPi0Table);
    histograms->AddTable("Table_PID","particle table (Row: PID, Col: MC)", kGCnElementsPidTable, kGCPidTable, kGCnElementsPidTable, kGCPidTable);
  }
    histograms->AddTable("Table_Reconstruction","Reconstruction",kGCnElementsElectronDalitzTable,kGCElectronDalitzTable);
    histograms->AddTable("Table_Cuts","Number of \\pi^{0}->e+e-\\gamma after cuts",kGCnElementsCutsTable,kGCCutsTable);
}


 /*
        // NOTE: We should use fESDpid from fV0Reader (right now it is private)
        // Better parameters for MonteCarlo from A. Kalweit 2010/01/8


              fESDpid->GetTPCResponse().SetBetheBlochParameters(  4.23232575531564326e+00/50.,
                                                                  8.68482806165147636e+00,
                                                                  1.34000000000000005e-05,
                                                                  2.30445734159456084e+00,
                                                                  2.25624744086878559e+00);//Trujillo

                
            fESDpid->GetTPCResponse().SetBetheBlochParameters(   2.49397e+00,
                                                                  1.46884e+01,
                                                                  7.44848e-04,
                                                                  2.24803e+00,
                                                                  4.85791e+00); //AliRoot v4-18-Rev-10 //Dalitz Generation*/


/*        fESDpid->GetTPCResponse().SetBetheBlochParameters(  2.15898e+00/50.,
                                                                  1.75295e+01,
                                                                  3.40030e-09,
                                                                  1.96178e+00,
                                                                  3.91720e+00 );

                //MC-CERN
        

    } 
        else { 
        // Better parameters for data from A. Kalweit 2010/01/8
        fESDpid->GetTPCResponse().SetBetheBlochParameters( 0.0283086,
                                                           2.63394e+01,
                                                           5.04114e-11,
                                                           2.12543e+00,
                                                           4.88663e+00 );
    */
