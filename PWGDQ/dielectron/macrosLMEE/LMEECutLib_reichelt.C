class LMEECutLib {
  
public:
	static  enum LMMECutSet {

    kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight,
    kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight,
    kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight,
    kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight,
    kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4,
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4,
    kPbPb2011_pidITSTPC_trkSPDfirst_3,            // (cutset w/o pairing)
    kPbPb2011_pidTPC_trkSPDfirst_3,               // (cutset w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose, // (cutset w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose, // (cutset w/o pairing)
    kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose,      // (cutset w/o pairing)
    kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose,      // (cutset w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1,
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1,       // Cutset for Technical Preliminaries for QM2014 (no prefilter used!)
    kPbPb2011_pidTPCTOF_trkSPDorSDD_1,
    kPbPb2011_pidTPCTOF_trkSPDfirst_1,
    kPbPb2011PID_ITSTPCTOFif2,        // (NO FULL CUTSET)
    kPbPb2011PID_TPCTOF3,             // (NO FULL CUTSET)
    kPbPb2011TRK_SDDfirstSPDnone,     // (NO FULL CUTSET) complimentary tracks, strictly without SPD, to be combined with others!
    kPbPb2011TRK_SPDfirst,            // (NO FULL CUTSET) main track selection, with SPD first
    kPbPb2011TRK_SDDfirstSPDnone4cls, // (NO FULL CUTSET) complimentary tracks, strictly without SPD, to be combined with others!
    kPbPb2011TRK_SPDfirst5cls,        // (NO FULL CUTSET) main track selection, with SPD first
    kPbPb2011_TPCITS_TOFif1,
    kPbPb2011_TPCTOF_Semi2, // changed PairCutsAna from PhiV to OpeningAngle. prefilter cuts renewed (if applicable)
    // following cutsets are not complete anymore!
    kPbPb2011_TPCTOF_Semi1, // old prefilter cuts (leg & pair), some are confusing
    kPbPb2011NoPID, // pairing disabled in config
    kPbPb2011TPCandTOF, // this was the final one activated by Christoph!
    kPbPb2011TPCandTOFHPT,
		kPbPb2011TPC, //TOF required, more relaxed cut on TPC
		kPbPb2011TPCandTOFwide, //TOF required, more relaxed cut on TPC
		kPbPb2011TPCorTOF,
		kpp2010TPCandTOF,
		kpp2010TPCorTOF,
		kCUTSETMAX
	};
  
	static  enum LMMECentSel {
	  kPbPb2011Central,
    kPbPb2011MidCentral,
		kPbPb2011SemiCentral,
		kPbPb2011Peripheral,
		kCENTSELMAX
	};
  
	//char* LMEECutNames[kCUTSETMAX] = { "PbPb2011TPCandTOF","PbPb2011TPCorTOF"};
  
  
	LMEECutLib() {}
  
	AliDielectronEventCuts*     GetEventCuts(Int_t cutSet);
	AliAnalysisCuts*            GetCentralityCuts(Int_t centSel);
	AliDielectronTrackRotator*  GetTrackRotator(Int_t cutSet);
	AliDielectronMixingHandler* GetMixingHandler(Int_t cutSet);
  
	AliAnalysisCuts* GetPairCutsAna(Int_t cutSet, Bool_t tooglePC=kFALSE);
	AliAnalysisCuts* GetPairCutsPre(Int_t cutSet);  
  
	AliAnalysisCuts* GetPIDCutsAna(Int_t cutSet);  
	AliAnalysisCuts* GetPIDCutsPre(Int_t cutSet);  
  
	AliAnalysisCuts* GetTrackCutsAna(Int_t cutSet);  
	AliAnalysisCuts* GetTrackCutsPre(Int_t cutSet);  
	AliAnalysisCuts* GetESDTrackCutsAna(Int_t cutSet);  
  
  
  // Note: event cuts are identical for all analysis 'cutDefinition's that run together!
  // the selection is hardcoded in the AddTask, currently to 'kPbPb2011_TPCTOF_Semi1'
	AliDielectronEventCuts* GetEventCuts(Int_t cutSet) {
	  AliDielectronEventCuts* eventCuts = 0x0;
	  switch (cutSet) {
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      case kPbPb2011_pidITSTPC_trkSPDfirst_3:
      case kPbPb2011_pidTPC_trkSPDfirst_3:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      case kPbPb2011_TPCITS_TOFif1:
      case kPbPb2011_TPCTOF_Semi2:
      case kPbPb2011_TPCTOF_Semi1:
      case kPbPb2011NoPID:
      case kPbPb2011TPCandTOF :
   		case kPbPb2011TPCandTOFHPT:
      case kPbPb2011TPC :
      case kPbPb2011TPCandTOFwide :
      case kPbPb2011TPCorTOF  :
      case kpp2010TPCandTOF :
      case kpp2010TPCorTOF  :
        //Basic Event Cuts for pp and Pb-Pb, additional cuts may be in the AddTask
        eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
        eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
        //eventCuts->SetVertexType(AliDielectronEventCuts::kVtxTPC); // AOD
        //           eventCuts->SetCentralityRange(0.0,80.0);
        eventCuts->SetRequireVertex();
        eventCuts->SetMinVtxContributors(1);
        eventCuts->SetVertexZ(-10.,10.);
        break;
      default: cout << "No Event Cut defined" << endl;
	  }
	  return eventCuts;
	}
  
  
	//Selection of relatively 'flat' centralities
	AliAnalysisCuts* GetCentralityCuts(Int_t centSel) {
	  AliDielectronVarCuts* centCuts = 0x0;
	  switch (centSel) {
      case kPbPb2011Central:
        centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb2011Central");
        centCuts->AddCut(AliDielectronVarManager::kCentrality,0.,10.);
        break;
      case kPbPb2011MidCentral:
        centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb2011MidCentral");
        centCuts->AddCut(AliDielectronVarManager::kCentrality,10.,20.);
        break;
      case kPbPb2011SemiCentral:
        centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb2011SemiCentral");
        centCuts->AddCut(AliDielectronVarManager::kCentrality,20.,50.);
        break;
      case kPbPb2011Peripheral:
        centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb2011Peripheral");
        centCuts->AddCut(AliDielectronVarManager::kCentrality,50.,90.);
        break;
			default: cout << "No Centrality selected" << endl;
	  }
	  return centCuts;
	}
  
  
	//Basic track rotator settings from J/Psi, more investigation needed
	AliDielectronTrackRotator* GetTrackRotator(Int_t cutSet) {
	  AliDielectronTrackRotator* trackRotator = 0x0;
	  switch (cutSet) {
      case kPbPb2011_pidITSTPC_trkSPDfirst_3:
      case kPbPb2011_pidTPC_trkSPDfirst_3:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      case kPbPb2011_TPCITS_TOFif1:
      case kPbPb2011_TPCTOF_Semi2:
      case kPbPb2011_TPCTOF_Semi1:
      case kPbPb2011NoPID:
      case kPbPb2011TPCandTOF :
   		case kPbPb2011TPCandTOFHPT:
      case kPbPb2011TPC :
      case kPbPb2011TPCandTOFwide :
      case kPbPb2011TPCorTOF  :
      case kpp2010TPCandTOF :
      case kpp2010TPCorTOF  :
        trackRotator = new AliDielectronTrackRotator();
        trackRotator->SetIterations(20);
        trackRotator->SetConeAnglePhi(TMath::Pi()/180*165);
        trackRotator->SetStartAnglePhi(TMath::Pi());
        break;
      default: cout << "No Rotator defined" << endl;
	  }
	  return trackRotator;
	}
  
  
	AliDielectronMixingHandler* GetMixingHandler(Int_t cutSet) {
	  AliDielectronMixingHandler* mixingHandler = 0x0;
	  switch (cutSet) {
      case kPbPb2011_pidITSTPC_trkSPDfirst_3:
      case kPbPb2011_pidTPC_trkSPDfirst_3:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      case kPbPb2011_TPCITS_TOFif1:
      case kPbPb2011_TPCTOF_Semi2:
      case kPbPb2011_TPCTOF_Semi1:
        mixingHandler = new AliDielectronMixingHandler;
        mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
        mixingHandler->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,30,50,80");
        // now using TPC event plane, uncorrected. (also, the old phi range was wrong, now same effective binning.)
        mixingHandler->AddVariable(AliDielectronVarManager::kTPCrpH2uc, 6, TMath::Pi()/-2., TMath::Pi()/2.);
        mixingHandler->SetDepth(15);
        mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
        break;
        //[...]
      default: cout << "No Mixer defined" << endl;
	  }
	  return mixingHandler;
	}
  
  
  
	//Pair Cuts for Analysis step - take care of logic - inverted compared to other PairCuts!!
  // cuts = SELECTION!!!
	AliAnalysisCuts* GetPairCutsAna(Int_t cutSet, Bool_t togglePC /*=kFALSE*/)  {
    cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
    AliAnalysisCuts* pairCuts=0x0;
    switch (cutSet) {
      case kPbPb2011_pidITSTPC_trkSPDfirst_3:
      case kPbPb2011_pidTPC_trkSPDfirst_3:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
        cout << "No Pair Cuts used - ok " << endl; // since 18.02.2014
        break;
        
      case kPbPb2011_TPCITS_TOFif1:
      case kPbPb2011_TPCTOF_Semi2:
        //        AliDielectronVarCuts* pairCutsPhivGood =new AliDielectronVarCuts("pairCutsPhivGood","pairCutsPhivGood");
        //        pairCutsPhivGood->AddCut(AliDielectronVarManager::kPhivPair, 0.0, 2.0); 
        AliDielectronVarCuts* pairCutsOpAngGood =new AliDielectronVarCuts("pairCutsOpAngGood","pairCutsOpAngGood");
        pairCutsOpAngGood->AddCut(AliDielectronVarManager::kOpeningAngle, 0.05, 999.); // in upgrade: 0.05
        AliDielectronVarCuts* pairCutsInvM =new AliDielectronVarCuts("pairCutsInvM","pairCutsInvM");
        pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.02); // in upgrade: 0.01
        AliDielectronVarCuts* pairCutsInvMgood =new AliDielectronVarCuts("pairCutsInvMgood","pairCutsInvMgood");
        pairCutsInvMgood->AddCut(AliDielectronVarManager::kM, 0.02, 99999.);
        
        AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
        pairCutsCG->AddCut(pairCutsInvM);
        pairCutsCG->AddCut(pairCutsOpAngGood);
        //        pairCutsCG->AddCut(pairCutsPhivGood);
        
        AliDielectronCutGroup* pairCutsCG2 =new AliDielectronCutGroup("pairCutsCG2","pairCutsCG2",AliDielectronCutGroup::kCompOR);
        pairCutsCG2->AddCut(pairCutsInvMgood);
        pairCutsCG2->AddCut(pairCutsCG);
        pairCuts = pairCutsCG2;
        break;
        
      case kPbPb2011_TPCTOF_Semi1:
        //[...] // PhiV and InvMass
		  default: cout << "No Pair Cuts defined " << endl;
    }
    return pairCuts;
	}
  
  
	//Pair Cuts for PREFILTER step
  // cuts = REJECTION!!!
	AliAnalysisCuts* GetPairCutsPre(Int_t cutSet)  {  
    cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
    AliAnalysisCuts* pairCuts=0x0;
    switch (cutSet) {
      case kPbPb2011_pidITSTPC_trkSPDfirst_3:
      case kPbPb2011_pidTPC_trkSPDfirst_3:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      case kPbPb2011_TPCITS_TOFif1:
      case kPbPb2011_TPCTOF_Semi2:
        AliDielectronVarCuts* pairCutsInvM =new AliDielectronVarCuts("pairCutsInvM","pairCutsInvM");
        pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.02); // in upgrade: 0.01
        AliDielectronVarCuts* pairCutsOpAng =new AliDielectronVarCuts("pairCutsOpAng","pairCutsOpAng");
        pairCutsOpAng->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.05); // in upgrade: 0.05
        
        AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
        pairCutsCG->AddCut(pairCutsInvM);
        pairCutsCG->AddCut(pairCutsOpAng);
        //pairCutsCG->AddCut(pairCutsPhiv);
        pairCuts = pairCutsCG;
        break;
        
      case kPbPb2011_TPCTOF_Semi1:
        //[...] // PhiV and InvMass
      default: cout << "No Prefilter Pair Cuts defined " << endl;
	  } 
	  return pairCuts;
	}
  
  
  
	AliAnalysisCuts* GetPIDCutsAna(Int_t cutSet) {
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
    //TPC: electron inclusion asymmetric
	  //     pion     exclusion 3sigma
	  //TOF: electron inclusion 3sigma in region where p,K cross electrons in TPC
	  AliDielectronPID *pidTPCTOF_Semi1 = new AliDielectronPID("pidTPCTOF_Semi1","pidTPCTOF_Semi1");
	  pidTPCTOF_Semi1->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3., 0. ,100., kFALSE);
	  pidTPCTOF_Semi1->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
	  pidTPCTOF_Semi1->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0. ,1.7 , kFALSE);
    //
    //
    // LOOSE PID TPC+TOF
    AliDielectronPID *pidTPCTOF_Semi_LOOSE = new AliDielectronPID("pidTPCTOF_Semi_LOOSE","pidTPCTOF_Semi_LOOSE");
	  pidTPCTOF_Semi_LOOSE->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
	  pidTPCTOF_Semi_LOOSE->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,1.7 , kFALSE);
    //
    //
    // PID TPC only
	  AliDielectronPID *pidTPC_3 = new AliDielectronPID("pidTPC_3","pidTPC_3");
	  pidTPC_3->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3., 0. ,100., kFALSE);
	  pidTPC_3->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
    
    
    //TPC: electron inclusion asymmetric
	  //     pion     exclusion 3sigma
	  //ITS: electron inclusion asymmetric in region where p,K cross electrons in TPC
	  //TOF: electron inclusion 3sigma in similar region - BUT ONLY IF AVAILABLE
	  AliDielectronPID *pidTPCITS_TOFif1 = new AliDielectronPID("pidTPCITS_TOFif1","pidTPCITS_TOFif1");
	  pidTPCITS_TOFif1->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3. , 0. ,100., kFALSE);
	  pidTPCITS_TOFif1->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
	  pidTPCITS_TOFif1->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 1. , 0. ,1.5 , kFALSE);
	  pidTPCITS_TOFif1->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,1.7 , kFALSE, AliDielectronPID::kIfAvailable);
    //
    //
    //TPC: electron inclusion asymmetric
	  //     pion     exclusion 3sigma
	  //ITS: electron inclusion asymmetric OVER FULL MOMENTUM RANGE
	  //TOF: electron inclusion 3sigma - BUT ONLY IF AVAILABLE
	  AliDielectronPID *pidTPCITS_TOFif2 = new AliDielectronPID("pidTPCITS_TOFif2","pidTPCITS_TOFif2");
	  pidTPCITS_TOFif2->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3. , 0. ,100., kFALSE);
	  pidTPCITS_TOFif2->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
	  pidTPCITS_TOFif2->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 1. , 0. ,100., kFALSE);
	  pidTPCITS_TOFif2->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
    //
    //
    // LOOSE PID ITS+TPC+TOFif
    AliDielectronPID *pidTPCITS_TOFif_LOOSE = new AliDielectronPID("pidTPCITS_TOFif_LOOSE","pidTPCITS_TOFif_LOOSE");
	  pidTPCITS_TOFif_LOOSE->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
	  pidTPCITS_TOFif_LOOSE->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-10. ,20. , 0. ,100., kFALSE);
	  pidTPCITS_TOFif_LOOSE->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
    //
    //
    // PID ITS+TPC
    AliDielectronPID *pidTPCITS_3 = new AliDielectronPID("pidTPCITS_3","pidTPCITS_3");
	  pidTPCITS_3->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3. , 0. ,100., kFALSE);
	  pidTPCITS_3->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
	  pidTPCITS_3->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 1. , 0. ,100., kFALSE);
    //
    //
    // tighter PID ITS+TPC+TOFif
    // ITS only up to momentum where proton contamination is seen in TPC signal
    AliDielectronPID *pidTPCITS_TOFif56 = new AliDielectronPID("pidTPCITS_TOFif56","pidTPCITS_TOFif56");
	  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 2.5, 0. ,100., kFALSE);
	  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
	  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 0.5, 0. ,  2., kFALSE);
	  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
    
    
	  //Apply ITS cuts (see Hongyan's talks):
	  //3 sigma inclusion of electrons in TPC
	  //3 sigma exclusion of pions in TPC
	  //3 sigma inclusion of electrons in ITS,TOF for p<1.5GeV, where p,K contamination
	  //TOF only IF available!
	  AliDielectronPID *pidTPCandITSTOF = new AliDielectronPID("pidTPCandITSTOF","pidTPCandITSTOF");//"TPC-TOF-ITS"
	  pidTPCandITSTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.,3., 0.,100., kFALSE);
	  pidTPCandITSTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3.,3., 0.,100., kTRUE);
	  pidTPCandITSTOF->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -6.,3., 0.,1.5, kFALSE );
	  pidTPCandITSTOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.,3., 0.,1.5, kFALSE,AliDielectronPID::kIfAvailable );
    
    // eta range:
    AliDielectronVarCuts *etaRange090 = new AliDielectronVarCuts("etaRange090","etaRange090");
    etaRange090->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
    AliDielectronVarCuts *etaRange084 = new AliDielectronVarCuts("etaRange084","etaRange084");
    etaRange084->AddCut(AliDielectronVarManager::kEta, -0.84, 0.84);
    AliDielectronVarCuts *etaRange076 = new AliDielectronVarCuts("etaRange076","etaRange076");
    etaRange076->AddCut(AliDielectronVarManager::kEta, -0.76, 0.76);
    // pt range:
    AliDielectronVarCuts *ptRange400to3500 = new AliDielectronVarCuts("ptRange400to3500","ptRange400to3500");
    ptRange400to3500->AddCut(AliDielectronVarManager::kPt, .4, 3.5);
    
    
	  //-----------------------------------------------
	  // Now see what Config actually loads and assemble final cuts
	  //-----------------------------------------------
    switch (cutSet) {
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight: // tighter "ITSTPCTOFif" PID
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
        AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
        cgPIDCutsAna->AddCut(etaRange076);
        cgPIDCutsAna->AddCut(ptRange400to3500);
        cgPIDCutsAna->AddCut(pidTPCITS_TOFif56);
        cgPIDCutsAna->AddCut(GetTrackCutsAna(cutSet));
        pidCuts = cgPIDCutsAna;
        break;
    
      case kPbPb2011_pidITSTPC_trkSPDfirst_3:
        AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
        cgPIDCutsAna->AddCut(etaRange076);
        cgPIDCutsAna->AddCut(ptRange400to3500);
        cgPIDCutsAna->AddCut(pidTPCITS_3);
        cgPIDCutsAna->AddCut(GetTrackCutsAna(cutSet));
        pidCuts = cgPIDCutsAna;
        break;
      case kPbPb2011_pidTPC_trkSPDfirst_3:
        AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
        cgPIDCutsAna->AddCut(etaRange090);
        cgPIDCutsAna->AddCut(ptRange400to3500);
        cgPIDCutsAna->AddCut(pidTPC_3);
        cgPIDCutsAna->AddCut(GetTrackCutsAna(cutSet));
        pidCuts = cgPIDCutsAna;
        break;
      
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose: // loose "ITSTPCTOFif" PID - for 2D contamination study
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
        AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
        cgPIDCutsAna->AddCut(etaRange076);
        cgPIDCutsAna->AddCut(ptRange400to3500);
        cgPIDCutsAna->AddCut(pidTPCITS_TOFif_LOOSE);
        cgPIDCutsAna->AddCut(GetTrackCutsAna(cutSet));
        pidCuts = cgPIDCutsAna;
        break;
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
        AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
        cgPIDCutsAna->AddCut(etaRange090);
        cgPIDCutsAna->AddCut(ptRange400to3500);
        cgPIDCutsAna->AddCut(pidTPCTOF_Semi_LOOSE);
        cgPIDCutsAna->AddCut(GetTrackCutsAna(cutSet));
        pidCuts = cgPIDCutsAna;
        break;
      
      case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4: // regular "ITSTPCTOFif" PID
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
        AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
        cgPIDCutsAna->AddCut(etaRange076);
        cgPIDCutsAna->AddCut(ptRange400to3500);
        cgPIDCutsAna->AddCut(pidTPCITS_TOFif2);
        cgPIDCutsAna->AddCut(GetTrackCutsAna(cutSet));
        pidCuts = cgPIDCutsAna;
        break;
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
        AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
        cgPIDCutsAna->AddCut(etaRange090);
        cgPIDCutsAna->AddCut(ptRange400to3500);
        cgPIDCutsAna->AddCut(pidTPCTOF_Semi1);
        cgPIDCutsAna->AddCut(GetTrackCutsAna(cutSet));
        pidCuts = cgPIDCutsAna;
        break;
        
      case kPbPb2011_TPCITS_TOFif1:
        AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
        cgPIDCutsAna->AddCut(etaRange084); // was 0.84 -> not ideal
        cgPIDCutsAna->AddCut(ptRange400to3500);
        cgPIDCutsAna->AddCut(pidTPCITS_TOFif1);
        cgPIDCutsAna->AddCut(GetTrackCutsAna(cutSet));
        pidCuts = cgPIDCutsAna;
        break;
      case kPbPb2011_TPCTOF_Semi2:
      case kPbPb2011_TPCTOF_Semi1:
        AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
        cgPIDCutsAna->AddCut(etaRange084);
        cgPIDCutsAna->AddCut(ptRange400to3500);
        cgPIDCutsAna->AddCut(pidTPCTOF_Semi1);
        cgPIDCutsAna->AddCut(GetTrackCutsAna(cutSet)); // for 'kPbPb2011_TPCTOF_Semi1', this was called in the Config
        pidCuts = cgPIDCutsAna;
        break;
        //[...]
   		case kPbPb2011TPCandTOFHPT:
        //test Hongyan's cut
        AliDielectronCutGroup* cgSecondTrackFilterPIDTPC1 = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);
        cgSecondTrackFilterPIDTPC1->AddCut(etaRange084);
        cgSecondTrackFilterPIDTPC1->AddCut(ptRange400to3500);
        cgSecondTrackFilterPIDTPC1->AddCut(pidTPCandITSTOF);
        pidCuts = cgSecondTrackFilterPIDTPC1;
        break;
        //[...]
      default: cout << "No Analysis PID Cut defined " << endl;
	  }
	  return pidCuts;
	}
  
  
	//Make/Tighten track Cuts that are *NOT* already
	//done in the AOD production
	//**IMPORTANT**: For AODs, select FilterBit
	//the method is ignored for ESDs
	
	AliAnalysisCuts* GetTrackCutsAna(Int_t cutSet) {
    cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
	  AliDielectronCutGroup* trackCuts=0x0;
	  switch (cutSet) {
        
      
        //----------
        // these MAIN settings have to combine different track selections:
        //----------
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
        // combine typical and new trackcuts with "kCompOR" condition:
        cgTrackCutsAnaSPDorSDD = new AliDielectronCutGroup("cgTrackCutsAnaSPDorSDD","cgTrackCutsAnaSPDorSDD",AliDielectronCutGroup::kCompOR);
        cgTrackCutsAnaSPDorSDD->AddCut(GetTrackCutsAna(kPbPb2011TRK_SPDfirst));         // typical trackcuts with requirement of SPD
        cgTrackCutsAnaSPDorSDD->AddCut(GetTrackCutsAna(kPbPb2011TRK_SDDfirstSPDnone)); // new additional trackcuts with SDD instead of SPD
        trackCuts = cgTrackCutsAnaSPDorSDD;
        break;
        
        //----------
        // these MAIN settings just load the main track selection directly below:
        //----------
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
      case kPbPb2011_pidITSTPC_trkSPDfirst_3:
      case kPbPb2011_pidTPC_trkSPDfirst_3:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
        //----------
      case kPbPb2011TRK_SPDfirst: // main track selection, now closer to what Hongyan does...
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 2 with PID
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
        break;
        
      case kPbPb2011TRK_SDDfirstSPDnone: // complimentary tracks, strictly without SPD, to be combined with others!
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0); // means at least 3 with PID
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<6); // GetStandardITSTPCTrackCuts2011(kTRUE), SPD none, SDD first
        
        cgTrackCutsAnaSDDnoSPD = new AliDielectronCutGroup("cgTrackCutsAnaSDDnoSPD","cgTrackCutsAnaSDDnoSPD",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsDiel);
        cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSDDnoSPD;
        break;
        
        //----------
        // MAIN settings - combined trackset - variation 1
        //----------
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
        // combine typical and new trackcuts with "kCompOR" condition:
        cgTrackCutsAnaSPDorSDD = new AliDielectronCutGroup("cgTrackCutsAnaSPDorSDD","cgTrackCutsAnaSPDorSDD",AliDielectronCutGroup::kCompOR);
        cgTrackCutsAnaSPDorSDD->AddCut(GetTrackCutsAna(kPbPb2011TRK_SPDfirst5cls));         // typical trackcuts with requirement of SPD
        cgTrackCutsAnaSPDorSDD->AddCut(GetTrackCutsAna(kPbPb2011TRK_SDDfirstSPDnone4cls)); // new additional trackcuts with SDD instead of SPD
        trackCuts = cgTrackCutsAnaSPDorSDD;
        break;
        
        //----------
        // MAIN settings - single trackset - variation 1
        //----------
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
        //----------
      case kPbPb2011TRK_SPDfirst5cls: // main track selection, 5+ ITS clusters
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     5.0, 100.0); // means at least 3 with PID
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
        break;
        
      case kPbPb2011TRK_SDDfirstSPDnone4cls: // complimentary tracks, 4+ ITS clusters, strictly without SPD, to be combined with others!
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 4 with PID
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<6); // GetStandardITSTPCTrackCuts2011(kTRUE), SPD none, SDD first
        
        cgTrackCutsAnaSDDnoSPD = new AliDielectronCutGroup("cgTrackCutsAnaSDDnoSPD","cgTrackCutsAnaSDDnoSPD",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsDiel);
        cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSDDnoSPD;
        break;
        
        // ==========
      case kPbPb2011_TPCITS_TOFif1:
      case kPbPb2011_TPCTOF_Semi2: // no pt and eta ranges in the trackcuts anymore!
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        // trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -0.84, 0.84); // eta commented out later. (05.02.2014)
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     110.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.0);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(16); //does nothing for ESDs // 16=2^4 -> filter bit 4!
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst); //function in AliDielectronTrackCuts
        
        cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAna","cgTrackCutsAna",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAna->AddCut(trackCutsDiel);
        cgTrackCutsAna->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAna;
        break;
        
      case kPbPb2011_TPCTOF_Semi1:
        //[...]
      default: cout << "No Analysis Track Cut defined " << endl;
	  }
	  return trackCuts;
	} 
  
  
  
	//Relaxed PID cuts for additional rejectin step, do not use blindly
	AliAnalysisCuts* GetPIDCutsPre(Int_t cutSet) {
    cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
	  AliAnalysisCuts* pidCuts=0x0;
	  switch (cutSet) {
      case kPbPb2011_pidITSTPC_trkSPDfirst_3:
      case kPbPb2011_pidTPC_trkSPDfirst_3:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      case kPbPb2011_TPCITS_TOFif1:
      case kPbPb2011_TPCTOF_Semi2:
        
        // eta range:
        AliDielectronVarCuts *etaRangePre1 = new AliDielectronVarCuts("etaRangePre1","etaRangePre1");
        etaRangePre1->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
        // pt range:
        AliDielectronVarCuts *ptRangePre1 = new AliDielectronVarCuts("ptRangePre1","ptRangePre1");
        ptRangePre1->AddCut(AliDielectronVarManager::kPt, .2, 3.5); // 0.2 is realistic. turnon at ~180MeV
        //AliDielectronVarCuts *ptRangePre2 = new AliDielectronVarCuts("ptRangePre2","ptRangePre2");
        //ptRangePre2->AddCut(AliDielectronVarManager::kPt, .4, 3.5);
        //AliDielectronVarCuts *ptRangePre3 = new AliDielectronVarCuts("ptRangePre3","ptRangePre3");
        //ptRangePre3->AddCut(AliDielectronVarManager::kPt, 0.05, 1.5);
        
        AliDielectronCutGroup* cgITSTPCTOFpre = new AliDielectronCutGroup("cgITSTPCTOFpre","cgITSTPCTOFpre",AliDielectronCutGroup::kCompAND);
        AliDielectronPID *pidITSTPCTOFpre = new AliDielectronPID("pidITSTPCTOFpre","pidITSTPCTOFpre");
        pidITSTPCTOFpre->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3., 0. ,100., kFALSE);
        pidITSTPCTOFpre->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
        // ITS will be used:
        //  pidITSTPCTOFpre->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3., 0. ,1.7 , kFALSE);
        // TOF will be used if available, and with pt instead of p:
        //  pidITSTPCTOFpre->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0.4,100., kFALSE, 
        //                          AliDielectronPID::kIfAvailable, AliDielectronVarManager::kPt);
        cgITSTPCTOFpre->AddCut(pidITSTPCTOFpre);
        cgITSTPCTOFpre->AddCut(etaRangePre1);
        cgITSTPCTOFpre->AddCut(ptRangePre1);
        cgITSTPCTOFpre->AddCut(GetTrackCutsAna(cutSet));
        
        //        AliDielectronCutGroup* cgTPCpre = new AliDielectronCutGroup("cgTPCpre","cgTPCpre",AliDielectronCutGroup::kCompAND);
        //        AliDielectronPID *pidTPCpre = new AliDielectronPID("pidTPCpre","pidTPCpre");
        //        pidTPCpre->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3., 0. ,100., kFALSE);
        //        pidTPCpre->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
        //        // TOF will be used if available, and with pt instead of p:
        //        pidTPCpre->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0.4,5.  , kFALSE, 
        //                          AliDielectronPID::kIfAvailable, AliDielectronVarManager::kPt);
        //        cgTPCpre->AddCut(pidTPCpre);
        //        cgITSTPCTOFpre->AddCut(etaRangePre1);
        //        cgTPCpre->AddCut(ptRangePre2);
        //        cgTPCpre->AddCut(GetTrackCutsAna(cutSet));
        
        //        AliDielectronCutGroup* cgITSSA = new AliDielectronCutGroup("cgITSSA","cgITSSA",AliDielectronCutGroup::kCompAND);
        //        AliDielectronPID *pidITSSA = new AliDielectronPID("pidITSSA","pidITSSA");
        //        pidITSSA->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 3.);
        //        // this means that very many pions will be used for rejection!
        //        cgITSSA->AddCut(pidITSSA);
        //        cgITSTPCTOFpre->AddCut(etaRangePre1);
        //        cgITSSA->AddCut(ptRangePre3);
        //        cgITSSA->AddCut(GetTrackCutsPre(cutSet));
        
        AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
        cgInitialTrackFilter->AddCut(GetPIDCutsAna(cutSet)); // in case the prefilter cuts do not include all needed global tracks.
        cgInitialTrackFilter->AddCut(cgITSTPCTOFpre);
        //cgInitialTrackFilter->AddCut(cgTPCpre);
        //cgInitialTrackFilter->AddCut(cgITSSA);
        pidCuts = cgInitialTrackFilter;   // kCompOR works!!! <- checked with 'SetNoPairing()' and commented out 'GetPIDCutsAna(selectedPID)'
        //cout << " ========== pidCuts prefilter: ========== " << endl;
        //pidCuts->Print();
        break;
        
      case kPbPb2011_TPCTOF_Semi1:
        //[...]
      default: cout << "No Prefilter PID Cut defined " << endl;
	  }
	  return pidCuts;
	}
  
  
	//Possibly different cut sets for Prefilter step
	//Not used at the moment
	AliAnalysisCuts* GetTrackCutsPre(Int_t cutSet) {
    cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
	  AliDielectronCutGroup* trackCuts=0x0;
	  switch (cutSet) {
      case kPbPb2011_pidITSTPC_trkSPDfirst_3:
      case kPbPb2011_pidTPC_trkSPDfirst_3:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kEta,-0.84,0.84);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1); //does nothing for ESDs, ITSSA(???) // maybe use FilterBit(2) instead!
        
        cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
        cgTrackCutsPre->AddCut(trackCutsDiel);
        cgTrackCutsPre->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsPre;
        break;
        
      case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1); //does nothing for ESDs, ITSSA(???) // maybe use FilterBit(2) instead!
        trackCutsDiel->SetRequireITSRefit(kTRUE); //function in AliDielectronTrackCuts
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst); //function in AliDielectronTrackCuts
        
        cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
        cgTrackCutsPre->AddCut(trackCutsDiel);
        cgTrackCutsPre->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsPre;
        break;
        
      case kPbPb2011_TPCITS_TOFif1:
      case kPbPb2011_TPCTOF_Semi2: // no pt ranges in the trackcuts anymore!
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kEta,-0.84,0.84);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1); //does nothing for ESDs, ITSSA(???) // maybe use FilterBit(2) instead!
        trackCutsDiel->SetRequireITSRefit(kTRUE); //function in AliDielectronTrackCuts
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst); //function in AliDielectronTrackCuts
        
        cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
        cgTrackCutsPre->AddCut(trackCutsDiel);
        cgTrackCutsPre->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsPre;
        break;
        
      case kPbPb2011_TPCTOF_Semi1:
        //[...]
      default: cout << "No Prefilter Track Cut defined " << endl;
	  }
	  return trackCuts;
	}
  
  
  
	//*******************************************************************************
	//*******************************************************************************
	//** ESD TRACK CUTS TUNED FOR AGREEMENT BETWEEN AODS AND ESDS  ******************
	//** NOT NECESSARILY 100% OPTIMIZED FOR DIEL-ANALYSIS          ******************
	//*******************************************************************************
	//*******************************************************************************
  
	//WHEN RUNNING ON ESDs: LOAD Default Cuts for AODs
	AliAnalysisCuts* GetESDTrackCutsAna(Int_t cutSet) {
    //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
    cout << " >>>>>>>>>>>>>>>>>>>>>>  GetESDTrackCutsAna()  >>>>>>>>>>>>>>>>>>>>>> " << endl;
    //cout << " >>>>>>>>>>>>>>>>>>>>>> Setting ESD Track Cuts >>>>>>>>>>>>>>>>>>>>>> " << endl;
    //cout << " >>>>>>>>>>>>>>>>>>>>>> ( do we run on ESD?! ) >>>>>>>>>>>>>>>>>>>>>> " << endl;
    //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
	  AliESDtrackCuts* esdTrackCutsH = 0x0;
	  switch (cutSet) {
      case kPbPb2011_pidITSTPC_trkSPDfirst_3:
      case kPbPb2011_pidTPC_trkSPDfirst_3:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
      case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
      case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
      case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
      case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      case kPbPb2011_TPCITS_TOFif1:
      case kPbPb2011_TPCTOF_Semi2:
      case kPbPb2011_TPCTOF_Semi1:
      case kPbPb2011NoPID:
      case kPbPb2011TPCandTOF :
   		case kPbPb2011TPCandTOFHPT:
      case kPbPb2011TPC :
      case kPbPb2011TPCandTOFwide :
      case kPbPb2011TPCorTOF  :
      case kpp2010TPCandTOF :
      case kpp2010TPCorTOF  :
        // standard cuts with very loose DCA: Bit4 (Int: 16), AOD095&115
        
        esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
        esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
        esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
        esdTrackCutsH->SetDCAToVertex2D(kTRUE);
        
        //The cuts below should be the onyl ones that are missing
        //explicitely in the TrackCutsAna method
        //To be sure, StandardITSTPCTrackCuts is loaded however
        /* 
         esdTrackCutsH = new AliESDtrackCuts();
         esdTrackCutsH->SetAcceptKinkDaughters(kFALSE);
         //Not done so far via dielectron cuts:
         */
        /*
         esdTrackCuts->SetDCAToVertex2D(kFALSE);
         esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
         esdTrackCuts->SetMaxChi2PerClusterITS(36);
         */
        
        break;
      default: cout << "No ESD Track Cut defined " << endl;
	  }
	  return esdTrackCutsH;
	} 
  
  
};
