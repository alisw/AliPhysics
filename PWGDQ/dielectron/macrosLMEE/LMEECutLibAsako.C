class LMEECutLibAsako {

public:
  static  enum LMMECutSet {
	kPbPb2011NoPID,
	kPbPb2011TPCandTOF,
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
	kPbPb2011SemiCentral1,
    kPbPb2011SemiCentral2,
	kPbPb2011Peripheral,
	kCENTSELMAX
  };

  static  enum LMEEPairCutSet{
	kPbPb2011RP,
	kPbPb2011Mag,
	//	kPbPb2011RPxy,
	//kPbPb2011RPxz,

	kPbPb2011MassLow,
	kPbPb2011MassMiddle,
	kPbPb2011MassHigh,
	kPbPb2011MassAll,



	kPAIRCUTSETMAX
  };

  //char* LMEECutNames[kCUTSETMAX] = { "PbPb2011TPCandTOF","PbPb2011TPCorTOF"};


  LMEECutLib() {}

  AliDielectronEventCuts*     GetEventCuts(Int_t cutSet);
  AliAnalysisCuts*            GetCentralityCuts(Int_t centSel);
  AliDielectronTrackRotator*  GetTrackRotator(Int_t cutSet);
  AliDielectronMixingHandler* GetMixingHandler(Int_t cutSet);

  AliAnalysisCuts* GetPIDCutsAna(Int_t cutSet);  
  AliAnalysisCuts* GetPIDCutsPre(Int_t cutSet);  

  AliAnalysisCuts* GetPairCuts2(Int_t cutSet,Bool_t tooglePC=kFALSE);
  AliAnalysisCuts* GetPairCuts(Int_t cutSet);  

  AliAnalysisCuts* GetPairCutsInOut(Int_t cutSet);
  AliAnalysisCuts* GetPairCuts4(Int_t cutSet);

  AliAnalysisCuts* GetTrackCutsAna(Int_t cutSet);  
  AliAnalysisCuts* GetTrackCutsPre(Int_t cutSet);  




  AliDielectronEventCuts* GetEventCuts(Int_t cutSet) {
	AliDielectronEventCuts* eventCuts = 0x0;
	switch (cutSet) {
	case kPbPb2011NoPID:
	case kPbPb2011TPCandTOF :
	case kPbPb2011TPCandTOFHPT:
	case kPbPb2011TPC :
	case kPbPb2011TPCandTOFwide :
	case kPbPb2011TPCorTOF  :
	case kpp2010TPCandTOF :
	case kpp2010TPCorTOF  :
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

  AliAnalysisCuts* GetCentralityCuts(Int_t centSel) {
	AliDielectronVarCuts* centCuts = 0x0;
	switch (centSel) {
	case kPbPb2011Central:
	  centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb2011Central");
	  centCuts->AddCut(AliDielectronVarManager::kCentrality,0.,10.);
	  break;
	case kPbPb2011SemiCentral1:
	  centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb2011SemiCentral1");
	  //Restrict to 50%, Trigger selction
	  centCuts->AddCut(AliDielectronVarManager::kCentrality,10.,30.);
	  break;
	case kPbPb2011SemiCentral2:
	  centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb2011SemiCentral2");
	  //Restrict to 50%, Trigger selction
	  centCuts->AddCut(AliDielectronVarManager::kCentrality,30.,50.);//
	  break;
	case kPbPb2011Peripheral:
	  centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb2011Peripheral");
	  centCuts->AddCut(AliDielectronVarManager::kCentrality,50.,90.);
	  break;
	default: cout << "No Centrality selected" << endl;
	}
	return centCuts;
  }


  AliDielectronTrackRotator* GetTrackRotator(Int_t cutSet) {
	AliDielectronTrackRotator* trackRotator = 0x0;
	switch (cutSet) {
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
	case kPbPb2011TPCorTOF  :
	  /*
		  mixingHandler = new AliDielectronMixingHandler;
		    mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
			  mixingHandler->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,40,80");
			    mixingHandler->SetDepth(25);
				  mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
				    break;
	  */
	case kPbPb2011NoPID:
	case kPbPb2011TPCandTOF :
	case kPbPb2011TPCandTOFHPT:
	case kPbPb2011TPC :
	case kPbPb2011TPCandTOFwide :
	  mixingHandler = new AliDielectronMixingHandler;
	  mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
	  mixingHandler->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,40,80");
	  mixingHandler->AddVariable(AliDielectronVarManager::kv0ACrpH2,"-6*(TMath::Pi()/6),-5*(TMath::Pi()/6),-4*(TMath::Pi()/6),-3*(TMath::Pi()/6),-2*(TMath::Pi()/6),-1*(TMath::Pi()/6),0,1*(TMath::Pi()/6),2*(TMath::Pi()/6),3*(TMath::Pi()/6),4*(TMath::Pi()/6),5*(TMath::Pi()/6),6*(TMath::Pi()/6)");
	  //mixingHandler->SetDepth(50);
	  mixingHandler->SetDepth(15);
	  mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
	  break;
	case kpp2010TPCandTOF :
	case kpp2010TPCorTOF  :
	  mixingHandler = new AliDielectronMixingHandler;
	  mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
	  mixingHandler->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
	  //might want to add multiplicity?
	  mixingHandler->SetDepth(50);
	  mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
	  break;
	default: cout << "No Rotator defined" << endl;
	}
	return mixingHandler;
  }


  AliAnalysisCuts* GetPIDCutsAna(Int_t cutSet) {
	AliAnalysisCuts* anaCuts=0x0;

	// + [2] added for Dec2010 Cut! 
	TF1 *lowerCut = new TF1("lowerCut", "[0] * TMath::Exp([1]*x) + [2]", 0, 20);
	/* until Nov2010
	    lowerCut->SetParameter(0, -2.7);
		 lowerCut->SetParameter(1, -0.4357);
	*/
	/* 18.01.2011 ALiHFEpid.cxx */
	lowerCut->SetParameter(0,-3.7);
	lowerCut->SetParameter(1,-0.8);
	lowerCut->SetParameter(2,-0.35);

	if (MCenabled) { //overwrite parameters
	  lowerCut->SetParameter(0,-2.5);
	  lowerCut->SetParameter(2,-2.2);
	}

	//---------------------------------------------
	AliDielectronPID *pidTPCTOFeOnly = new AliDielectronPID("TPC-TOF","TPC-TOF");
	pidTPCTOFeOnly->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3,3.,0.0,100.,kFALSE);
	pidTPCTOFeOnly->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.0 , 100., kFALSE );


	AliDielectronPID *pidTPCandTOF = new AliDielectronPID("TPC-TOF-HFE","TPC-TOF-HFE");
	pidTPCandTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.,0.0,100.,kFALSE);
	pidTPCandTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.,3.,0.,100.,kTRUE);
	pidTPCandTOF->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.0 , 100., kFALSE );

	AliDielectronPID *pidTPChardTOF = new AliDielectronPID("TPC-TOF-HFE","TPC-TOF-HFE");
	pidTPChardTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,lowerCut,3.,0.0,100.,kFALSE);
	pidTPChardTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3,3.,0.,100.,kTRUE);
	pidTPChardTOF->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.,3.,0.,100.,kTRUE);
	pidTPChardTOF->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3.,3.,0.,100.,kTRUE);
	pidTPChardTOF->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.0 , 100., kFALSE );
	//___________________________________________
	AliDielectronPID *pidTT = new AliDielectronPID("TPC-TOF","TPC-TOF");
	pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.5,3.,0.2,0.4,kFALSE);
	pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3,3.,0.4,100.,kFALSE);
	pidTT->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.2 , 100., kFALSE );

	pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3,3.,0.,100.,kTRUE);
	pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.,3.,0.,100.,kTRUE);
	pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3.,3.,0.,100.,kTRUE);
	//___________________________________________
	AliDielectronVarCuts *pidTPCsignal = new AliDielectronVarCuts("pidTPCsignal","cut on the TPC signal");
	if (MCenabled) {
	  pidTPCsignal->AddCut(AliDielectronVarManager::kTPCsignal,65.,85.); 
	}
	else {
	  pidTPCsignal->AddCut(AliDielectronVarManager::kTPCsignal,75.,90.); 
	}
	//___________________________________________

	AliDielectronVarCuts *pidTPCsignalWide = new AliDielectronVarCuts("pidTPCsignalWide","cut on the TPC signal");
	pidTPCsignalWide->AddCut(AliDielectronVarManager::kTPCsignal,75.,90.);


	AliDielectronVarCuts *pTPCHPT = new AliDielectronVarCuts("P>.4hpt","P>.4hpt");
	pTPCHPT->AddCut(AliDielectronVarManager::kPt,.4,3.0);

	AliDielectronVarCuts *pTPC = new AliDielectronVarCuts("P>.4","P>.4");
	pTPC->AddCut(AliDielectronVarManager::kPt,.4,2.0);
	  
	AliDielectronVarCuts *pMin = new AliDielectronVarCuts("P>.2","P>.2");
	pMin->AddCut(AliDielectronVarManager::kPt,.2,2.5);

	switch (cutSet) {
	case kPbPb2011NoPID:
	  AliDielectronCutGroup* cgSecondTrackFilterNoPID = new AliDielectronCutGroup("cgNoPID","cgNoPID",AliDielectronCutGroup::kCompAND);
	  cgSecondTrackFilterNoPID->AddCut(pTPC);
	  anaCuts= cgSecondTrackFilterNoPID;
	  break;
	case kPbPb2011TPCandTOFHPT:
	  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC1 = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);
	  cgSecondTrackFilterPIDTPC1->AddCut(pTPCHPT);
	  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCTOFeOnly);
	  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCsignal);
	  anaCuts = cgSecondTrackFilterPIDTPC1;

	case kPbPb2011TPCandTOF :
	  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC1 = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);
	  cgSecondTrackFilterPIDTPC1->AddCut(pTPC);
	  //cgSecondTrackFilterPIDTPC1->AddCut(pidTPChardTOF);
	  // cgSecondTrackFilterPIDTPC1->AddCut(pidTPCTOFeOnly);
	  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCandTOF);
	  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCsignal);
	  anaCuts = cgSecondTrackFilterPIDTPC1;
	  break;

	case kPbPb2011TPC :
	  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC1 = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);
	  cgSecondTrackFilterPIDTPC1->AddCut(pTPC);
	  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCsignalWide);
	  anaCuts = cgSecondTrackFilterPIDTPC1;
	  break;

	case kPbPb2011TPCandTOFwide :
	  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC1 = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);
	  cgSecondTrackFilterPIDTPC1->AddCut(pTPC);
	  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCandTOF);
	  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCsignalWide);
	  anaCuts = cgSecondTrackFilterPIDTPC1;
	  break;

	case kPbPb2011TPCorTOF  :
	  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC2 = new AliDielectronCutGroup("cgPIDTPC2","cgPIDTPC2",AliDielectronCutGroup::kCompAND);
	  cgSecondTrackFilterPIDTPC2->AddCut(pMin);
	  cgSecondTrackFilterPIDTPC2->AddCut(pidTT);
	  cgSecondTrackFilterPIDTPC2->AddCut(pidTPCsignal);
	  anaCuts = cgSecondTrackFilterPIDTPC2;
	  break;

	case kpp2010TPCandTOF :
	  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC = new AliDielectronCutGroup("cgPIDTPC","cgPIDTPC",AliDielectronCutGroup::kCompAND);
	  cgSecondTrackFilterPIDTPC->AddCut(pTPC);
	  //cgSecondTrackFilterPIDTPC->AddCut(pidTPChardTOF);
	  cgSecondTrackFilterPIDTPC->AddCut(pidTPCTOFeOnly);
	  anaCuts = cgSecondTrackFilterPIDTPC;
	  break;
	case kpp2010TPCorTOF  :
	  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC = new AliDielectronCutGroup("cgPIDTPC","cgPIDTPC",AliDielectronCutGroup::kCompAND);
	  cgSecondTrackFilterPIDTPC->AddCut(pMin);
	  cgSecondTrackFilterPIDTPC->AddCut(pidTT);
	  anaCuts = cgSecondTrackFilterPIDTPC;
	  break;
	default: cout << "No Analysis PID Cut defined " << endl;
	}
	return anaCuts;
  }

  AliAnalysisCuts* GetPIDCutsPre(Int_t cutSet) {
	AliAnalysisCuts* anaCuts=0x0;
	switch (cutSet) {
	case kPbPb2011NoPID:
	case kPbPb2011TPCandTOF :
	case kPbPb2011TPCandTOFHPT:
	case kPbPb2011TPC :
	case kPbPb2011TPCandTOFwide :
	case kPbPb2011TPCorTOF  :
	case kpp2010TPCandTOF :
	case kpp2010TPCorTOF  :
	  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC = new AliDielectronCutGroup("cgPIDTPC","cgPIDTPC",AliDielectronCutGroup::kCompAND);
	  AliDielectronCutGroup* cgITSTPC = new AliDielectronCutGroup("cgITSTPC","cgITSTPC",AliDielectronCutGroup::kCompAND);
	  AliDielectronPID *pidITSTPC = new AliDielectronPID("TPCpre","TPCpre");
	  pidITSTPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
	  cgITSTPC->AddCut(pidITSTPC);


	  AliDielectronCutGroup* cgITSSA = new AliDielectronCutGroup("cgITSSA","cgITSSA",AliDielectronCutGroup::kCompAND);
	  AliDielectronPID *pidITSSA = new  AliDielectronPID("pidITSSA","pidITSSA");
	  pidITSSA->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-3.,3.);
	  cgITSSA->AddCut(pidITSSA);
	  //  cgITSSA->AddCut(GetTrackCutsPre(cutSet));

	  AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
	  cgInitialTrackFilter->AddCut(cgITSTPC);
	  cgInitialTrackFilter->AddCut(cgITSSA);
	  anaCuts = cgInitialTrackFilter;
	  break;
	default: cout << "No Pre-PID Cut defined " << endl;
	}
	return anaCuts;
  }



  AliAnalysisCuts* GetPairCuts2(Int_t cutSet, Bool_t togglePC /*=kFALSE*/) {
	AliAnalysisCuts* pairCuts=0x0;
	switch (cutSet) {
	case kPbPb2011TPCorTOF  :
	case kpp2010TPCandTOF :
	  AliDielectronVarCuts* pairCutsV= new AliDielectronVarCuts("InvMass","InvMass > 150 MeV");
	  pairCutsV->AddCut(AliDielectronVarManager::kM,0.15,100.,kTRUE);
	  pairCuts = pairCutsV;
	  break;
	case kPbPb2011NoPID:
	case kPbPb2011TPCandTOF :
	case kPbPb2011TPC :
	case kPbPb2011TPCandTOFHPT:
	case kPbPb2011TPCandTOFwide :
	case kpp2010TPCorTOF  :
	  if (!togglePC) {

		AliDielectronCutGroup* pairCutsCG2 =new AliDielectronCutGroup("pairCutsCG2","pairCutsCG2",AliDielectronCutGroup::kCompOR);
		AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
		pairCutsPhiv =new AliDielectronVarCuts("Phiv Cuts","Phiv<2.0rad");
		pairCutsPhiv->AddCut(AliDielectronVarManager::kPhivPair, 0.0, 2.0); 
		pairCutsInvM =new AliDielectronVarCuts("InvM Cuts","InvM<0.3");
		pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.3); 
		pairCutsInvMgood =new AliDielectronVarCuts("InvM Cuts good","InvM>0.3");
		pairCutsInvMgood->AddCut(AliDielectronVarManager::kM, 0.3, 99999.); 
		pairCutsCG->AddCut(pairCutsPhiv);
		pairCutsCG->AddCut(pairCutsInvM);
		pairCutsCG2->AddCut(pairCutsInvMgood);
		pairCutsCG2->AddCut(pairCutsCG);
		pairCuts = pairCutsCG2;
	  }
	  else {
		AliDielectronVarCuts* pairCutsV =new AliDielectronVarCuts("OpeningAngle","Opening angle > .035rad");
		pairCutsV->AddCut(AliDielectronVarManager::kOpeningAngle, 0. , 0.035,kTRUE);
		pairCuts = pairCutsV;
	  }
	  break;
	default: cout << "No Pair Cuts defined " << endl;
	}
	return pairCuts;
  }


  AliAnalysisCuts* GetPairCuts(Int_t cutSet) {  
	AliAnalysisCuts* pairCuts=0x0;
	switch (cutSet) {
	case kPbPb2011TPCorTOF  :
	case kpp2010TPCandTOF :
	  AliDielectronVarCuts* pairCutsM=0x0;
	  pairCutsM = new AliDielectronVarCuts("InvMass","InvMass > 150 MeV");
	  pairCutsM->AddCut(AliDielectronVarManager::kM,0.15,100.,kTRUE);
	  pairCuts = pairCutsM;
	  break;
	case kPbPb2011NoPID:
	case kPbPb2011TPCandTOF :
	case kPbPb2011TPC :
	case kPbPb2011TPCandTOFwide :
	case kPbPb2011TPCandTOFHPT:
	case kpp2010TPCorTOF  :

	  AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
	  AliDielectronVarCuts* pairCutsPhiv =new AliDielectronVarCuts("Phiv Cuts","Phiv<2.0rad");
	  pairCutsPhiv->AddCut(AliDielectronVarManager::kPhivPair, 2.0, 3.2); 
	  AliDielectronVarCuts* pairCutsInvM =new AliDielectronVarCuts("InvM Cuts","InvM<0.3");
	  pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.3); 
	  pairCutsCG->AddCut(pairCutsPhiv);
	  pairCutsCG->AddCut(pairCutsInvM);

	  pairCuts = pairCutsCG;


	  //pairCuts =new AliDielectronVarCuts("OpeningAngle","Opening angle > .035rad");
	  //pairCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0. , 0.035);
	  break;
	default: cout << "No Pair Cuts defined " << endl;
	} 
	return pairCuts;
  }

  AliAnalysisCuts* GetPairCutsInOut(Int_t cutSet){
	AliAnalysisCuts* pairCut=0x0;
	switch (cutSet) {
	case kPbPb2011RP:
	  AliDielectronCutGroup* pairCutsPhiRP =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompOR);
      AliDielectronVarCuts* pairCutsPhi1 =new AliDielectronVarCuts("Phi Cuts","-pi/4<Phi<pi/4");
      pairCutsPhi1->AddCut(AliDielectronVarManager::kDeltaPhiv0CrpH2, (-1.0)*TMath::Pi()/4., TMath::Pi()/4.);
      AliDielectronVarCuts* pairCutsPhi2 =new AliDielectronVarCuts("Phi Cuts","3*pi/4<Phi");
      pairCutsPhi2->AddCut(AliDielectronVarManager::kDeltaPhiv0CrpH2, 3.0*TMath::Pi()/4., TMath::Pi());
      AliDielectronVarCuts* pairCutsPhi3 =new AliDielectronVarCuts("Phi Cuts","-3*pi/4<Phi");
      pairCutsPhi3->AddCut(AliDielectronVarManager::kDeltaPhiv0CrpH2,  (-1.0)*TMath::Pi(),(-3.0)*TMath::Pi()/4.);

      pairCutsPhiRP->AddCut(pairCutsPhi1);
	  pairCutsPhiRP->AddCut(pairCutsPhi2);
	  pairCutsPhiRP->AddCut(pairCutsPhi3);

      pairCuts = pairCutsPhiRP;

	  break;
	case kPbPb2011Mag:
	  AliDielectronCutGroup* pairCutsPhiMag=new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompOR);
      AliDielectronVarCuts* pairCutsPhi4 =new AliDielectronVarCuts("Phi Cuts","-3*pi/4<Phi");
      pairCutsPhi4->AddCut(AliDielectronVarManager::kDeltaPhiv0CrpH2, TMath::Pi()/4.,3.0*TMath::Pi()/4.);
      AliDielectronVarCuts* pairCutsPhi5 =new AliDielectronVarCuts("Phi Cuts","-3*pi/4<Phi");
      pairCutsPhi5->AddCut(AliDielectronVarManager::kDeltaPhiv0CrpH2,(-3.0)*TMath::Pi()/4.,(-1.0)*TMath::Pi()/4.);

	pairCutsPhiMag->AddCut(pairCutsPhi4);
	pairCutsPhiMag->AddCut(pairCutsPhi5);
	
	pairCuts = pairCutsPhiMag;
      //pairCuts =new AliDielectronVarCuts("OpeningAngle","Opening angle > .035rad");
      //pairCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0. , 0.035);
      break;
    default: cout << "No Pair Cuts defined " << endl;
    }
    return pairCuts;
  }


  AliAnalysisCuts* GetPairCuts4(Int_t cutSet) {
    AliDielectronVarCuts* pairCuts=0x0;
    switch (cutSet) {
    case kPbPb2011MassLow :
      pairCuts = new AliDielectronVarCuts("InvMass","0 MeV  < InvMass < 30 MeV");
      pairCuts->AddCut(AliDielectronVarManager::kM, 0.,0.03);
      break;
    case kPbPb2011MassMiddle :
      pairCuts = new AliDielectronVarCuts("InvMass","120 MeV  < InvMass < 300 MeV");
      pairCuts->AddCut(AliDielectronVarManager::kM, 0.12,0.30);
      break;
    case kPbPb2011MassHigh :
      pairCuts = new AliDielectronVarCuts("InvMass","300 MeV  < InvMass < 500 MeV");
      pairCuts->AddCut(AliDielectronVarManager::kM, 0.30,0.50);
      break;
    case kPbPb2011MassAll :
      pairCuts = new AliDielectronVarCuts("InvMass","0 GeV  < InvMass < 10 GeV");
      pairCuts->AddCut(AliDielectronVarManager::kM, 0.0,10.0);
      break;

    default: cout << "No Pair Cuts defined " << endl;
    }
    return pairCuts;
  }






  AliAnalysisCuts* GetESDTrackCutsAna(Int_t cutSet) {
	AliESDtrackCuts* esdTrackCutsH = 0x0;
	switch (cutSet) {
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
	default: cout << "No Analysis Track Cut defined " << endl;
	}
	return esdTrackCutsH;
  } 

  AliAnalysisCuts* GetTrackCutsAna(Int_t cutSet) {
	AliDielectronCutGroup* trackCuts=0x0;
	switch (cutSet) {
	case kPbPb2011NoPID:
	case kPbPb2011TPCandTOF :
	case kPbPb2011TPCandTOFHPT:
	case kPbPb2011TPC :
	case kPbPb2011TPCandTOFwide :
	case kPbPb2011TPCorTOF  :
	case kpp2010TPCandTOF :
	case kpp2010TPCorTOF  :
	  trackCuts = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);

	  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
	  trackCutsAOD->AddCut(AliDielectronVarManager::kPt,0.05,6.);
	  trackCutsAOD->AddCut(AliDielectronVarManager::kEta,-0.84,0.84);
	  //DCA Cut
	  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
	  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
	  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
	  trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.5);
	  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
	  trackCutsDiel->SetAODFilterBit(16); //does nothing for ESDs
	  trackCutsDiel->SetRequireITSRefit(kTRUE);
	  trackCutsDiel->SetRequireTPCRefit(kTRUE);

	  trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
	  //trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,     80., 140.0);
	  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     110.0, 160.0);
	  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.0);//tighter than before,
	  //due to AOD production
	  trackCuts->AddCut(trackCutsDiel);
	  trackCuts->AddCut(trackCutsAOD);
	  break;
	default: cout << "No Analysis Track Cut defined " << endl;
	}
	return trackCuts;
  } 

  AliAnalysisCuts* GetTrackCutsPre(Int_t cutSet) {
	AliESDtrackCuts* trackCuts=0x0;
	switch (cutSet) {
	case kPbPb2011NoPID:
	case kPbPb2011TPCandTOF :
	case kPbPb2011TPCandTOFHPT:
	case kPbPb2011TPCorTOF  :
	case kpp2010TPCandTOF :
	case kpp2010TPCorTOF  :
	  /*
		trackCuts = new AliESDtrackCuts();
		trackCuts->SetDCAToVertex2D(kTRUE);
		trackCuts->SetMaxDCAToVertexZ(3.0);
		trackCuts->SetMaxDCAToVertexXY(1.0);
		  trackCuts->SetEtaRange( -0.84 , 0.84 );
		    trackCuts->SetPtRange(  0.05 , 0.5);
			  trackCuts->SetAcceptKinkDaughters(kFALSE);
			    trackCuts->SetRequireITSRefit(kTRUE);
				  trackCuts->SetRequireITSStandAlone(kTRUE);
				    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
					  trackCuts->SetMinNClustersITS(3); //PhotonGroup-Dalitz: 2?!
	  */
	  cout << "No Pre-Track Cut defined for AODs at the moment " << endl;
	  break;
	default: cout << "No Pre-Track Cut defined " << endl;
	}
	return trackCuts;
  }

};