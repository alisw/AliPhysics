class LMEECutLib {

  public:
	static  enum LMMECutSet {
   		kPbPb2011TPCandTOF,
		kPbPb2011TPCorTOF,
		kpp2010TPCandTOF,
		kpp2010TPCorTOF,
		kCUTSETMAX
	};

	static  enum LMMECentSel {
	  kPbPb2011Central,
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

	AliAnalysisCuts* GetPIDCutsAna(Int_t cutSet);  
	AliAnalysisCuts* GetPIDCutsPre(Int_t cutSet);  

	AliAnalysisCuts* GetPairCuts(Int_t cutSet);  

	AliAnalysisCuts* GetTrackCutsAna(Int_t cutSet);  
	AliAnalysisCuts* GetTrackCutsPre(Int_t cutSet);  




	AliDielectronEventCuts* GetEventCuts(Int_t cutSet) {
	  AliDielectronEventCuts* eventCuts = 0x0;
	  switch (cutSet) {
		case kPbPb2011TPCandTOF :
		case kPbPb2011TPCorTOF  :
		case kpp2010TPCandTOF :
		case kpp2010TPCorTOF  :
		  eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
		  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny); // AOD
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
		case kPbPb2011SemiCentral:
		  centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb2011SemiCentral");
		  //Restrict to 50%, Trigger selction
		  centCuts->AddCut(AliDielectronVarManager::kCentrality,10.,50.);
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
		case kPbPb2011TPCandTOF :
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
		  mixingHandler = new AliDielectronMixingHandler;
		  mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
		  mixingHandler->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,40,80");
		  mixingHandler->SetDepth(25);
		  mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
		case kPbPb2011TPCandTOF :
		  mixingHandler = new AliDielectronMixingHandler;
		  mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
		  mixingHandler->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,40,80");
		   mixingHandler->AddVariable(AliDielectronVarManager::kv0ACrpH2,"-6*(TMath::Pi()/6),-5*(TMath::Pi()/6),-4*(TMath::Pi()/6),-3*(TMath::Pi()/6),-2*(TMath::Pi()/6),-1*(TMath::Pi()/6),0,1*(TMath::Pi()/6),2*(TMath::Pi()/6),3*(TMath::Pi()/6),4*(TMath::Pi()/6),5*(TMath::Pi()/6),6*(TMath::Pi()/6)");
		  //mixingHandler->SetDepth(50);
		  mixingHandler->SetDepth(20);
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
	  pidTT->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.4 , 100., kFALSE );

	  pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3,3.,0.,100.,kTRUE);
	  pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.,3.,0.,100.,kTRUE);
	  pidTT->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3.,3.,0.,100.,kTRUE);
	  //___________________________________________
	  AliDielectronVarCuts *pidTPCsignal = new AliDielectronVarCuts("pidTPCsignal","cut on the TPC signal");
	  pidTPCsignal->AddCut(AliDielectronVarManager::kTPCsignal,75.,90.); 
	  //___________________________________________

	  AliDielectronVarCuts *pTPC = new AliDielectronVarCuts("P>.4","P>.4");
	  pTPC->AddCut(AliDielectronVarManager::kP,.4,1e30);
	  
	  AliDielectronVarCuts *pMin = new AliDielectronVarCuts("P>.2","P>.2");
	  pMin->AddCut(AliDielectronVarManager::kP,.2,1e30);

	  switch (cutSet) {
		case kPbPb2011TPCandTOF :
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC1 = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);
		  cgSecondTrackFilterPIDTPC1->AddCut(pTPC);
		  //cgSecondTrackFilterPIDTPC1->AddCut(pidTPChardTOF);
		  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCTOFeOnly);
		  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCsignal);
		  cgSecondTrackFilterPIDTPC1->AddCut(GetTrackCutsAna(cutSet));
		  anaCuts = cgSecondTrackFilterPIDTPC1;
		  break;
		case kPbPb2011TPCorTOF  :
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC2 = new AliDielectronCutGroup("cgPIDTPC2","cgPIDTPC2",AliDielectronCutGroup::kCompAND);
		  cgSecondTrackFilterPIDTPC2->AddCut(pMin);
		  cgSecondTrackFilterPIDTPC2->AddCut(pidTT);
		  cgSecondTrackFilterPIDTPC2->AddCut(pidTPCsignal);
		  cgSecondTrackFilterPIDTPC2->AddCut(GetTrackCutsAna(cutSet));
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
		case kPbPb2011TPCandTOF :
		case kPbPb2011TPCorTOF  :
		case kpp2010TPCandTOF :
		case kpp2010TPCorTOF  :
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC = new AliDielectronCutGroup("cgPIDTPC","cgPIDTPC",AliDielectronCutGroup::kCompAND);
		  AliDielectronCutGroup* cgITSTPC = new AliDielectronCutGroup("cgITSTPC","cgITSTPC",AliDielectronCutGroup::kCompAND);
		  AliDielectronPID *pidITSTPC = new AliDielectronPID("TPCpre","TPCpre");
		  pidITSTPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
		  cgITSTPC->AddCut(pidITSTPC);
		  //	  cgITSTPC->AddCut(GetTrackCutsAna(cutSet));

		  AliDielectronCutGroup* cgITSSA = new AliDielectronCutGroup("cgITSSA","cgITSSA",AliDielectronCutGroup::kCompAND);
		  AliDielectronPID *pidITSSA = new  AliDielectronPID("pidITSSA","pidITSSA");
		  pidITSSA->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-3.,3.);
		  cgITSSA->AddCut(pidITSSA);
		  //	  cgITSSA->AddCut(GetTrackCutsPre(cutSet));

		  AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
		  cgInitialTrackFilter->AddCut(cgITSTPC);
		  cgInitialTrackFilter->AddCut(cgITSSA);
		  anaCuts = cgInitialTrackFilter;
		  break;
		default: cout << "No Pre-PID Cut defined " << endl;
	  }
	  return anaCuts;
	}


	AliAnalysisCuts* GetPairCuts(Int_t cutSet)  {  
	  AliDielectronVarCuts* pairCuts=0x0;
	  switch (cutSet) {
		case kPbPb2011TPCorTOF  :
		case kpp2010TPCandTOF :
		  pairCuts = new AliDielectronVarCuts("InvMass","InvMass > 150 MeV");
		  pairCuts->AddCut(AliDielectronVarManager::kM,0.15,100.,kTRUE);
		  break;
		case kPbPb2011TPCandTOF :
		case kpp2010TPCorTOF  :
		pairCuts =new AliDielectronVarCuts("OpeningAngle","Opening angle > .035rad");
		pairCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0. , 0.035);
		break;
		default: cout << "No Pair Cuts defined " << endl;
	  } 
	  return pairCuts;
	}

	AliAnalysisCuts* GetTrackCutsAna(Int_t cutSet) {
	  AliDielectronCutGroup* trackCuts=0x0;
	  switch (cutSet) {
		case kPbPb2011TPCandTOF :
		case kPbPb2011TPCorTOF  :
		case kpp2010TPCandTOF :
		case kpp2010TPCorTOF  :
			trackCuts = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);

			trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
 
		        //Not necessary for AOD?
		        //-AOD-trackCuts->SetDCAToVertex2D(kTRUE);
			
		  
			//Not possible in AOD(?) 
			//-AOD-trackCuts->SetMinNCrossedRowsTPC(110);
			//-AOD-trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
		  
			trackCutsAOD->AddCut(AliDielectronVarManager::kPt,0.05,200.);
			trackCutsAOD->AddCut(AliDielectronVarManager::kEta,-0.84,0.84);

			trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
			trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
			//Not working for AODs
		//	trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.5);
			
			//trackCuts->SetMinNClustersTPC(60);
			//Legacy cut: Use Crossed Rows in ESD, in AOD ASAP
			trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,     60.0, 160.0);
			//-AOD-trackCuts->SetMinNClustersITS(3);
			trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
 		        
			//trackCuts->SetAcceptKinkDaughters(kFALSE);  
			trackCutsAOD->AddCut(AliDielectronVarManager::kKinkIndex0,   0.0);
			



			AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
			trackCutsDiel->SetRequireITSRefit(kTRUE);
			trackCutsDiel->SetRequireTPCRefit(kTRUE);
			trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst)
			// trackCutsAOD->AddCut(AliDielectronVarManager::kITSLayerFirstCls,-0.01,0.5); //ITS(0) = SPDfirst
			
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
		case kPbPb2011TPCandTOF :
		case kPbPb2011TPCorTOF  :
		case kpp2010TPCandTOF :
		case kpp2010TPCorTOF  :
/*
	//	  trackCuts = new AliESDtrackCuts();
	//	  trackCuts->SetDCAToVertex2D(kTRUE);
	//	  trackCuts->SetMaxDCAToVertexZ(3.0);
	//	  trackCuts->SetMaxDCAToVertexXY(1.0);
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
