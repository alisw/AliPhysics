class LMEECutLib {

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

	AliAnalysisCuts* GetPairCuts2(Int_t cutSet,Bool_t tooglePC=kFALSE);
	AliAnalysisCuts* GetPairCuts(Int_t cutSet);  

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
		case kPbPb2011SemiCentral:
		  centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb2011SemiCentral");
		  //Restrict to 50%, Trigger selction
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
		  //ATTENTION: Trivial 1 Bin Variable on Nacc needed: Not understood bug, mixing breaks
		  //when just adding one variable *****************!!! 
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

	  //-----------------------------------------------
	  //Define different PID Cuts, that are used later
	  //-----------------------------------------------
	  
      //TPC: UPPER HALF inclusion of electron
	  //     3sigma exclusion of Pions
	  //TOF: 3sigma inclusion of electrons
	  AliDielectronPID *pidTPCTOFeOnly = new AliDielectronPID("TPC-TOF","TPC-TOF");
	  pidTPCTOFeOnly->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-0.,3.,0.0,100.,kFALSE);
	  pidTPCTOFeOnly->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.,3.,0.,100.,kTRUE);
	  pidTPCTOFeOnly->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.0 , 1.5, kFALSE );

      //TPC: 3sigma inclusion of electron
	  //     3sigma exclusion of Pions
	  //TOF: 3sigma inclusion of electrons in region where p,K cross electrons in TPC
	  AliDielectronPID *pidTPCandTOF = new AliDielectronPID("TPC-TOF-HFE","TPC-TOF-HFE");
	  pidTPCandTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.,0.0,100.,kFALSE);
	  pidTPCandTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.,3.,0.,100.,kTRUE);
	  pidTPCandTOF->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.0 , 1.5, kFALSE );

	  //Apply ITS cuts (see Hongyan's talks):
	  //3 sigma inclusion of electrons in TPC
	  //3 sigma exclusion of pions in TPC
	  //3 sigma inclusion of electrons in ITS,TOF for p<1.5GeV, where p,K contamination
	  AliDielectronPID *pidTPCandITSandTOF = new AliDielectronPID("TPC-TOFANDITS","TPC-TOFANDITS");
	  pidTPCandITSandTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.,0.0,100.,kFALSE);
	  pidTPCandITSandTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.,3.,0.,100.,kTRUE);
	  pidTPCandITSandTOF->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.0 , 1.5, kFALSE);
	  pidTPCandITSandTOF->AddCut(AliDielectronPID::kITS ,AliPID::kElectron , -6. , 3. , 0.0 , 1.5, kFALSE );

	  //Apply ITS cuts (see Hongyan's talks):
	  //3 sigma inclusion of electrons in TPC
	  //3 sigma exclusion of pions in TPC
	  //3 sigma inclusion of electrons in ITS,TOF for p<1.5GeV, where p,K contamination
	  //TOF only IF available!
	  AliDielectronPID *pidTPCandITSTOF = new AliDielectronPID("TPC-TOF-ITS","TPC-TOF-ITS");
	  pidTPCandITSTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.,0.0,100.,kFALSE);
	  pidTPCandITSTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3.,3.,0.,100.,kTRUE);
	  pidTPCandITSTOF->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.0 , 1.5, kFALSE,AliDielectronPID::kIfAvailable );
	  pidTPCandITSTOF->AddCut(AliDielectronPID::kITS ,AliPID::kElectron , -6. , 3. , 0.0 , 1.5, kFALSE );

	  //___________________________________________
	  //Direct cuts on TPC signal used for QM12
	  AliDielectronVarCuts *pidTPCsignal = new AliDielectronVarCuts("pidTPCsignal","cut on the TPC signal");
	  if (MCenabled) {
	  pidTPCsignal->AddCut(AliDielectronVarManager::kTPCsignal,65.,85.); 
	  }	
	  else {
	  pidTPCsignal->AddCut(AliDielectronVarManager::kTPCsignal,75.,90.); 
	  }
	  //___________________________________________

	  //low pT cut-off 0.4 - Pb-Pb
	  AliDielectronVarCuts *pTPC = new AliDielectronVarCuts("P>.4","P>.4");
	  pTPC->AddCut(AliDielectronVarManager::kPt,.4,3.5);
	  
	  //low pT cut-off 0.4 - pp
	  AliDielectronVarCuts *pMin = new AliDielectronVarCuts("P>.2","P>.2");
	  pMin->AddCut(AliDielectronVarManager::kPt,.2,2.5);


	  //Now see what Config actually loads and assemble final cuts
	  switch (cutSet) {
		case kPbPb2011NoPID:
		  AliDielectronCutGroup* cgSecondTrackFilterNoPID = new AliDielectronCutGroup("cgNoPID","cgNoPID",AliDielectronCutGroup::kCompAND);
		  cgSecondTrackFilterNoPID->AddCut(pTPC);
		  anaCuts= cgSecondTrackFilterNoPID;
		  break;
   		case kPbPb2011TPCandTOFHPT:
		  //test Hongyan's cut
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC1 = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);
		  cgSecondTrackFilterPIDTPC1->AddCut(pTPC);
		  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCandITSTOF);
		  anaCuts = cgSecondTrackFilterPIDTPC1;
		  break;
		case kPbPb2011TPCandTOF :
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC1 = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);
		  cgSecondTrackFilterPIDTPC1->AddCut(pTPC);
		  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCandTOF);
		  anaCuts = cgSecondTrackFilterPIDTPC1;
		  break;
		case kPbPb2011TPC :
		  //Old, QM12
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC1 = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);
		  cgSecondTrackFilterPIDTPC1->AddCut(pTPC);
		  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCsignal);
		  anaCuts = cgSecondTrackFilterPIDTPC1;
		  break;

		case kPbPb2011TPCandTOFwide :
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC1 = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);
		  cgSecondTrackFilterPIDTPC1->AddCut(pTPC);
		  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCandTOF);
		  cgSecondTrackFilterPIDTPC1->AddCut(pidTPCsignal);
		  anaCuts = cgSecondTrackFilterPIDTPC1;
		  break;

		case kPbPb2011TPCorTOF  :
		  //unused
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC2 = new AliDielectronCutGroup("cgPIDTPC2","cgPIDTPC2",AliDielectronCutGroup::kCompAND);
		  cgSecondTrackFilterPIDTPC2->AddCut(pTPC);
//		  cgSecondTrackFilterPIDTPC2->AddCut(pidTT);
		  anaCuts = cgSecondTrackFilterPIDTPC2;
		  break;
		case kpp2010TPCandTOF :
		  //unused
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC = new AliDielectronCutGroup("cgPIDTPC","cgPIDTPC",AliDielectronCutGroup::kCompAND);
		  cgSecondTrackFilterPIDTPC->AddCut(pTPC);
		  cgSecondTrackFilterPIDTPC->AddCut(pidTPCandTOF);
		  anaCuts = cgSecondTrackFilterPIDTPC;
		  break;
		case kpp2010TPCorTOF  :
		  //unused
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC = new AliDielectronCutGroup("cgPIDTPC","cgPIDTPC",AliDielectronCutGroup::kCompAND);
		  cgSecondTrackFilterPIDTPC->AddCut(pTPC);
		  ///		  cgSecondTrackFilterPIDTPC->AddCut(pidTT);
		  anaCuts = cgSecondTrackFilterPIDTPC;
		  break;
		default: cout << "No Analysis PID Cut defined " << endl;
	  }
	  return anaCuts;
	}


	//Relaxed PID cuts for additional rejectin step, do not use blindly
	AliAnalysisCuts* GetPIDCutsPre(Int_t cutSet) {
	  AliAnalysisCuts* anaCuts=0x0;
	  switch (cutSet) {
		case kPbPb2011NoPID:
		case kPbPb2011TPCandTOF :
   		case kPbPb2011TPCandTOFHPT:
		case kPbPb2011TPC :
		case kPbPb2011TPCandTOFwide :
		case kpp2010TPCandTOF :
		case kpp2010TPCorTOF  :
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC = new AliDielectronCutGroup("cgPIDTPC","cgPIDTPC",AliDielectronCutGroup::kCompAND);
		  AliDielectronCutGroup* cgITSTPC = new AliDielectronCutGroup("cgITSTPC","cgITSTPC",AliDielectronCutGroup::kCompAND);
		  AliDielectronPID *pidITSTPC = new AliDielectronPID("TPCpre","TPCpre");

		  pidITSTPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
		  pidITSTPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3,3.,0.,100.,kTRUE);
		  pidITSTPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.,3.,0.,0.4,kTRUE);
		  pidITSTPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3.,3.,0.,0.4,kTRUE);
		  pidITSTPC->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.4 , 100., kFALSE );
		  cgITSTPC->AddCut(pidITSTPC);

		  AliDielectronVarCuts *pTPCHPT = new AliDielectronVarCuts("P>.4hpt","P>.4hpt");
		  pTPCHPT->AddCut(AliDielectronVarManager::kPt,.2,3.5);
		  cgITSTPC->AddCut(pTPCHPT);

		  cgITSTPC->AddCut(GetTrackCutsAna(cutSet));


		  AliDielectronCutGroup* cgITSSA = new AliDielectronCutGroup("cgITSSA","cgITSSA",AliDielectronCutGroup::kCompAND);
		  AliDielectronPID *pidITSSA = new  AliDielectronPID("pidITSSA","pidITSSA");
		  pidITSSA->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-3.,3.);
		  cgITSSA->AddCut(pidITSSA);
		  AliDielectronVarCuts *pITSPT = new AliDielectronVarCuts("P>.4hpt","P>.4hpt");
		  pITSPT->AddCut(AliDielectronVarManager::kPt,0.0,0.8);
		  cgITSSA->AddCut(pITSPT);
		  cgITSSA->AddCut(GetTrackCutsPre(cutSet));

		  AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
		  cgInitialTrackFilter->AddCut(cgITSTPC);
		  cgInitialTrackFilter->AddCut(cgITSSA);
		  anaCuts = cgInitialTrackFilter;
		  break;


		case kPbPb2011TPCorTOF  :
		  AliDielectronCutGroup* cgSecondTrackFilterPIDTPC = new AliDielectronCutGroup("cgPIDTPC","cgPIDTPC",AliDielectronCutGroup::kCompAND);
		  AliDielectronCutGroup* cgITSTPCalone = new AliDielectronCutGroup("cgITSTPCalone","cgITSTPCalone",AliDielectronCutGroup::kCompAND);
		  AliDielectronPID *pidITSTPCalone = new AliDielectronPID("TPCpre","TPCpre");

		  pidITSTPCalone->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
		  pidITSTPCalone->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-3,3.,0.,100.,kTRUE);
		  pidITSTPCalone->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-3.,3.,0.,0.4,kTRUE);
		  pidITSTPCalone->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3.,3.,0.,0.4,kTRUE);
		  pidITSTPCalone->AddCut(AliDielectronPID::kTOF ,AliPID::kElectron , -3. , 3. , 0.4 , 100., kFALSE );
		  cgITSTPCalone->AddCut(pidITSTPCalone);

		  AliDielectronVarCuts *pTPCHPT = new AliDielectronVarCuts("P>.4hpt","P>.4hpt");
		  pTPCHPT->AddCut(AliDielectronVarManager::kPt,.2,3.5);
		  cgITSTPCalone->AddCut(pTPCHPT);

		  cgITSTPCalone->AddCut(GetTrackCutsAna(cutSet));



		  anaCuts = cgITSTPCalone;
		  break;
		default: cout << "No Pre-PID Cut defined " << endl;
	  }
	  return anaCuts;
	}




	//Pair Cuts for Analysis step - take care of logic - inverted compared to other PairCuts!!
	AliAnalysisCuts* GetPairCuts2(Int_t cutSet, Bool_t togglePC /*=kFALSE*/)  {
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
				pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.05); 
				pairCutsInvMgood =new AliDielectronVarCuts("InvM Cuts good","InvM>0.3");
				pairCutsInvMgood->AddCut(AliDielectronVarManager::kM, 0.05, 99999.); 
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



	//Pair Cuts for PREFILTER step
	AliAnalysisCuts* GetPairCuts(Int_t cutSet)  {  
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
/*		case kpp2010TPCorTOF  :

		  AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
		  //AliDielectronVarCuts* pairCutsPhiv =new AliDielectronVarCuts("Phiv Cuts","Phiv<2.0rad");
		  //pairCutsPhiv->AddCut(AliDielectronVarManager::kPhivPair, 2.0, 3.2); 
		  AliDielectronVarCuts* pairCutsPhiv =new AliDielectronVarCuts("Phiv Cuts","Phiv<2.0rad");
		  pairCutsPhiv->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.05); 
		  AliDielectronVarCuts* pairCutsInvM =new AliDielectronVarCuts("InvM Cuts","InvM<0.3");
		  pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.01); 
		  pairCutsCG->AddCut(pairCutsPhiv);
		  pairCutsCG->AddCut(pairCutsInvM);
*/
		  AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
		  AliDielectronVarCuts* pairCutsPhiv =new AliDielectronVarCuts("Phiv Cuts","Phiv<2.0rad");
		  pairCutsPhiv->AddCut(AliDielectronVarManager::kPhivPair, 2.0, 3.2); 
		  //AliDielectronVarCuts* pairCutsPhiv =new AliDielectronVarCuts("Phiv Cuts","Phiv<2.0rad");
		  //pairCutsPhiv->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.05); 
		  AliDielectronVarCuts* pairCutsInvM =new AliDielectronVarCuts("InvM Cuts","InvM<0.3");
		  pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.05); 
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





	//*******************************************************************************
	//*******************************************************************************
	//** ESD TRACK CUTS TUNED FOR AGREEMENT BETWEEN AODS AND ESDS  ******************
	//** NOT NECESSARILY 100% OPTIMIZED FOR DIEL-ANALYSIS          ******************
	//*******************************************************************************
	//*******************************************************************************

	//WHEN RUNNING ON ESDs: LOAD Default Cuts for AODs
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
		default: cout << "No Analysis Track Cut defined " << endl;
	  }
	  return esdTrackCutsH;
	} 


	//Make/Tighten track Cuts that are *NOT* already
	//done in the AOD production
	//**IMPORTANT**: For AODs, select FilterBit
	//the method is ignored for ESDs
	
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
//			trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,     80., 140.0);
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


	//Possibly different cut sets for Prefilter step
	//Not used at the moment
	AliAnalysisCuts* GetTrackCutsPre(Int_t cutSet) {
	  AliDielectronCutGroup* trackCuts=0x0;
	  switch (cutSet) {
		case kPbPb2011NoPID:
		case kPbPb2011TPCandTOF :
   		case kPbPb2011TPCandTOFHPT:
		case kPbPb2011TPCorTOF  :
		case kpp2010TPCandTOF :
		case kpp2010TPCorTOF  :
		   trackCuts = new AliDielectronCutGroup("cgPIDTPC1","cgPIDTPC1",AliDielectronCutGroup::kCompAND);

		   AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
		   trackCutsAOD->AddCut(AliDielectronVarManager::kPt,0.05,0.2);
		   trackCutsAOD->AddCut(AliDielectronVarManager::kEta,-0.84,0.84);
		   //DCA Cut
		   trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
		   trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
		   trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
		   AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
		   trackCutsDiel->SetAODFilterBit(1); //does nothing for ESDs, ITSSA
		   trackCutsDiel->SetRequireITSRefit(kTRUE);

		   trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
		   trackCuts->AddCut(trackCutsDiel);
		   trackCuts->AddCut(trackCutsAOD);
		   //		  cout << "No Pre-Track Cut defined for AODs at the moment " << endl;
		  break;
		default: cout << "No Pre-Track Cut defined " << endl;
	  }
	  return trackCuts;
	}

};
