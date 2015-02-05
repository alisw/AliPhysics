AliAnalysisTask* AddTaskKinksFilimon(const char* const taskName, const Bool_t useMC, const Short_t lCollidingSystems=0 /*0 = pp, 1 = PbPb, 2 = pPb, 3 = AA*/, const AliVEvent::EOfflineTriggerTypes offlineTriggerType=AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral) {
#if 0
 const AliAnalysisTaskKinksFilimon::ECollisionType collisionType = AliAnalysisTaskKinksFilimon::kPP;
 switch (lCollidingSystems) {
	 case 0:
		 collisionType = AliAnalysisTaskKinksFilimon::kPP;
		 break;
	 case 1:
		 collisionType = AliAnalysisTaskKinksFilimon::kPbPb;
		 break;
	 case 2:
		 collisionType = AliAnalysisTaskKinksFilimon::kPPb;
		 break;
	 case 3:
		 collisionType = AliAnalysisTaskKinksFilimon::kAA;
		 break;
	 default:
		 break;
 };
#endif
 Bool_t fillCutHist = kTRUE;
 Double_t minPt=0.2;
 Double_t absY07=0.7;
 Double_t absY05=0.5;
 Double_t absEta=1.1;
 /*AliCFTrackKineCuts* commonKineTrackCuts = new AliCFTrackKineCuts("commonKineTrackCuts", "commonKineTrackCuts"); // MSS compatible (to be checked)
 commonKineTrackCuts->SetPtRange(minPt);
 commonKineTrackCuts->SetRapidityRange(-absY07, +absY07);
 commonKineTrackCuts->SetRapidityRange(-absEta, +absEta);*/
 AliESDtrackCuts* esdTrackCutsKinkMother = new AliESDtrackCuts("esdTrackCutsKinkMother", "esdTrackCutsKinkMother"); // MSS compatible (to be checked)
 esdTrackCutsKinkMother->SetPtRange(minPt, 1e10);
 //esdTrackCutsKinkMother->SetEtaRange(-absEta, +absEta);
 esdTrackCutsKinkMother->SetRapRange(-absY07, +absY07);
        //fCutsMul=new AliESDtrackCuts("Mul","Mul");
        //fCutsMul->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
        //                                      AliESDtrackCuts::kAny);
 AliESDtrackCuts* esdTrackCutsKinkDaughter = esdTrackCutsKinkMother->Clone("esdTrackCutsKinkDaughter"); // MSS compatible (to be checked)
 AliESDtrackCuts* esdTrackCutsKinkPartner = esdTrackCutsKinkMother->Clone("esdTrackCutsKinkPartner"); // MSS compatible (to be checked)
 // Kink mother specific ESD cuts
 esdTrackCutsKinkMother->SetMinNClustersTPC(20);
 esdTrackCutsKinkMother->SetMaxChi2PerClusterTPC(4);
 esdTrackCutsKinkMother->SetRequireITSRefit(kTRUE);
 esdTrackCutsKinkMother->SetRequireTPCRefit(kTRUE);
 esdTrackCutsKinkMother->SetAcceptKinkDaughters(kFALSE);
 esdTrackCutsKinkMother->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
 esdTrackCutsKinkMother->SetMaxDCAToVertexZ(2);
 esdTrackCutsKinkMother->SetMaxChi2TPCConstrainedGlobal(36);
 esdTrackCutsKinkMother->SetHistogramsOn(kTRUE);
 esdTrackCutsKinkMother->SaveHistograms();
 //fCutsMul->SetDCAToVertex2D(kFALSE);
 //esdTrackCutsKinkMother->SetRequireSigmaToVertex(kFALSE);
 // Kink daughter specific ESD cuts
 esdTrackCutsKinkDaughter->SetMinNClustersTPC(10);
 //esdTrackCutsKinkDaughter->SetMaxChi2PerClusterTPC(4);
 //esdTrackCutsKinkDaughter->SetRequireITSRefit(kTRUE);
 //esdTrackCutsKinkDaughter->SetRequireTPCRefit(kTRUE);
 esdTrackCutsKinkDaughter->SetAcceptKinkDaughters(kTRUE);
 //esdTrackCutsKinkDaughter->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
 //esdTrackCutsKinkDaughter->SetMaxDCAToVertexZ(2);
 //esdTrackCutsKinkDaughter->SetMaxChi2TPCConstrainedGlobal(36);
 esdTrackCutsKinkDaughter->SetHistogramsOn(kTRUE);
 esdTrackCutsKinkDaughter->SaveHistograms();
 // Resonance partner specific ESD cuts
 esdTrackCutsKinkPartner->SetMinNClustersTPC(60);
 esdTrackCutsKinkPartner->SetMaxChi2PerClusterTPC(4);
 esdTrackCutsKinkPartner->SetRequireITSRefit(kTRUE);
 esdTrackCutsKinkPartner->SetRequireTPCRefit(kTRUE);
 esdTrackCutsKinkPartner->SetAcceptKinkDaughters(kFALSE);
 esdTrackCutsKinkPartner->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
 esdTrackCutsKinkPartner->SetMaxDCAToVertexZ(2);
 esdTrackCutsKinkPartner->SetMaxChi2TPCConstrainedGlobal(36);
 esdTrackCutsKinkPartner->SetHistogramsOn(kTRUE);
 esdTrackCutsKinkPartner->SaveHistograms();
 
 AliESDtrackCuts* esdTrackCutsKinkMother_mothertpcnclu50 = esdTrackCutsKinkMother->Clone("esdTrackCutsKinkMother_mothertpcnclu50"); // MSS compatible (to be checked)
 esdTrackCutsKinkMother_mothertpcnclu50->SetMinNClustersTPC(50);
 //esdTrackCutsKinkMother->Dump();
 //esdTrackCutsKinkDaughter->Dump();
 //esdTrackCutsKinkPartner->Dump();

//
 AliESDtrackCuts* esdTrackCutsKinkMother_y05 = new AliESDtrackCuts("esdTrackCutsKinkMother_y05", "esdTrackCutsKinkMother_y05"); // MSS compatible (to be checked)
 esdTrackCutsKinkMother_y05->SetPtRange(minPt, 1e10);
 //esdTrackCutsKinkMother_y05->SetEtaRange(-absEta, +absEta);
 esdTrackCutsKinkMother_y05->SetRapRange(-absY05, +absY05);
        //fCutsMul=new AliESDtrackCuts("Mul","Mul");
        //fCutsMul->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
        //                                      AliESDtrackCuts::kAny);
 AliESDtrackCuts* esdTrackCutsKinkDaughter_y05 = esdTrackCutsKinkMother_y05->Clone("esdTrackCutsKinkDaughter_y05"); // MSS compatible (to be checked)
 AliESDtrackCuts* esdTrackCutsKinkPartner_y05 = esdTrackCutsKinkMother_y05->Clone("esdTrackCutsKinkPartner_y05"); // MSS compatible (to be checked)
 // Kink mother specific ESD cuts
 esdTrackCutsKinkMother_y05->SetMinNClustersTPC(20);
 esdTrackCutsKinkMother_y05->SetMaxChi2PerClusterTPC(4);
 esdTrackCutsKinkMother_y05->SetRequireITSRefit(kTRUE);
 esdTrackCutsKinkMother_y05->SetRequireTPCRefit(kTRUE);
 esdTrackCutsKinkMother_y05->SetAcceptKinkDaughters(kFALSE);
 esdTrackCutsKinkMother_y05->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
 esdTrackCutsKinkMother_y05->SetMaxDCAToVertexZ(2);
 esdTrackCutsKinkMother_y05->SetMaxChi2TPCConstrainedGlobal(36);
 esdTrackCutsKinkMother_y05->SetHistogramsOn(kTRUE);
 esdTrackCutsKinkMother_y05->SaveHistograms();
 //fCutsMul->SetDCAToVertex2D(kFALSE);
 //esdTrackCutsKinkMother_y05->SetRequireSigmaToVertex(kFALSE);
 // Kink daughter specific ESD cuts
 esdTrackCutsKinkDaughter_y05->SetMinNClustersTPC(10);
 //esdTrackCutsKinkDaughter_y05->SetMaxChi2PerClusterTPC(4);
 //esdTrackCutsKinkDaughter_y05->SetRequireITSRefit(kTRUE);
 //esdTrackCutsKinkDaughter_y05->SetRequireTPCRefit(kTRUE);
 esdTrackCutsKinkDaughter_y05->SetAcceptKinkDaughters(kTRUE);
 //esdTrackCutsKinkDaughter_y05->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
 //esdTrackCutsKinkDaughter_y05->SetMaxDCAToVertexZ(2);
 //esdTrackCutsKinkDaughter_y05->SetMaxChi2TPCConstrainedGlobal(36);
 esdTrackCutsKinkDaughter_y05->SetHistogramsOn(kTRUE);
 esdTrackCutsKinkDaughter_y05->SaveHistograms();
 // Resonance partner specific ESD cuts
 esdTrackCutsKinkPartner_y05->SetMinNClustersTPC(60);
 esdTrackCutsKinkPartner_y05->SetMaxChi2PerClusterTPC(4);
 esdTrackCutsKinkPartner_y05->SetRequireITSRefit(kTRUE);
 esdTrackCutsKinkPartner_y05->SetRequireTPCRefit(kTRUE);
 esdTrackCutsKinkPartner_y05->SetAcceptKinkDaughters(kFALSE);
 esdTrackCutsKinkPartner_y05->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
 esdTrackCutsKinkPartner_y05->SetMaxDCAToVertexZ(2);
 esdTrackCutsKinkPartner_y05->SetMaxChi2TPCConstrainedGlobal(36);
 esdTrackCutsKinkPartner_y05->SetHistogramsOn(kTRUE);
 esdTrackCutsKinkPartner_y05->SaveHistograms();
 
 //esdTrackCutsKinkMother_y05->Dump();
 //esdTrackCutsKinkDaughter_y05->Dump();
 //esdTrackCutsKinkPartner_y05->Dump();

//
 
                   //fMaxDCAtoVtxCut=new AliESDtrackCuts("fMaxDCA", "fMaxDCA");
       //fMaxDCAtoVtxCut->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

		const Int_t aliPIDnSpecies = AliPID::kSPECIES;
		/*AliESDpidCuts* esdPIDcutsKinks=new AliESDpidCuts("esdPIDcutsKinks", "esdPIDcutsKinks");
		esdPIDcutsKinks->SetTPCnSigmaCut(AliPID::kKaon, 3.5);
		esdPIDcutsKinks->SetTPCnSigmaCut(AliPID::kPion, 3.5);
		esdPIDcutsKinks->SetTPCnSigmaCut(AliPID::kMuon, 3.5);*/
		TArrayF esdPIDcutsKinksArray(aliPIDnSpecies);
		esdPIDcutsKinksArray[AliPID::kKaon]=3.5;
		esdPIDcutsKinksArray[AliPID::kPion]=3.5;
		esdPIDcutsKinksArray[AliPID::kMuon]=3.5;

		TArrayF esdPIDcutsKinksArray_nodedx(aliPIDnSpecies);
		esdPIDcutsKinksArray_nodedx[AliPID::kKaon]=15.0;
		esdPIDcutsKinksArray_nodedx[AliPID::kPion]=15.0;
		esdPIDcutsKinksArray_nodedx[AliPID::kMuon]=15.0;
				
		/*AliESDpidCuts* esdPIDcutsResonances=new AliESDpidCuts("esdPIDcutsResonances", "esdPIDcutsResonances");
		esdPIDcutsResonances->SetTPCnSigmaCut(AliPID::kKaon, 15.0);
		esdPIDcutsResonances->SetTPCnSigmaCut(AliPID::kPion, 15.0);
		esdPIDcutsResonances->SetTPCnSigmaCut(AliPID::kMuon, 15.0);
		esdPIDcutsResonances->SetTPCnSigmaCut(AliPID::kProton, 4.0);*/
		TArrayF esdPIDcutsResonancesArray(aliPIDnSpecies);
		esdPIDcutsResonancesArray[AliPID::kKaon]=4.0;
		esdPIDcutsResonancesArray[AliPID::kPion]=4.0;
		esdPIDcutsResonancesArray[AliPID::kMuon]=4.0;
		esdPIDcutsResonancesArray[AliPID::kProton]=4.0;

		TArrayF esdPIDcutsResonancesArray_nodedx(aliPIDnSpecies);
		esdPIDcutsResonancesArray_nodedx[AliPID::kKaon]=15.0;
		esdPIDcutsResonancesArray_nodedx[AliPID::kPion]=15.0;
		esdPIDcutsResonancesArray_nodedx[AliPID::kMuon]=15.0;
		esdPIDcutsResonancesArray_nodedx[AliPID::kProton]=15.0;
				
	Float_t hLoPt = 0, hHiPt = 6, hLoY = -1.0, hHiY = 1.0, hLoEta = -1.5, hHiEta = 1.5;//, fLowX=0.6, fHighX=1.8;
	Int_t nBinsPt = TMath::Nint(5*(hHiPt-hLoPt)), nBinsPtInvMassSlice = nBinsPt, nBinsY = TMath::Nint(10*(hHiY-hLoY)), nBinsEta = TMath::Nint(10*(hHiEta-hLoEta));
	
  Double_t gPt7K0[] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
                         2.2, 2.4, 2.6, 2.8,  3.0,   3.3, 3.6, 3.9,
                         4.2, 4.6,5.0, 5.4, 5.9,  6.5,   7.0,7.5, 8.0,8.5,  9.2, 10., 11., 12., 13.5,15.0 };  // David K0
  Int_t nPt7K0 = sizeof(gPt7K0)/sizeof(Double_t)-1;
  Double_t gPt7TOF[] = { 0.2,0.25, 0.3,0.35,  0.4,0.45,  0.5,0.55,  0.6,0.65,  0.7,0.75,  0.8, 0.85, 0.9, 0.95, 1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
                         2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,2.8, 2.9, 3.0,
                         3.2, 3.4, 3.6, 3.8, 4.0, /*4.2, 4.4, 4.6, 4.8,*/4.5,5.0,/* 5.5,*/ 6.0 };  //  Barbara TOF Kch
  Int_t nPt7TOF = sizeof(gPt7TOF)/sizeof(Double_t)-1;
	Double_t gPt7Comb[] = {
		0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,
		1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 
		3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0, 5.5, 6 };  // 25/11/2013 from Francesco
  Int_t nPt7Comb = sizeof(gPt7Comb)/sizeof(Double_t)-1;
                         
	TH2* histPtYTemplate = new TH2F("histPtYTemplate", "histPtYTemplate", nPt7Comb/*nBinsPt*/, gPt7Comb,/*hLoPt, hHiPt,*/ nBinsY, hLoY, hHiY);
	histPtYTemplate->SetOption("colz");
/*	 
	 Int_t bins[3] = {250, 40, 300};
   Double_t xmin[3] = {0., -4.0, -300.};
   Double_t xmax[3] = {500., 4.0, 300.};
   THnSparseF* histRPhiZTemplate = new THnSparseF("histRPhiZTemplate", "histRPhiZTemplate", 3, bins, xmin, xmax);
*/

  //input hander
  AliAnalysisManager* man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) AliFatal("Input handler needed");
  
  //pid response object
  //inputHandler->CreatePIDResponse(fIsMC);
  fPIDResponse=inputHandler->GetPIDResponse();
  //if (!fPIDResponse) AliFatal("PIDResponse object was not created");

	//AliESDkinkCuts(const char* name, const char* title, AliESDtrackCuts* esdTrackCutsKinkMother, AliESDtrackCuts* esdTrackCutsKinkDaughter/*, AliKineCuts* fiducialDecayCut*/, const AliPIDResponse* pidResponse, /*AliESDpidCuts**/ const TArrayF esdPIDcutsKinks, const Bool_t fillCutHist, const Double_t minPt, const Double_t maxAbsY, const Double_t minQt, const Double_t maxQt, const Double_t minQtPion, const Double_t maxQtPion, const Double_t minKinkR=120.0, const Double_t maxKinkR=210.0, const Double_t minAbsKinkZ=0.1, const Double_t maxAbsKinkZ=225.0, const Double_t minAngleKaon=2.0, const Double_t minAnglePion=1.0)
  AliESDkinkCuts* fRecKinkCutsKaon = new AliESDkinkCuts("fRecKinkCutsKaon", "Rec Kaon Kink Cuts", esdTrackCutsKinkMother, esdTrackCutsKinkDaughter, fPIDResponse, esdPIDcutsKinksArray, fillCutHist, minPt, absY07, 0.12/*minQt*/, 0.30/*maxQt*/, 0.020/*minQtPion*/, 0.040/*maxQtPion*/); // Need to have this init here for the moment as AliESDkinkCuts is a "local" class
  AliESDkinkCuts* fRecKinkCutsKaon_y05 = new AliESDkinkCuts("fRecKinkCutsKaon_y05", "Rec Kaon Kink Cuts", esdTrackCutsKinkMother_y05, esdTrackCutsKinkDaughter_y05, fPIDResponse, esdPIDcutsKinksArray, fillCutHist, minPt, absY05, 0.12/*minQt*/, 0.30/*maxQt*/, 0.020/*minQtPion*/, 0.040/*maxQtPion*/); // Need to have this init here for the moment as AliESDkinkCuts is a "local" class
  AliESDkinkCuts* fRecKinkCutsKaon_qt40 = new AliESDkinkCuts("fRecKinkCutsKaon_qt40", "Rec Kaon Kink Cuts", esdTrackCutsKinkMother, esdTrackCutsKinkDaughter, fPIDResponse, esdPIDcutsKinksArray, fillCutHist, minPt, absY07, 0.040/*minQt*/, 0.30/*maxQt*/, 0.020/*minQtPion*/, 0.040/*maxQtPion*/); // Need to have this init here for the moment as AliESDkinkCuts is a "local" class
  AliESDkinkCuts* fRecKinkCutsKaon_r130_200 = new AliESDkinkCuts("fRecKinkCutsKaon_r130_200", "Rec Kaon Kink Cuts", esdTrackCutsKinkMother, esdTrackCutsKinkDaughter, fPIDResponse, esdPIDcutsKinksArray, fillCutHist, minPt, absY07, 0.12/*minQt*/, 0.30/*maxQt*/, 0.020/*minQtPion*/, 0.040/*maxQtPion*/, 130.0, 200.0); // Need to have this init here for the moment as AliESDkinkCuts is a "local" class
  AliESDkinkCuts* fRecKinkCutsKaon_nodedx = new AliESDkinkCuts("fRecKinkCutsKaon_nodedx", "Rec Kaon Kink Cuts", esdTrackCutsKinkMother, esdTrackCutsKinkDaughter, fPIDResponse, esdPIDcutsKinksArray_nodedx, fillCutHist, minPt, absY07, 0.12/*minQt*/, 0.30/*maxQt*/, 0.020/*minQtPion*/, 0.040/*maxQtPion*/); // Need to have this init here for the moment as AliESDkinkCuts is a "local" class
  AliESDkinkCuts* fRecKinkCutsKaon_mothertpcnclu50 = new AliESDkinkCuts("fRecKinkCutsKaon_mothertpcnclu50", "Rec Kaon Kink Cuts", esdTrackCutsKinkMother_mothertpcnclu50, esdTrackCutsKinkDaughter, fPIDResponse, esdPIDcutsKinksArray, fillCutHist, minPt, absY07, 0.12/*minQt*/, 0.30/*maxQt*/, 0.020/*minQtPion*/, 0.040/*maxQtPion*/); // Need to have this init here for the moment as AliESDkinkCuts is a "local" class
  AliESDkinkCuts* fRecKinkCutsKaon_pionqt = new AliESDkinkCuts("fRecKinkCutsKaon_pionqt", "Rec Kaon Kink Cuts", esdTrackCutsKinkMother, esdTrackCutsKinkDaughter, fPIDResponse, esdPIDcutsKinksArray, fillCutHist, minPt, absY07, 0.12/*minQt*/, 0.30/*maxQt*/, 0.010/*minQtPion*/, 0.050/*maxQtPion*/); // Need to have this init here for the moment as AliESDkinkCuts is a "local" class

  //AliAnalysisTaskKinksFilimon(const char *name, const AliVEvent::EOfflineTriggerTypes offlineTriggerType, AliESDkinkCuts* esdKinkCuts, AliESDtrackCuts* esdTrackCutsKinkMother, AliESDtrackCuts* esdTrackCutsKinkDaughter, AliESDtrackCuts* esdTrackCutsPartner, TArrayF esdPIDcutsResonances, /*const AliCFParticleGenCuts* genTrackCuts, const TH1* histPtTemplate=0,*/ const Bool_t useMC=kTRUE, const Bool_t fillCutHist=kFALSE, const TH2* histPtYTemplate=0, const THnSparseF* histRPhiZTemplate=0, const ECollisionType collisionType=kPP, TClass* outputContClass=TList::Class());
  AliAnalysisTaskKinksFilimon* taskKinksFilimon=new AliAnalysisTaskKinksFilimon(taskName, offlineTriggerType, fRecKinkCutsKaon, esdTrackCutsKinkMother, 0/*esdTrackCutsKinkDaughter*/, esdTrackCutsKinkPartner, esdPIDcutsResonancesArray, useMC, fillCutHist, histPtYTemplate, 0/*histRPhiZTemplate*/, AliAnalysisTaskKinksFilimon::kPP);
	//taskKinksFilimon->SelectCollisionCandidates(offlineTriggerType); // Better use only for AOD
	man->AddTask(taskKinksFilimon);
  AliAnalysisTaskKinksFilimon* taskKinksFilimon_y05=new AliAnalysisTaskKinksFilimon(Form("%s_y05", taskName), offlineTriggerType, fRecKinkCutsKaon_y05, esdTrackCutsKinkMother_y05, 0/*esdTrackCutsKinkDaughter_y05*/, esdTrackCutsKinkPartner_y05, esdPIDcutsResonancesArray, useMC, fillCutHist, histPtYTemplate, 0/*histRPhiZTemplate*/, AliAnalysisTaskKinksFilimon::kPP);
	//taskKinksFilimon->SelectCollisionCandidates(offlineTriggerType); // Better use only for AOD
	man->AddTask(taskKinksFilimon_y05);
  AliAnalysisTaskKinksFilimon* taskKinksFilimon_qt40=new AliAnalysisTaskKinksFilimon(Form("%s_qt40", taskName), offlineTriggerType, fRecKinkCutsKaon_qt40, esdTrackCutsKinkMother, 0/*esdTrackCutsKinkDaughter*/, esdTrackCutsKinkPartner, esdPIDcutsResonancesArray, useMC, fillCutHist, histPtYTemplate, 0/*histRPhiZTemplate*/, AliAnalysisTaskKinksFilimon::kPP);
	//taskKinksFilimon->SelectCollisionCandidates(offlineTriggerType); // Better use only for AOD
	man->AddTask(taskKinksFilimon_qt40);
  AliAnalysisTaskKinksFilimon* taskKinksFilimon_r130_200=new AliAnalysisTaskKinksFilimon(Form("%s_r130_200", taskName), offlineTriggerType, fRecKinkCutsKaon_r130_200, esdTrackCutsKinkMother, 0/*esdTrackCutsKinkDaughter*/, esdTrackCutsKinkPartner, esdPIDcutsResonancesArray, useMC, fillCutHist, histPtYTemplate, 0/*histRPhiZTemplate*/, AliAnalysisTaskKinksFilimon::kPP);
	//taskKinksFilimon->SelectCollisionCandidates(offlineTriggerType); // Better use only for AOD
	man->AddTask(taskKinksFilimon_r130_200);
  AliAnalysisTaskKinksFilimon* taskKinksFilimon_nodedx=new AliAnalysisTaskKinksFilimon(Form("%s_nodedx", taskName), offlineTriggerType, fRecKinkCutsKaon_nodedx, esdTrackCutsKinkMother, 0/*esdTrackCutsKinkDaughter*/, esdTrackCutsKinkPartner, esdPIDcutsResonancesArray, useMC, fillCutHist, histPtYTemplate, 0/*histRPhiZTemplate*/, AliAnalysisTaskKinksFilimon::kPP);
	//taskKinksFilimon->SelectCollisionCandidates(offlineTriggerType); // Better use only for AOD
	man->AddTask(taskKinksFilimon_nodedx);
  AliAnalysisTaskKinksFilimon* taskKinksFilimon_mothertpcnclu50=new AliAnalysisTaskKinksFilimon(Form("%s_mothertpcnclu50", taskName), offlineTriggerType, fRecKinkCutsKaon_mothertpcnclu50, esdTrackCutsKinkMother_mothertpcnclu50, 0/*esdTrackCutsKinkDaughter*/, esdTrackCutsKinkPartner, esdPIDcutsResonancesArray, useMC, fillCutHist, histPtYTemplate, 0/*histRPhiZTemplate*/, AliAnalysisTaskKinksFilimon::kPP);
	//taskKinksFilimon->SelectCollisionCandidates(offlineTriggerType); // Better use only for AOD
	man->AddTask(taskKinksFilimon_mothertpcnclu50);
  AliAnalysisTaskKinksFilimon* taskKinksFilimon_pionqt=new AliAnalysisTaskKinksFilimon(Form("%s_pionqt", taskName), offlineTriggerType, fRecKinkCutsKaon_pionqt, esdTrackCutsKinkMother, 0/*esdTrackCutsKinkDaughter*/, esdTrackCutsKinkPartner, esdPIDcutsResonancesArray, useMC, fillCutHist, histPtYTemplate, 0/*histRPhiZTemplate*/, AliAnalysisTaskKinksFilimon::kPP);
	//taskKinksFilimon->SelectCollisionCandidates(offlineTriggerType); // Better use only for AOD
	man->AddTask(taskKinksFilimon_pionqt);
  AliAnalysisTaskKinksFilimon* taskKinksFilimon_resnodedx=new AliAnalysisTaskKinksFilimon(Form("%s_resnodedx", taskName), offlineTriggerType, fRecKinkCutsKaon, esdTrackCutsKinkMother, 0/*esdTrackCutsKinkDaughter*/, esdTrackCutsKinkPartner, esdPIDcutsResonancesArray_nodedx, useMC, fillCutHist, histPtYTemplate, 0/*histRPhiZTemplate*/, AliAnalysisTaskKinksFilimon::kPP);
	//taskKinksFilimon->SelectCollisionCandidates(offlineTriggerType); // Better use only for AOD
	man->AddTask(taskKinksFilimon_resnodedx);
	//for ( Int_t iTask = 0; iTask < 5; ++iTask ) {
		//man->AddTask(new AliAnalysisTaskKinksFilimon(Form("%s%d", taskName, iTask), offlineTriggerType, fRecKinkCutsKaon, esdTrackCutsKinkMother, 0/*esdTrackCutsKinkDaughter*/, esdTrackCutsKinkPartner, esdPIDcutsResonancesArray, useMC, fillCutHist, histPtYTemplate, 0/*histRPhiZTemplate*/, AliAnalysisTaskKinksFilimon::kPP));
	//}
  //kinkResonanceObjectPESD->InitOutputHistograms(70, 0.99, 1.088, 36, -0.9, 0.9, 100, 0.0, 10.0);
  //kinkResonanceObjectPESD->SetPDGCodes(kKPlus, kKPlus, AliResonanceKink::kPhi);
  //AddTaskKinkResEviCommonConfig(kinkResonanceObjectPESD);

  //AliAnalysisTaskKinkResonance* taskresonancePhiESD = new AliAnalysisTaskKinkResonance(taskName);
  //taskresonancePhiESD->SetAnalysisKinkObject(kinkResonanceObjectPESD);
  return(0);

}
