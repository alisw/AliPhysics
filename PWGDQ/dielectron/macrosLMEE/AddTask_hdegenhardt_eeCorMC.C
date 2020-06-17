
AliAnalysisTask* AddTask_hdegenhardt_eeCorMC(
			char *name = 			"Output"
			,char *period = 		"17"
			,Bool_t recabPID = 		kTRUE
			,Bool_t smearing = 		kTRUE
			,Bool_t fillRecMaps = 	kTRUE
			,Bool_t sysUnc = 		kFALSE
			,Bool_t dcaSmr = 		kTRUE
			,char* dcaSmrMethod = 	"parDir" //math,maps,parDir,parGaus
			,Double_t ptMin = 		0.4
			,TString outName =	 	"LMEE.root"
			,char *smrMapsFrom = 	"17"	//16 (2016), 17 (2017), 18 (2018)
			,char *smrDCAMapsFrom = "1678" 	//16 (2016), 17 (2017), 18 (2018), 1678 (all)
			,Bool_t dcaMapsFromMC = kFALSE 	//Data in [cm] and MC in [m]
){
    
	Int_t nConfigs = 1;
	if (sysUnc) nConfigs = 25;
	
	AliAnalysisTaskeeCor *tasklmee;
	for (int config = 0; config < nConfigs; config++){
		AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

		if (!mgr) {
			::Error("AliAnalysisTaskeeCor", "No analysis manager to connect to.");
			return NULL;
		}

		if (!mgr->GetInputEventHandler()) {
			::Error("AliAnalysisTaskeeCor", "This task requires an input event handler");
			return NULL;
		}
		
		tasklmee = new AliAnalysisTaskeeCor(ConfigNames[config]);
		// tasklmee->SelectCollisionCandidates(AliVEvent::kINT7);
		
		printf("----------------------------------\n");
		printf("Track cuts: %s\n",trackCutsVar[config]);
		Double_t impactParXY = 3.;
		if (trackCutsVar[config][0] == '0') impactParXY = 5.;
		else if (trackCutsVar[config][0] == '2') impactParXY = 2.;
		Double_t impactParZ = 1.;
		if (trackCutsVar[config][1] == '0') impactParZ = 2.;
		else if (trackCutsVar[config][1] == '2') impactParZ = 0.7;
		Double_t nClusITS = 3.;
		if (trackCutsVar[config][2] == '0') nClusITS = 2.;
		else if (trackCutsVar[config][2] == '2') nClusITS = 4.;
		Double_t chi2ITS = 4.5;
		if (trackCutsVar[config][3] == '0') chi2ITS = 6.;
		else if (trackCutsVar[config][3] == '2') chi2ITS = 3.5;
		Double_t nClusTPC = 80;
		Double_t nFclsTPCr = 100;
		if (trackCutsVar[config][4] == '0'){
			nClusTPC = 60;
			nFclsTPCr = 80;
		}
		else if (trackCutsVar[config][4] == '2'){
			nClusTPC = 100;
			nFclsTPCr = 120;
		}
		Double_t nFclsTPCfCross = 0.8;
		if (trackCutsVar[config][5] == '0') nFclsTPCfCross = 0.6;
		else if (trackCutsVar[config][5] == '2') nFclsTPCfCross = 0.9;
		Double_t nClsSFracTPC = 0.4;
		if (trackCutsVar[config][6] == '0') nClsSFracTPC = 0.2;
		else if (trackCutsVar[config][6] == '2') nClsSFracTPC = 1.0;
		Double_t chi2TPC = 4.;
		if (trackCutsVar[config][7] == '0') chi2TPC = 6.;
		else if (trackCutsVar[config][7] == '2') chi2TPC = 3.;	
		
		//----- add track cuts ---------------------------------------------
		AliESDtrackCuts* trackCuts = new AliESDtrackCuts();

		//----- Fill Calibration Maps? -------------------------------------
		tasklmee->FillPIDrecMaps(fillRecMaps);
		tasklmee->FillSmrMaps(fillRecMaps);

		//----- pT and eta -------------------------------------------------
		Double_t cutPtMin = ptMin;

		trackCuts->SetPtRange(cutPtMin, 1e30);
		trackCuts->SetEtaRange(-0.8, 0.8);

		//----- TPC --------------------------------------------------------
		trackCuts->SetRequireTPCRefit(kTRUE);
		trackCuts->SetMinNCrossedRowsTPC(nFclsTPCr);
		trackCuts->SetMinNClustersTPC(nClusTPC);
		trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(nFclsTPCfCross);
		trackCuts->SetMaxChi2PerClusterTPC(chi2TPC);
		trackCuts->SetMaxFractionSharedTPCClusters(nClsSFracTPC);

		//----- ITS --------------------------------------------------------
		trackCuts->SetRequireITSRefit(kTRUE);
		trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
		trackCuts->SetMinNClustersITS(nClusITS);
		trackCuts->SetMaxChi2PerClusterITS(chi2ITS);
		     if (trackCutsVar[config][8] == '0') trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
		else if (trackCutsVar[config][8] == '1') trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
		else if (trackCutsVar[config][8] == '2') trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);


		//----- primary selection ------------------------------------------
		trackCuts->SetAcceptKinkDaughters(kFALSE);
		trackCuts->SetDCAToVertex2D(kFALSE);
		trackCuts->SetRequireSigmaToVertex(kFALSE);
		trackCuts->SetMaxDCAToVertexZ(impactParZ);
		trackCuts->SetMaxDCAToVertexXY(impactParXY);

		//----- set track cuts to the task ---------------------------------
		tasklmee->SetTrackCuts(trackCuts);
		tasklmee->SetMinPtCut(cutPtMin);

		//----- PID cuts ---------------------------------------------------
		Double_t TPCv = 0.;
		Double_t TOFv = 0.;
		if (config == 1) TPCv = -0.5;	//TPC tight
		if (config == 2) TPCv =  0.5;	//TPC loose
		if (config == 3) TOFv = -0.5;	//TOF tight
		if (config == 4) TOFv =  0.5;	//TOF loose

		tasklmee->SetEleTPCcuts(-3.-TPCv,3.+TPCv);
		tasklmee->SetEleTOFcuts(-3.-TOFv,3.+TOFv);
		tasklmee->SetPionTPCrej(-100.   ,4.-TPCv);
		tasklmee->SetProtTPCrej(-4.+TOFv,4.-TPCv);
		tasklmee->SetKaonTPCrej(-4.+TOFv,4.-TPCv);

		//------------------- INDEPENDENT OF THE CONFIG --------------------
		if (period[1] == '6') tasklmee->SetPyHeader(1); //Set position of the header to find information about CC/BB production
		else tasklmee->SetPyHeader(0);

		//----- PID Recalibration ------------------------------------------
		tasklmee->SetRecabPID(recabPID);
		if (recabPID){
			printf("##### PID Post-Calibration Enabled!\n");
			TFile *f = TFile::Open(Form("calMCmaps1%c.root",period[1]));
			if (!f){
				char rootFile[100];
				sprintf(rootFile,"alien://alice/cern.ch/user/h/hfranzde/corMC/calMCmaps1%c.root",period[1]);
				TGrid::Connect("alien://");
				gSystem->Exec(Form("alien_cp %s .",rootFile));
				f = TFile::Open(Form("calMCmaps1%c.root",period[1]));
				if (!f) return NULL;
			}
			TH2F *TPCm; f->GetObject("TPCm",TPCm);
			TH2F *TPCw; f->GetObject("TPCw",TPCw);
			TH2F *TOFm; f->GetObject("TOFm",TOFm);
			TH2F *TOFw; f->GetObject("TOFw",TOFw);
			for (Int_t i = 0; i <= TPCm->GetNbinsX()+1; i++){
				for (Int_t k = 0; k <= TPCm->GetNbinsY()+1; k++){
					if ( (i == 0) || (k == 0) || (i > TPCm->GetNbinsX()) || (k > TPCm->GetNbinsY())) { // under/overflows
						TPCm->SetBinContent(i, k, 0.0 );
						TPCw->SetBinContent(i, k, 1.0 );
					}
				}
			}
			for (Int_t i = 0; i <= TOFm->GetNbinsX()+1; i++){
				for (Int_t k = 0; k <= TOFm->GetNbinsY()+1; k++){
					if ( (i == 0) || (k == 0) || (i > TOFm->GetNbinsX()) || (k > TOFm->GetNbinsY())) { // under/overflows
						TOFm->SetBinContent(i, k, 0.0 );
						TOFw->SetBinContent(i, k, 1.0 );
					}
				}
			}

			if (TPCm){
				tasklmee->SetTPCmeanPIDmap(TPCm);
				printf(">>>>> TPCm Post-calibration added!\n");
			}
			if (TPCw){
				tasklmee->SetTPCwidthPIDmap(TPCw);
				printf(">>>>> TPCw Post-calibration added!\n");
			}
			if (TOFm){
				tasklmee->SetTOFmeanPIDmap(TOFm);
				printf(">>>>> TOFm Post-calibration added!\n");
			}
			if (TOFw){
				tasklmee->SetTOFwidthPIDmap(TOFw);
				printf(">>>>> TOFw Post-calibration added!\n\n");
			}
		}
		//----- (Pt,Eta,Phi) Smearing --------------------------------------
		tasklmee->SetSmearing(smearing);
		if (smearing){
			printf("##### Smearing (Pt,Eta,Phi) Enabled!\n");
			TFile *f2 = TFile::Open(Form("smrMaps%s.root",smrMapsFrom));
			if (!f2){
				char rootFile2[100];
				sprintf(rootFile2,"alien://alice/cern.ch/user/h/hfranzde/corMC/smrMaps%s.root",smrMapsFrom);
				TGrid::Connect("alien://");
				gSystem->Exec(Form("alien_cp %s .",rootFile2));
				f2 = TFile::Open(Form("smrMaps%s.root",smrMapsFrom));
				if (!f2) return NULL;
			}
			TH2F *smrPt; f2->GetObject("pt",smrPt);
			TH2F *smrEta; f2->GetObject("eta",smrEta);
			TH2F *smrPhiEle; f2->GetObject("phiEle",smrPhiEle);
			TH2F *smrPhiPos; f2->GetObject("phiPos",smrPhiPos);
			tasklmee->SetPtSmrMap(smrPt);
			printf(">>>>> Pt Smearing Map added!\n");
			tasklmee->SetEtaSmrMap(smrEta);
			printf(">>>>> Eta Smearing Map added!\n");
			tasklmee->SetPhiEleSmrMap(smrPhiEle);
			printf(">>>>> Phi(e-) Smearing Map added!\n");
			tasklmee->SetPhiPosSmrMap(smrPhiPos);
			printf(">>>>> Phi(e+) Smearing Map added!\n");
		}
		//----- (DCAee) Smearing --------------------------------------
		if (dcaSmr){
			if (dcaSmrMethod == "math")	 		tasklmee->SetDCASmearingByMath(kTRUE);
			else if (dcaSmrMethod == "maps") 	tasklmee->SetDCASmearingByMaps(kTRUE);
			else if (dcaSmrMethod == "parDir") 	tasklmee->SetDCASmearingByPars(kTRUE,kFALSE);//most probable value set (kFALSE - center)
			else if (dcaSmrMethod == "parGaus") tasklmee->SetDCASmearingByPars(kTRUE,kTRUE);//for gaussian distribution set (kTRUE)
			
			tasklmee->SetDCAmapsFromMC(dcaMapsFromMC);
			
			printf("\n##### Adding DCA maps...\n");
			TFile *f3 = TFile::Open(Form("smrDCA%s.root",smrDCAMapsFrom));
			if (!f3){
				char rootFile3[100];
				sprintf(rootFile3,"alien://alice/cern.ch/user/h/hfranzde/corMC/smrDCA%s.root",smrDCAMapsFrom);
				TGrid::Connect("alien://");
				gSystem->Exec(Form("alien_cp %s .",rootFile3));
				f3 = TFile::Open(Form("smrDCA%s.root",smrDCAMapsFrom));
				if (!f3) return NULL;
			}
			TH2F *smrDCApt0; f3->GetObject("smrDCApt0",smrDCApt0);
			TH2F *smrDCApt1; f3->GetObject("smrDCApt1",smrDCApt1);
			tasklmee->SetDCASmrMap0(smrDCApt0);
			printf(">>>>> DCA Smearing Map for pt0 added!\n");
			tasklmee->SetDCASmrMap1(smrDCApt1);
			printf(">>>>> DCA Smearing Map for pt1 added!\n");
			printf("##### Adding DCA pars...\n");
			TFile *f4 = TFile::Open(Form("smrDCApar%s.root",smrDCAMapsFrom));
			if (!f4){
				char rootFile4[100];
				sprintf(rootFile4,"alien://alice/cern.ch/user/h/hfranzde/corMC/smrDCApar%s.root",smrDCAMapsFrom);
				TGrid::Connect("alien://");
				gSystem->Exec(Form("alien_cp %s .",rootFile4));
				f4 = TFile::Open(Form("smrDCApar%s.root",smrDCAMapsFrom));
				if (!f4) return NULL;
			}
			TH1D *smrDCAcen; f4->GetObject("pCen",smrDCAcen);
			TH1D *smrDCAsig; f4->GetObject("pRes",smrDCAsig);
			TH1D *smrDCAmax; f4->GetObject("pMax1",smrDCAmax);
			tasklmee->SetDCASmrParCen(smrDCAcen);
			printf(">>>>> DCA Smearing Center added!\n");
			tasklmee->SetDCASmrParSig(smrDCAsig);
			printf(">>>>> DCA Smearing Sigma added!\n");
			tasklmee->SetDCASmrParMax(smrDCAmax);
			printf(">>>>> DCA Smearing Maximum added!\n\n");
		}
		//-----------------------------------------------------------
		
		mgr->AddTask(tasklmee);

		// connect the manager to the task
		mgr->ConnectInput(tasklmee,0,mgr->GetCommonInputContainer());
		mgr->ConnectOutput(tasklmee,1,mgr->CreateContainer(Form("%s",ConfigNames[config]), TList::Class(), AliAnalysisManager::kOutputContainer, outName.Data()));

		//return tasklmee;
	}
	return tasklmee;
}

char *trackCutsVar[25];
trackCutsVar[0] = "1111111111"; //Default

trackCutsVar[1] = "1111111111"; //TPC loose
trackCutsVar[2] = "1111111111"; //TPC tight
trackCutsVar[3] = "1111111111"; //TOF loose
trackCutsVar[4] = "1111111111"; //TOF tight

trackCutsVar[5] = "0120002021"; //track01
trackCutsVar[6] = "1212202100";
trackCutsVar[7] = "2121200002";
trackCutsVar[8] = "2012011220";
trackCutsVar[9] = "0112021112"; //track05
trackCutsVar[10] = "2022210010";
trackCutsVar[11] = "0202221220";
trackCutsVar[12] = "2111001020";
trackCutsVar[13] = "0200221122";
trackCutsVar[14] = "2200010100"; //track10
trackCutsVar[15] = "1020110110";
trackCutsVar[16] = "0001120212";
trackCutsVar[17] = "1222122110";
trackCutsVar[18] = "2201101100";
trackCutsVar[19] = "2222212002"; //track15
trackCutsVar[20] = "2212102120";
trackCutsVar[21] = "1212100100";
trackCutsVar[22] = "2101212101";
trackCutsVar[23] = "2002120100";
trackCutsVar[24] = "2201210011"; //track20

char *ConfigNames[25];
ConfigNames[0] = "pt200_TPCTOFcombITSshared";
ConfigNames[1] = "TPCtight";
ConfigNames[2] = "TPCloose";
ConfigNames[3] = "TOFtight";
ConfigNames[4] = "TOFloose";
ConfigNames[5] = "track1";
ConfigNames[6] = "track2";
ConfigNames[7] = "track3";
ConfigNames[8] = "track4";
ConfigNames[9] = "track5";
ConfigNames[10] = "track6";
ConfigNames[11] = "track7";
ConfigNames[12] = "track8";
ConfigNames[13] = "track9";
ConfigNames[14] = "track10";
ConfigNames[15] = "track11";
ConfigNames[16] = "track12";
ConfigNames[17] = "track13";
ConfigNames[18] = "track14";
ConfigNames[19] = "track15";
ConfigNames[20] = "track16";
ConfigNames[21] = "track17";
ConfigNames[22] = "track18";
ConfigNames[23] = "track19";
ConfigNames[24] = "track20";
