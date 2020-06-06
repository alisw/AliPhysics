AliAnalysisTaskBJetTC *AddTaskBJetTC(
	const char *ntracks = "PicoTracks",
	const char *nclusters = "",
	const char *njets = "Jets",
	const char *nrho = "",
	Double_t jetradius = 0.4,
	Bool_t isMC = kFALSE,
	const char *type = "TPC",
	const char *taskname = "AliAnalysisTaskBJetTC",
	const char *njetsMC = "Jets",
	const char *nrhoMC = "RhoMC",
	Bool_t DoSVAnalysis = kFALSE,
	Bool_t DoPtRelAna = kFALSE,
	Bool_t DoJetProb = kFALSE,
	TString pathToResolFunc = "",
	TString pathToResolFuncb = "",
	TString pathToResolFuncc = "",
	TString pathToResolFunclf = "",
	TString pathToCorrFuncPscat = "",
	TString pathToCorrFuncNvtxContrib = "",
	Bool_t V0PhotonRejection = kFALSE,
	TString cutnumberAODBranch = "", // cutnumber for AOD branch
	Int_t ptHardBin = -999,
	const char *suffix = "")
{

	// Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
		return NULL;
	}

	// Check the analysis type using the event handlers connected to the analysis manager.
	//==============================================================================
	if (!mgr->GetInputEventHandler())
	{
		::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
		return NULL;
	}

	TString V0ReaderName = "";
	if (V0PhotonRejection)
	{

		//=========  Set Cutnumber for V0Reader ================================
		TString cutnumberPhoton = "10000029200000003220400000";
		TString cutnumberEvent = "80010103";

		Bool_t doEtaShift = kFALSE;
		AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
		//========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
		V0ReaderName = Form("V0ReaderV1_%s_%s", cutnumberEvent.Data(), cutnumberPhoton.Data());
		if (!(AliV0ReaderV1 *)mgr->GetTask(V0ReaderName.Data()))
		{
			AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
			fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
			fV0ReaderV1->SetCreateAODs(kFALSE); // AOD Output
			fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
			if (!mgr)
			{
				Error("AddTask_V0ReaderV1", "No analysis manager found.");
				return NULL;
			}

			AliConvEventCuts *fEventCuts = NULL;
			if (cutnumberEvent != "")
			{
				fEventCuts = new AliConvEventCuts(cutnumberEvent.Data(), cutnumberEvent.Data());
				fEventCuts->SetPreSelectionCutFlag(kTRUE);
				fEventCuts->SetPeriodEnumExplicit(AliConvEventCuts::kLHC16qt);
				fEventCuts->SetV0ReaderName(V0ReaderName);
				fEventCuts->SetLightOutput(kTRUE);
				if (fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data()))
				{
					fEventCuts->DoEtaShift(doEtaShift);
					fV0ReaderV1->SetEventCuts(fEventCuts);
					fEventCuts->SetFillCutHistograms("", kTRUE);
				}
			}

			// Set AnalysisCut Number
			AliConversionPhotonCuts *fCuts = NULL;
			if (cutnumberPhoton != "")
			{
				fCuts = new AliConversionPhotonCuts(cutnumberPhoton.Data(), cutnumberPhoton.Data());
				fCuts->SetPreSelectionCutFlag(kTRUE);
				fCuts->SetIsHeavyIon(kTRUE);
				fCuts->SetV0ReaderName(V0ReaderName);
				fCuts->SetLightOutput(kTRUE);
				fCuts->SetProcessAODCheck(kTRUE);
				if (fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data()))
				{
					cout << " Photon conversion is working \n";
					fV0ReaderV1->SetConversionCuts(fCuts);
					fCuts->SetFillCutHistograms("", kTRUE);
				}
			}
			if ((mgr->GetInputEventHandler())->IsA() == AliAODInputHandler::Class())
			{
				// AOD mode
				cout << "AOD handler: adding " << cutnumberAODBranch.Data() << " as conversion branch" << endl;
				fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma", cutnumberAODBranch.Data()));
			}
			fV0ReaderV1->Init();

			AliLog::SetGlobalLogLevel(AliLog::kInfo);

			//connect input V0Reader
			mgr->AddTask(fV0ReaderV1);
			mgr->ConnectInput(fV0ReaderV1, 0, cinput);
		}
	}

	TString name(taskname);

	TString combinedName;
	combinedName.Form("%s%s", name.Data(), suffix);

	AliAnalysisTaskBJetTC *jetTask = new AliAnalysisTaskBJetTC(combinedName);

	jetTask->SetPtHardBin(ptHardBin);

	jetTask->UseGammaV0Rejection(V0PhotonRejection);

	jetTask->DoJetProbabilityAnalysis(DoJetProb);

	jetTask->DoPtRelAnalysis(DoPtRelAna);

	AliTrackContainer *trackCont = jetTask->AddTrackContainer(ntracks);
	trackCont->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter);
	trackCont->SetAODFilterBits((1 << 4) | (1 << 9));

	AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);

	TString strType(type);

	AliJetContainer *jetCont = jetTask->AddJetContainer(njets, strType, jetradius);

	if (jetCont)
	{
		jetCont->SetRhoName(nrho);

		jetCont->ConnectParticleContainer(trackCont);

		jetCont->ConnectClusterContainer(clusterCont);

		jetCont->SetJetEtaLimits(-0.5, 0.5);

		jetCont->SetJetPtCut(0.0);

		jetCont->SetMaxTrackPt(1000);

		jetCont->SetPercAreaCut(0.6);
	}

	if (isMC)
	{
		AliJetContainer *jetContMC = jetTask->AddJetContainer(njetsMC, strType, jetradius);

		if (jetContMC)
		{

			jetContMC->SetRhoName(nrhoMC);

			jetContMC->SetIsParticleLevel(kTRUE);

			jetContMC->SetJetEtaLimits(-0.5, 0.5);

			jetContMC->SetJetPtCut(0.0);

			jetContMC->SetMaxTrackPt(1000);

			jetContMC->SetPercAreaCut(0.6);
		}
	}
	//-------------------------------------------------------
	//  Configure analysis task
	//-------------------------------------------------------

	jetTask->SetIsPythia(isMC);

	jetTask->SetDoSVAnalysis(DoSVAnalysis);
	jetTask->SetDoTCAnalysis(kTRUE);

	if (V0PhotonRejection)
		jetTask->SetV0ReaderName(V0ReaderName);

	if (DoJetProb && !pathToResolFunc.IsNull())
	{
		TFile *file = TFile::Open(pathToResolFunc.Data());
		for (int i = 0; i < 7; i++)
		{
			jetTask->SetResFunction((TF1 *)file->Get(Form("QualityClass%i", i)), i);
		}

		if (isMC)
		{
			TFile *fileb = TFile::Open(pathToResolFuncb.Data());
			TFile *filec = TFile::Open(pathToResolFuncc.Data());
			TFile *filelf = TFile::Open(pathToResolFunclf.Data());

			for (int i = 0; i < 7; i++)
			{
				jetTask->SetResFunctionb((TF1 *)fileb->Get(Form("QualityClass%i", i)), i);
				jetTask->SetResFunctionc((TF1 *)filec->Get(Form("QualityClass%i", i)), i);
				jetTask->SetResFunctionlf((TF1 *)filelf->Get(Form("QualityClass%i", i)), i);
			}
		}
	}

	if (!pathToCorrFuncPscat.IsNull())
	{
		TFile *file = TFile::Open(pathToCorrFuncPscat.Data());
		for (int i = 0; i < 5; i++)
		{
			if ((TH1D *)file->Get(Form("fno_%i", i)))
			{
				TF1 *CorrectionFunction = (TF1 *)file->Get(Form("fno_%i", i));
				CorrectionFunction->SetName(Form("IPsVsPscat_%i", i));
				jetTask->SetCorrectionFunctionPscat(CorrectionFunction, i);
			}
		}
	}

	if (!pathToCorrFuncNvtxContrib.IsNull())
	{
		TFile *file = TFile::Open(pathToCorrFuncNvtxContrib.Data());
		for (int i = 0; i < 5; i++)
		{
			if ((TH1D *)file->Get(Form("fno_%i", i)))
			{
				TF1 *CorrectionFunction = (TF1 *)file->Get(Form("fno_%i", i));
				CorrectionFunction->SetName(Form("IPsVsNvtxContrib_%i", i));
				jetTask->SetCorrectionFunctionNvtxContrib(CorrectionFunction, i);
			}
		}
	}

	//-------------------------------------------------------
	// Final settings, pass to manager and set the containers
	//-------------------------------------------------------

	mgr->AddTask(jetTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

	TString contname(combinedName);
	contname += "_histos";
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
															  TList::Class(), AliAnalysisManager::kOutputContainer,
															  Form("%s", AliAnalysisManager::GetCommonFileName()));

	mgr->ConnectInput(jetTask, 0, cinput1);

	mgr->ConnectOutput(jetTask, 1, coutput1);

	return jetTask;
}