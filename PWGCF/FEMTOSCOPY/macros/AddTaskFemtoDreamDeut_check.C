#include "TROOT.h"
#include "TSystem.h"

AliAnalysisTaskSE *AddTaskFemtoDreamDeut_check(bool isMC = false,						 // 1
																							 TString trigger = "kINT7",		 // 2
																							 bool fullBlastQA = true,			 // 3
																							 bool SystematicLowpT = false, // 4
																							 TString suffix = "")
{
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		printf("No analysis manager to connect to!\n");
		return nullptr;
	}
	if (!mgr->GetInputEventHandler())
	{
		printf("This task requires an input event handler!\n");
		return nullptr;
	}

	AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
	evtCuts->CleanUpMult(false, false, false, true);

	AliFemtoDreamTrackCuts *TrackCutsDeuteronDCA = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
			isMC, true, false, false);
	TrackCutsDeuteronDCA->SetCutCharge(1);
	TrackCutsDeuteronDCA->SetPtRange(0.8, 1.4);
	TrackCutsDeuteronDCA->SetFilterBit(128);
	TrackCutsDeuteronDCA->SetPIDkd(false, false, 3, 3);

	AliFemtoDreamTrackCuts *TrackCutsDeuteronMass = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
			isMC, true, false, false);
	TrackCutsDeuteronMass->SetCutCharge(1);
	TrackCutsDeuteronMass->SetPID(AliPID::kDeuteron, 999.0);

	AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronDCA = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
			isMC, true, false, false);
	TrackCutsAntiDeuteronDCA->SetCutCharge(-1);
	TrackCutsAntiDeuteronDCA->SetPtRange(0.8, 1.4);
	TrackCutsAntiDeuteronDCA->SetFilterBit(128);
	TrackCutsAntiDeuteronDCA->SetPIDkd(false, false, 3, 3);

	AliFemtoDreamTrackCuts *TrackCutsAntiDeuteronMass = AliFemtoDreamTrackCuts::PrimDeuteronCuts(
			isMC, true, false, false);
	TrackCutsAntiDeuteronMass->SetCutCharge(-1);
	TrackCutsAntiDeuteronMass->SetPID(AliPID::kDeuteron, 999.0);

	AliFemtoDreamTrackCuts *TrackCutsProtonDCA = AliFemtoDreamTrackCuts::PrimProtonCuts(
			isMC, true, false, false);
	TrackCutsProtonDCA->SetCutCharge(1);

	AliFemtoDreamTrackCuts *TrackCutsAntiProtonDCA = AliFemtoDreamTrackCuts::PrimProtonCuts(
			isMC, true, false, false);
	TrackCutsAntiProtonDCA->SetCutCharge(-1);

	std::vector<int> PDGParticles = {2212, 2212, 1000010020, 1000010020};

	std::vector<bool> closeRejection;
	std::vector<float> mTBins = {1.14, 1.26, 999.};
	std::vector<int> pairQA;
	// pairs:
	//  pp             0
	//  p bar p        1
	//  p d            2
	//  p bar d        3
	//  bar p bar p    4
	//  bar p d        5
	//  bar p bar d    6
	//  d d            7
	//  d bar d        8
	//  bar d bar d    9
	const int nPairs = 10;
	for (int i = 0; i < nPairs; ++i)
	{
		closeRejection.push_back(false);
		pairQA.push_back(0);
	}
	closeRejection[0] = true; // pp
	closeRejection[2] = true; // pd
	closeRejection[4] = true; // barp barp
	closeRejection[6] = true; // barp bard
	closeRejection[7] = true; // dd
	closeRejection[9] = true; // bard bar
	pairQA[0] = 11;						// pp
	pairQA[2] = 11;						// pd
	pairQA[4] = 11;						// barp barp
	pairQA[6] = 11;						// barp bard
	pairQA[7] = 11;						// dd
	pairQA[9] = 11;						// bard bard
	// We need to set the ZVtx bins
	std::vector<float> ZVtxBin = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};
	// The Multiplicity bins are set here
	std::vector<int> MultBins = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80};

	std::vector<int> NBins = {750, 750, 750, 750, 750, 750, 750, 750, 750, 750};

	// minimum k* value
	std::vector<float> kMin = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
	// maximum k* value
	std::vector<float> kMax = {3., 3., 3., 3., 3., 3., 3., 3., 3., 3.};

	AliFemtoDreamCollConfig *config = new AliFemtoDreamCollConfig("Femto", "Femto", false);
	config->SetZBins(ZVtxBins);
	config->SetMultBins(MultBins);
	config->SetMultBinning(true);
	config->SetPDGCodes(PDGParticles);
	config->SetNBinsHist(NBins);
	config->SetMinKRel(kMin);
	config->SetMaxKRel(kMax);
	config->SetClosePairRejection(closeRejection);
	config->SetExtendedQAPairs(pairQA);
	config->SetDeltaEtaMax(0.017); // and here you set the actual values
	config->SetDeltaPhiMax(0.017); // and here you set the actual values
	config->SetMixingDepth(10);
	config->SetmTBins(mTBins);
	config->SetDomTMultBinning(true);
	config->SetmTBinning(true);
	config->SetMultiplicityEstimator(AliFemtoDreamEvent::kRef08);

	if (isMC)
	{
		config->SetMomentumResolution(true);
	}
	if (fullBlastQA)
	{
		// config->SetkTBinning(true);
		config->SetPtQA(true);
	}

	if (!fullBlastQA)
	{
		evtCuts->SetMinimalBooking(true);
		TrackCutsProtonDCA->SetMinimalBooking(true);
		TrackCutsAntiProtonDCA->SetMinimalBooking(true);
		TrackCutsDeuteronDCA->SetMinimalBooking(true);
		TrackCutsAntiDeuteronDCA->SetMinimalBooking(true);
		config->SetMinimalBookingME(true);
		config->SetMinimalBookingSample(true);
	}

	AliAnalysisTaskFemtoDreamDeuteron *task =
			new AliAnalysisTaskFemtoDreamDeuteron("FemtoDreamDefault", isMC);
	if (trigger == "kINT7")
	{
		task->SelectCollisionCandidates(AliVEvent::kINT7);
		std::cout << "Added kINT7 Trigger \n";
	}
	else if (trigger == "kHM")
	{
		task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
		std::cout << "Added kHighMult Trigger \n";
	}
	else
	{
		std::cout
				<< "====================================================================="
				<< std::endl;
		std::cout
				<< "====================================================================="
				<< std::endl;
		std::cout
				<< "Centrality Estimator not set, fix it else your Results will be empty!"
				<< std::endl;
		std::cout
				<< "====================================================================="
				<< std::endl;
		std::cout
				<< "====================================================================="
				<< std::endl;
	}
	if (!fullBlastQA)
	{
		task->SetRunTaskLightWeight(true);
	}
	task->SetEventCuts(evtCuts);
	task->SetTrackCutsDeuteronDCA(TrackCutsDeuteronDCA);
	task->SetTrackCutsAntiDeuteronDCA(TrackCutsAntiDeuteronDCA);
	task->SetTrackCutsProtonDCA(TrackCutsProtonDCA);
	task->SetTrackCutsAntiProtonDCA(TrackCutsAntiProtonDCA);
	task->SetCollectionConfig(config);

	mgr->AddTask(task);

	TString addon = "";

	if (trigger == "kINT7")
	{
		addon += "MB";
	}
	else if (trigger == "kHM")
	{
		addon += "HM";
	}

	TString file = AliAnalysisManager::GetCommonFileName();
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	mgr->ConnectInput(task, 0, cinput);

	TString EvtCutsName = Form("%sEvtCuts%s", addon.Data(), suffix.Data());
	AliAnalysisDataContainer *coutputEvtCuts = mgr->CreateContainer(
			EvtCutsName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), EvtCutsName.Data()));
	mgr->ConnectOutput(task, 1, coutputEvtCuts);

	TString TrackCutsName = Form("%sProtonDCA%s", addon.Data(), suffix.Data());
	AliAnalysisDataContainer *coutputTrkCuts = mgr->CreateContainer(
			TrackCutsName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), TrackCutsName.Data()));
	mgr->ConnectOutput(task, 2, coutputTrkCuts);

	TString AntiTrackCutsName = Form("%sAntiProtonDCA%s", addon.Data(), suffix.Data());
	AliAnalysisDataContainer *coutputAntiTrkCuts = mgr->CreateContainer(
			AntiTrackCutsName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
	mgr->ConnectOutput(task, 3, coutputAntiTrkCuts);

	TString TrackCutsDeuteronName = Form("%sDeuteronDCA%s", addon.Data(), suffix.Data());
	AliAnalysisDataContainer *coutputTrkCutsDeuteron = mgr->CreateContainer(
			TrackCutsDeuteronName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), TrackCutsDeuteronName.Data()));
	mgr->ConnectOutput(task, 4, coutputTrkCutsDeuteron);

	TString AntiTrackCutsDeuteronName = Form("%sAntiDeuteronDCA%s", addon.Data(), suffix.Data());
	AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron = mgr->CreateContainer(
			AntiTrackCutsDeuteronName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), AntiTrackCutsDeuteronName.Data()));
	mgr->ConnectOutput(task, 5, coutputAntiTrkCutsDeuteron);

	TString TrackCutsDeuteronNoTOFName = Form("%sDeuteronMass%s", addon.Data(), suffix.Data());
	AliAnalysisDataContainer *coutputTrkCutsDeuteronNoTOF = mgr->CreateContainer(
			TrackCutsDeuteronNoTOFName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), TrackCutsDeuteronNoTOFName.Data()));
	mgr->ConnectOutput(task, 6, coutputTrkCutsDeuteronNoTOF);

	TString AntiTrackCutsDeuteronNoTOFName = Form("%sAntiDeuteronMass%s", addon.Data(), suffix.Data());
	AliAnalysisDataContainer *coutputAntiTrkCutsDeuteronNoTOF = mgr->CreateContainer(
			AntiTrackCutsDeuteronNoTOFName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), AntiTrackCutsDeuteronNoTOFName.Data()));
	mgr->ConnectOutput(task, 7, coutputAntiTrkCutsDeuteronNoTOF);

	AliAnalysisDataContainer *coutputResults;
	TString ResultsName = Form("%sResults%s", addon.Data(), suffix.Data());
	coutputResults = mgr->CreateContainer(ResultsName.Data(),
																				TList::Class(), AliAnalysisManager::kOutputContainer,
																				Form("%s:%s", file.Data(), ResultsName.Data()));
	mgr->ConnectOutput(task, 8, coutputResults);

	AliAnalysisDataContainer *coutputResultsQA;
	TString ResultsQAName = Form("%sResultsQA%s", addon.Data(), suffix.Data());
	coutputResultsQA = mgr->CreateContainer(
			//@suppress("Invalid arguments") it works ffs
			ResultsQAName.Data(),
			TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), ResultsQAName.Data()));
	mgr->ConnectOutput(task, 9, coutputResultsQA);

	if (isMC)
	{
		AliAnalysisDataContainer *coutputTrkCutsMC;
		TString TrkCutsMCName = Form("%sProtonMC%s", addon.Data(), suffix.Data());
		coutputTrkCutsMC = mgr->CreateContainer(
				TrkCutsMCName.Data(),
				TList::Class(),
				AliAnalysisManager::kOutputContainer,
				Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
		mgr->ConnectOutput(task, 10, coutputTrkCutsMC);

		AliAnalysisDataContainer *coutputAntiTrkCutsMC;
		TString AntiTrkCutsMCName = Form("%sAntiProtonMC%s", addon.Data(),
																		 suffix.Data());
		coutputAntiTrkCutsMC = mgr->CreateContainer(
				//@suppress("Invalid arguments") it works ffs
				AntiTrkCutsMCName.Data(),
				TList::Class(),
				AliAnalysisManager::kOutputContainer,
				Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
		mgr->ConnectOutput(task, 11, coutputAntiTrkCutsMC);

		AliAnalysisDataContainer *coutputv0CutsMC;
		TString v0CutsMCName = Form("%sDeuteronMC%s", addon.Data(), suffix.Data());
		coutputv0CutsMC = mgr->CreateContainer(
				//@suppress("Invalid arguments") it works ffs
				v0CutsMCName.Data(),
				TList::Class(),
				AliAnalysisManager::kOutputContainer,
				Form("%s:%s", file.Data(), v0CutsMCName.Data()));
		mgr->ConnectOutput(task, 12, coutputv0CutsMC);

		AliAnalysisDataContainer *coutputAntiv0CutsMC;
		TString Antiv0CutsMCName = Form("%sAntiDeuteronMC%s", addon.Data(),
																		suffix.Data());
		coutputAntiv0CutsMC = mgr->CreateContainer(
				//@suppress("Invalid arguments") it works ffs
				Antiv0CutsMCName.Data(),
				TList::Class(),
				AliAnalysisManager::kOutputContainer,
				Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
		mgr->ConnectOutput(task, 13, coutputAntiv0CutsMC);
	}
	return task;
}
