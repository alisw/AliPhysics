#include "TROOT.h"
#include "TSystem.h"
AliAnalysisTaskSE* AddTaskFemtoDream(
    bool isMC=false,
    TString CentEst="kInt7",
    bool notpp=true,
    bool DCAPlots=false,
    bool CPAPlots=false,
    bool MomReso=false,
    bool etaPhiPlots=false,
    bool CombSigma=false,
    bool PileUpRej=true,
    bool ContributionSplitting=false,
    bool ContributionSplittingDaug=false)
{
	// the manager is static, so get the existing manager via the static method
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	if (!mgr)
	{
		printf("No analysis manager to connect to!\n");
		return nullptr;
	}
	// just to see if all went well, check if the input event handler has been
	// connected
	if (!mgr->GetInputEventHandler()) {
		printf("This task requires an input event handler!\n");
		return nullptr;
	}

	if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
	  if (isMC) {
	    // IMPORTANT - SET WHEN USING DIFFERENT PASS
	    AliAnalysisTaskPIDResponse *pidResponse =
	        reinterpret_cast<AliAnalysisTaskPIDResponse *>(
	            gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/"
	                "AddTaskPIDResponse.C (kTRUE, kTRUE, "
	                "kTRUE, \"1\")"));
	  } else {
	    AliAnalysisTaskPIDResponse *pidResponse =
	        reinterpret_cast<AliAnalysisTaskPIDResponse *>(
	            gInterpreter->ExecuteMacro(
	                "$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C)"));
	  }
	}


	AliFemtoDreamEventCuts *evtCuts=
			AliFemtoDreamEventCuts::StandardCutsRun2();
	evtCuts->CleanUpMult(false,false,false,true);

	//Track Cuts
	AliFemtoDreamTrackCuts *TrackCuts=
			AliFemtoDreamTrackCuts::PrimProtonCuts(
					isMC,DCAPlots,CombSigma,ContributionSplitting);
	TrackCuts->SetCutCharge(1);
	AliFemtoDreamTrackCuts *AntiTrackCuts=
			AliFemtoDreamTrackCuts::PrimProtonCuts(
					isMC,DCAPlots,CombSigma,ContributionSplitting);
	AntiTrackCuts->SetCutCharge(-1);
	AliFemtoDreamv0Cuts *v0Cuts;
	AliFemtoDreamv0Cuts *Antiv0Cuts;
	AliFemtoDreamCascadeCuts *CascadeCuts;
	AliFemtoDreamCascadeCuts *AntiCascadeCuts;

	//Lambda Cuts
	v0Cuts=
			AliFemtoDreamv0Cuts::LambdaCuts(isMC,CPAPlots,ContributionSplitting);
	AliFemtoDreamTrackCuts *Posv0Daug=AliFemtoDreamTrackCuts::DecayProtonCuts(
			isMC,PileUpRej,false);
	AliFemtoDreamTrackCuts *Negv0Daug=AliFemtoDreamTrackCuts::DecayPionCuts(
			isMC,PileUpRej,false);
	v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
	v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
	v0Cuts->SetPDGCodePosDaug(2212);//Proton
	v0Cuts->SetPDGCodeNegDaug(211);//Pion
	v0Cuts->SetPDGCodev0(3122);//Lambda
	Antiv0Cuts=
			AliFemtoDreamv0Cuts::LambdaCuts(isMC,CPAPlots,ContributionSplitting);
	AliFemtoDreamTrackCuts *PosAntiv0Daug=AliFemtoDreamTrackCuts::DecayPionCuts(
			isMC,PileUpRej,false);
	PosAntiv0Daug->SetCutCharge(1);
	AliFemtoDreamTrackCuts *NegAntiv0Daug=AliFemtoDreamTrackCuts::DecayProtonCuts(
			isMC,PileUpRej,false);
	NegAntiv0Daug->SetCutCharge(-1);
	Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
	Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
	Antiv0Cuts->SetPDGCodePosDaug(211);//Pion
	Antiv0Cuts->SetPDGCodeNegDaug(2212);//Proton
	Antiv0Cuts->SetPDGCodev0(-3122);//Lambda

	//Cascade Cuts
	CascadeCuts=
			AliFemtoDreamCascadeCuts::XiCuts(isMC,ContributionSplitting);
	CascadeCuts->SetXiCharge(-1);
	AliFemtoDreamTrackCuts *XiNegCuts=
			AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC,PileUpRej,false);
	AliFemtoDreamTrackCuts *XiPosCuts=
			AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC,PileUpRej,false);
	AliFemtoDreamTrackCuts *XiBachCuts=
			AliFemtoDreamTrackCuts::XiBachPionCuts(isMC,PileUpRej,false);
	CascadeCuts->Setv0Negcuts(XiNegCuts);
	CascadeCuts->Setv0PosCuts(XiPosCuts);
	CascadeCuts->SetBachCuts(XiBachCuts);
	CascadeCuts->SetPDGCodeCasc(-3312);
	CascadeCuts->SetPDGCodev0(3122);
	CascadeCuts->SetPDGCodePosDaug(2212);
	CascadeCuts->SetPDGCodeNegDaug(-211);
	CascadeCuts->SetPDGCodeBach(-211);

	AntiCascadeCuts=
			AliFemtoDreamCascadeCuts::XiCuts(isMC,ContributionSplitting);
	AntiCascadeCuts->SetXiCharge(1);
	AliFemtoDreamTrackCuts *AntiXiNegCuts=
			AliFemtoDreamTrackCuts::Xiv0ProtonCuts(isMC,PileUpRej,false);
	AntiXiNegCuts->SetCutCharge(-1);
	AliFemtoDreamTrackCuts *AntiXiPosCuts=
			AliFemtoDreamTrackCuts::Xiv0PionCuts(isMC,PileUpRej,false);
	AntiXiPosCuts->SetCutCharge(1);
	AliFemtoDreamTrackCuts *AntiXiBachCuts=
			AliFemtoDreamTrackCuts::XiBachPionCuts(isMC,PileUpRej,false);
	AntiXiBachCuts->SetCutCharge(1);
	AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
	AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
	AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
	AntiCascadeCuts->SetPDGCodeCasc(3312);
	AntiCascadeCuts->SetPDGCodev0(-3122);
	AntiCascadeCuts->SetPDGCodePosDaug(211);
	AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
	AntiCascadeCuts->SetPDGCodeBach(-211);

  //Thanks, CINT - will not compile due to an illegal constructor
  //std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3312);
  PDGParticles.push_back(3312);
  //std::vector<double> ZVtxBins = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
  std::vector<float> ZVtxBins;
  ZVtxBins.push_back(-10);
  ZVtxBins.push_back(-8);
  ZVtxBins.push_back(-6);
  ZVtxBins.push_back(-4);
  ZVtxBins.push_back(-2);
  ZVtxBins.push_back(0);
  ZVtxBins.push_back(2);
  ZVtxBins.push_back(4);
  ZVtxBins.push_back(6);
  ZVtxBins.push_back(8);
  ZVtxBins.push_back(10);
  //std::vector<int> NBins= {750,750,150,150,150,150,750,150,150,150,150,150,150,150,150,150,150,150,150,150,150};
  std::vector<int> NBins;
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(750);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  //std::vector<double> kMin= {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  std::vector<float> kMin;
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  //std::vector<double> kMax= {3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.};
  std::vector<float> kMax;
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  AliFemtoDreamCollConfig *config=new AliFemtoDreamCollConfig("Femto","Femto");
  if (notpp) {
    //std::vector<int> MultBins = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,80};
    std::vector<int> MultBins;
    MultBins.push_back(0);
    MultBins.push_back(4);
    MultBins.push_back(8);
    MultBins.push_back(12);
    MultBins.push_back(16);
    MultBins.push_back(20);
    MultBins.push_back(24);
    MultBins.push_back(28);
    MultBins.push_back(32);
    MultBins.push_back(36);
    MultBins.push_back(40);
    MultBins.push_back(44);
    MultBins.push_back(48);
    MultBins.push_back(52);
    MultBins.push_back(56);
    MultBins.push_back(60);
    MultBins.push_back(64);
    MultBins.push_back(68);
    MultBins.push_back(72);
    MultBins.push_back(76);
    MultBins.push_back(80);
    config->SetMultBins(MultBins);
  } else {
    //std::vector<int> MultBins = {0,4,8,12,16,20,24,28,32,36,40,60,80};
    std::vector<int> MultBins;
    MultBins.push_back(0);
    MultBins.push_back(4);
    MultBins.push_back(8);
    MultBins.push_back(12);
    MultBins.push_back(16);
    MultBins.push_back(20);
    MultBins.push_back(24);
    MultBins.push_back(28);
    MultBins.push_back(32);
    MultBins.push_back(36);
    MultBins.push_back(40);
    MultBins.push_back(60);
    MultBins.push_back(80);
    config->SetMultBins(MultBins);
  }
	config->SetMultBinning(true);
	config->SetZBins(ZVtxBins);
	if (MomReso) {
	  if (isMC) {
	    config->SetMomentumResolution(true);
	  } else {
	    std::cout << "You are trying to request the Momentum Resolution without MC Info; fix it wont work! \n";
	  }
 	}
	if (etaPhiPlots) {
	  if (isMC) {
	    config->SetPhiEtaBinnign(true);
	  } else {
	    std::cout << "You are trying to request the Eta Phi Plots without MC Info; fix it wont work! \n";
	  }
 	}
//	config->SetMultBins(MultBins);
	config->SetPDGCodes(PDGParticles);
	config->SetNBinsHist(NBins);
	config->SetMinKRel(kMin);
	config->SetMaxKRel(kMax);
	config->SetMixingDepth(10);

	AliAnalysisTaskFemtoDream *task=
	    new AliAnalysisTaskFemtoDream("FemtoDreamDefault",isMC,false);
	if(CentEst == "kInt7"){
		task->SelectCollisionCandidates(AliVEvent::kINT7);
		std::cout << "Added kINT7 Trigger \n";
		task->SetMVPileUp(kTRUE);
	}else if(CentEst == "kMB"){
		task->SelectCollisionCandidates(AliVEvent::kMB);
		std::cout << "Added kMB Trigger \n";
		task->SetMVPileUp(kFALSE);
	} else if (CentEst == "kHM") {
		task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
		std::cout << "Added kHighMultV0 Trigger \n";
    task->SetMVPileUp(kFALSE);
	}else{
		std::cout << "=====================================================================" << std::endl;
		std::cout << "=====================================================================" << std::endl;
		std::cout << "Centrality Estimator not set, fix it else your Results will be empty!" << std::endl;
		std::cout << "=====================================================================" << std::endl;
		std::cout << "=====================================================================" << std::endl;
	}
	task->SetDebugLevel(0);
	task->SetEvtCutQA(true);
	task->SetTrackBufferSize(10000);
	task->SetEventCuts(evtCuts);
	task->SetTrackCuts(TrackCuts);
	task->SetAntiTrackCuts(AntiTrackCuts);
	task->Setv0Cuts(v0Cuts);
	task->SetAntiv0Cuts(Antiv0Cuts);
	task->SetCascadeCuts(CascadeCuts);
	task->SetAntiCascadeCuts(AntiCascadeCuts);
	task->SetCollectionConfig(config);
	mgr->AddTask(task);
	//  TString QAOutput="AnalysisQA.root";
	//  TString TrackOutput="AnalysisQATracks.root";
	//  TString v0Output="AnalysisQAv0.root";
	//  TString CascadeOutput="AnalysisQACascade.root";

	TString file = AliAnalysisManager::GetCommonFileName();

	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

	mgr->ConnectInput(task, 0, cinput);

	AliAnalysisDataContainer *coutputQA;
	TString addon="";
	if (CentEst=="kInt7") {
	  addon+="MB";
	} else if (CentEst=="kHM") {
	  addon+="HM";
	}
	TString QAName = Form("%sQA",addon.Data());
	coutputQA = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			QAName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), QAName.Data()));
	mgr->ConnectOutput(task, 1, coutputQA);

	AliAnalysisDataContainer *coutputEvtCuts;
	TString EvtCutsName = Form("%sEvtCuts",addon.Data());
	coutputEvtCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			EvtCutsName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), EvtCutsName.Data()));
	mgr->ConnectOutput(task, 2, coutputEvtCuts);

	AliAnalysisDataContainer *couputTrkCuts;
	TString TrackCutsName = Form("%sTrackCuts",addon.Data());
	couputTrkCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			TrackCutsName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), TrackCutsName.Data()));
	mgr->ConnectOutput(task, 3, couputTrkCuts);

	AliAnalysisDataContainer *coutputAntiTrkCuts;
	TString AntiTrackCutsName = Form("%sAntiTrackCuts",addon.Data());
	coutputAntiTrkCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			AntiTrackCutsName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), AntiTrackCutsName.Data()));
	mgr->ConnectOutput(task, 4, coutputAntiTrkCuts);

	AliAnalysisDataContainer *coutputv0Cuts;
	TString v0CutsName = Form("%sv0Cuts",addon.Data());
	coutputv0Cuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			v0CutsName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), v0CutsName.Data()));
	mgr->ConnectOutput(task, 5, coutputv0Cuts);

	AliAnalysisDataContainer *coutputAntiv0Cuts;
	TString Antiv0CutsName = Form("%sAntiv0Cuts",addon.Data());
	coutputAntiv0Cuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			Antiv0CutsName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), Antiv0CutsName.Data()));
	mgr->ConnectOutput(task, 6, coutputAntiv0Cuts);

	AliAnalysisDataContainer *coutputCascadeCuts;
	TString CascadeCutsName = Form("%sCascadeCuts",addon.Data());
	coutputCascadeCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			CascadeCutsName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), CascadeCutsName.Data()));
	mgr->ConnectOutput(task, 7, coutputCascadeCuts);

	AliAnalysisDataContainer *coutputAntiCascadeCuts;
	TString AntiCascadeCutsName = Form("%sAntiCascadeCuts",addon.Data());
	coutputAntiCascadeCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			AntiCascadeCutsName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), AntiCascadeCutsName.Data()));
	mgr->ConnectOutput(task, 8, coutputAntiCascadeCuts);

	AliAnalysisDataContainer *coutputResults;
	TString ResultsName = Form("%sResults",addon.Data());
	coutputResults = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			ResultsName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), ResultsName.Data()));
	mgr->ConnectOutput(task, 9, coutputResults);

	AliAnalysisDataContainer *coutputResultQA;
	TString ResultQAName = Form("%sResultQA",addon.Data());
	coutputResultQA = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			ResultQAName.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), ResultQAName.Data()));
	mgr->ConnectOutput(task, 10, coutputResultQA);
	if (isMC) {
		AliAnalysisDataContainer *coutputTrkCutsMC;
		TString TrkCutsMCName = Form("%sTrkCutsMC",addon.Data());
		coutputTrkCutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
				TrkCutsMCName.Data(), TList::Class(),
				AliAnalysisManager::kOutputContainer,
				Form("%s:%s", file.Data(), TrkCutsMCName.Data()));
		mgr->ConnectOutput(task, 11, coutputTrkCutsMC);

		AliAnalysisDataContainer *coutputAntiTrkCutsMC;
		TString AntiTrkCutsMCName = Form("%sAntiTrkCutsMC",addon.Data());
		coutputAntiTrkCutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
				AntiTrkCutsMCName.Data(), TList::Class(),
				AliAnalysisManager::kOutputContainer,
				Form("%s:%s", file.Data(), AntiTrkCutsMCName.Data()));
		mgr->ConnectOutput(task, 12, coutputAntiTrkCutsMC);

		AliAnalysisDataContainer *coutputv0CutsMC;
		TString v0CutsMCName = Form("%sv0CutsMC",addon.Data());
		coutputv0CutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
				v0CutsMCName.Data(), TList::Class(),
				AliAnalysisManager::kOutputContainer,
				Form("%s:%s", file.Data(), v0CutsMCName.Data()));
		mgr->ConnectOutput(task, 13, coutputv0CutsMC);

		AliAnalysisDataContainer *coutputAntiv0CutsMC;
		TString Antiv0CutsMCName = Form("%sAntiv0CutsMC",addon.Data());
		coutputAntiv0CutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
				Antiv0CutsMCName.Data(), TList::Class(),
				AliAnalysisManager::kOutputContainer,
				Form("%s:%s", file.Data(), Antiv0CutsMCName.Data()));
		mgr->ConnectOutput(task, 14, coutputAntiv0CutsMC);

		AliAnalysisDataContainer *coutputXiCutsMC;
		TString XiCutsMCName = Form("%sXiCutsMC",addon.Data());
		coutputXiCutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
				XiCutsMCName.Data(), TList::Class(),
				AliAnalysisManager::kOutputContainer,
				Form("%s:%s", file.Data(), XiCutsMCName.Data()));
		mgr->ConnectOutput(task, 15, coutputXiCutsMC);

		AliAnalysisDataContainer *coutputAntiXiCutsMC;
		TString AntiXiCutsMCName = Form("%sAntiXiCutsMC",addon.Data());
		coutputAntiXiCutsMC = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
				AntiXiCutsMCName.Data(), TList::Class(),
				AliAnalysisManager::kOutputContainer,
				Form("%s:%s", file.Data(), AntiXiCutsMCName.Data()));
		mgr->ConnectOutput(task, 16, coutputAntiXiCutsMC);
	}
//	if (!mgr->InitAnalysis()) {
//		return nullptr;
//	}
	return task;
}


//  if (ContributionSplitting) {
//    v0Cuts=new AliFemtoDreamv0Cuts();
//    v0Cuts->SetAxisInvMassPlots(400,1.0, 1.2);
//    v0Cuts->SetIsMonteCarlo(isMC);
//    v0Cuts->SetPlotCPADist(false);
//    v0Cuts->SetPlotContrib(true);
//    AliFemtoDreamTrackCuts *Posv0Daug=
//        new AliFemtoDreamTrackCuts();
//    Posv0Daug->SetIsMonteCarlo(isMC);
//    Posv0Daug->SetPlotContrib(true);
//    AliFemtoDreamTrackCuts *Negv0Daug=
//        new AliFemtoDreamTrackCuts();
//    Negv0Daug->SetIsMonteCarlo(isMC);
//    Negv0Daug->SetPlotContrib(true);
//    v0Cuts->SetPosDaugterTrackCuts(Posv0Daug);
//    v0Cuts->SetNegDaugterTrackCuts(Negv0Daug);
//    v0Cuts->SetPDGCodePosDaug(2212);//Proton
//    v0Cuts->SetPDGCodeNegDaug(211);//Pion
//    v0Cuts->SetPDGCodev0(3122);//Lambda
//
//    Antiv0Cuts=new AliFemtoDreamv0Cuts();
//    Antiv0Cuts->SetAxisInvMassPlots(400,1.0, 1.2);
//    Antiv0Cuts->SetIsMonteCarlo(isMC);
//    Antiv0Cuts->SetPlotCPADist(false);
//    Antiv0Cuts->SetPlotContrib(true);
//    AliFemtoDreamTrackCuts *PosAntiv0Daug=
//        new AliFemtoDreamTrackCuts();
//    PosAntiv0Daug->SetCutCharge(1);
//    PosAntiv0Daug->SetIsMonteCarlo(isMC);
//    PosAntiv0Daug->SetPlotContrib(true);
//    AliFemtoDreamTrackCuts *NegAntiv0Daug=
//        new AliFemtoDreamTrackCuts();
//    NegAntiv0Daug->SetCutCharge(-1);
//    NegAntiv0Daug->SetIsMonteCarlo(isMC);
//    NegAntiv0Daug->SetPlotContrib(true);
//    Antiv0Cuts->SetPosDaugterTrackCuts(PosAntiv0Daug);
//    Antiv0Cuts->SetNegDaugterTrackCuts(NegAntiv0Daug);
//    Antiv0Cuts->SetPDGCodePosDaug(211);//Pion
//    Antiv0Cuts->SetPDGCodeNegDaug(2212);//Proton
//    Antiv0Cuts->SetPDGCodev0(-3122);//Lambda
//
//    CascadeCuts=
//        new AliFemtoDreamCascadeCuts();
//    CascadeCuts->SetIsMonteCarlo(isMC);
//    CascadeCuts->SetContributionSplitting(true);
//    CascadeCuts->SetXiCharge(-1);
//    CascadeCuts->SetXiMassRange(1.31486,5000);
//    AliFemtoDreamTrackCuts *XiNegCuts=
//        new AliFemtoDreamTrackCuts();
//    XiNegCuts->SetCutCharge(-1);
//    XiNegCuts->SetIsMonteCarlo(isMC);
//    XiNegCuts->SetPlotContrib(true);
//    XiNegCuts->SetCutCharge(-1);
//    AliFemtoDreamTrackCuts *XiPosCuts=
//        new AliFemtoDreamTrackCuts();
//    XiPosCuts->SetCutCharge(-1);
//    XiPosCuts->SetIsMonteCarlo(isMC);
//    XiPosCuts->SetPlotContrib(true);
//    XiPosCuts->SetCutCharge(1);
//    AliFemtoDreamTrackCuts *XiBachCuts=
//        new AliFemtoDreamTrackCuts();
//    XiBachCuts->SetCutCharge(-1);
//    XiBachCuts->SetIsMonteCarlo(isMC);
//    XiBachCuts->SetPlotContrib(true);
//    XiBachCuts->SetCutCharge(-1);
//    CascadeCuts->Setv0Negcuts(XiNegCuts);
//    CascadeCuts->Setv0PosCuts(XiPosCuts);
//    CascadeCuts->SetBachCuts(XiBachCuts);
//    CascadeCuts->SetPDGCodeCasc(-3312);
//    CascadeCuts->SetPDGCodev0(3122);
//    CascadeCuts->SetPDGCodePosDaug(2212);
//    CascadeCuts->SetPDGCodeNegDaug(-211);
//    CascadeCuts->SetPDGCodeBach(-211);
//
//    AntiCascadeCuts=
//        new AliFemtoDreamCascadeCuts();
//    AntiCascadeCuts->SetIsMonteCarlo(isMC);
//    AntiCascadeCuts->SetContributionSplitting(true);
//    AntiCascadeCuts->SetXiCharge(1);
//    AntiCascadeCuts->SetXiMassRange(1.31486,5000);
//    AliFemtoDreamTrackCuts *AntiXiNegCuts=
//        new AliFemtoDreamTrackCuts();
//    AntiXiNegCuts->SetCutCharge(-1);
//    AntiXiNegCuts->SetIsMonteCarlo(isMC);
//    AntiXiNegCuts->SetPlotContrib(true);
//    AntiXiNegCuts->SetCutCharge(-1);
//    AliFemtoDreamTrackCuts *AntiXiPosCuts=
//        new AliFemtoDreamTrackCuts();
//    AntiXiPosCuts->SetCutCharge(-1);
//    AntiXiPosCuts->SetIsMonteCarlo(isMC);
//    AntiXiPosCuts->SetPlotContrib(true);
//    AntiXiPosCuts->SetCutCharge(1);
//    AliFemtoDreamTrackCuts *AntiXiBachCuts=
//        new AliFemtoDreamTrackCuts();
//    AntiXiBachCuts->SetCutCharge(-1);
//    AntiXiBachCuts->SetIsMonteCarlo(isMC);
//    AntiXiBachCuts->SetPlotContrib(true);
//    AntiXiBachCuts->SetCutCharge(1);
//    AntiCascadeCuts->Setv0Negcuts(AntiXiNegCuts);
//    AntiCascadeCuts->Setv0PosCuts(AntiXiPosCuts);
//    AntiCascadeCuts->SetBachCuts(AntiXiBachCuts);
//    AntiCascadeCuts->SetPDGCodeCasc(3312);
//    AntiCascadeCuts->SetPDGCodev0(-3122);
//    AntiCascadeCuts->SetPDGCodePosDaug(211);
//    AntiCascadeCuts->SetPDGCodeNegDaug(-2212);
//    AntiCascadeCuts->SetPDGCodeBach(211);
//  } else {
