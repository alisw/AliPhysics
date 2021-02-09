//#ifndef ADDTASKLIGHTN_C
//#define ADDTASKLIGHTN_C
//#include "AliLightNEventCuts.h"
//#include "AliLightNTrackCuts.h"
//#include "AliAnalysisTaskLightN.h"
//#include "TROOT.h"

#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskLightN.h"
#include "AliLightNTrackCuts.h"
#include "AliLightNEventCuts.h"
#include <TString.h>
#include <TList.h>
#include "TROOT.h"
#include "TSystem.h"
#endif


AliAnalysisTaskLightN* AddTaskLightNGrid(Bool_t isMC, TString CentEst, Bool_t DoSystematics) {

	gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
	gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
	gInterpreter->ProcessLine(".include $ROOTSYS/include");

	// the manager is static, so get the existing manager via the static method
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	if (!mgr){
		printf("No analysis manager to connect to!\n");
		return nullptr;
	}
	// just to see if all went well, check if the input event handler has been
	// connected
	if (!mgr->GetInputEventHandler()) {
		printf("This task requires an input event handler!\n");
		return nullptr;
	}
	AliLightNEventCuts *evtCutsParticle = AliLightNEventCuts::StandardCutsRun1();
	bool DCAPlots=true;
	bool CPAPlots=false;
	bool CombSigma=false;
	bool ContributionSplitting=false;
	bool ContributionSplittingDaug=false;
	if(isMC){
		ContributionSplitting=false;
		ContributionSplittingDaug=false;
	}else{
		ContributionSplitting=false;
		ContributionSplittingDaug=false;
	}

	//Track Cuts Proton
	AliLightNTrackCuts *ProtonTrackCuts = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	ProtonTrackCuts->SetCutCharge(1);
	AliLightNTrackCuts *AntiProtonTrackCuts = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	AntiProtonTrackCuts->SetCutCharge(-1);
	
	//Track Cuts Deuteron
	AliLightNTrackCuts *DeuteronTrackCuts = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	DeuteronTrackCuts->SetCutCharge(1);
	AliLightNTrackCuts *AntiDeuteronTrackCuts = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	AntiDeuteronTrackCuts->SetCutCharge(-1);
    
    //Systematic: track cuts (proton)
    AliLightNTrackCuts *ProtonTrackCuts_sysTrackT = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    ProtonTrackCuts_sysTrackT->SetCutCharge(1);
    ProtonTrackCuts_sysTrackT->SetNClsTPC(120);
    ProtonTrackCuts_sysTrackT->SetDCAVtxZ(0.1);
    ProtonTrackCuts_sysTrackT->SetDCAVtxXY(0.05);
    AliLightNTrackCuts *AntiProtonTrackCuts_sysTrackT = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    AntiProtonTrackCuts_sysTrackT->SetCutCharge(-1);
    AntiProtonTrackCuts_sysTrackT->SetNClsTPC(120);
    AntiProtonTrackCuts_sysTrackT->SetDCAVtxZ(0.1);
    AntiProtonTrackCuts_sysTrackT->SetDCAVtxXY(0.05);
    AliLightNTrackCuts *ProtonTrackCuts_sysTrackL = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    ProtonTrackCuts_sysTrackL->SetCutCharge(1);
    ProtonTrackCuts_sysTrackL->SetNClsTPC(50);
    ProtonTrackCuts_sysTrackL->SetDCAVtxZ(1.0);
    ProtonTrackCuts_sysTrackL->SetDCAVtxXY(0.5);
    AliLightNTrackCuts *AntiProtonTrackCuts_sysTrackL = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    AntiProtonTrackCuts_sysTrackL->SetCutCharge(-1);
    AntiProtonTrackCuts_sysTrackL->SetNClsTPC(50);
    AntiProtonTrackCuts_sysTrackL->SetDCAVtxZ(1.0);
    AntiProtonTrackCuts_sysTrackL->SetDCAVtxXY(0.5);
    

	//Systematic: track cuts (deuteron)
	AliLightNTrackCuts *DeuteronTrackCuts_sysTrackT = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	DeuteronTrackCuts_sysTrackT->SetCutCharge(1);
	DeuteronTrackCuts_sysTrackT->SetNClsTPC(120);
	DeuteronTrackCuts_sysTrackT->SetDCAVtxZ(0.2);
	DeuteronTrackCuts_sysTrackT->SetDCAVtxXY(0.1);
	AliLightNTrackCuts *AntiDeuteronTrackCuts_sysTrackT = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	AntiDeuteronTrackCuts_sysTrackT->SetCutCharge(-1);
	AntiDeuteronTrackCuts_sysTrackT->SetNClsTPC(120);
	AntiDeuteronTrackCuts_sysTrackT->SetDCAVtxZ(0.2);
	AntiDeuteronTrackCuts_sysTrackT->SetDCAVtxXY(0.1);
	AliLightNTrackCuts *DeuteronTrackCuts_sysTrackL = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	DeuteronTrackCuts_sysTrackL->SetCutCharge(1);
	DeuteronTrackCuts_sysTrackL->SetNClsTPC(50);
	DeuteronTrackCuts_sysTrackL->SetDCAVtxZ(2.0);
	DeuteronTrackCuts_sysTrackL->SetDCAVtxXY(0.5);
	AliLightNTrackCuts *AntiDeuteronTrackCuts_sysTrackL = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	AntiDeuteronTrackCuts_sysTrackL->SetCutCharge(-1);
	AntiDeuteronTrackCuts_sysTrackL->SetNClsTPC(50);
	AntiDeuteronTrackCuts_sysTrackL->SetDCAVtxZ(2.0);
	AntiDeuteronTrackCuts_sysTrackL->SetDCAVtxXY(0.5);
    
    //Systematic: PID cuts (proton)
    AliLightNTrackCuts *ProtonTrackCuts_sysPIDT = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    ProtonTrackCuts_sysPIDT->SetCutCharge(1);
    ProtonTrackCuts_sysPIDT->SetPID(AliPID::kProton, 0.7, 2.0,1e30);
    AliLightNTrackCuts *AntiProtonTrackCuts_sysPIDT = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    AntiProtonTrackCuts_sysPIDT->SetCutCharge(-1);
    AntiProtonTrackCuts_sysPIDT->SetPID(AliPID::kProton, 0.7, 2.0,1e30);
    AliLightNTrackCuts *ProtonTrackCuts_sysPIDL = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    ProtonTrackCuts_sysPIDL->SetCutCharge(1);
    ProtonTrackCuts_sysPIDL->SetPID(AliPID::kProton, 0.7, 5.0,1e30);
    AliLightNTrackCuts *AntiProtonTrackCuts_sysPIDL = AliLightNTrackCuts::PrimProtonCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
    AntiProtonTrackCuts_sysPIDL->SetCutCharge(-1);
    AntiProtonTrackCuts_sysPIDL->SetPID(AliPID::kProton, 0.7, 5.0,1e30);
    
	
	//Systematic: PID cuts (deuteron)
	AliLightNTrackCuts *DeuteronTrackCuts_sysPIDT = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	DeuteronTrackCuts_sysPIDT->SetCutCharge(1);
	DeuteronTrackCuts_sysPIDT->SetPID(AliPID::kDeuteron, 1.4, 2.0,1e30);
	AliLightNTrackCuts *AntiDeuteronTrackCuts_sysPIDT = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	AntiDeuteronTrackCuts_sysPIDT->SetCutCharge(-1);
	AntiDeuteronTrackCuts_sysPIDT->SetPID(AliPID::kDeuteron, 1.4, 2.0,1e30);
	AliLightNTrackCuts *DeuteronTrackCuts_sysPIDL = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	DeuteronTrackCuts_sysPIDL->SetCutCharge(1);
	DeuteronTrackCuts_sysPIDL->SetPID(AliPID::kDeuteron, 1.4, 5.0,1e30);
	AliLightNTrackCuts *AntiDeuteronTrackCuts_sysPIDL = AliLightNTrackCuts::PrimDeuteronCuts(isMC,DCAPlots,CombSigma,ContributionSplitting);
	AntiDeuteronTrackCuts_sysPIDL->SetCutCharge(-1);
	AntiDeuteronTrackCuts_sysPIDL->SetPID(AliPID::kDeuteron, 1.4, 5.0,1e30);


    AliAnalysisTaskLightN *task=new AliAnalysisTaskLightN("LightN",isMC);
    AliAnalysisTaskLightN *task1=new AliAnalysisTaskLightN("systematics1",isMC);
    AliAnalysisTaskLightN *task2=new AliAnalysisTaskLightN("systematics2",isMC);
    AliAnalysisTaskLightN *task3=new AliAnalysisTaskLightN("systematics3",isMC);
    AliAnalysisTaskLightN *task4=new AliAnalysisTaskLightN("systematics4",isMC);
    
    
    
	if(CentEst == "kInt7"){
		task->SelectCollisionCandidates(AliVEvent::kINT7);
        task->SetMVPileUp(kTRUE);
        if(DoSystematics){
            task1->SelectCollisionCandidates(AliVEvent::kINT7);
            task2->SelectCollisionCandidates(AliVEvent::kINT7);
            task3->SelectCollisionCandidates(AliVEvent::kINT7);
            task4->SelectCollisionCandidates(AliVEvent::kINT7);
            task1->SetMVPileUp(kTRUE);
            task2->SetMVPileUp(kTRUE);
            task3->SetMVPileUp(kTRUE);
            task4->SetMVPileUp(kTRUE);
        }
	}else if(CentEst == "kMB"){
		task->SelectCollisionCandidates(AliVEvent::kMB);
		task->SetMVPileUp(kFALSE);
    }else if(CentEst == "kHighMultV0"){
        task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
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
	task->SetTrackBufferSize(2500);
	task->SetEventCutsParticle(evtCutsParticle);
	task->SetTrackCutsProton(ProtonTrackCuts);
	task->SetAntiTrackCutsProton(AntiProtonTrackCuts);
	task->SetTrackCutsDeuteron(DeuteronTrackCuts);
	task->SetAntiTrackCutsDeuteron(AntiDeuteronTrackCuts);
	//Systematics
    if(DoSystematics){
        task1->SetDebugLevel(0);
        task1->SetEvtCutQA(true);
        task1->SetTrackBufferSize(2500);
        task1->SetEventCutsParticle(evtCutsParticle);
        task1->SetTrackCutsProton(ProtonTrackCuts);
        task1->SetAntiTrackCutsProton(AntiProtonTrackCuts);
        task1->SetTrackCutsDeuteron(DeuteronTrackCuts);
        task1->SetAntiTrackCutsDeuteron(AntiDeuteronTrackCuts);
        
        task2->SetDebugLevel(0);
        task2->SetEvtCutQA(true);
        task2->SetTrackBufferSize(2500);
        task2->SetEventCutsParticle(evtCutsParticle);
        task2->SetTrackCutsProton(ProtonTrackCuts_sysTrackL);
        task2->SetAntiTrackCutsProton(AntiProtonTrackCuts_sysTrackL);
        task2->SetTrackCutsDeuteron(DeuteronTrackCuts_sysTrackL);
        task2->SetAntiTrackCutsDeuteron(AntiDeuteronTrackCuts_sysTrackL);
        
        task3->SetDebugLevel(0);
        task3->SetEvtCutQA(true);
        task3->SetTrackBufferSize(2500);
        task3->SetEventCutsParticle(evtCutsParticle);
        task3->SetTrackCutsProton(ProtonTrackCuts_sysPIDT);
        task3->SetAntiTrackCutsProton(AntiProtonTrackCuts_sysPIDT);
        task3->SetTrackCutsDeuteron(DeuteronTrackCuts_sysPIDT);
        task3->SetAntiTrackCutsDeuteron(AntiDeuteronTrackCuts_sysPIDT);
        
        task4->SetDebugLevel(0);
        task4->SetEvtCutQA(true);
        task4->SetTrackBufferSize(2500);
        task4->SetEventCutsParticle(evtCutsParticle);
        task4->SetTrackCutsProton(ProtonTrackCuts_sysPIDL);
        task4->SetAntiTrackCutsProton(AntiProtonTrackCuts_sysPIDL);
        task4->SetTrackCutsDeuteron(DeuteronTrackCuts_sysPIDL);
        task4->SetAntiTrackCutsDeuteron(AntiDeuteronTrackCuts_sysPIDL);
    }

	mgr->AddTask(task);
    if(DoSystematics){
        mgr->AddTask(task1);
        mgr->AddTask(task2);
        mgr->AddTask(task3);
        mgr->AddTask(task4);
    }

	//  TString QAOutput="AnalysisQA.root";
	//  TString TrackOutput="AnalysisQATracks.root";
	//  TString v0Output="AnalysisQAv0.root";
	//  TString CascadeOutput="AnalysisQACascade.root";

	TString file = AliAnalysisManager::GetCommonFileName();

	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

	mgr->ConnectInput(task, 0, cinput);
    if(DoSystematics){
        mgr->ConnectInput(task1, 0, cinput);
        mgr->ConnectInput(task2, 0, cinput);
        mgr->ConnectInput(task3, 0, cinput);
        mgr->ConnectInput(task4, 0, cinput);
    }

	AliAnalysisDataContainer *coutputQA;
	TString QAName = Form("QA");
	coutputQA = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
		QAName.Data(), TList::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:%s", file.Data(), QAName.Data()));
	mgr->ConnectOutput(task, 1, coutputQA);

	AliAnalysisDataContainer *coutputEvtCuts;
	TString EvtCutsName = Form("EvtCuts");
	coutputEvtCuts = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
		EvtCutsName.Data(), TList::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:%s", file.Data(), EvtCutsName.Data()));
	mgr->ConnectOutput(task, 2, coutputEvtCuts);

	AliAnalysisDataContainer *couputTrkCutsProton;
	TString TrackCutsNameProton = Form("ProtonTrackCuts");
	couputTrkCutsProton = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
		TrackCutsNameProton.Data(), TList::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:%s", file.Data(), TrackCutsNameProton.Data()));
	mgr->ConnectOutput(task, 3, couputTrkCutsProton);

	AliAnalysisDataContainer *coutputAntiTrkCutsProton;
	TString AntiTrackCutsNameProton = Form("AntiProtonTrackCuts");
	coutputAntiTrkCutsProton = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
		AntiTrackCutsNameProton.Data(), TList::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:%s", file.Data(), AntiTrackCutsNameProton.Data()));
	mgr->ConnectOutput(task, 4, coutputAntiTrkCutsProton);

	AliAnalysisDataContainer *couputTrkCutsDeuteron;
	TString TrackCutsNameDeuteron = Form("DeuteronTrackCuts");
	couputTrkCutsDeuteron = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
		TrackCutsNameDeuteron.Data(), TList::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:%s", file.Data(), TrackCutsNameDeuteron.Data()));
	mgr->ConnectOutput(task, 5, couputTrkCutsDeuteron);

	AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron;
	TString AntiTrackCutsNameDeuteron = Form("AntiDeuteronTrackCuts");
	coutputAntiTrkCutsDeuteron = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
		AntiTrackCutsNameDeuteron.Data(), TList::Class(),
		AliAnalysisManager::kOutputContainer,
		Form("%s:%s", file.Data(), AntiTrackCutsNameDeuteron.Data()));
	mgr->ConnectOutput(task, 6, coutputAntiTrkCutsDeuteron);
    
    if(DoSystematics){
        AliAnalysisDataContainer *coutputTrkCutsProton_sysTrackT;
        TString TrackCutsNameProton_sysTrackT = Form("ProtonTrackCuts_sysTrackT");
        coutputTrkCutsProton_sysTrackT = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                              TrackCutsNameProton_sysTrackT.Data(), TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("%s:%s", file.Data(), TrackCutsNameProton_sysTrackT.Data()));
        mgr->ConnectOutput(task1, 7, coutputTrkCutsProton_sysTrackT);
        
        AliAnalysisDataContainer *coutputAntiTrkCutsProton_sysTrackT;
        TString AntiTrackCutsNameProton_sysTrackT = Form("AntiProtonTrackCuts_sysTrackT");
        coutputAntiTrkCutsProton_sysTrackT = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                                  AntiTrackCutsNameProton_sysTrackT.Data(), TList::Class(),
                                                                  AliAnalysisManager::kOutputContainer,
                                                                  Form("%s:%s", file.Data(), AntiTrackCutsNameProton_sysTrackT.Data()));
        mgr->ConnectOutput(task1, 8, coutputAntiTrkCutsProton_sysTrackT);
        
        AliAnalysisDataContainer *coutputTrkCutsDeuteron_sysTrackT;
        TString TrackCutsNameDeuteron_sysTrackT = Form("DeuteronTrackCuts_sysTrackT");
        coutputTrkCutsDeuteron_sysTrackT = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                                TrackCutsNameDeuteron_sysTrackT.Data(), TList::Class(),
                                                                AliAnalysisManager::kOutputContainer,
                                                                Form("%s:%s", file.Data(), TrackCutsNameDeuteron_sysTrackT.Data()));
        mgr->ConnectOutput(task1, 9, coutputTrkCutsDeuteron_sysTrackT);
        
        AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron_sysTrackT;
        TString AntiTrackCutsNameDeuteron_sysTrackT = Form("AntiDeuteronTrackCuts_sysTrackT");
        coutputAntiTrkCutsDeuteron_sysTrackT = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                                    AntiTrackCutsNameDeuteron_sysTrackT.Data(), TList::Class(),
                                                                    AliAnalysisManager::kOutputContainer,
                                                                    Form("%s:%s", file.Data(), AntiTrackCutsNameDeuteron_sysTrackT.Data()));
        mgr->ConnectOutput(task1, 10, coutputAntiTrkCutsDeuteron_sysTrackT);
        
        AliAnalysisDataContainer *coutputTrkCutsProton_sysTrackL;
        TString TrackCutsNameProton_sysTrackL = Form("ProtonTrackCuts_sysTrackL");
        coutputTrkCutsProton_sysTrackL = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                              TrackCutsNameProton_sysTrackL.Data(), TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("%s:%s", file.Data(), TrackCutsNameProton_sysTrackL.Data()));
        mgr->ConnectOutput(task2, 11, coutputTrkCutsProton_sysTrackL);
        
        AliAnalysisDataContainer *coutputAntiTrkCutsProton_sysTrackL;
        TString AntiTrackCutsNameProton_sysTrackL = Form("AntiProtonTrackCuts_sysTrackL");
        coutputAntiTrkCutsProton_sysTrackL = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                                  AntiTrackCutsNameProton_sysTrackL.Data(), TList::Class(),
                                                                  AliAnalysisManager::kOutputContainer,
                                                                  Form("%s:%s", file.Data(), AntiTrackCutsNameProton_sysTrackL.Data()));
        mgr->ConnectOutput(task2, 12, coutputAntiTrkCutsProton_sysTrackL);
        
        AliAnalysisDataContainer *coutputTrkCutsDeuteron_sysTrackL;
        TString TrackCutsNameDeuteron_sysTrackL = Form("DeuteronTrackCuts_sysTrackL");
        coutputTrkCutsDeuteron_sysTrackL = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                                TrackCutsNameDeuteron_sysTrackL.Data(), TList::Class(),
                                                                AliAnalysisManager::kOutputContainer,
                                                                Form("%s:%s", file.Data(), TrackCutsNameDeuteron_sysTrackL.Data()));
        mgr->ConnectOutput(task2, 13, coutputTrkCutsDeuteron_sysTrackL);
        
        AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron_sysTrackL;
        TString AntiTrackCutsNameDeuteron_sysTrackL = Form("AntiDeuteronTrackCuts_sysTrackL");
        coutputAntiTrkCutsDeuteron_sysTrackL = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                                    AntiTrackCutsNameDeuteron_sysTrackL.Data(), TList::Class(),
                                                                    AliAnalysisManager::kOutputContainer,
                                                                    Form("%s:%s", file.Data(), AntiTrackCutsNameDeuteron_sysTrackL.Data()));
        mgr->ConnectOutput(task2, 14, coutputAntiTrkCutsDeuteron_sysTrackL);
        
        AliAnalysisDataContainer *coutputTrkCutsProton_sysPIDT;
        TString TrackCutsNameProton_sysPIDT = Form("ProtonTrackCuts_sysPIDT");
        coutputTrkCutsProton_sysPIDT = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                            TrackCutsNameProton_sysPIDT.Data(), TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("%s:%s", file.Data(), TrackCutsNameProton_sysPIDT.Data()));
        mgr->ConnectOutput(task3, 15, coutputTrkCutsProton_sysPIDT);
        
        AliAnalysisDataContainer *coutputAntiTrkCutsProton_sysPIDT;
        TString AntiTrackCutsNameProton_sysPIDT = Form("AntiProtonTrackCuts_sysPIDT");
        coutputAntiTrkCutsProton_sysPIDT = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                                AntiTrackCutsNameProton_sysPIDT.Data(), TList::Class(),
                                                                AliAnalysisManager::kOutputContainer,
                                                                Form("%s:%s", file.Data(), AntiTrackCutsNameProton_sysPIDT.Data()));
        mgr->ConnectOutput(task3, 16, coutputAntiTrkCutsProton_sysPIDT);
        
        AliAnalysisDataContainer *coutputTrkCutsDeuteron_sysPIDT;
        TString TrackCutsNameDeuteron_sysPIDT = Form("DeuteronTrackCuts_sysPIDT");
        coutputTrkCutsDeuteron_sysPIDT = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                              TrackCutsNameDeuteron_sysPIDT.Data(), TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("%s:%s", file.Data(), TrackCutsNameDeuteron_sysPIDT.Data()));
        mgr->ConnectOutput(task3, 17, coutputTrkCutsDeuteron_sysPIDT);
        
        AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron_sysPIDT;
        TString AntiTrackCutsNameDeuteron_sysPIDT = Form("AntiDeuteronTrackCuts_sysPIDT");
        coutputAntiTrkCutsDeuteron_sysPIDT = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                                  AntiTrackCutsNameDeuteron_sysPIDT.Data(), TList::Class(),
                                                                  AliAnalysisManager::kOutputContainer,
                                                                  Form("%s:%s", file.Data(), AntiTrackCutsNameDeuteron_sysPIDT.Data()));
        mgr->ConnectOutput(task3, 18, coutputAntiTrkCutsDeuteron_sysPIDT);
        
        AliAnalysisDataContainer *coutputTrkCutsProton_sysPIDL;
        TString TrackCutsNameProton_sysPIDL = Form("ProtonTrackCuts_sysPIDL");
        coutputTrkCutsProton_sysPIDL = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                            TrackCutsNameProton_sysPIDL.Data(), TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("%s:%s", file.Data(), TrackCutsNameProton_sysPIDL.Data()));
        mgr->ConnectOutput(task4, 19, coutputTrkCutsProton_sysPIDL);
        
        AliAnalysisDataContainer *coutputAntiTrkCutsProton_sysPIDL;
        TString AntiTrackCutsNameProton_sysPIDL = Form("AntiProtonTrackCuts_sysPIDL");
        coutputAntiTrkCutsProton_sysPIDL = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                                AntiTrackCutsNameProton_sysPIDL.Data(), TList::Class(),
                                                                AliAnalysisManager::kOutputContainer,
                                                                Form("%s:%s", file.Data(), AntiTrackCutsNameProton_sysPIDL.Data()));
        mgr->ConnectOutput(task4, 20, coutputAntiTrkCutsProton_sysPIDL);
        
        AliAnalysisDataContainer *coutputTrkCutsDeuteron_sysPIDL;
        TString TrackCutsNameDeuteron_sysPIDL = Form("DeuteronTrackCuts_sysPIDL");
        coutputTrkCutsDeuteron_sysPIDL = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                              TrackCutsNameDeuteron_sysPIDL.Data(), TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("%s:%s", file.Data(), TrackCutsNameDeuteron_sysPIDL.Data()));
        mgr->ConnectOutput(task4, 21, coutputTrkCutsDeuteron_sysPIDL);
        
        AliAnalysisDataContainer *coutputAntiTrkCutsDeuteron_sysPIDL;
        TString AntiTrackCutsNameDeuteron_sysPIDL = Form("AntiDeuteronTrackCuts_sysPIDL");
        coutputAntiTrkCutsDeuteron_sysPIDL = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
                                                                  AntiTrackCutsNameDeuteron_sysPIDL.Data(), TList::Class(),
                                                                  AliAnalysisManager::kOutputContainer,
                                                                  Form("%s:%s", file.Data(), AntiTrackCutsNameDeuteron_sysPIDL.Data()));
        mgr->ConnectOutput(task4, 22, coutputAntiTrkCutsDeuteron_sysPIDL);
    }
    
    if (isMC) {
		AliAnalysisDataContainer *coutputTrkCutsMCProton;
		TString TrkCutsMCNameProton = Form("ProtonTrkCutsMC");
		coutputTrkCutsMCProton = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			TrkCutsMCNameProton.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), TrkCutsMCNameProton.Data()));
		mgr->ConnectOutput(task, 23, coutputTrkCutsMCProton);

		AliAnalysisDataContainer *coutputAntiTrkCutsMCProton;
		TString AntiTrkCutsMCNameProton = Form("AntiProtonTrkCutsMC");
		coutputAntiTrkCutsMCProton = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			AntiTrkCutsMCNameProton.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), AntiTrkCutsMCNameProton.Data()));
		mgr->ConnectOutput(task, 24, coutputAntiTrkCutsMCProton);

		AliAnalysisDataContainer *coutputTrkCutsMCDeuteron;
		TString TrkCutsMCNameDeuteron = Form("DeuteronTrkCutsMC");
		coutputTrkCutsMCDeuteron = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			TrkCutsMCNameDeuteron.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), TrkCutsMCNameDeuteron.Data()));
		mgr->ConnectOutput(task, 25, coutputTrkCutsMCDeuteron);

		AliAnalysisDataContainer *coutputAntiTrkCutsMCDeuteron;
		TString AntiTrkCutsMCNameDeuteron = Form("AntiDeuteronTrkCutsMC");
		coutputAntiTrkCutsMCDeuteron = mgr->CreateContainer(//@suppress("Invalid arguments") it works ffs
			AntiTrkCutsMCNameDeuteron.Data(), TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s:%s", file.Data(), AntiTrkCutsMCNameDeuteron.Data()));
		mgr->ConnectOutput(task, 26, coutputAntiTrkCutsMCDeuteron);

	
	}
	if (!mgr->InitAnalysis()) {
		return nullptr;
	}
	return task;
}
//#endif
