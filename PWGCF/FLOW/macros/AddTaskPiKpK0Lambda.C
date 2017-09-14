AliAnalysisTaskPiKpK0Lamba* AddTaskPiKpK0Lamba(UInt_t filterbit = 768,
        Double_t etaCut = 0.8,
        Double_t vtxCut = 10.,
        Double_t minPt = 0.2,
        Double_t maxPt = 10.0,
        Int_t noclus = 70,
        Short_t nCentFl = 0,
        Double_t nHarm = 2.,
        Double_t nSigCut = 3.,
        Double_t minPiPtHCut = 1.,
        Double_t maxPiPtHCut = 7.,
        Double_t minPPtHCut = -25.,
        Double_t maxPPtHCut = -15.,
        Bool_t IsQPos = kFALSE,
        Bool_t IsQNeg = kFALSE,
        Bool_t IsPileUp = kTRUE,
        Bool_t IsPileUpTOF = kFALSE,
        Bool_t IsPhiCut = kFALSE,
        Bool_t IsQA = kFALSE,
        Bool_t IsExclPID = kFALSE,
        Int_t nocluspid = 70,
        Bool_t IsQAV0 = kFALSE,
        Bool_t IsQAOut = kFALSE,
        Bool_t IsBin1Cen = kFALSE,
        Bool_t IsRemch46V0A = kFALSE,
        Short_t flEtaRange = 2,
        Bool_t IsRemPhiReg = kFALSE,
        Float_t valCutMultESDdif = 15000.,
        Bool_t isCrRowShFrCls = kFALSE,
        Bool_t isFlagPsi24a = kFALSE,
        Int_t nPtb = 17,
        Float_t noCrFind = 0.8,
        Float_t dcaDghtPV = 0.1,
        Float_t maxDCADght = 0.5,
        Float_t cosPA = 0.998,
        Float_t minRad = 5.,
        Float_t maxRad = 100.,
        Bool_t isArmPodo = kFALSE,
        Bool_t minPtDght = kFALSE,
        TString uniqueID = "def")
{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskPiKpK0Lamba", "No analysis manager to connect to.");
        return NULL;
    }


    if (!mgr->GetInputEventHandler()) {
        Error("AddTaskPiKpK0Lamba", "This task requires an input event handler");
        return NULL;
    }


    Double_t ptBins[] = {0.2, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0};


    AliAnalysisTaskPiKpK0Lamba* taskD = new AliAnalysisTaskPiKpK0Lamba(Form("taskFlow_%s", uniqueID.Data()));

    taskD->SetDebugLevel(3);
    taskD->SetFilterbit(filterbit);
    taskD->SetNoClus(noclus);
    taskD->SetEtaCut(etaCut);
    taskD->SetVtxCut(vtxCut);
    taskD->SetMinPt(minPt);
    taskD->SetMaxPt(maxPt);
    taskD->SetFlagPileUp(IsPileUp);
    taskD->SetCentFlag(nCentFl);
    taskD->SetNHarmonic(nHarm);
    taskD->SetNsigCut(nSigCut);
    taskD->SetFlagQA(IsQA);
    taskD->SetFlagPileUpTOF(IsPileUpTOF);
    taskD->SetFlagPhiCut(IsPhiCut);
    taskD->SetMinPiPtHCut(minPiPtHCut);
    taskD->SetMaxPiPtHCut(maxPiPtHCut);
    taskD->SetMinPPtHCut(minPPtHCut);
    taskD->SetMaxPPtHCut(maxPPtHCut);
    taskD->SetFlagQPos(IsQPos);
    taskD->SetFlagQNeg(IsQNeg);
    taskD->SetFlagExclPID(IsExclPID);
    taskD->SetNoClusPID(nocluspid);
    taskD->SetFlagQAV0(IsQAV0);
    taskD->SetFlagQAOutl(IsQAOut);
    taskD->SetFlagBin1Cent(IsBin1Cen);
    taskD->SetFlagRemoveCh46V0A(IsRemch46V0A);
    taskD->SetFlagEtaRange(flEtaRange);
    taskD->SetRemovePhiReg(IsRemPhiReg);
    taskD->SetCutMultESDdif(valCutMultESDdif);
    taskD->SetCrsRowsFrcShClsCuts(isCrRowShFrCls);
    taskD->SetFlagPsi42A(isFlagPsi24a);
    taskD->SetNPtBins(nPtb);
    taskD->SetPtBins(ptBins);
    taskD->SetNcrFind(noCrFind);
    taskD->SetDCADghtPV(dcaDghtPV);
    taskD->SetMaxDCADght(maxDCADght);
    taskD->SetCosPA(cosPA);
    taskD->SetMinRad(minRad);
    taskD->SetMaxRad(maxRad);
    taskD->SetFlagArmPodCut(isArmPodo);
    taskD->SetFlagMinPtCutDght(minPtDght);
    mgr->AddTask(taskD);


    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer* cout = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("vn_SPV0_run2.root:%s", uniqueID.Data()));
    mgr->ConnectInput (taskD, 0, cinput);
    mgr->ConnectOutput(taskD, 1, cout);

    // Return task pointer at the end
    return taskD;

}
