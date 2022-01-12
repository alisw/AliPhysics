// TODO LIST
// TODO: Set fName as last par per Train set-up indications
// TODO: Set Trigger Mask to be set from the AddAnalysisTask argument list

AliAnalysisTaskPhiCount* AddAnalysisTaskPhiCount( Bool_t MCFlag, Bool_t PhiFlag, Bool_t KaonFlag, TString fName = "PhiCount_STD", Int_t kFilterBit = 5, Float_t kVertexCut = 10., Float_t kDCAZcut = 2., Int_t kDCAXYcut = 7, Float_t kMinTPCclst = 70., Float_t kChi2TPCclst = 4., Float_t kChi2TPCglob = 36., Float_t kChi2ITSclst = 36.,  Float_t kTPCClsOverFndbl = .8, Float_t kSgTPC_Alone = 3., Float_t kSgTPC_TOFVt = 5., Float_t kSgTOF_Veto = 3. )
{
    // Analysis Manager
    AliAnalysisManager         *fAliAnlManager      =   AliAnalysisManager::GetAnalysisManager();
    
    // Filename
    TString fileName                                =   AliAnalysisManager::GetCommonFileName();
    
    // Task initialisation
    AliAnalysisTaskPhiCount    *fAliAnlTask         =   new AliAnalysisTaskPhiCount(fName.Data());
    
    // Checks
    if (!fAliAnlManager)                            return 0x0;
    if (!fAliAnlManager->GetInputEventHandler())    return 0x0;
    if (!fAliAnlTask)                               return 0x0;
    
    // Task options
    fAliAnlTask ->  SetMCFlag(MCFlag);
    fAliAnlTask ->  SetPhiFlag(PhiFlag);
    fAliAnlTask ->  SetKaonFlag(KaonFlag);
    fAliAnlTask ->  SetfRunName(fName);
    
    //  Standard Analysis Options
    
    //  -   //  Event Selection
    fAliAnlTask ->  SelectCollisionCandidates(AliVEvent::kAny);
    fAliAnlTask ->  SetkTriggerMask(AliVEvent::kAnyINT); // | AliVEvent::kHighMultV0
    fAliAnlTask ->  SetVertexCut(kVertexCut);
    fAliAnlTask ->  SetSPCompute(kFALSE);
    fAliAnlTask ->  SetRTCompute(kFALSE);
    
    //  -   //  Track Selection
    //  -   //  -   //  General
    fAliAnlTask ->  SetFilterBit(kFilterBit);               //  5   is Standard but it is custom made w/ differences
    fAliAnlTask ->  SetDCAzCut(kDCAZcut);                   //  2   is Standard
    fAliAnlTask ->  SetNSigmaPtDepXYDCA(kDCAXYcut);         //  7   is Standard
    
    //  -   //  -   //  TPC
    fAliAnlTask ->  SetMinTPCclusters(kMinTPCclst);         //  70  is Standard
    fAliAnlTask ->  SetChi2TPCcluster(kChi2TPCclst);        //  4   is Standard
    fAliAnlTask ->  SetChi2TPCGlobal(kChi2TPCglob);         //  36  is Standard
    fAliAnlTask ->  SetTPCClsOverFndbl(kTPCClsOverFndbl);   //  .8  is Standard
    
    //  -   //  -   //  ITS
    fAliAnlTask ->  SetChi2ITScluster(kChi2ITSclst);        //  36  is Standard
    
    //  -   //  -   //  PID
    fAliAnlTask ->  SetkSgTPC_Alone(kSgTPC_Alone);
    fAliAnlTask ->  SetkSgTPC_TOFVt(kSgTPC_TOFVt);
    fAliAnlTask ->  SetkSgTOF_Veto(kSgTOF_Veto);
    
    // Input
    fAliAnlManager->ConnectInput(fAliAnlTask,0,fAliAnlManager->GetCommonInputContainer());
    
    // Output
    // - // TLists
    
    fAliAnlManager->ConnectOutput(fAliAnlTask,1,fAliAnlManager->CreateContainer(Form("fAnalysisOutputList_%s",fName.Data()),                    TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    fAliAnlManager->ConnectOutput(fAliAnlTask,2,fAliAnlManager->CreateContainer(Form("fQCOutputList_%s",fName.Data()),                          TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    // - // TTrees
    if ( PhiFlag )          fAliAnlManager->ConnectOutput(fAliAnlTask,3,fAliAnlManager->CreateContainer(Form("PhiCandidate_%s",fName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    if ( KaonFlag )         fAliAnlManager->ConnectOutput(fAliAnlTask,4,fAliAnlManager->CreateContainer(Form("KaonCandidate_%s",fName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    if ( MCFlag && PhiFlag )fAliAnlManager->ConnectOutput(fAliAnlTask,5,fAliAnlManager->CreateContainer(Form("PhiEfficiency_%s",fName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    if ( MCFlag && KaonFlag)fAliAnlManager->ConnectOutput(fAliAnlTask,6,fAliAnlManager->CreateContainer(Form("KaonEfficiency_%s",fName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    // Add task
    fAliAnlManager->AddTask(fAliAnlTask);
    
    return fAliAnlTask;
}
