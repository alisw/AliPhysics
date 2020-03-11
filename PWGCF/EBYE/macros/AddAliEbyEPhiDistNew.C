//=========================================================================//
//             AliEbyE Analysis for Net-Particle study                     //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                            Surya Prakash Pathak                         //
//                       surya.prakash.pathak@cern.ch                      //
//                         (Last Modified 2020/02/27)                      //
//                                                                         //
//Some parts of the code are taken from J. Thaeder/ M. Weber NetParticle analysis code//
//=========================================================================//

TString fileNameBase="AnalysisResults.root";

AliAnalysisTask *AddAliEbyEPhiDistNew(
                                      TString runName = "LHC10h",
                                      Bool_t isModeAOD = 0,
                                      Int_t aodFilterBit = 768,
                                      Bool_t IsMC  = 1,
                                      Bool_t IsQA = 0,
                                      Int_t pidtype = 1,//pidtype:charge-0, pion-1, kaon-2, proton-3
                                      Double_t DCAxy = 2.,
                                      Double_t DCAz = 2.,
                                      Double_t Vx = 0.3,
                                      Double_t Vy = 0.3,
                                      Double_t Vz = 10.,
                                      Double_t Eta = 0.8,
                                      Double_t Philow = 0.0,
                                      Double_t Phihigh = 1.04,
                                      Int_t TPCCrossRow = 80,
                                      Double_t Chi2NDF = 4.,
                                      const char* CentEstimator = "V0M",
                                      Int_t pidstrategy = 0,
                                      Float_t nSigmaITS = 2.,
                                      Float_t nSigmaTPC = 2.,
                                      Float_t nSigmaTOF = 3.,
                                      TString taskname = "TestHM") {
    
    Double_t ptl, pth;
    
    //Set the track cuts here ---
    if( pidtype == 0){
        ptl = 0.18;
        pth = 2.1;
    }
    else{
        ptl = 0.4;
        pth = 1.5;
    }
    
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddAliEbyEPhiDistNew", "No analysis manager to connect to.");
        return NULL;
    }
    
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddAliEbyEPhiDistNew", "This task requires an input event handler");
        return NULL;
    }
    
    const Char_t *pidname[4] = {"Nch","Npi","Nka","Npr"};
    
    const Char_t *ctsk = Form("%sNET%s",pidname[pidtype], taskname);
    
    AliEbyEPhiDistNew *task = new AliEbyEPhiDistNew(ctsk);
    
    
    task->SetRunPeriod(runName);
    
    if(isModeAOD) {
        task->SetIsAOD(isModeAOD);
        task->SetAODtrackCutBit(aodFilterBit);
    }
    
    task->SetIsMC(IsMC);
    task->RunQA(IsQA);
    
    task->SetDca( DCAxy,DCAz );
    task->SetVertexDiamond( Vx, Vy, Vz );
    task->SetKinematicsCuts( ptl, pth, Eta );
    task->SetPhi(Philow,Phihigh);
    if( pidtype == 0){
        task->SetNumberOfPtBins( 16 );
    }
    else{
        task->SetNumberOfPtBins( 8 );
    }
    task->SetTPCTrackQualityCuts( TPCCrossRow, Chi2NDF );
    task->SetCentralityEstimator( CentEstimator );
    
    //--------- PID -----
    task->SetPidType(pidtype);
    task->SetPidStrategy(pidstrategy);
    task->SetNSigmaMaxITS(nSigmaITS);
    task->SetNSigmaMaxTPC(nSigmaTPC);
    task->SetNSigmaMaxTOF(nSigmaTOF);
    
    if( runName == "LHC10h") task->SelectCollisionCandidates(AliVEvent::kMB);
    else if ( runName == "LHC15o") task->SelectCollisionCandidates(AliVEvent::kINT7);
    
    mgr->AddTask(task);
    
    
    //TString basefilename = AliAnalysisManager::GetCommonFileName();
    //TString basefilename = Form("AnalysisResult.%s.root",taskname.Data());
    
    AliAnalysisDataContainer *coutt   = mgr->CreateContainer(Form("Net%s%s",pidname[pidtype], taskname.Data()),
                                                             TList::Class(),AliAnalysisManager::kOutputContainer,
                                                             fileNameBase.Data() );
    mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutt);
    
    return task;
}


