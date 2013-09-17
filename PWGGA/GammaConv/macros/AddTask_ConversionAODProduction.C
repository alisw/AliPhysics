AliAnalysisTask *AddTask_ConversionAODProduction(Int_t dataset=0, Bool_t isMC = kFALSE){

   // Before doing anything, we load the needed library
    gSystem->Load("libPWGGAGammaConv.so");
    // dataset 0: pp
    // dataset 1: PbPb
    // dataset 2: pPb

    //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return 0;
    }

   //========= Add PID Reponse to ANALYSIS manager ====
    if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
       gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
       AddTaskPIDResponse(isMC);
    }
    
    TString analysiscut;
    TString analysiscutB;

    if(dataset == 1){
     // Old cut string, no longer compatible with AliConversionCuts
     // analysiscut="900177009350113211200001000000000";
     // New cut string as of April 2013
       analysiscut= "100000006008400000150000000";
       analysiscutB="100000016008400000150000000";
    } else if (dataset == 2){
       analysiscut= "800000006008400000150000000";
       analysiscutB="800000016008400000150000000";
    } else{
      // analysiscut="
       analysiscut ="000000006008400100150000000";
       analysiscutB="000000016008400100150000000";
    }

    //========= Add V0 Reader to  ANALYSIS manager =====

    AliV0ReaderV1 *fV0Reader=new AliV0ReaderV1("ConvGammaAODProduction");
    fV0Reader->SetCreateAODs(kTRUE);
    fV0Reader->SetUseOwnXYZCalculation(kTRUE);
    fV0Reader->SetUseAODConversionPhoton(kTRUE);
//     fV0Reader->CheckAODConsistency();

    AliV0ReaderV1 *fV0ReaderB=new AliV0ReaderV1("ConvGammaAODProductionB");
    fV0ReaderB->SetCreateAODs(kTRUE);
    fV0ReaderB->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderB->SetUseAODConversionPhoton(kTRUE);
//     fV0ReaderB->CheckAODConsistency();

    // Set AnalysisCut Number
    AliConversionCuts *fCuts= new AliConversionCuts(analysiscut.Data(),analysiscut.Data());
    AliConversionCuts *fCutsB= new AliConversionCuts(analysiscutB.Data(),analysiscutB.Data());
    if(fCuts->InitializeCutsFromCutString(analysiscut.Data())){
      fV0Reader->SetConversionCuts(fCuts);
    }
    fV0Reader->Init();

    if(fCutsB->InitializeCutsFromCutString(analysiscutB.Data())){
      fV0ReaderB->SetConversionCuts(fCutsB);
    }
    fV0ReaderB->Init();

    
    AliLog::SetGlobalLogLevel(AliLog::kInfo);

    //================================================
    //              data containers
    //================================================
    //            find input container
    //below the trunk version
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

    // connect input V0Reader
    //fV0Reader->SelectCollisionCandidates(AliVEvent::kAny);
    mgr->AddTask(fV0Reader);
    mgr->ConnectInput (fV0Reader,0,cinput);

    //fV0ReaderB->SelectCollisionCandidates(AliVEvent::kAny);
    mgr->AddTask(fV0ReaderB);
    mgr->ConnectInput (fV0ReaderB,0,cinput);

    return fV0Reader;
}
