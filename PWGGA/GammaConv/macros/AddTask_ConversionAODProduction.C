AliAnalysisTask *AddTask_ConversionAODProduction(Int_t dataset=1){

   // Before doing anything, we load the needed library
    gSystem->Load("libPWGGAGammaConv.so");
    // dataset 0: pp
    // dataset 1: PbPb

    Bool_t IsHeavyIon=(dataset==1);

    //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return 0;
    }

  
    TString analysiscut;

    if(IsHeavyIon){
     // Old cut string, no longer compatible with AliConversionCuts
     // analysiscut="900177009350113211200001000000000";
     // New cut string as of April 2013
      analysiscut="1000000032091071001000000";
    }
    else{
      // analysiscut="900397209450304221200000002000000";
      analysiscut="0000000002084000002200000";
    }

    //========= Add V0 Reader to  ANALYSIS manager =====

    AliV0ReaderV1 *fV0Reader=new AliV0ReaderV1("ConvGammaAODProduction");
    fV0Reader->SetCreateAODs(kTRUE);
    fV0Reader->SetUseOwnXYZCalculation(kTRUE);
    // Set AnalysisCut Number
    AliConversionCuts *fCuts= new AliConversionCuts(analysiscut.Data(),analysiscut.Data());
    if(fCuts->InitializeCutsFromCutString(analysiscut.Data())){
      fV0Reader->SetConversionCuts(fCuts);
    }
    fV0Reader->Init();
    mgr->AddTask(fV0Reader);

    AliLog::SetGlobalLogLevel(AliLog::kInfo);

    //================================================
    //              data containers
    //================================================
    //            find input container
    //below the trunk version
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

    // connect input V0Reader
    fV0Reader->SelectCollisionCandidates(AliVEvent::kAny);
    mgr->ConnectInput (fV0Reader,0,cinput);

    return fV0Reader;
}
