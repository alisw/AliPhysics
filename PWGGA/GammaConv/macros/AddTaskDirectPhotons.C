void AddTaskDirectPhotons(){

    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return ;

    AliAnalysisTaskDirectPhotons* task = new AliAnalysisTaskDirectPhotons();

    AliAnalysisDataContainer *ic = mgr->GetCommonInputContainer(); //input container

    Int_t numberOfCuts = 1;

    TString EventCutnumber = "00000013"; //pp collision
    TString ConvCutnumber = "00000021810020000500000770";
    Bool_t enableV0findingEffi = kFALSE;
    TString periodNameV0Reader = "";

    TList* eventCutList = new TList();
    TList* conversionCutList = new TList();
    TList* neutralPionCutList = new TList();

    if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");

    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);

    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

    AliConvEventCuts* fCutsEv=NULL;

    if( EventCutnumber !=""){
       fCutsEv= new AliConvEventCuts(EventCutnumber.Data(),EventCutnumber.Data());
       fCutsEv->SetPreSelectionCutFlag(kFALSE);
       if(fCutsEv->InitializeCutsFromCutString(EventCutnumber.Data())){
       fV0ReaderV1->SetEventCuts(fCutsEv);
       }

       conversionCutList->Add(fCutsEv);
    }

    AliConversionPhotonCuts* fCuts=NULL;

    if( ConvCutnumber !=""){
       fCuts= new AliConversionPhotonCuts(ConvCutnumber.Data(),ConvCutnumber.Data());
       fCuts->SetPreSelectionCutFlag(kFALSE);
       if(fCuts->InitializeCutsFromCutString(ConvCutnumber.Data())){
       fV0ReaderV1->SetConversionCuts(fCuts);
       }

       eventCutList->Add(fCuts);
    }

    fV0ReaderV1->Init();

    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,ic);
    task->SetV0Reader(fV0ReaderV1);
    task->SetEventCutList(numberOfCuts, eventCutList);
    task->SetConversionCutList(conversionCutList);

    }

    AliAnalysisDataContainer *oc = mgr->CreateContainer("fOutputList", TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, ic);
    mgr->ConnectOutput(task, 1, oc);

    return task;

}


