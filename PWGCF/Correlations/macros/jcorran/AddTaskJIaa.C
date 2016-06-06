//_____________________________________________________________________
AliAnalysisTask *AddTaskJIaa(TString taskName, TString cardName, TString jtrigg, TString jassoc, TString cardSetting, TString inclusFileName=""){

    // Load Custom Configuration and parameters
    // override values with parameters

    cout<<"### DEGUG Input is "<< cardName <<"\t"<<jtrigg<<"\t"<<jassoc<<"\t"<<inclusFileName<<"\t"<<"#########"<<endl;
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();


    //==== JCORRAN Efficiency TASK
    AliJDiHadronIaaTask *jiaatask = new AliJDiHadronIaaTask(taskName.Data(),"JOD");
    jiaatask->SetDebugLevel(5);
    jiaatask->SetFilterTaskName("PWGCFJCORRANTask");
    cout << jiaatask->GetName() << endl;


    // === Create AliJCORRAN ====
    AliJCard *card = new AliJCard(cardName.Data());
    card->PrintOut();
    card->ReadLine( cardSetting.Data() );
    card->ReCompile();
    card->PrintOut();

    AliJIaaAnalysis *fJIaa;
    fJIaa = new AliJIaaAnalysis( kFALSE );

    fJIaa->SetCard( card );
    fJIaa->SetTrigger( jtrigg.Data() );
    fJIaa->SetAssoc( jassoc.Data() );
    if( inclusFileName ) fJIaa->SetInclusiveFile(inclusFileName.Data());

    jiaatask->SetIaaAnalysis( fJIaa );

    mgr->AddTask((AliAnalysisTask*) jiaatask);

    // Create containers for input/output
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

    // Connect input/output
    mgr->ConnectInput(jiaatask, 0, cinput);
    AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",jiaatask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), jiaatask->GetName()));
    mgr->ConnectOutput(jiaatask, 1, jHist );

    return jiaatask;
}

