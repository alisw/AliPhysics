///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskMcKnoUe Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////

//AliAnalysisTaskGenMcKnoUe* AddTaskGenMcKnoUe(const Char_t* taskname="McKnoUe",Bool_t isPythia=kFALSE,Bool_t isEPOS=kFALSE,Bool_t isAnyGen=kFALSE, Double_t minpT=0.5)

AliAnalysisTask *AddTaskGenUeVsMult(Double_t minpT=0.5,TString suffixName ="")
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function

    AliAnalysisTaskGenUeVsMult* taskUE = new AliAnalysisTaskGenUeVsMult("taskKno");

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Printf("AliAnalysisTaskSimSpectraLF: No analysis manager to connect to.");
        return 0x0;
    }
    // get the input event handler this handler is part of the managing system and feeds events to your task
    if (!mgr->GetMCtruthEventHandler()) {
        Printf("AliAnalysisTaskSimSpectraLF: This task requires an input MC event handler.");
        return 0x0;
    }

    // now you create an instance of your task
    
    if(!taskUE) return 0x0;
   
    // add your task to the manager
    //taskUE -> SetGenerator("Pythia8");
    taskUE->SetPtMin(minpT);
    
    mgr->AddTask(taskUE);

    
    // Create containers for input/output

    TString finDirname    = "";
    TString inname    = "cinput";
    TString outBasic    = "cList";

    finDirname    += suffixName.Data();
    inname    += finDirname.Data();
    outBasic    += finDirname.Data();


    // Input and Output Slots
    //===========================================================================

    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += Form(":PWGMM_UeVsMult_%1.2f",minpT);

    AliAnalysisDataContainer *coutSim = mgr->CreateContainer(outBasic,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

    mgr->ConnectInput (taskUE, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskUE, 1, coutSim);

    return taskUE;
}
