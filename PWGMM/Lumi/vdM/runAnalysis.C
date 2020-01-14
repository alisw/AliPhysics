void runAnalysis() {
    // header location
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
    
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisMyTask");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
    
    // compile the class (locally) with debug symbols
    gInterpreter->LoadMacro("AliAnalysisTaskVdmStability.cxx++g");
    
    // load the addtask macro and create the task
    AliAnalysisTaskVdmStability *task = reinterpret_cast<AliAnalysisTaskVdmStability*>(gInterpreter->ExecuteMacro("AddTask_ilofnes_Vdm.C"));
    
    // if you want to run locally, we need to define some input
    TChain* chain = new TChain("aodTree");
    chain->Add("$HOME/ServiceTask/AliAOD.root");
    
    // start the analysis locally
    mgr->StartAnalysis("local", chain);
}
