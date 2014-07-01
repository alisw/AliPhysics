void runTTreeFilter() {
    // author: Redmer Alexander Bertens, Utrecht University
    // rbertens@cern.ch , rbertens@nikhef.nl , r.a.bertens@uu.nl
    //
    // example which converts input data (in this case local aod's put into a chain)
    // to a tree which holds
    // - AliFlowTTreeEvent : event object
    // - AliFlowTTreeTrack : track objects
    // see source of these classes for more details

    // load libraries
    gSystem->Load("libCore.so");        
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");

    // create the analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager("MyManager");

    // create a tchain which will point to an aod tree
    TChain* chain = new TChain("aodTree");
    // add a few files to the chain (change this so that your local files are added)
    chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0003/AliAOD.root");
    chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0003/AliAOD.root");
    chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0004/AliAOD.root");
    chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0005/AliAOD.root");
    chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0006/AliAOD.root");
    chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0007/AliAOD.root");
    chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0008/AliAOD.root");
    chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0009/AliAOD.root");
    chain->Add("/home/rbertens/Documents/CERN/ALICE_DATA/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0010/AliAOD.root");
    // create an input handler
    AliVEventHandler* inputH = new AliAODInputHandler();
    // and connect it to the manager
    mgr->SetInputEventHandler(inputH);

     // the manager is static, so get the existing manager via the static method
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        printf("No analysis manager to connect to!\n");
        return NULL;
    }
        
    // just to see if all went well, check if the input event handler has been connected
    if (!mgr->GetInputEventHandler()) {
        printf("This task requires an input event handler!\n");
        return NULL;
      }

    // compile the relevant classes
      // include paths, necessary for compilation
    gSystem->AddIncludePath("-Wno-deprecated");
    gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/EMCAL");
    gSystem->AddIncludePath("-I$ALICE_ROOT/PWGDQ/dielectron -I$ALICE_ROOT/PWGHF/hfe");
    gSystem->AddIncludePath("-I$ALICE_ROOT/JETAN -I$ALICE_ROOT/JETAN/fastjet");

    gROOT->LoadMacro("../objects/AliFlowTTreeEvent.cxx+");
    gROOT->LoadMacro("../objects/AliFlowTTreeTrack.cxx+");
    gROOT->LoadMacro("../objects/AliAnalysisTaskTTreeFilter.cxx+");

    // load the addtask
    gROOT->LoadMacro("AddTaskTTreeFilter.C");

    // launch the task
    AddTaskTTreeFilter();

    // check if we can initialize the manager
    if(!mgr->InitAnalysis()) return;   
    // print the status of the manager to screen 
    mgr->PrintStatus();
    // print to screen how the analysis is progressing
    mgr->SetUseProgressBar(1, 25);
    // start the analysis locally, reading the events from the tchain
    mgr->StartAnalysis("local", chain);
}
