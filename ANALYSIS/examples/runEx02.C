// run.C
//
// Template run macro for producer/consumer tasks
//
// Author: Andrei Gheata
//
void CreateAlienHandler();

//______________________________________________________________________________
void runEx02()
{
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
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");
  
    // add aliroot indlude path
    gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
    gROOT->SetStyle("Plain");
        
    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager("example for using producer/consumer tasks");
    mgr->SetCommonFileName("output.root");

    // create the alien handler and attach it to the manager
    AliAnalysisGrid *plugin = CreateAlienHandler(); 
    mgr->SetGridHandler(plugin);
    
    // create the alien handler and attach it to the manager
    
    AliVEventHandler* iH = new AliESDInputHandler();
    mgr->SetInputEventHandler(iH);        
                
    // create 2 producer tasks
    gROOT->LoadMacro("TaskExchange.cxx+g");
    gROOT->LoadMacro("AddTaskProducer.C");
    TaskProducer *prod1 = AddTaskProducer("producer1");
    TaskProducer *prod2 = AddTaskProducer("producer2");
    if (!prod1 || !prod2) return;
    
    // Create 1 consumer task
    gROOT->LoadMacro("AddTaskConsumer.C");
    TaskConsumer *task = AddTaskConsumer("consumer", "producer1", "producer2");
        
    // enable debug printouts
    mgr->SetDebugLevel(2);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
  
    // start analysis
    Printf("Starting analysis example for producer/consumers...");
    mgr->StartAnalysis("local",3);
}

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler()
{
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode("test");
    // Plugin test mode works only providing a file containing test file locations, used in "local" mode also
    plugin->SetFileForTestMode("files.txt"); // file should contain path name to a local directory containg *ESDs.root etc
    return plugin;
}
