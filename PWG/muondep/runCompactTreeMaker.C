
//______________________________________________________________________________
TChain* CreateLocalChain(const char* filelist, const char* treeName)
{
  TChain* c = new TChain(treeName);
  
  char line[1024];
  
  ifstream in(filelist);
  while ( in.getline(line,1024,'\n') )
  {
    c->Add(line);
  }
  return c;
}

void runCompactTreeMaker(const char* inputFileName="list.esd.txt", const char* outputFileName="compacttreemaker.root")
{
    AliAnalysisManager mgr("COMPACTESDMUON");

    mgr.SetInputEventHandler(new AliESDInputHandler);

   mgr.SetMCtruthEventHandler(new AliMCEventHandler);

    //AliMuonCompactTreeMaker task("raw://");
    AliMuonCompactTreeMaker task("local:///alice/data/2015/OCDB");

    AliAnalysisDataContainer* output = mgr.CreateContainer("compactevents",TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);

    mgr.AddTask(&task);

    mgr.ConnectInput(&task,0,mgr.GetCommonInputContainer());
    mgr.ConnectOutput(&task,1,output);

    if (!mgr.InitAnalysis()) 
    {
        std::cout << "Could not InitAnalysis" << std::endl;
        return 1;
    }

    mgr.Print();

    TChain* chain = CreateLocalChain(inputFileName,"esdTree");

    TStopwatch timer;
    mgr.StartAnalysis("local",chain);
    timer.Print();
}
