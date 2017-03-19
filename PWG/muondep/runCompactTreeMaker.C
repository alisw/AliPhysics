
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

void runCompactTreeMaker(const char* inputFileName="list.esd.txt", const char* outputFileName="compacttreemaker.root", const char* ocdbPath="raw://")
{
    AliAnalysisManager* mgr = new AliAnalysisManager("COMPACTESDMUON");

    mgr->SetInputEventHandler(new AliESDInputHandler);

    mgr->SetMCtruthEventHandler(new AliMCEventHandler);

    gROOT->LoadMacro("AddTaskCompactTreeMaker.C");

    AddTaskCompactTreeMaker(outputFileName,ocdbPath);

    if (!mgr->InitAnalysis()) 
    {
        std::cout << "Could not InitAnalysis" << std::endl;
        return 1;
    }

    mgr->Print();

    TChain* chain = CreateLocalChain(inputFileName,"esdTree");

    TStopwatch timer;
    mgr->StartAnalysis("local",chain);
    timer.Print();
}
