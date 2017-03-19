#include "FILTER_COMPACTESDMUON.h"

#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMuonCompactTreeMaker.h"
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include <iostream>

namespace AAF {

int FILTER_COMPACTESDMUON(const char* from, const char* to)
{
    AliAnalysisManager mgr("COMPACTESDMUON");

    mgr.SetInputEventHandler(new AliESDInputHandler);

   mgr.SetMCtruthEventHandler(new AliMCEventHandler);

    AliMuonCompactTreeMaker task("raw://");

    TString destination(gSystem->BaseName(to));

    gSystem->ChangeDirectory(gSystem->DirName(to));

    AliAnalysisDataContainer* output = mgr.CreateContainer("compactevents",TTree::Class(),AliAnalysisManager::kOutputContainer,destination.Data());
    mgr.AddTask(&task);

    mgr.ConnectInput(&task,0,mgr.GetCommonInputContainer());
    mgr.ConnectOutput(&task,1,output);

    if (!mgr.InitAnalysis()) 
    {
        std::cout << "Could not InitAnalysis" << std::endl;
        return 1;
    }

    mgr.Print();

    TChain chain("esdTree");
    chain.Add(from);

    mgr.StartAnalysis("local",&chain);

    return 0;
}

}
