#include "FILTER_COMPACTESDMUON.h"

#include <iostream>
#include "AliAnalysisManager.h"
#include "AliMuonCompactTreeMaker.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisDataContainer.h"
#include "TChain.h"
#include "AliMCEventHandler.h"

namespace AAF {

int FILTER_COMPACTESDMUON(const char* from, const char* to)
{
    AliAnalysisManager mgr("COMPACTESDMUON");

    mgr.SetInputEventHandler(new AliESDInputHandler);

   mgr.SetMCtruthEventHandler(new AliMCEventHandler);

    AliMuonCompactTreeMaker task("raw://");

    AliAnalysisDataContainer* output = mgr.CreateContainer("compactevents",TTree::Class(),AliAnalysisManager::kOutputContainer,to);

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
