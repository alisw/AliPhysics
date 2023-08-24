#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSEOmegacZero2XiPifromKFP.h"
#include "AliRDHFCutsKFP.h"
#include <TString.h>
#include <TList.h>
#endif

AliAnalysisTaskSEOmegacZero2XiPifromKFP* AddTaskOmegacZero2XiPiFromKFParticle(TString finname="", Bool_t IsMC=kFALSE, TString cuttype="")
{
    //Settings
    TString particle = "Omegac0";

    // get the manager via the static access member. Since it's static, you don't need
    // an instance of the class to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddMyTask", "No analysis manager to connect to.");
        return NULL;
    }

    // get the input event handler, again via a static method.
    // This handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return NULL;
    }

    // check input cut object
    Bool_t stdCuts = kFALSE;
    TFile* fileCuts;
    if ( finname.EqualTo("") ) {
      stdCuts = kTRUE;
    } else {
      fileCuts = TFile::Open(finname.Data());
      if( !fileCuts || (fileCuts && !fileCuts->IsOpen()) ) {
        cout << "Input file not found : check your cut object" << endl;
        return NULL;
      }else{
        cout << "Apply CutObject selection" <<endl;
      }
    }

    AliRDHFCutsKFP *RDHFCutsKFP = new AliRDHFCutsKFP();
    if (stdCuts) RDHFCutsKFP->SetStandardCutsPP2010();
    else {
      RDHFCutsKFP = (AliRDHFCutsKFP*)fileCuts->Get("Omegac0AnaCuts");
    }
      RDHFCutsKFP->SetName("Omegac0AnaCuts");

    if (!RDHFCutsKFP) {
      cout << "Specific AliRDHFCutsKFP not found" << endl;
      return NULL;
    }

    printf("CREATE TASK\n");

    // by default, a file is open for writing. Here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":PWGHF_D2H_KFP_";      // create a subfolder in the file
    fileName += cuttype.Data();

    // now we create an instance of your task
    AliAnalysisTaskSEOmegacZero2XiPifromKFP* task = new AliAnalysisTaskSEOmegacZero2XiPifromKFP("AliAnalysisTaskSEOmegacZero2XiPifromKFP", RDHFCutsKFP, IsMC);

    if(!task) return NULL;

    task->SetDebugLevel(1);

    //add your task to the manager
    mgr->AddTask(task);

    //your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

    //same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("CutsObj_%s", cuttype.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("Counter_%s", cuttype.Data()), AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("Tree_%s_%s", particle.Data(),cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,4,mgr->CreateContainer("HistoList", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    //in the end, this macro returns a pointer to your task. This will be convenient later on
    //when you will run your analysis in an analysis train on grid
    return task;
}

//if everything worked fine you should get the file AnalysisResults.root
