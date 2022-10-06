#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSELc2pKs0fromKFP.h"
#include <TString.h>
#include <TList.h>
#endif

AliAnalysisTaskSELc2pKs0fromKFP* AddTaskLc2pKs0FromKFParticle(TString finname="", Bool_t IsMC=kTRUE, TString cuttype="", Bool_t writeQATree=kTRUE, Bool_t IsAnaLc2Lpi=kFALSE, Bool_t useWeights = kFALSE, Bool_t keepOnlyMCSignal = kTRUE, Bool_t useMultiplicity = kFALSE, TString multProfiles = "", Int_t analysisType = AliAnalysisTaskSELc2pKs0fromKFP::kpPb2016, Double_t refMult = 29.2, Bool_t keepAllVariables = kFALSE, Bool_t useOnFlyV0 = kFALSE)
{
    Bool_t writeLcRecTree = kTRUE;
    Bool_t writeLcMCGenTree = kFALSE;
    if (IsMC) writeLcMCGenTree = kTRUE;

    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddMyTask", "No analysis manager to connect to.");
        return NULL;
    }

    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
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
      }
    }

    AliRDHFCutsKFP *RDHFCutsKFP = new AliRDHFCutsKFP();
    if (stdCuts) RDHFCutsKFP->SetStandardCutsPP2010();
    else {
      if (!IsAnaLc2Lpi) {
        RDHFCutsKFP = (AliRDHFCutsKFP*)fileCuts->Get("Lc2pKs0AnaCuts");
        RDHFCutsKFP->SetName("Lc2pKs0AnaCuts");
      }
      if (IsAnaLc2Lpi) {
        RDHFCutsKFP = (AliRDHFCutsKFP*)fileCuts->Get("Lc2LpiAnaCuts");
        RDHFCutsKFP->SetName("Lc2LpiAnaCuts");
      }
    }

    if (!RDHFCutsKFP) {
      cout << "Specific AliRDHFCutsKFP not found" << endl;
      return NULL;
    }

    printf("CREATE TASK\n");

    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":PWGHF_D2H_KFP_";      // create a subfolder in the file
    fileName += cuttype.Data();
    // now we create an instance of your task
    AliAnalysisTaskSELc2pKs0fromKFP* task = new AliAnalysisTaskSELc2pKs0fromKFP("AliAnalysisTaskSELc2pKs0fromKFP", RDHFCutsKFP);

    if(!task) return NULL;
    task->SetMC(IsMC);
    task->SetAnaLc2Lpi(IsAnaLc2Lpi);
    task->SetUseWeights(useWeights);
    task->SetDebugLevel(1);
    task->SetWriteLcMCGenTree(writeLcMCGenTree);
    task->SetWriteLcTree(writeLcRecTree);
    task->SetWriteLcQATree(writeQATree);
    task->SetKeepAllVariables(keepAllVariables);
    task->SetUseOnTheFlyV0(useOnFlyV0);
    // weight
    TF1 *weight = new TF1("weight", "expo", 0., 50.);
    weight->SetParameter(0, 0.853544);
    weight->SetParameter(1, -0.325586);
    task->SetWeightFunction(weight);


    if (useMultiplicity) {
      task->SetUseMult(useMultiplicity);
      TFile *fileEstimator=TFile::Open(multProfiles.Data());
      if (!fileEstimator) {
             Printf("FATAL: File with mult estimator not found");
             return NULL;
            }
      task->SetReferenceMultiplicity(refMult);
      const Char_t* profilebasename="SPDmult10";

      if (analysisType == AliAnalysisTaskSELc2pKs0fromKFP::kpPb2016) {
	    const Char_t* periodNames[4] = {"LHC16q_265499to265525_265309to265387", "LHC16q_265435","LHC16q_265388to265427","LHC16t_267163to267166"};
	    TProfile* multEstimatorAvg[4];
	    for(Int_t ip=0; ip<4; ip++) {
	    cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
	    multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
	    if (!multEstimatorAvg[ip]) {
	      Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
	      return NULL;
	    }
	  }
	  task->SetMultVsZProfileLHC16qt1stBunch(multEstimatorAvg[0]);
	  task->SetMultVsZProfileLHC16qt2ndBunch(multEstimatorAvg[1]);
	  task->SetMultVsZProfileLHC16qt3rdBunch(multEstimatorAvg[2]);
	  task->SetMultVsZProfileLHC16qt4thBunch(multEstimatorAvg[3]);
     task->SetAnalysisType(analysisType);

     }
    }

    // select type of event
//    task->SelectCollisionCandidates(AliVEvent::kAnyINT); // kAnyINT = kMB | kINT7 | kINT5 | kINT8 | kSPI7
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("CutsObj_%s", cuttype.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("Counter_%s", cuttype.Data()), AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("tree_event_char_%s", cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,4,mgr->CreateContainer(Form("tree_Lc_%s", cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,5,mgr->CreateContainer(Form("tree_Lc_MCGen_%s", cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,6,mgr->CreateContainer(Form("weight_%s", cuttype.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,7,mgr->CreateContainer(Form("tree_Lc_QA_%s", cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
