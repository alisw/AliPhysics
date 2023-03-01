//#if !defined (__CINT__) || defined (__CLING__)
//#include "AliAnalysisManager.h"
//#include "AliHFMLXicZeroToXiPifromKFP.h"
//#include <TString.h>
//#include <TList.h>
//#endif

AliHFMLXicZeroToXiPifromKFP* AddTaskXicZeroToXiPiFromKFParticleForPbPb (TString cutsfile="", // Cut object
                                                                        TString confFileML="", // ML configuration file
                                                                        Bool_t IsMC=kFALSE, // kFALSE: data; kTRUE: MC
                                                                        TString cuttype="", // Cut type: "std", "loose", "tight", "veryloose", "verytight", "veryverytight"
                                                                        Bool_t IsAnaOmegac0=kFALSE, // kFALSE: Xic0; kTRUE: Omegac0
                                                                        Bool_t IsPbPb=kTRUE, // kTRUE: Pb-Pb (enable ML output tree); kFALSE: pp and p-Pb
                                                                        Bool_t IsStoreOnlyMLoutput=kTRUE, // kTRUE: store only ML output tree; kFALSE: store tree with all variables and ML output tree
                                                                        Bool_t IsStoreLS=kFALSE, // kTRUE: store (Pi+ Xi+), (Pi- Xi-), (Pi+ Xi-) and (Pi- Xi+); kFALSE: store (Pi+ Xi-) and (Pi- Xi+)
                                                                        Int_t centmin=0, Int_t centmax=100, // centrality
                                                                        Int_t pTmin=0, Int_t pTmax=100) // pT
{

    Bool_t writeXic0RecTree = kTRUE;
    Bool_t writeXic0MCGenTree = kFALSE;
    if (IsMC) writeXic0MCGenTree = kTRUE;

    TString particle = "Xic0";
    if (IsAnaOmegac0) particle = "Omegac0";

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
    if ( cutsfile.EqualTo("") ) {
      stdCuts = kTRUE;
    } else {
      fileCuts = TFile::Open(cutsfile.Data());
      if( !fileCuts || (fileCuts && !fileCuts->IsOpen()) ) {
        cout << "Input file not found : check your cut object" << endl;
        return NULL;
      }
    }

    AliRDHFCutsKFP *RDHFCutsKFP = new AliRDHFCutsKFP();
    if (stdCuts) RDHFCutsKFP->SetStandardCutsPP2010();
    else {
      if (!IsAnaOmegac0) RDHFCutsKFP = (AliRDHFCutsKFP*)fileCuts->Get("Xic0AnaCuts");
      if (IsAnaOmegac0)  RDHFCutsKFP = (AliRDHFCutsKFP*)fileCuts->Get("Omegac0AnaCuts");
    }
    if (!IsAnaOmegac0) RDHFCutsKFP->SetName("Xic0AnaCuts");
    if (IsAnaOmegac0)  RDHFCutsKFP->SetName("Omegac0AnaCuts");

    if (!RDHFCutsKFP) {
      cout << "Specific AliRDHFCutsKFP not found" << endl;
      return NULL;
    }

    printf("CREATE TASK\n");

    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":PWGHF_D2H_KFP_";      // create a subfolder in the file
    fileName += Form("cent_%d_%d_pT_%d_%d_", centmin, centmax, pTmin, pTmax);
    fileName += cuttype.Data();
    // now we create an instance of your task
    AliHFMLXicZeroToXiPifromKFP* task = new AliHFMLXicZeroToXiPifromKFP("AliHFMLXicZeroToXiPifromKFP", RDHFCutsKFP);

    if(!task) return NULL;
    task->SetMC(IsMC);
    task->SetAnaOmegac0(IsAnaOmegac0);
    task->SetDebugLevel(1);
    task->SetWriteXic0MCGenTree(writeXic0MCGenTree);
    task->SetWriteXic0Tree(writeXic0RecTree);

    task->SetAnaPbPb(IsPbPb);
    task->SetMLConfigFile(confFileML);
    task->SetStoreOnlyMLoutput(IsStoreOnlyMLoutput);
    task->SetStoreLikeSign(IsStoreLS);

    /*
    // weight
    TF1 *weight = new TF1("weight", "expo", 0., 50.);
    TF1 *weight_up = new TF1("weight_up", "expo", 0., 50.);
    TF1 *weight_dw = new TF1("weight_dw", "expo", 0., 50.);
    // === PYTHIA 6 ===
//    weight->SetParameter(0, 0.853544);
//    weight->SetParameter(1, -0.325586);
    // === PYTHIA 8 + WeakDecayFinder ===
    if (!IsAnaOmegac0) {
      weight->SetParameter(0, 1.05904);
      weight->SetParameter(1, -0.380048);
      weight_up->SetParameter(0, 1.47631);
      weight_up->SetParameter(1, -0.474415);
      weight_dw->SetParameter(0, 0.570731);
      weight_dw->SetParameter(1, -0.279403);

    }
    if (IsAnaOmegac0) {
      weight->SetParameter(0, 1.67885);
      weight->SetParameter(1, -0.232748);
      weight_up->SetParameter(0, 2.06839);
      weight_up->SetParameter(1, -0.304407);
      weight_dw->SetParameter(0, 1.27061);
      weight_dw->SetParameter(1, -0.162865);
    }
    task->SetWeightFunction(weight);
    task->SetWeightFunctionUp(weight_up);
    task->SetWeightFunctionDw(weight_dw);
    */

    // select type of event
//    task->SelectCollisionCandidates(AliVEvent::kAnyINT); // kAnyINT = kMB | kINT7 | kINT5 | kINT8 | kSPI7
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("CutsObj_cent_%d_%d_pT_%d_%d_%s", centmin, centmax, pTmin, pTmax, cuttype.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("Counter_cent_%d_%d_pT_%d_%d_%s", centmin, centmax, pTmin, pTmax, cuttype.Data()), AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("tree_%s_cent_%d_%d_pT_%d_%d_%s", particle.Data(), centmin, centmax, pTmin, pTmax, cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,4,mgr->CreateContainer(Form("tree_%s_MCGen_cent_%d_%d_pT_%d_%d_%s", particle.Data(), centmin, centmax, pTmin, pTmax, cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,5,mgr->CreateContainer(Form("hist_%s_cent_%d_%d_pT_%d_%d_%s", particle.Data(), centmin, centmax, pTmin, pTmax, cuttype.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    if (IsPbPb) mgr->ConnectOutput(task,6,mgr->CreateContainer(Form("tree_ML_%s_cent_%d_%d_pT_%d_%d_%s", particle.Data(), centmin, centmax, pTmin, pTmax, cuttype.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
